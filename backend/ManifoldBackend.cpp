// ═══════════════════════════════════════════════════════════════════════════════
// EMBER ManifoldBackend — Boolean Operations via Manifold Library
// ═══════════════════════════════════════════════════════════════════════════════

#include "ManifoldBackend.h"
#include "../ember/MeshImport.h"
#include <manifold/manifold.h>
#include <glm/glm.hpp>
#include <vector>
#include <unordered_map>
#include <chrono>

namespace ember {
namespace {

//=============================================================================
// File-local types and helpers (not exposed in the header)
//=============================================================================

struct ManifoldMesh {
    std::vector<glm::vec3>  vertices;
    std::vector<glm::ivec3> triangles;
};

struct VHash {
    size_t operator()(const std::array<int32_t,3>& v) const {
        size_t h = 14695981039346656037ULL;
        h ^= static_cast<size_t>(v[0]); h *= 1099511628211ULL;
        h ^= static_cast<size_t>(v[1]); h *= 1099511628211ULL;
        h ^= static_cast<size_t>(v[2]); h *= 1099511628211ULL;
        return h;
    }
};
struct VEq {
    bool operator()(const std::array<int32_t,3>& a,
                    const std::array<int32_t,3>& b) const {
        return a[0]==b[0] && a[1]==b[1] && a[2]==b[2];
    }
};

ManifoldMesh toManifoldFormat(const PolygonSoup& soup, uint32_t mesh_id)
{
    ManifoldMesh result;
    std::unordered_map<std::array<int32_t,3>, uint32_t, VHash, VEq> vmap;

    for (const auto& tri : soup.triangles) {
        if (static_cast<uint32_t>(tri.mesh_id) != mesh_id) continue;

        glm::ivec3 indices;
        for (int v = 0; v < 3; ++v) {
            std::array<int32_t,3> key = { tri.iv[v][0], tri.iv[v][1], tri.iv[v][2] };
            auto it = vmap.find(key);
            if (it == vmap.end()) {
                uint32_t idx = static_cast<uint32_t>(result.vertices.size());
                vmap[key]  = idx;
                indices[v] = static_cast<int>(idx);
                result.vertices.emplace_back(
                    static_cast<float>(tri.iv[v][0]),
                    static_cast<float>(tri.iv[v][1]),
                    static_cast<float>(tri.iv[v][2]));
            } else {
                indices[v] = static_cast<int>(it->second);
            }
        }
        result.triangles.push_back(indices);
    }
    return result;
}

BackendResult fromManifoldResult(const manifold::Mesh& mesh)
{
    BackendResult result;
    result.success       = true;
    result.num_vertices  = static_cast<uint32_t>(mesh.vertPos.size());
    result.num_triangles = static_cast<uint32_t>(mesh.triVerts.size());

    result.output_soup.int_vertices.reserve(mesh.vertPos.size());
    for (const auto& pos : mesh.vertPos) {
        result.output_soup.int_vertices.emplace_back(
            static_cast<int32_t>(pos.x),
            static_cast<int32_t>(pos.y),
            static_cast<int32_t>(pos.z));
    }

    result.output_soup.output_polygons.reserve(mesh.triVerts.size());
    for (const auto& tri : mesh.triVerts) {
        OutputPolygon poly;
        poly.vertex_indices = {
            static_cast<uint32_t>(tri[0]),
            static_cast<uint32_t>(tri[1]),
            static_cast<uint32_t>(tri[2])
        };
        poly.mesh_id        = 0;
        poly.winding_number = 1;
        poly.is_interior    = false;
        poly.src_prim_idx   = 0;
        result.output_soup.output_polygons.push_back(std::move(poly));
    }
    return result;
}

} // anonymous namespace

//=============================================================================
// ManifoldBackend::execute
//=============================================================================

BackendResult ManifoldBackend::execute(const PolygonSoup& soup, BooleanOp op)
{
    BackendResult result;
    auto t0 = std::chrono::high_resolution_clock::now();

    ManifoldMesh dataA = toManifoldFormat(soup, 0);
    ManifoldMesh dataB = toManifoldFormat(soup, 1);

    if (dataA.vertices.empty() || dataB.vertices.empty()) {
        result.success       = false;
        result.error_message = "Empty input mesh (A or B)";
        return result;
    }

    manifold::Mesh meshA, meshB;
    meshA.vertPos  = std::move(dataA.vertices);
    meshA.triVerts = std::move(dataA.triangles);
    meshB.vertPos  = std::move(dataB.vertices);
    meshB.triVerts = std::move(dataB.triangles);

    manifold::Manifold mA, mB;
    try {
        mA = manifold::Manifold(meshA);
        mB = manifold::Manifold(meshB);
    } catch (const std::exception& e) {
        result.success       = false;
        result.error_message = std::string("Failed to create Manifold: ") + e.what();
        return result;
    }

    if (mA.Status() != manifold::Manifold::Error::NO_ERROR) {
        result.success = false; result.error_message = "Mesh A is not a valid manifold"; return result;
    }
    if (mB.Status() != manifold::Manifold::Error::NO_ERROR) {
        result.success = false; result.error_message = "Mesh B is not a valid manifold"; return result;
    }

    manifold::Manifold out;
    try {
        switch (op) {
            case BooleanOp::Union:        out = mA + mB;                       break;
            case BooleanOp::Intersection: out = mA ^ mB;                       break;
            case BooleanOp::DiffAB:       out = mA - mB;                       break;
            case BooleanOp::DiffBA:       out = mB - mA;                       break;
            case BooleanOp::Xor:          out = mA ^ mB;                       break;
            case BooleanOp::Shatter:      out = (mA^mB) + (mA-mB) + (mB-mA);  break;
            default:
                result.success = false;
                result.error_message = "Unsupported Boolean operation";
                return result;
        }
    } catch (const std::exception& e) {
        result.success       = false;
        result.error_message = std::string("Boolean operation failed: ") + e.what();
        return result;
    }

    if (out.Status() != manifold::Manifold::Error::NO_ERROR) {
        result.success = false; result.error_message = "Boolean produced invalid result"; return result;
    }

    manifold::Mesh outMesh = out.ToMesh();
    if (outMesh.vertPos.empty()) {
        result.success = true;
        result.num_vertices = result.num_triangles = 0;
    } else {
        result = fromManifoldResult(outMesh);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    result.execution_time_ms = std::chrono::duration<double,std::milli>(t1-t0).count();
    return result;
}

} // namespace ember
