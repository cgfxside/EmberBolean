// ═══════════════════════════════════════════════════════════════════════════════
// EMBER KigumiBackend — Boolean Operations via kigumi Library
// ═══════════════════════════════════════════════════════════════════════════════

#include "KigumiBackend.h"
#include "../ember/MeshImport.h"
#include <kigumi/mesh.h>
#include <kigumi/boolean_operation.h>
#include <vector>
#include <unordered_map>
#include <chrono>

namespace ember {
namespace {

//=============================================================================
// File-local types and helpers (not exposed in the header)
//=============================================================================

struct KVHash {
    size_t operator()(const std::array<int32_t,3>& v) const {
        size_t h = 14695981039346656037ULL;
        h ^= static_cast<size_t>(v[0]); h *= 1099511628211ULL;
        h ^= static_cast<size_t>(v[1]); h *= 1099511628211ULL;
        h ^= static_cast<size_t>(v[2]); h *= 1099511628211ULL;
        return h;
    }
};
struct KVEq {
    bool operator()(const std::array<int32_t,3>& a,
                    const std::array<int32_t,3>& b) const {
        return a[0]==b[0] && a[1]==b[1] && a[2]==b[2];
    }
};

void toKigumiFormat(const PolygonSoup& soup,
                    kigumi::Mesh& meshA,
                    kigumi::Mesh& meshB)
{
    std::unordered_map<std::array<int32_t,3>, uint32_t, KVHash, KVEq> vmapA, vmapB;
    std::vector<kigumi::Point3>           vertsA, vertsB;
    std::vector<std::array<uint32_t,3>>  facesA, facesB;

    for (const auto& tri : soup.triangles) {
        const bool isA = (tri.mesh_id == 0);
        auto& verts = isA ? vertsA : vertsB;
        auto& faces = isA ? facesA : facesB;
        auto& vmap  = isA ? vmapA  : vmapB;

        std::array<uint32_t,3> indices;
        for (int v = 0; v < 3; ++v) {
            std::array<int32_t,3> key = { tri.iv[v][0], tri.iv[v][1], tri.iv[v][2] };
            auto it = vmap.find(key);
            if (it == vmap.end()) {
                uint32_t idx = static_cast<uint32_t>(verts.size());
                vmap[key]   = idx;
                indices[v]  = idx;
                verts.emplace_back(
                    static_cast<double>(tri.iv[v][0]),
                    static_cast<double>(tri.iv[v][1]),
                    static_cast<double>(tri.iv[v][2]));
            } else {
                indices[v] = it->second;
            }
        }
        faces.push_back(indices);
    }

    meshA = kigumi::Mesh(vertsA, facesA);
    meshB = kigumi::Mesh(vertsB, facesB);
}

BackendResult fromKigumiResult(const kigumi::Mesh& mesh)
{
    BackendResult result;
    result.success       = true;
    result.num_vertices  = static_cast<uint32_t>(mesh.vertices().size());
    result.num_triangles = static_cast<uint32_t>(mesh.faces().size());

    result.output_soup.int_vertices.reserve(mesh.vertices().size());
    for (const auto& pos : mesh.vertices()) {
        result.output_soup.int_vertices.emplace_back(
            static_cast<int32_t>(pos.x()),
            static_cast<int32_t>(pos.y()),
            static_cast<int32_t>(pos.z()));
    }

    result.output_soup.output_polygons.reserve(mesh.faces().size());
    for (const auto& face : mesh.faces()) {
        OutputPolygon poly;
        poly.vertex_indices = {
            static_cast<uint32_t>(face[0]),
            static_cast<uint32_t>(face[1]),
            static_cast<uint32_t>(face[2])
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
// KigumiBackend::execute
//=============================================================================

BackendResult KigumiBackend::execute(const PolygonSoup& soup, BooleanOp op)
{
    BackendResult result;
    auto t0 = std::chrono::high_resolution_clock::now();

    kigumi::Mesh meshA, meshB;
    toKigumiFormat(soup, meshA, meshB);

    if (meshA.vertices().empty() || meshB.vertices().empty()) {
        result.success       = false;
        result.error_message = "Empty input mesh (A or B)";
        return result;
    }

    kigumi::Boolean_operation kigumi_op;
    bool reverse = false;

    switch (op) {
        case BooleanOp::Union:        kigumi_op = kigumi::Boolean_operation::UNION;               break;
        case BooleanOp::Intersection: kigumi_op = kigumi::Boolean_operation::INTERSECTION;        break;
        case BooleanOp::DiffAB:       kigumi_op = kigumi::Boolean_operation::DIFFERENCE;          break;
        case BooleanOp::DiffBA:       kigumi_op = kigumi::Boolean_operation::DIFFERENCE; reverse = true; break;
        case BooleanOp::Xor:          kigumi_op = kigumi::Boolean_operation::SYMMETRIC_DIFFERENCE; break;
        case BooleanOp::Shatter:      kigumi_op = kigumi::Boolean_operation::SHATTER;             break;
        default:
            result.success       = false;
            result.error_message = "Unsupported Boolean operation";
            return result;
    }

    try {
        kigumi::Mesh out = reverse
            ? kigumi::apply_boolean_operation(meshB, meshA, kigumi_op)
            : kigumi::apply_boolean_operation(meshA, meshB, kigumi_op);

        if (out.faces().empty()) {
            result.success       = true;
            result.num_vertices  = 0;
            result.num_triangles = 0;
        } else {
            result = fromKigumiResult(out);
        }
    } catch (const std::exception& e) {
        result.success       = false;
        result.error_message = std::string("kigumi Boolean failed: ") + e.what();
        return result;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    result.execution_time_ms = std::chrono::duration<double,std::milli>(t1-t0).count();
    return result;
}

} // namespace ember
