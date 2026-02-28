// ═══════════════════════════════════════════════════════════════════════════════
// EMBER MCUTBackend — Boolean Operations via MCUT Library
// ═══════════════════════════════════════════════════════════════════════════════

#include "MCUTBackend.h"
#include "../ember/MeshImport.h"
#include <mcut/mcut.h>
#include <vector>
#include <unordered_map>
#include <chrono>

namespace ember {

//=============================================================================
// Error helpers
//=============================================================================

static std::string mcutError(McResult r) {
    switch (r) {
        case MC_NO_ERROR:              return "no error";
        case MC_INVALID_VALUE:         return "invalid value";
        case MC_INVALID_OPERATION:     return "invalid operation";
        case MC_OUT_OF_MEMORY:         return "out of memory";
        case MC_NO_CONTEXT:            return "no context";
        case MC_INVALID_MESH_DEF:      return "invalid mesh definition";
        case MC_INVALID_VERTEX_ARRAY:  return "invalid vertex array";
        case MC_INVALID_FACE_ARRAY:    return "invalid face array";
        default:                       return "unknown error (" + std::to_string(r) + ")";
    }
}

static McFlags mcutDispatchFlags(BooleanOp op) {
    switch (op) {
        case BooleanOp::Union:        return MC_DISPATCH_FILTER_FRAGMENT_SEALING_INSIDE;
        case BooleanOp::Intersection: return MC_DISPATCH_FILTER_FRAGMENT_SEALING_OUTSIDE;
        case BooleanOp::DiffAB:       return MC_DISPATCH_FILTER_FRAGMENT_SEALING_OUTSIDE_EXHAUSTIVE;
        case BooleanOp::DiffBA:       return MC_DISPATCH_FILTER_FRAGMENT_SEALING_INSIDE_EXHAUSTIVE;
        case BooleanOp::Shatter:      return MC_DISPATCH_FILTER_ALL;
        default:                      return MC_DISPATCH_FILTER_FRAGMENT_SEALING_INSIDE;
    }
}

//=============================================================================
// PolygonSoup → flat MCUT arrays
//    Vertex array: [ x0,y0,z0, x1,y1,z1, ... ]  (double)
//    Face array:   [ v0,v1,v2, v0,v1,v2, ... ]   (uint32 triangle soup)
//    Face sizes:   [ 3,3,3, ... ]                 (McUint32)
//=============================================================================

void MCUTBackend::convertToMCUTFormat(
    const PolygonSoup& soup,
    std::vector<double>&   vA, std::vector<uint32_t>& fA,
    std::vector<double>&   vB, std::vector<uint32_t>& fB)
{
    // Vertex dedup maps: quantized key → flat index
    struct VHash {
        size_t operator()(const std::array<int32_t,3>& k) const {
            size_t h = 14695981039346656037ULL;
            h ^= static_cast<size_t>(k[0]); h *= 1099511628211ULL;
            h ^= static_cast<size_t>(k[1]); h *= 1099511628211ULL;
            h ^= static_cast<size_t>(k[2]); h *= 1099511628211ULL;
            return h;
        }
    };
    struct VEq {
        bool operator()(const std::array<int32_t,3>& a,
                        const std::array<int32_t,3>& b) const {
            return a[0]==b[0] && a[1]==b[1] && a[2]==b[2];
        }
    };

    std::unordered_map<std::array<int32_t,3>, uint32_t, VHash, VEq> mapA, mapB;

    for (const auto& tri : soup.triangles) {
        const bool isA = (tri.mesh_id == 0);
        auto& verts = isA ? vA : vB;
        auto& faces = isA ? fA : fB;
        auto& vmap  = isA ? mapA : mapB;

        for (int v = 0; v < 3; ++v) {
            std::array<int32_t,3> key = { tri.iv[v][0], tri.iv[v][1], tri.iv[v][2] };
            auto it = vmap.find(key);
            if (it == vmap.end()) {
                uint32_t idx = static_cast<uint32_t>(verts.size() / 3);
                vmap[key] = idx;
                verts.push_back(static_cast<double>(tri.iv[v][0]));
                verts.push_back(static_cast<double>(tri.iv[v][1]));
                verts.push_back(static_cast<double>(tri.iv[v][2]));
                faces.push_back(idx);
            } else {
                faces.push_back(it->second);
            }
        }
    }
}

//=============================================================================
// MCUT flat arrays → BackendResult (populates output_soup)
//=============================================================================

BackendResult MCUTBackend::convertFromMCUTResult(
    const std::vector<double>&   vertices,
    const std::vector<uint32_t>& faces,
    BooleanOp /*op*/)
{
    BackendResult result;
    result.success = true;

    // Populate int_vertices
    const size_t nv = vertices.size() / 3;
    result.output_soup.int_vertices.reserve(nv);
    for (size_t i = 0; i < nv; ++i) {
        result.output_soup.int_vertices.emplace_back(
            static_cast<int32_t>(vertices[i*3+0]),
            static_cast<int32_t>(vertices[i*3+1]),
            static_cast<int32_t>(vertices[i*3+2]));
    }
    result.num_vertices = static_cast<uint32_t>(nv);

    // faces is a flat triangle index array (3 indices per tri)
    const size_t nt = faces.size() / 3;
    result.output_soup.output_polygons.reserve(nt);
    for (size_t i = 0; i < nt; ++i) {
        OutputPolygon poly;
        poly.vertex_indices = {
            faces[i*3+0],
            faces[i*3+1],
            faces[i*3+2]
        };
        poly.mesh_id        = 0;
        poly.winding_number = 1;
        poly.is_interior    = false;
        poly.src_prim_idx   = 0;
        result.output_soup.output_polygons.push_back(std::move(poly));
    }
    result.num_triangles = static_cast<uint32_t>(nt);

    return result;
}

//=============================================================================
// Main execution
//=============================================================================

BackendResult MCUTBackend::execute(const PolygonSoup& soup, BooleanOp op)
{
    BackendResult result;
    auto t0 = std::chrono::high_resolution_clock::now();

    // Convert to MCUT format
    std::vector<double>   vA, vB;
    std::vector<uint32_t> fA, fB;
    convertToMCUTFormat(soup, vA, fA, vB, fB);

    if (vA.empty() || vB.empty()) {
        result.success = false;
        result.error_message = "Empty input mesh (A or B)";
        return result;
    }

    // Build face-size arrays (all triangles → all 3s)
    const McUint32 nfA = static_cast<McUint32>(fA.size() / 3);
    const McUint32 nfB = static_cast<McUint32>(fB.size() / 3);
    std::vector<McUint32> fsA(nfA, 3);
    std::vector<McUint32> fsB(nfB, 3);

    // Create context
    McContext ctx = MC_NULL_HANDLE;
    McResult  rc  = mcCreateContext(&ctx, MC_NULL_HANDLE);
    if (rc != MC_NO_ERROR) {
        result.success = false;
        result.error_message = "mcCreateContext failed: " + mcutError(rc);
        return result;
    }

    // Dispatch
    rc = mcDispatch(
        ctx,
        MC_DISPATCH_VERTEX_ARRAY_DOUBLE | MC_DISPATCH_ENFORCE_GENERAL_POSITION |
            mcutDispatchFlags(op),
        vA.data(), fA.data(), fsA.data(),
        static_cast<McUint32>(vA.size() / 3), nfA,
        vB.data(), fB.data(), fsB.data(),
        static_cast<McUint32>(vB.size() / 3), nfB);

    if (rc != MC_NO_ERROR) {
        result.success = false;
        result.error_message = "mcDispatch failed: " + mcutError(rc);
        mcReleaseContext(ctx);
        return result;
    }

    // Query connected components
    McUint32 numCC = 0;
    mcGetConnectedComponents(ctx, MC_CONNECTED_COMPONENT_TYPE_ALL, 0, nullptr, &numCC);

    if (numCC == 0) {
        // Valid empty result
        result.success = true;
        mcReleaseContext(ctx);
        auto t1 = std::chrono::high_resolution_clock::now();
        result.execution_time_ms =
            std::chrono::duration<double,std::milli>(t1-t0).count();
        return result;
    }

    std::vector<McConnectedComponent> ccs(numCC);
    mcGetConnectedComponents(ctx, MC_CONNECTED_COMPONENT_TYPE_ALL,
                             numCC, ccs.data(), nullptr);

    // Accumulate all CC geometry
    std::vector<double>   allVerts;
    std::vector<uint32_t> allFaces;

    for (McUint32 i = 0; i < numCC; ++i) {
        // Vertices
        uint64_t nBytes = 0;
        mcGetConnectedComponentData(ctx, ccs[i],
            MC_CONNECTED_COMPONENT_DATA_VERTEX_DOUBLE, 0, nullptr, &nBytes);
        if (nBytes == 0) continue;

        const uint32_t nv = static_cast<uint32_t>(nBytes / (3 * sizeof(double)));
        const uint32_t vOffset = static_cast<uint32_t>(allVerts.size() / 3);

        allVerts.resize(allVerts.size() + nv * 3);
        mcGetConnectedComponentData(ctx, ccs[i],
            MC_CONNECTED_COMPONENT_DATA_VERTEX_DOUBLE,
            nBytes, allVerts.data() + vOffset * 3, nullptr);

        // Triangulated faces
        nBytes = 0;
        mcGetConnectedComponentData(ctx, ccs[i],
            MC_CONNECTED_COMPONENT_DATA_FACE_TRIANGULATION, 0, nullptr, &nBytes);
        if (nBytes == 0) continue;

        const uint32_t ni = static_cast<uint32_t>(nBytes / sizeof(uint32_t));
        const size_t   fBase = allFaces.size();
        allFaces.resize(fBase + ni);
        mcGetConnectedComponentData(ctx, ccs[i],
            MC_CONNECTED_COMPONENT_DATA_FACE_TRIANGULATION,
            nBytes, allFaces.data() + fBase, nullptr);

        // Remap vertex indices for this CC to global offset
        for (size_t j = fBase; j < allFaces.size(); ++j) {
            allFaces[j] += vOffset;
        }
    }

    for (auto& cc : ccs) mcReleaseConnectedComponent(ctx, cc);
    mcReleaseContext(ctx);

    result = convertFromMCUTResult(allVerts, allFaces, op);

    auto t1 = std::chrono::high_resolution_clock::now();
    result.execution_time_ms =
        std::chrono::duration<double,std::milli>(t1-t0).count();
    return result;
}

} // namespace ember
