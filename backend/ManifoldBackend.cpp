// ═══════════════════════════════════════════════════════════════════════════════
// EMBER ManifoldBackend.cpp — Boolean Operations via elalish/manifold v3.x
// ═══════════════════════════════════════════════════════════════════════════════
//
// v3 API: MeshGL (no GLM dependency).
//
// WINDING FIX: Houdini's polygon vertex order produces CW triangles (viewed
// from outside), but Manifold expects CCW. We reverse the winding on INPUT
// (swap v1↔v2 in triVerts) so Manifold sees correct outward-pointing normals,
// then reverse the OUTPUT winding back to Houdini convention.
//
// ═══════════════════════════════════════════════════════════════════════════════

#include "ManifoldBackend.h"

#include <manifold/manifold.h>

#include <vector>
#include <unordered_map>
#include <array>
#include <cstring>
#include <chrono>

namespace ember {

// ─────────────────────────────────────────────────────────────────────────────
// Name / description / capabilities
// ─────────────────────────────────────────────────────────────────────────────

std::string ManifoldBackend::name()        const { return "manifold"; }
std::string ManifoldBackend::description() const { return "Manifold (guaranteed watertight)"; }
bool ManifoldBackend::supportsOpenMesh()   const { return false; }
bool ManifoldBackend::providesExactTopology() const { return false; }

namespace {

// ─────────────────────────────────────────────────────────────────────────────
// Constants
// ─────────────────────────────────────────────────────────────────────────────

static constexpr uint32_t CAP_FLAG = 0x80000000u;

// ─────────────────────────────────────────────────────────────────────────────
// §1  PolygonSoup → MeshGL  (per mesh_id)
// ─────────────────────────────────────────────────────────────────────────────
//
// Deduplicates vertices by exact float-bit equality (FNV-1a on raw bits).
// Uses tri.v[][] (float) NOT tri.iv[][] (quantized int).
//
// WINDING: Pushes triVerts as [v0, v2, v1] to reverse Houdini's CW winding
// to the CCW winding that Manifold requires for outward-facing normals.

struct FloatVtxHash {
    size_t operator()(const std::array<float,3>& p) const {
        size_t h = 14695981039346656037ULL;
        const unsigned char* raw = reinterpret_cast<const unsigned char*>(p.data());
        for (size_t i = 0; i < sizeof(float)*3; ++i) {
            h ^= raw[i];
            h *= 1099511628211ULL;
        }
        return h;
    }
};

struct FloatVtxEq {
    bool operator()(const std::array<float,3>& a,
                    const std::array<float,3>& b) const {
        return std::memcmp(a.data(), b.data(), sizeof(float)*3) == 0;
    }
};

manifold::MeshGL soupToMeshGL(const PolygonSoup& soup,
                               uint32_t mesh_id,
                               uint32_t originalID)
{
    manifold::MeshGL gl;
    gl.numProp = 3;  // x, y, z

    std::unordered_map<std::array<float,3>, uint32_t,
                       FloatVtxHash, FloatVtxEq> vmap;

    uint32_t numTris = 0;

    for (const auto& tri : soup.triangles) {
        if (static_cast<uint32_t>(tri.mesh_id) != mesh_id) continue;

        // Collect the 3 vertex indices for this triangle
        uint32_t vidx[3];
        for (int v = 0; v < 3; ++v) {
            std::array<float,3> key = { tri.v[v][0], tri.v[v][1], tri.v[v][2] };

            auto it = vmap.find(key);
            if (it == vmap.end()) {
                uint32_t idx = static_cast<uint32_t>(vmap.size());
                vmap[key] = idx;
                gl.vertProperties.push_back(key[0]);
                gl.vertProperties.push_back(key[1]);
                gl.vertProperties.push_back(key[2]);
                vidx[v] = idx;
            } else {
                vidx[v] = it->second;
            }
        }

        // WINDING FIX: push as [v0, v2, v1] to reverse CW → CCW
        gl.triVerts.push_back(vidx[0]);
        gl.triVerts.push_back(vidx[2]);
        gl.triVerts.push_back(vidx[1]);

        numTris++;
    }

    // Face tracking: single run covering all triangles, tagged with originalID
    if (numTris > 0) {
        gl.runIndex.push_back(0);
        gl.runIndex.push_back(numTris);
        gl.runOriginalID.push_back(originalID);
    }

    return gl;
}

// ─────────────────────────────────────────────────────────────────────────────
// §2  MeshGL → Triangle array  (output extraction)
// ─────────────────────────────────────────────────────────────────────────────
//
// WINDING FIX: Reads Manifold's CCW output and reverses back to Houdini's
// CW convention by swapping v1↔v2 in the output triangle.

void appendMeshGLToTriangles(const manifold::MeshGL& gl,
                              std::vector<Triangle>& out,
                              int piece_id,
                              uint32_t idB)
{
    const uint32_t numTris = static_cast<uint32_t>(gl.triVerts.size() / 3);
    const uint32_t numProp = gl.numProp;

    auto getOriginalID = [&](uint32_t triIdx) -> uint32_t {
        for (size_t r = 0; r + 1 < gl.runIndex.size(); ++r) {
            if (triIdx >= gl.runIndex[r] && triIdx < gl.runIndex[r + 1]) {
                return (r < gl.runOriginalID.size()) ? gl.runOriginalID[r] : 0;
            }
        }
        return 0;
    };

    for (uint32_t f = 0; f < numTris; ++f) {
        Triangle tri;
        tri.mesh_id = piece_id;

        // Read Manifold's vertex order [v0, v1, v2] (CCW)
        // and store as [v0, v2, v1] to reverse back to Houdini's CW convention
        static const int remap[3] = {0, 2, 1};

        for (int i = 0; i < 3; ++i) {
            int v = remap[i];
            uint32_t vi = gl.triVerts[f * 3 + v];
            tri.v[i][0] = gl.vertProperties[vi * numProp + 0];
            tri.v[i][1] = gl.vertProperties[vi * numProp + 1];
            tri.v[i][2] = gl.vertProperties[vi * numProp + 2];
            tri.iv[i] = {0, 0, 0};
        }

        uint32_t faceOrigin = getOriginalID(f);
        bool is_cap = (faceOrigin == idB);
        tri.src_prim_idx = is_cap ? CAP_FLAG : 0;

        out.push_back(tri);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// §3  OpType mapping
// ─────────────────────────────────────────────────────────────────────────────

manifold::OpType toManifoldOp(BooleanOp op)
{
    switch (op) {
        case BooleanOp::Union:        return manifold::OpType::Add;
        case BooleanOp::Intersection: return manifold::OpType::Intersect;
        case BooleanOp::DiffAB:       return manifold::OpType::Subtract;
        default:                      return manifold::OpType::Add;
    }
}

} // anonymous namespace

// ═══════════════════════════════════════════════════════════════════════════════
// §4  ManifoldBackend::execute
// ═══════════════════════════════════════════════════════════════════════════════

BackendResult ManifoldBackend::execute(const PolygonSoup& soup, BooleanOp op)
{
    BackendResult result;
    auto t0 = std::chrono::high_resolution_clock::now();

    // ── Reserve unique IDs for face tracking ─────────────────────────────
    uint32_t idA = manifold::Manifold::ReserveIDs(1);
    uint32_t idB = manifold::Manifold::ReserveIDs(1);

    // ── Convert inputs to MeshGL (with winding reversal) ─────────────────
    manifold::MeshGL glA = soupToMeshGL(soup, 0, idA);
    manifold::MeshGL glB = soupToMeshGL(soup, 1, idB);

    if (glA.triVerts.empty() || glB.triVerts.empty()) {
        result.success       = false;
        result.error_message = "Empty input mesh (A or B has no triangles)";
        return result;
    }

    // ── Create Manifold objects ──────────────────────────────────────────
    manifold::Manifold mA, mB;
    try {
        mA = manifold::Manifold(glA);
        mB = manifold::Manifold(glB);
    } catch (const std::exception& e) {
        result.success       = false;
        result.error_message = std::string("Failed to create Manifold: ") + e.what();
        return result;
    }

    if (mA.Status() != manifold::Manifold::Error::NoError) {
        result.success       = false;
        result.error_message = "Mesh A is not a valid manifold";
        return result;
    }
    if (mB.Status() != manifold::Manifold::Error::NoError) {
        result.success       = false;
        result.error_message = "Mesh B is not a valid manifold";
        return result;
    }

    // ── Execute operation ────────────────────────────────────────────────
    try {
        if (op == BooleanOp::Shatter) {
            // ── SHATTER: Two Boolean ops → two watertight pieces of A ────
            //
            //   Piece 0: A - B  (part of A outside B)
            //   Piece 1: A ∩ B  (part of A inside B)
            //
            // Houdini convention: B is a cutter only. B's own geometry
            // does NOT appear in the output. Only pieces of A are emitted.

            manifold::Manifold diffAB = mA.Boolean(mB, manifold::OpType::Subtract);
            manifold::Manifold isect  = mA.Boolean(mB, manifold::OpType::Intersect);

            if (diffAB.Status() != manifold::Manifold::Error::NoError) {
                result.success = false;
                result.error_message = "Shatter A-B failed";
                return result;
            }
            if (isect.Status() != manifold::Manifold::Error::NoError) {
                result.success = false;
                result.error_message = "Shatter A∩B failed";
                return result;
            }

            manifold::MeshGL out0 = diffAB.GetMeshGL();
            manifold::MeshGL out1 = isect.GetMeshGL();

            appendMeshGLToTriangles(out0, result.output_soup.triangles, 0, idB);
            appendMeshGLToTriangles(out1, result.output_soup.triangles, 1, idB);

        } else if (op == BooleanOp::DiffBA) {
            manifold::Manifold out = mB.Boolean(mA, manifold::OpType::Subtract);

            if (out.Status() != manifold::Manifold::Error::NoError) {
                result.success = false;
                result.error_message = "Boolean B-A failed";
                return result;
            }

            manifold::MeshGL outGL = out.GetMeshGL();
            appendMeshGLToTriangles(outGL, result.output_soup.triangles, 0, idB);

        } else if (op == BooleanOp::Seam) {
            result.success       = false;
            result.error_message = "Manifold backend does not support Seam operation. "
                                   "Use Cherchi or MCUT backend for seam extraction.";
            return result;

        } else {
            // Standard Boolean: Union, Intersection, DiffAB
            manifold::OpType mop = toManifoldOp(op);
            manifold::Manifold out = mA.Boolean(mB, mop);

            if (out.Status() != manifold::Manifold::Error::NoError) {
                result.success = false;
                result.error_message = "Boolean operation produced invalid result";
                return result;
            }

            manifold::MeshGL outGL = out.GetMeshGL();
            appendMeshGLToTriangles(outGL, result.output_soup.triangles, 0, idB);
        }

    } catch (const std::exception& e) {
        result.success       = false;
        result.error_message = std::string("Boolean operation failed: ") + e.what();
        return result;
    }

    // ── Stats ────────────────────────────────────────────────────────────
    result.success        = true;
    result.num_triangles  = static_cast<uint32_t>(result.output_soup.triangles.size());
    result.num_vertices   = 0;

    std::memcpy(result.output_soup.quantization_scale,
                soup.quantization_scale, sizeof(double) * 3);
    std::memcpy(result.output_soup.quantization_inv_scale,
                soup.quantization_inv_scale, sizeof(double) * 3);
    std::memcpy(result.output_soup.bbox_min, soup.bbox_min, sizeof(double) * 3);
    std::memcpy(result.output_soup.bbox_max, soup.bbox_max, sizeof(double) * 3);

    auto t1 = std::chrono::high_resolution_clock::now();
    result.execution_time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    return result;
}

} // namespace ember
