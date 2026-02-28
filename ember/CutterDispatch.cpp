#include "CutterDispatch.h"
#include "Diagnostics.h"
#include <cmath>
#include <numeric>      // for std::iota
#include <unordered_map> // for std::unordered_map

namespace ember {

// ═══════════════════════════════════════════════════════════════
// Internal: Connected Component Extraction
// ═══════════════════════════════════════════════════════════════

struct ComponentInfo {
    std::vector<uint32_t> triangle_indices;
    uint32_t tri_start;
    uint32_t tri_count;
};

// Extract connected components from triangles belonging to a
// specific mesh_id. Uses edge-adjacency via union-find.
// Triangles are adjacent if they share an edge (2 vertices with same coordinates).
static std::vector<ComponentInfo> extractConnectedComponents(
    const PolygonSoup& soup,
    uint32_t mesh_id_filter)
{
    // Collect all triangles belonging to the target mesh
    std::vector<uint32_t> target_tris;
    for (uint32_t i = 0; i < soup.triangles.size(); ++i) {
        if (soup.triangles[i].mesh_id == static_cast<int>(mesh_id_filter)) {
            target_tris.push_back(i);
        }
    }

    if (target_tris.empty()) return {};

    // Union-Find for edge-adjacency grouping
    std::vector<uint32_t> parent(target_tris.size());
    std::iota(parent.begin(), parent.end(), 0);

    auto find = [&](uint32_t x) -> uint32_t {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    auto unite = [&](uint32_t a, uint32_t b) {
        a = find(a); b = find(b);
        if (a != b) parent[a] = b;
    };

    // Build edge -> triangle index map for adjacency
    // Edge key: hash of two vertex coordinates
    struct EdgeKey {
        int32_t v0[3], v1[3];  // Store actual coordinates, not indices
        bool operator==(const EdgeKey& o) const {
            return v0[0] == o.v0[0] && v0[1] == o.v0[1] && v0[2] == o.v0[2] &&
                   v1[0] == o.v1[0] && v1[1] == o.v1[1] && v1[2] == o.v1[2];
        }
    };
    struct EdgeHash {
        size_t operator()(const EdgeKey& e) const {
            // Simple hash combining coordinates
            size_t h = 0;
            for (int i = 0; i < 3; ++i) {
                h = h * 31 + static_cast<size_t>(e.v0[i]);
                h = h * 31 + static_cast<size_t>(e.v1[i]);
            }
            return h;
        }
    };

    std::unordered_map<EdgeKey, uint32_t, EdgeHash> edge_to_tri;

    for (uint32_t local = 0; local < target_tris.size(); ++local) {
        const auto& tri = soup.triangles[target_tris[local]];
        for (int e = 0; e < 3; ++e) {
            // Get integer vertices for this edge
            const auto& a = tri.iv[e];
            const auto& b = tri.iv[(e + 1) % 3];
            
            // Create canonical edge key (sorted by coordinates)
            EdgeKey key;
            bool a_less_b = (a[0] < b[0]) || 
                           (a[0] == b[0] && a[1] < b[1]) ||
                           (a[0] == b[0] && a[1] == b[1] && a[2] < b[2]);
            const auto& first = a_less_b ? a : b;
            const auto& second = a_less_b ? b : a;
            for (int i = 0; i < 3; ++i) {
                key.v0[i] = first[i];
                key.v1[i] = second[i];
            }

            auto it = edge_to_tri.find(key);
            if (it != edge_to_tri.end()) {
                unite(local, it->second);
            } else {
                edge_to_tri[key] = local;
            }
        }
    }

    // Group triangles by root
    std::unordered_map<uint32_t, uint32_t> root_to_component;
    std::vector<ComponentInfo> components;

    for (uint32_t local = 0; local < target_tris.size(); ++local) {
        uint32_t root = find(local);
        auto it = root_to_component.find(root);
        uint32_t comp_idx;
        if (it == root_to_component.end()) {
            comp_idx = components.size();
            root_to_component[root] = comp_idx;
            components.emplace_back();
        } else {
            comp_idx = it->second;
        }
        components[comp_idx].triangle_indices.push_back(target_tris[local]);
    }

    // Set tri_start/tri_count for each component
    for (auto& comp : components) {
        comp.tri_start = comp.triangle_indices.front();
        comp.tri_count = comp.triangle_indices.size();
    }

    return components;
}

// ═══════════════════════════════════════════════════════════════
// Internal: Exact Coplanarity Test
// ═══════════════════════════════════════════════════════════════

// Test if all triangles in a component share the same plane.
// Uses exact cross-product for the normal from the first
// non-degenerate triangle, then exact dot-product to verify
// every vertex satisfies the plane equation.
// NO EPSILONS. dot != 0 means not coplanar.

static bool isComponentPlanar(
    const PolygonSoup& soup,
    const std::vector<uint32_t>& component_tris,
    IntPlane& out_plane)
{
    if (component_tris.empty()) return false;

    // Find first non-degenerate triangle to establish the plane
    bool found_plane = false;
    for (uint32_t ti : component_tris) {
        const auto& tri = soup.triangles[ti];
        // Use integer vertices directly from triangle
        const auto& v0 = tri.iv[0];
        const auto& v1 = tri.iv[1];
        const auto& v2 = tri.iv[2];

        // Exact cross product (26-bit inputs → 53-bit result)
        int64_t nx = int64_t(v1[1] - v0[1]) * int64_t(v2[2] - v0[2])
                   - int64_t(v1[2] - v0[2]) * int64_t(v2[1] - v0[1]);
        int64_t ny = int64_t(v1[2] - v0[2]) * int64_t(v2[0] - v0[0])
                   - int64_t(v1[0] - v0[0]) * int64_t(v2[2] - v0[2]);
        int64_t nz = int64_t(v1[0] - v0[0]) * int64_t(v2[1] - v0[1])
                   - int64_t(v1[1] - v0[1]) * int64_t(v2[0] - v0[0]);

        if (nx == 0 && ny == 0 && nz == 0) continue; // degenerate tri

        // d = -(n · v0)  (53-bit * 26-bit = 80-bit, fits in Int128)
        Int128 d = -(Int128(nx) * v0[0]
                   + Int128(ny) * v0[1]
                   + Int128(nz) * v0[2]);

        out_plane = IntPlane{nx, ny, nz, d};
        found_plane = true;
        break;
    }

    if (!found_plane) return false;

    // Verify ALL vertices satisfy the plane equation exactly
    for (uint32_t ti : component_tris) {
        const auto& tri = soup.triangles[ti];
        for (int v = 0; v < 3; ++v) {
            const auto& pt = tri.iv[v];
            Int128 dot = Int128(out_plane.a) * pt[0]
                       + Int128(out_plane.b) * pt[1]
                       + Int128(out_plane.c) * pt[2]
                       + out_plane.d;
            if (dot != Int128::zero()) return false;  // NOT coplanar — exact test
        }
    }
    return true;
}

// ═══════════════════════════════════════════════════════════════
// Grout Plane Computation
// ═══════════════════════════════════════════════════════════════

void computeGroutPlanes(
    const IntPlane& base,
    float grout_width,
    double quantization_scale,
    IntPlane& out_pos,
    IntPlane& out_neg)
{
    // Normal magnitude squared (exact, 53-bit inputs)
    Int128 n_sq = Int128(base.a) * base.a
               + Int128(base.b) * base.b
               + Int128(base.c) * base.c;

    // Integer offset: half_grout / quantization_scale * |n|
    double half_grout = grout_width * 0.5;
    double n_len = std::sqrt(static_cast<double>(n_sq));
    Int128 d_offset = Int128(static_cast<int64_t>(std::llround(
        half_grout * n_len / quantization_scale
    )));

    // Positive side: d + offset
    out_pos = base;
    out_pos.d = base.d + d_offset;

    // Negative side: d - offset
    out_neg = base;
    out_neg.d = base.d - d_offset;
}

// ═══════════════════════════════════════════════════════════════
// Main Analysis Entry Point
// ═══════════════════════════════════════════════════════════════

DispatchResult analyzeCutters(
    const PolygonSoup& soup,
    const Parameters& params)
{
    DispatchResult result;
    result.force_exact  = params.force_exact;
    result.all_planar   = true;
    result.any_grout    = false;
    result.any_noised   = false;
    result.tier1_count  = 0;
    result.tier2_count  = 0;

    // Extract connected components from all non-target meshes
    // For N-way: iterate mesh_id 1..N
    // For 2-way: only mesh_id 1
    uint32_t max_mesh_id = soup.mesh_count > 0 ? soup.mesh_count - 1 : 0;

    for (uint32_t mid = 1; mid <= max_mesh_id; ++mid) {
        auto components = extractConnectedComponents(soup, mid);

        for (uint32_t ci = 0; ci < components.size(); ++ci) {
            CutterDescriptor desc{};
            desc.component_id = ci;
            desc.mesh_id      = mid;
            desc.grout_width  = params.grout_width;
            desc.noise_seed   = params.noise_seed;

            // Force Exact override — bypass all analysis
            if (params.force_exact) {
                desc.type      = CutterType::GeneralMesh;
                desc.tri_start = components[ci].tri_start;
                desc.tri_count = components[ci].tri_count;
                result.all_planar = false;
                result.tier2_count++;
                result.cutters.push_back(desc);
                continue;
            }

            // Test coplanarity (exact integer arithmetic)
            IntPlane plane;
            bool planar = isComponentPlanar(
                soup, components[ci].triangle_indices, plane
            );

            if (!planar) {
                // General mesh → Tier 2 (exact mesh processing)
                desc.type = CutterType::GeneralMesh;
                desc.tri_start = components[ci].tri_start;
                desc.tri_count = components[ci].tri_count;
                result.all_planar = false;
                result.tier2_count++;
            } else {
                // Planar — check for grout/noise
                bool has_grout = (params.grout_width > 0.0f);
                bool has_noise = (params.noise_amplitude > 0.0f);

                if (has_noise) {
                    // Noised planar → Tier 2 (tessellation required)
                    desc.type = CutterType::NoisedMesh;
                    desc.tri_start = components[ci].tri_start;
                    desc.tri_count = components[ci].tri_count;
                    result.any_noised = true;
                    result.tier2_count++;
                } else if (has_grout) {
                    // Planar with grout → Tier 1G (dual-plane fast-path)
                    desc.type = CutterType::PlanarGrout;
                    desc.plane = plane;
                    computeGroutPlanes(
                        plane,
                        params.grout_width,
                        soup.quantization_scale[0],  // Use X scale as reference
                        desc.grout_plane_pos,
                        desc.grout_plane_neg
                    );
                    result.any_grout = true;
                    result.tier1_count++;
                } else {
                    // Simple planar → Tier 1 (single-plane fast-path)
                    desc.type = CutterType::Planar;
                    desc.plane = plane;
                    result.tier1_count++;
                }
            }

            result.cutters.push_back(desc);
        }
    }

    return result;
}

} // namespace ember
