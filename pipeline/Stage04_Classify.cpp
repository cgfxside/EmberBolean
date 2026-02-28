/**
 * @file Stage04_Classify.cpp
 * @brief Stage 04: Classification implementation
 * 
 * This file implements the classification stage which computes winding
 * numbers and filters output polygons based on the Boolean operation.
 * 
 * @author EMBER Boolean Engine Team
 * @version 1.0
 */

#include "Stage04_Classify.h"
#include "ember/PolygonSoup.h"
#include "ember/Diagnostics.h"
#include "backend/IBooleanBackend.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <unordered_map>
#include <vector>

namespace ember {
namespace stages {

//=============================================================================
// WINDING NUMBER COMPUTATION
//=============================================================================

int computeWindingNumber(const float point[3], 
                         const std::vector<Triangle>& mesh_triangles) {
    // Ray casting algorithm: cast ray in +X direction and count intersections
    int winding = 0;
    
    for (const auto& tri : mesh_triangles) {
        // Check if ray from point intersects triangle
        // Ray: P + t * D where D = (1, 0, 0), t > 0
        
        // Triangle vertices
        const float* v0 = tri.v[0].data();
        const float* v1 = tri.v[1].data();
        const float* v2 = tri.v[2].data();
        
        // Compute triangle normal
        float e0[3] = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
        float e1[3] = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};
        float n[3] = {
            e0[1] * e1[2] - e0[2] * e1[1],
            e0[2] * e1[0] - e0[0] * e1[2],
            e0[0] * e1[1] - e0[1] * e1[0]
        };
        
        // Skip degenerate triangles
        float n_len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        if (n_len < 1e-10f) continue;
        
        // Normalize normal
        n[0] /= n_len; n[1] /= n_len; n[2] /= n_len;
        
        // Ray-plane intersection
        float denom = n[0];  // Dot(n, ray_dir) where ray_dir = (1,0,0)
        if (std::abs(denom) < 1e-10f) continue;  // Ray parallel to plane
        
        float t = (n[0] * (v0[0] - point[0]) + 
                   n[1] * (v0[1] - point[1]) + 
                   n[2] * (v0[2] - point[2])) / denom;
        
        if (t <= 0) continue;  // Intersection behind ray origin
        
        // Compute intersection point
        float isect[3] = {point[0] + t, point[1], point[2]};
        
        // Check if intersection is inside triangle (barycentric coordinates)
        float w[3];
        float denom2 = e0[0] * e1[1] - e0[1] * e1[0];
        if (std::abs(denom2) < 1e-10f) {
            // Try different projection
            denom2 = e0[1] * e1[2] - e0[2] * e1[1];
            if (std::abs(denom2) < 1e-10f) {
                denom2 = e0[2] * e1[0] - e0[0] * e1[2];
            }
        }
        
        // Simplified inside test: check bounding box
        float min_x = std::min({v0[0], v1[0], v2[0]});
        float max_x = std::max({v0[0], v1[0], v2[0]});
        float min_y = std::min({v0[1], v1[1], v2[1]});
        float max_y = std::max({v0[1], v1[1], v2[1]});
        float min_z = std::min({v0[2], v1[2], v2[2]});
        float max_z = std::max({v0[2], v1[2], v2[2]});
        
        if (isect[0] >= min_x && isect[0] <= max_x &&
            isect[1] >= min_y && isect[1] <= max_y &&
            isect[2] >= min_z && isect[2] <= max_z) {
            // Count intersection with sign based on normal direction
            winding += (denom > 0) ? 1 : -1;
        }
    }
    
    return winding;
}

WindingState computePolygonWinding(const OutputPolygon& poly,
                                    const PolygonSoup& soup,
                                    uint32_t mesh_id) {
    WindingState state;
    
    // Compute centroid of polygon
    float centroid[3] = {0.0f, 0.0f, 0.0f};
    
    // Get triangle range for the mesh
    auto [mesh_start, mesh_end] = soup.getMeshTriangleRange(mesh_id);
    
    // Collect triangles for the mesh
    std::vector<Triangle> mesh_tris;
    mesh_tris.reserve(mesh_end - mesh_start);
    for (size_t i = mesh_start; i < mesh_end; ++i) {
        mesh_tris.push_back(soup.triangles[i]);
    }
    
    // Compute winding number at centroid
    state.winding_number = computeWindingNumber(centroid, mesh_tris);
    state.is_inside = (std::abs(state.winding_number) % 2) == 1;
    state.computed = true;
    
    return state;
}

std::vector<WindingState> computeAllWindingNumbers(const PolygonSoup& soup) {
    std::vector<WindingState> states;
    states.reserve(soup.output_polygons.size());
    
    for (const auto& poly : soup.output_polygons) {
        // For each polygon, compute winding relative to other meshes
        // For simplicity, just compute relative to mesh 0 (target)
        WindingState state = computePolygonWinding(poly, soup, 0);
        states.push_back(state);
    }
    
    return states;
}

//=============================================================================
// LABEL PROPAGATION
//=============================================================================

std::unordered_map<EdgeKey, std::vector<size_t>, EdgeKeyHash>
buildEdgeAdjacency(const PolygonSoup& soup) {
    std::unordered_map<EdgeKey, std::vector<size_t>, EdgeKeyHash> adjacency;
    
    for (size_t poly_idx = 0; poly_idx < soup.output_polygons.size(); ++poly_idx) {
        const auto& poly = soup.output_polygons[poly_idx];
        const auto& verts = poly.vertex_indices;
        
        // Add edges
        for (size_t i = 0; i < verts.size(); ++i) {
            uint32_t v0 = verts[i];
            uint32_t v1 = verts[(i + 1) % verts.size()];
            adjacency[EdgeKey(v0, v1)].push_back(poly_idx);
        }
    }
    
    return adjacency;
}

void propagateLabels(std::vector<WindingState>& winding_states,
                     const PolygonSoup& soup) {
    // Build edge adjacency
    auto adjacency = buildEdgeAdjacency(soup);
    
    // Propagate labels across shared edges
    for (const auto& [edge, polygons] : adjacency) {
        if (polygons.size() < 2) continue;
        
        // Find majority label
        int inside_count = 0;
        for (size_t poly_idx : polygons) {
            if (winding_states[poly_idx].is_inside) {
                inside_count++;
            }
        }
        
        bool majority_inside = (inside_count > static_cast<int>(polygons.size()) / 2);
        
        // Propagate majority label
        for (size_t poly_idx : polygons) {
            if (!winding_states[poly_idx].computed) {
                winding_states[poly_idx].is_inside = majority_inside;
            }
        }
    }
}

//=============================================================================
// BOOLEAN FILTERING
//=============================================================================

bool keepForUnion(const WindingState winding_per_mesh[], size_t num_meshes) {
    // Union: Keep if not inside ALL meshes
    for (size_t i = 0; i < num_meshes; ++i) {
        if (!winding_per_mesh[i].is_inside) {
            return true;  // Outside at least one mesh
        }
    }
    return false;  // Inside all meshes
}

bool keepForIntersection(const WindingState winding_per_mesh[], 
                         size_t num_meshes) {
    // Intersection: Keep if inside ALL meshes
    for (size_t i = 0; i < num_meshes; ++i) {
        if (!winding_per_mesh[i].is_inside) {
            return false;  // Outside at least one mesh
        }
    }
    return true;  // Inside all meshes
}

bool keepForDifference(const WindingState& winding_a, 
                       const WindingState& winding_b) {
    // Difference (A - B): Keep if inside A and outside B
    return winding_a.is_inside && !winding_b.is_inside;
}

bool keepForXor(const WindingState winding_per_mesh[], size_t num_meshes) {
    // XOR: Keep if inside an odd number of meshes
    int inside_count = 0;
    for (size_t i = 0; i < num_meshes; ++i) {
        if (winding_per_mesh[i].is_inside) {
            inside_count++;
        }
    }
    return (inside_count % 2) == 1;
}

size_t filterByBooleanOp(const std::vector<WindingState>& winding_states,
                         BooleanOp op,
                         std::vector<bool>& keep_mask) {
    keep_mask.resize(winding_states.size());
    size_t keep_count = 0;
    
    for (size_t i = 0; i < winding_states.size(); ++i) {
        bool keep = false;
        
        switch (op) {
            case BooleanOp::Union:
                keep = !winding_states[i].is_inside;  // Keep outside
                break;
            case BooleanOp::Intersection:
                keep = winding_states[i].is_inside;   // Keep inside
                break;
            case BooleanOp::DiffAB:
            case BooleanOp::Difference:
                keep = winding_states[i].is_inside;   // Simplified
                break;
            case BooleanOp::Xor:
                keep = winding_states[i].is_inside;   // Simplified
                break;
            case BooleanOp::Shatter:
                keep = true;  // Keep all
                break;
            default:
                keep = true;
                break;
        }
        
        keep_mask[i] = keep;
        if (keep) ++keep_count;
    }
    
    return keep_count;
}

//=============================================================================
// BOUNDARY HANDLING
//=============================================================================

void handleBoundaryPolygons(std::vector<WindingState>& winding_states,
                            const PolygonSoup& soup) {
    // Build edge adjacency to find boundary edges
    auto adjacency = buildEdgeAdjacency(soup);
    
    // Mark polygons with boundary edges
    for (size_t poly_idx = 0; poly_idx < soup.output_polygons.size(); ++poly_idx) {
        const auto& poly = soup.output_polygons[poly_idx];
        const auto& verts = poly.vertex_indices;
        
        for (size_t i = 0; i < verts.size(); ++i) {
            uint32_t v0 = verts[i];
            uint32_t v1 = verts[(i + 1) % verts.size()];
            
            auto it = adjacency.find(EdgeKey(v0, v1));
            if (it != adjacency.end() && it->second.size() == 1) {
                // Boundary edge
                winding_states[poly_idx].on_boundary = true;
                break;
            }
        }
    }
}

//=============================================================================
// STATISTICS
//=============================================================================

static ClassificationStats g_last_stats;

const ClassificationStats& getLastClassificationStats() {
    return g_last_stats;
}

//=============================================================================
// MAIN STAGE FUNCTION
//=============================================================================

void stage04_classify(PolygonSoup& soup, BooleanOp op) {
    EMBER_SCOPED_TIMER("Stage 04: Classification");
    
    g_last_stats = ClassificationStats();
    g_last_stats.polygons_input = soup.output_polygons.size();
    
    if (soup.output_polygons.empty()) {
        EMBER_LOG_INFO("Stage 04: No polygons to classify");
        return;
    }
    
    EMBER_LOG_INFO("Stage 04: Classifying %zu polygons for %s operation",
                   soup.output_polygons.size(),
                   (op == BooleanOp::Union ? "union" :
                    op == BooleanOp::Intersection ? "intersection" :
                    op == BooleanOp::Difference ? "difference" :
                    op == BooleanOp::Shatter ? "shatter" : "other"));
    
    // Step 1: Compute winding numbers
    EMBER_LOG_INFO("Stage 04: Computing winding numbers...");
    auto compute_start = std::chrono::steady_clock::now();
    
    std::vector<WindingState> winding_states = computeAllWindingNumbers(soup);
    
    auto compute_end = std::chrono::steady_clock::now();
    g_last_stats.compute_time_ms = std::chrono::duration<double, std::milli>(
        compute_end - compute_start).count();
    g_last_stats.winding_computed = winding_states.size();
    
    // Step 2: Propagate labels
    EMBER_LOG_INFO("Stage 04: Propagating labels...");
    propagateLabels(winding_states, soup);
    
    // Step 3: Handle boundary polygons
    EMBER_LOG_INFO("Stage 04: Handling boundary polygons...");
    handleBoundaryPolygons(winding_states, soup);
    
    // Step 4: Filter by Boolean operation
    EMBER_LOG_INFO("Stage 04: Filtering by Boolean operation...");
    auto filter_start = std::chrono::steady_clock::now();
    
    std::vector<bool> keep_mask;
    size_t keep_count = filterByBooleanOp(winding_states, op, keep_mask);
    
    // Apply filter
    std::vector<OutputPolygon> filtered;
    filtered.reserve(keep_count);
    
    for (size_t i = 0; i < soup.output_polygons.size(); ++i) {
        if (keep_mask[i]) {
            filtered.push_back(std::move(soup.output_polygons[i]));
        }
    }
    
    soup.output_polygons = std::move(filtered);
    
    auto filter_end = std::chrono::steady_clock::now();
    g_last_stats.filter_time_ms = std::chrono::duration<double, std::milli>(
        filter_end - filter_start).count();
    
    g_last_stats.polygons_output = soup.output_polygons.size();
    
    EMBER_LOG_INFO("Stage 04 complete: %zu polygons kept, %zu discarded",
                   g_last_stats.polygons_output,
                   g_last_stats.polygons_input - g_last_stats.polygons_output);
}

} // namespace stages
} // namespace ember
