/**
 * @file Stage01_Diagnostic.cpp
 * @brief Stage 01: Diagnostic Gatekeeper implementation
 * 
 * This file implements the diagnostic stage which validates and prepares
 * input meshes for Boolean operations in the EMBER pipeline.
 * 
 * @author EMBER Boolean Engine Team
 * @version 1.0
 */

#include "Stage01_Diagnostic.h"
#include "ember/PolygonSoup.h"
#include "ember/Diagnostics.h"
#include "ember/Plane.h"
#include "ember/IntegerTypes.h"

#include "../upload/CutterDispatch.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>

namespace ember {
namespace stages {

//=============================================================================
// TRIANGULATION IMPLEMENTATION
//=============================================================================

size_t fanTriangulate(const std::vector<uint32_t>& vertices,
                      int mesh_id,
                      uint32_t src_prim_idx,
                      PolygonSoup& soup) {
    const size_t n = vertices.size();
    if (n < 3) {
        return 0;  // Need at least 3 vertices for a triangle
    }
    
    // Fan triangulation from vertex 0
    // Creates triangles: (0,1,2), (0,2,3), (0,3,4), ...
    size_t triangles_created = 0;
    for (size_t i = 1; i + 1 < n; ++i) {
        Triangle tri;
        
        // Initialize vertices to zero
        tri.v[0] = {0.0f, 0.0f, 0.0f};
        tri.v[1] = {0.0f, 0.0f, 0.0f};
        tri.v[2] = {0.0f, 0.0f, 0.0f};
        
        // Set integer vertex indices (will be filled later from actual geometry)
        // For now, store polygon vertex indices in a custom way
        // In actual implementation, these would reference global vertex buffer
        tri.iv[0][0] = static_cast<int32_t>(vertices[0]);
        tri.iv[1][0] = static_cast<int32_t>(vertices[i]);
        tri.iv[2][0] = static_cast<int32_t>(vertices[i + 1]);
        
        tri.mesh_id = mesh_id;
        tri.src_prim_idx = src_prim_idx;
        
        soup.triangles.push_back(tri);
        ++triangles_created;
    }
    
    return triangles_created;
}

size_t triangulateAllPolygons(PolygonSoup& soup, DiagnosticLog& log) {
    // For EMBER, we assume input is already triangulated or handled
    // by the Houdini import. This function is a placeholder for
    // future n-gon support.
    
    // Count triangles per mesh
    soup.mesh_tri_count.clear();
    soup.mesh_tri_count.resize(soup.mesh_count, 0);
    
    for (const auto& tri : soup.triangles) {
        if (tri.mesh_id >= 0 && static_cast<uint32_t>(tri.mesh_id) < soup.mesh_count) {
            soup.mesh_tri_count[tri.mesh_id]++;
        }
    }
    
    log.info("Triangulation: " + std::to_string(soup.triangles.size()) + " triangles");
    
    return soup.triangles.size();
}

//=============================================================================
// DEGENERATE TRIANGLE REMOVAL
//=============================================================================

int64_t triangleSignedArea2x(const int32_t v0[3], const int32_t v1[3], 
                              const int32_t v2[3]) {
    // Compute 2D projected area (use XY plane by default)
    // Area = (v1.x - v0.x) * (v2.y - v0.y) - (v2.x - v0.x) * (v1.y - v0.y)
    int64_t e0x = static_cast<int64_t>(v1[0]) - v0[0];
    int64_t e0y = static_cast<int64_t>(v1[1]) - v0[1];
    int64_t e1x = static_cast<int64_t>(v2[0]) - v0[0];
    int64_t e1y = static_cast<int64_t>(v2[1]) - v0[1];
    
    return e0x * e1y - e1x * e0y;
}

bool isTriangleDegenerate(const Triangle& tri) {
    // Check for duplicate vertices
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            if (tri.iv[i][0] == tri.iv[j][0] &&
                tri.iv[i][1] == tri.iv[j][1] &&
                tri.iv[i][2] == tri.iv[j][2]) {
                return true;  // Duplicate vertices
            }
        }
    }
    
    // Compute 3D cross product magnitude squared
    // Edge vectors
    int64_t e0x = static_cast<int64_t>(tri.iv[1][0]) - tri.iv[0][0];
    int64_t e0y = static_cast<int64_t>(tri.iv[1][1]) - tri.iv[0][1];
    int64_t e0z = static_cast<int64_t>(tri.iv[1][2]) - tri.iv[0][2];
    
    int64_t e1x = static_cast<int64_t>(tri.iv[2][0]) - tri.iv[0][0];
    int64_t e1y = static_cast<int64_t>(tri.iv[2][1]) - tri.iv[0][1];
    int64_t e1z = static_cast<int64_t>(tri.iv[2][2]) - tri.iv[0][2];
    
    // Cross product
    int64_t cx = e0y * e1z - e0z * e1y;
    int64_t cy = e0z * e1x - e0x * e1z;
    int64_t cz = e0x * e1y - e0y * e1x;
    
    // Check if cross product is zero (zero area)
    return (cx == 0 && cy == 0 && cz == 0);
}

size_t removeDegenerateTriangles(PolygonSoup& soup, DiagnosticLog& log) {
    size_t removed = 0;
    std::vector<Triangle> valid_triangles;
    valid_triangles.reserve(soup.triangles.size());
    
    for (const auto& tri : soup.triangles) {
        if (isTriangleDegenerate(tri)) {
            ++removed;
        } else {
            valid_triangles.push_back(tri);
        }
    }
    
    if (removed > 0) {
        log.warn(DiagCategory::Degenerate, 
                 "Removed " + std::to_string(removed) + " degenerate triangles");
    }
    
    soup.triangles = std::move(valid_triangles);
    return removed;
}

//=============================================================================
// WINDING NORMALIZATION
//=============================================================================

bool computeExpectedWinding(const PolygonSoup& soup, uint32_t mesh_id) {
    // Default to CCW (counter-clockwise) winding
    // In a proper implementation, this would analyze the mesh
    // to determine the dominant winding order
    return true;  // CCW
}

size_t normalizeWinding(PolygonSoup& soup, DiagnosticLog& log) {
    size_t flipped = 0;
    
    // For each mesh, compute expected winding and flip triangles as needed
    for (uint32_t mesh_id = 0; mesh_id < soup.mesh_count; ++mesh_id) {
        bool expect_ccw = computeExpectedWinding(soup, mesh_id);
        
        auto [start, end] = soup.getMeshTriangleRange(mesh_id);
        
        for (size_t i = start; i < end; ++i) {
            Triangle& tri = soup.triangles[i];
            
            // Compute actual winding using cross product Z component
            int64_t e0x = static_cast<int64_t>(tri.iv[1][0]) - tri.iv[0][0];
            int64_t e0y = static_cast<int64_t>(tri.iv[1][1]) - tri.iv[0][1];
            int64_t e1x = static_cast<int64_t>(tri.iv[2][0]) - tri.iv[0][0];
            int64_t e1y = static_cast<int64_t>(tri.iv[2][1]) - tri.iv[0][1];
            
            int64_t cross_z = e0x * e1y - e1x * e0y;
            
            bool is_ccw = cross_z > 0;
            
            if (is_ccw != expect_ccw) {
                // Flip triangle by swapping vertices 1 and 2
                std::swap(tri.v[1], tri.v[2]);
                std::swap(tri.iv[1][0], tri.iv[2][0]);
                std::swap(tri.iv[1][1], tri.iv[2][1]);
                std::swap(tri.iv[1][2], tri.iv[2][2]);
                ++flipped;
            }
        }
    }
    
    if (flipped > 0) {
        log.info("Normalized winding: flipped " + std::to_string(flipped) + " triangles");
    }
    
    return flipped;
}

//=============================================================================
// WATERTIGHTNESS CHECK
//=============================================================================

bool checkWatertightness(const PolygonSoup& soup, uint32_t mesh_id, 
                         DiagnosticLog& log) {
    auto [start, end] = soup.getMeshTriangleRange(mesh_id);
    
    // Count edge occurrences
    std::unordered_map<Edge, int, EdgeHash> edge_counts;
    
    for (size_t i = start; i < end; ++i) {
        const Triangle& tri = soup.triangles[i];
        
        // Get global vertex indices for this triangle
        // In actual implementation, these would come from a global vertex buffer
        // For now, use the integer coordinates as unique identifiers
        uint32_t v0 = static_cast<uint32_t>(tri.iv[0][0] + tri.iv[0][1] * 1000 + tri.iv[0][2] * 1000000);
        uint32_t v1 = static_cast<uint32_t>(tri.iv[1][0] + tri.iv[1][1] * 1000 + tri.iv[1][2] * 1000000);
        uint32_t v2 = static_cast<uint32_t>(tri.iv[2][0] + tri.iv[2][1] * 1000 + tri.iv[2][2] * 1000000);
        
        Edge e01(v0, v1);
        Edge e12(v1, v2);
        Edge e20(v2, v0);
        
        edge_counts[e01]++;
        edge_counts[e12]++;
        edge_counts[e20]++;
    }
    
    // Check for boundary edges (should appear exactly twice for watertight)
    size_t boundary_edges = 0;
    for (const auto& [edge, count] : edge_counts) {
        if (count != 2) {
            ++boundary_edges;
        }
    }
    
    if (boundary_edges > 0) {
        log.warn(DiagCategory::NotWatertight,
                 "Mesh " + std::to_string(mesh_id) + " has " + 
                 std::to_string(boundary_edges) + " boundary edges (not watertight)");
        return false;
    }
    
    return true;
}

//=============================================================================
// PLANE COMPUTATION
//=============================================================================

void computeTrianglePlanes(PolygonSoup& soup) {
    for (auto& tri : soup.triangles) {
        tri.plane = planeFromTriangle(tri.iv[0], tri.iv[1], tri.iv[2]);
    }
}

//=============================================================================
// CONNECTED COMPONENTS
//=============================================================================

std::vector<uint32_t> computeConnectedComponents(const PolygonSoup& soup, 
                                                  uint32_t mesh_id) {
    auto [start, end] = soup.getMeshTriangleRange(mesh_id);
    size_t num_tris = end - start;
    
    if (num_tris == 0) {
        return {};
    }
    
    // Build adjacency graph
    // Two triangles are adjacent if they share an edge
    std::vector<std::vector<size_t>> adjacency(num_tris);
    
    // Map edges to triangles
    std::unordered_map<Edge, std::vector<size_t>, EdgeHash> edge_to_tris;
    
    for (size_t i = 0; i < num_tris; ++i) {
        const Triangle& tri = soup.triangles[start + i];
        
        uint32_t v0 = static_cast<uint32_t>(tri.iv[0][0] + tri.iv[0][1] * 1000 + tri.iv[0][2] * 1000000);
        uint32_t v1 = static_cast<uint32_t>(tri.iv[1][0] + tri.iv[1][1] * 1000 + tri.iv[1][2] * 1000000);
        uint32_t v2 = static_cast<uint32_t>(tri.iv[2][0] + tri.iv[2][1] * 1000 + tri.iv[2][2] * 1000000);
        
        edge_to_tris[Edge(v0, v1)].push_back(i);
        edge_to_tris[Edge(v1, v2)].push_back(i);
        edge_to_tris[Edge(v2, v0)].push_back(i);
    }
    
    // Build adjacency from shared edges
    for (const auto& [edge, tris] : edge_to_tris) {
        for (size_t i = 0; i < tris.size(); ++i) {
            for (size_t j = i + 1; j < tris.size(); ++j) {
                adjacency[tris[i]].push_back(tris[j]);
                adjacency[tris[j]].push_back(tris[i]);
            }
        }
    }
    
    // BFS to find connected components
    std::vector<uint32_t> component(num_tris, UINT32_MAX);
    uint32_t current_component = 0;
    
    for (size_t i = 0; i < num_tris; ++i) {
        if (component[i] != UINT32_MAX) continue;
        
        // BFS from this triangle
        std::queue<size_t> queue;
        queue.push(i);
        component[i] = current_component;
        
        while (!queue.empty()) {
            size_t cur = queue.front();
            queue.pop();
            
            for (size_t neighbor : adjacency[cur]) {
                if (component[neighbor] == UINT32_MAX) {
                    component[neighbor] = current_component;
                    queue.push(neighbor);
                }
            }
        }
        
        ++current_component;
    }
    
    return component;
}

//=============================================================================
// CUTTER TYPE DETECTION
//=============================================================================

bool isMeshCoplanar(const PolygonSoup& soup, uint32_t mesh_id, Plane& plane) {
    auto [start, end] = soup.getMeshTriangleRange(mesh_id);
    
    if (start >= end) {
        return false;
    }
    
    // Use first triangle's plane as reference
    plane = soup.triangles[start].plane;
    
    // Check if all other triangles are coplanar
    for (size_t i = start + 1; i < end; ++i) {
        const Plane& tri_plane = soup.triangles[i].plane;
        
        // Check if planes are parallel (cross product of normals = 0)
        int64_t cx = static_cast<int64_t>(plane.b) * tri_plane.c - 
                     static_cast<int64_t>(plane.c) * tri_plane.b;
        int64_t cy = static_cast<int64_t>(plane.c) * tri_plane.a - 
                     static_cast<int64_t>(plane.a) * tri_plane.c;
        int64_t cz = static_cast<int64_t>(plane.a) * tri_plane.b - 
                     static_cast<int64_t>(plane.b) * tri_plane.a;
        
        if (cx != 0 || cy != 0 || cz != 0) {
            return false;  // Not parallel, so not coplanar
        }
        
        // Check if d values are consistent (planes are the same, not just parallel)
        // Normalize and compare
        int64_t n1_len_sq = static_cast<int64_t>(plane.a) * plane.a + 
                            static_cast<int64_t>(plane.b) * plane.b + 
                            static_cast<int64_t>(plane.c) * plane.c;
        int64_t n2_len_sq = static_cast<int64_t>(tri_plane.a) * tri_plane.a + 
                            static_cast<int64_t>(tri_plane.b) * tri_plane.b + 
                            static_cast<int64_t>(tri_plane.c) * tri_plane.c;
        
        if (n1_len_sq == 0 || n2_len_sq == 0) {
            continue;  // Degenerate plane
        }
        
        // Check if d1/|n1| == d2/|n2|
        // Use cross multiplication to avoid division
        int64_t d1_scaled = static_cast<int64_t>(plane.d) * 
                            static_cast<int64_t>(std::sqrt(n2_len_sq));
        int64_t d2_scaled = static_cast<int64_t>(tri_plane.d) * 
                            static_cast<int64_t>(std::sqrt(n1_len_sq));
        
        // Allow small tolerance for integer rounding
        if (std::abs(d1_scaled - d2_scaled) > 1000) {
            return false;
        }
    }
    
    return true;
}

void detectCutterTypes(const PolygonSoup& soup, 
                       const PipelineConfig& config,
                       DispatchResult& dispatch) {
    dispatch.cutters.clear();
    dispatch.all_planar = true;
    dispatch.any_grout = false;
    dispatch.any_noised = false;
    dispatch.tier1_count = 0;
    dispatch.tier2_count = 0;
    dispatch.force_exact = config.force_exact;
    
    // Analyze each cutter mesh (mesh_id >= 1)
    for (uint32_t mesh_id = 1; mesh_id < soup.mesh_count; ++mesh_id) {
        auto [start, end] = soup.getMeshTriangleRange(mesh_id);
        
        if (start >= end) {
            continue;  // Empty mesh
        }
        
        // Compute connected components
        std::vector<uint32_t> components = computeConnectedComponents(soup, mesh_id);
        
        // Analyze each connected component
        std::unordered_set<uint32_t> unique_components(components.begin(), components.end());
        
        for (uint32_t comp_id : unique_components) {
            CutterDescriptor cutter;
            cutter.component_id = comp_id;
            cutter.mesh_id = mesh_id;
            cutter.grout_width = config.grout_width;
            cutter.noise_seed = config.noise_seed;
            
            // Determine if component is planar
            Plane plane;
            bool is_coplanar = isMeshCoplanar(soup, mesh_id, plane);
            
            // Check for noise
            bool has_noise = config.hasNoise();
            
            if (has_noise) {
                // Noise promotes planar to Tier 2
                cutter.type = CutterType::NoisedMesh;
                cutter.tri_start = static_cast<uint32_t>(start);
                cutter.tri_count = static_cast<uint32_t>(end - start);
                dispatch.any_noised = true;
                dispatch.tier2_count++;
            } else if (is_coplanar && !config.force_exact) {
                // Tier 1: Planar cutter
                if (config.hasGrout()) {
                    cutter.type = CutterType::PlanarGrout;
                    dispatch.any_grout = true;
                } else {
                    cutter.type = CutterType::Planar;
                }
                
                // Convert Plane to IntPlane
                cutter.plane.a = plane.a;
                cutter.plane.b = plane.b;
                cutter.plane.c = plane.c;
                cutter.plane.d = plane.d;
                
                // Compute grout planes if needed
                if (cutter.type == CutterType::PlanarGrout) {
                    double scale = soup.quantization_scale[0];
                    computeGroutPlanes(cutter.plane, config.grout_width, scale,
                                       cutter.grout_plane_pos, cutter.grout_plane_neg);
                }
                
                dispatch.tier1_count++;
                dispatch.plane_equations.push_back(cutter.plane);
            } else {
                // Tier 2: General mesh
                cutter.type = CutterType::GeneralMesh;
                cutter.tri_start = static_cast<uint32_t>(start);
                cutter.tri_count = static_cast<uint32_t>(end - start);
                dispatch.any_noised = true;
                dispatch.all_planar = false;
                dispatch.tier2_count++;
            }
            
            dispatch.cutters.push_back(cutter);
        }
    }
    
    // Update aggregate flags
    dispatch.all_planar = (dispatch.tier2_count == 0);
}

//=============================================================================
// MAIN STAGE FUNCTION
//=============================================================================

DiagnosticLog stage01_diagnostic(PolygonSoup& soup, 
                                  const PipelineConfig& config,
                                  DispatchResult& dispatch) {
    DiagnosticLog log;
    EMBER_SCOPED_TIMER("Stage 01: Diagnostic");
    
    log.info("=== Stage 01: Diagnostic Gatekeeper ===");
    
    // Step 1: Triangulate all polygons
    log.info("Step 1: Triangulating polygons...");
    size_t tri_count = triangulateAllPolygons(soup, log);
    log.info("  Created " + std::to_string(tri_count) + " triangles");
    
    // Step 2: Remove degenerate triangles
    log.info("Step 2: Removing degenerate triangles...");
    size_t removed = removeDegenerateTriangles(soup, log);
    log.info("  Removed " + std::to_string(removed) + " degenerate triangles");
    
    // Step 3: Normalize winding order
    log.info("Step 3: Normalizing winding order...");
    size_t flipped = normalizeWinding(soup, log);
    log.info("  Flipped " + std::to_string(flipped) + " triangles");
    
    // Step 4: Compute plane equations
    log.info("Step 4: Computing plane equations...");
    computeTrianglePlanes(soup);
    log.info("  Computed planes for all triangles");
    
    // Step 5: Check watertightness
    log.info("Step 5: Checking watertightness...");
    for (uint32_t mesh_id = 0; mesh_id < soup.mesh_count; ++mesh_id) {
        bool watertight = checkWatertightness(soup, mesh_id, log);
        log.info("  Mesh " + std::to_string(mesh_id) + ": " + 
                 (watertight ? "watertight" : "not watertight"));
    }
    
    // Step 6: Detect cutter types
    log.info("Step 6: Detecting cutter types...");
    detectCutterTypes(soup, config, dispatch);
    
    log.info("  Tier 1 (planar) cutters: " + std::to_string(dispatch.tier1_count));
    log.info("  Tier 2 (general) cutters: " + std::to_string(dispatch.tier2_count));
    log.info("  All planar: " + std::string(dispatch.all_planar ? "yes" : "no"));
    log.info("  Fast path eligible: " + std::string(dispatch.isFastPathOnly() ? "yes" : "no"));
    
    // Final validation
    if (soup.triangles.empty()) {
        log.fatal(DiagCategory::Error, "No valid triangles after diagnostic stage");
    }
    
    log.info("=== Stage 01 Complete ===");
    
    return log;
}

} // namespace stages
} // namespace ember
