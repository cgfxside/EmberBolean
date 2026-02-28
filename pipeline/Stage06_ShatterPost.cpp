/**
 * @file Stage06_ShatterPost.cpp
 * @brief Stage 06: Shatter Post-Processing implementation
 * 
 * This file implements the shatter post-processing stage which extracts
 * connected components and generates RBD metadata.
 * 
 * @author EMBER Boolean Engine Team
 * @version 1.0
 */

#include "Stage06_ShatterPost.h"
#include "ember/PolygonSoup.h"
#include "ember/Diagnostics.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace ember {
namespace stages {

//=============================================================================
// TRIANGLE ADJACENCY
//=============================================================================

std::vector<std::vector<uint32_t>> buildTriangleAdjacency(const PolygonSoup& soup) {
    const size_t num_tris = soup.triangles.size();
    std::vector<std::vector<uint32_t>> adjacency(num_tris);
    
    // Map edges to triangles
    std::unordered_map<ShatterEdgeKey, std::vector<uint32_t>, ShatterEdgeKeyHash> edge_to_tris;
    
    for (size_t tri_idx = 0; tri_idx < num_tris; ++tri_idx) {
        const Triangle& tri = soup.triangles[tri_idx];
        
        // Get vertex indices (using integer coordinates as unique IDs)
        uint32_t v0 = static_cast<uint32_t>(tri.iv[0][0] + tri.iv[0][1] * 1000 + tri.iv[0][2] * 1000000);
        uint32_t v1 = static_cast<uint32_t>(tri.iv[1][0] + tri.iv[1][1] * 1000 + tri.iv[1][2] * 1000000);
        uint32_t v2 = static_cast<uint32_t>(tri.iv[2][0] + tri.iv[2][1] * 1000 + tri.iv[2][2] * 1000000);
        
        edge_to_tris[ShatterEdgeKey(v0, v1)].push_back(static_cast<uint32_t>(tri_idx));
        edge_to_tris[ShatterEdgeKey(v1, v2)].push_back(static_cast<uint32_t>(tri_idx));
        edge_to_tris[ShatterEdgeKey(v2, v0)].push_back(static_cast<uint32_t>(tri_idx));
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
    
    return adjacency;
}

//=============================================================================
// CONNECTED COMPONENT EXTRACTION
//=============================================================================

std::vector<ConnectedComponent> extractConnectedComponents(const PolygonSoup& soup) {
    const size_t num_tris = soup.triangles.size();
    if (num_tris == 0) {
        return {};
    }
    
    // Build adjacency
    auto adjacency = buildTriangleAdjacency(soup);
    
    // BFS to find connected components
    std::vector<bool> visited(num_tris, false);
    std::vector<ConnectedComponent> components;
    
    for (size_t start_idx = 0; start_idx < num_tris; ++start_idx) {
        if (visited[start_idx]) continue;
        
        // BFS from this triangle
        ConnectedComponent comp;
        comp.id = static_cast<uint32_t>(components.size());
        
        std::queue<uint32_t> queue;
        queue.push(static_cast<uint32_t>(start_idx));
        visited[start_idx] = true;
        comp.triangle_indices.push_back(static_cast<uint32_t>(start_idx));
        
        while (!queue.empty()) {
            uint32_t cur = queue.front();
            queue.pop();
            
            for (uint32_t neighbor : adjacency[cur]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    queue.push(neighbor);
                    comp.triangle_indices.push_back(neighbor);
                }
            }
        }
        
        components.push_back(std::move(comp));
    }
    
    return components;
}

//=============================================================================
// COMPONENT ANALYSIS
//=============================================================================

void computeComponentCentroid(ConnectedComponent& comp, const PolygonSoup& soup) {
    float cx = 0.0f, cy = 0.0f, cz = 0.0f;
    
    for (uint32_t tri_idx : comp.triangle_indices) {
        const Triangle& tri = soup.triangles[tri_idx];
        for (int v = 0; v < 3; ++v) {
            cx += tri.v[v][0];
            cy += tri.v[v][1];
            cz += tri.v[v][2];
        }
    }
    
    float num_verts = static_cast<float>(comp.triangle_indices.size() * 3);
    comp.centroid = {cx / num_verts, cy / num_verts, cz / num_verts};
}

void computeComponentBBox(ConnectedComponent& comp, const PolygonSoup& soup) {
    if (comp.triangle_indices.empty()) {
        comp.bbox_min = {0.0f, 0.0f, 0.0f};
        comp.bbox_max = {0.0f, 0.0f, 0.0f};
        return;
    }
    
    // Initialize with first vertex
    const Triangle& first_tri = soup.triangles[comp.triangle_indices[0]];
    comp.bbox_min = {first_tri.v[0][0], first_tri.v[0][1], first_tri.v[0][2]};
    comp.bbox_max = comp.bbox_min;
    
    // Expand with all vertices
    for (uint32_t tri_idx : comp.triangle_indices) {
        const Triangle& tri = soup.triangles[tri_idx];
        for (int v = 0; v < 3; ++v) {
            for (int i = 0; i < 3; ++i) {
                comp.bbox_min[i] = std::min(comp.bbox_min[i], tri.v[v][i]);
                comp.bbox_max[i] = std::max(comp.bbox_max[i], tri.v[v][i]);
            }
        }
    }
}

float computeComponentVolume(const ConnectedComponent& comp, const PolygonSoup& soup) {
    // Use signed volume method (tetrahedralization from origin)
    float volume = 0.0f;
    
    for (uint32_t tri_idx : comp.triangle_indices) {
        const Triangle& tri = soup.triangles[tri_idx];
        
        // Signed volume of tetrahedron (origin, v0, v1, v2)
        // V = (1/6) * dot(v0, cross(v1, v2))
        float v0[3] = {tri.v[0][0], tri.v[0][1], tri.v[0][2]};
        float v1[3] = {tri.v[1][0], tri.v[1][1], tri.v[1][2]};
        float v2[3] = {tri.v[2][0], tri.v[2][1], tri.v[2][2]};
        
        float cross[3] = {
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0]
        };
        
        float dot = v0[0] * cross[0] + v0[1] * cross[1] + v0[2] * cross[2];
        volume += dot / 6.0f;
    }
    
    return std::abs(volume);
}

float computeComponentSurfaceArea(const ConnectedComponent& comp, 
                                   const PolygonSoup& soup) {
    float area = 0.0f;
    
    for (uint32_t tri_idx : comp.triangle_indices) {
        const Triangle& tri = soup.triangles[tri_idx];
        
        // Cross product magnitude / 2
        float e0[3] = {
            tri.v[1][0] - tri.v[0][0],
            tri.v[1][1] - tri.v[0][1],
            tri.v[1][2] - tri.v[0][2]
        };
        float e1[3] = {
            tri.v[2][0] - tri.v[0][0],
            tri.v[2][1] - tri.v[0][1],
            tri.v[2][2] - tri.v[0][2]
        };
        
        float cross[3] = {
            e0[1] * e1[2] - e0[2] * e1[1],
            e0[2] * e1[0] - e0[0] * e1[2],
            e0[0] * e1[1] - e0[1] * e1[0]
        };
        
        float cross_len = std::sqrt(cross[0]*cross[0] + 
                                     cross[1]*cross[1] + 
                                     cross[2]*cross[2]);
        area += cross_len * 0.5f;
    }
    
    return area;
}

//=============================================================================
// METADATA GENERATION
//=============================================================================

std::vector<std::string> generatePieceNames(
    const std::vector<ConnectedComponent>& components) {
    std::vector<std::string> names;
    names.reserve(components.size());
    
    for (size_t i = 0; i < components.size(); ++i) {
        std::ostringstream oss;
        oss << "piece_" << std::setw(3) << std::setfill('0') << (i + 1);
        names.push_back(oss.str());
    }
    
    return names;
}

std::vector<std::array<float, 3>> generatePivotPoints(
    const std::vector<ConnectedComponent>& components,
    const PolygonSoup& soup) {
    (void)soup;
    
    std::vector<std::array<float, 3>> pivots;
    pivots.reserve(components.size());
    
    for (const auto& comp : components) {
        // Use centroid as pivot point
        pivots.push_back(comp.centroid);
    }
    
    return pivots;
}

std::vector<RBDMetadata> generateRBDMetadata(
    const std::vector<ConnectedComponent>& components,
    const PolygonSoup& soup) {
    (void)soup;
    
    std::vector<RBDMetadata> metadata;
    metadata.reserve(components.size());
    
    for (const auto& comp : components) {
        RBDMetadata meta;
        
        // Mass proportional to volume
        meta.mass = comp.volume;
        
        // Simple inertia approximation (box inertia)
        float dx = comp.bbox_max[0] - comp.bbox_min[0];
        float dy = comp.bbox_max[1] - comp.bbox_min[1];
        float dz = comp.bbox_max[2] - comp.bbox_min[2];
        
        meta.inertia[0] = (dy*dy + dz*dz) * meta.mass / 12.0f;
        meta.inertia[1] = (dx*dx + dz*dz) * meta.mass / 12.0f;
        meta.inertia[2] = (dx*dx + dy*dy) * meta.mass / 12.0f;
        
        metadata.push_back(meta);
    }
    
    return metadata;
}

//=============================================================================
// PIECE ASSIGNMENT
//=============================================================================

void assignPieceIds(PolygonSoup& soup, 
                    const std::vector<ConnectedComponent>& components) {
    // Initialize piece_ids array
    soup.piece_ids.clear();
    soup.piece_ids.resize(soup.output_polygons.size(), -1);
    
    // Map triangles to pieces
    // In a full implementation, this would use the triangle indices
    // from each component to assign piece IDs to output polygons
    
    // For now, assign piece IDs based on component ID
    for (size_t i = 0; i < soup.output_polygons.size() && i < components.size(); ++i) {
        soup.piece_ids[i] = static_cast<int>(components[i].id);
    }
    
    soup.piece_count = static_cast<int>(components.size());
}

//=============================================================================
// STATISTICS
//=============================================================================

static ShatterStats g_last_stats;

const ShatterStats& getLastShatterStats() {
    return g_last_stats;
}

//=============================================================================
// MAIN STAGE FUNCTION
//=============================================================================

void stage06_shatterPost(PolygonSoup& soup) {
    EMBER_SCOPED_TIMER("Stage 06: Shatter Post-Processing");
    
    g_last_stats = ShatterStats();
    auto start_time = std::chrono::steady_clock::now();
    
    if (soup.triangles.empty()) {
        EMBER_LOG_INFO("Stage 06: No triangles to process");
        return;
    }
    
    EMBER_LOG_INFO("Stage 06: Extracting connected components from %zu triangles",
                   soup.triangles.size());
    
    // Step 1: Extract connected components
    EMBER_LOG_INFO("Stage 06: Extracting connected components...");
    auto extract_start = std::chrono::steady_clock::now();
    
    std::vector<ConnectedComponent> components = extractConnectedComponents(soup);
    
    auto extract_end = std::chrono::steady_clock::now();
    g_last_stats.extraction_time_ms = std::chrono::duration<double, std::milli>(
        extract_end - extract_start).count();
    
    g_last_stats.pieces_total = components.size();
    
    EMBER_LOG_INFO("Stage 06: Found %zu connected components", components.size());
    
    // Step 2: Compute component properties
    EMBER_LOG_INFO("Stage 06: Computing component properties...");
    
    for (auto& comp : components) {
        computeComponentCentroid(comp, soup);
        computeComponentBBox(comp, soup);
        comp.volume = computeComponentVolume(comp, soup);
        comp.surface_area = computeComponentSurfaceArea(comp, soup);
        
        g_last_stats.total_volume += comp.volume;
        
        // Count small pieces (volume < 1% of average)
        // Threshold computed after we know total
    }
    
    if (!components.empty()) {
        g_last_stats.avg_piece_volume = g_last_stats.total_volume / components.size();
        
        // Count small pieces
        float small_threshold = g_last_stats.avg_piece_volume * 0.01f;
        for (const auto& comp : components) {
            if (comp.volume < small_threshold) {
                ++g_last_stats.pieces_small;
            }
        }
    }
    
    // Step 3: Generate metadata
    EMBER_LOG_INFO("Stage 06: Generating piece metadata...");
    auto meta_start = std::chrono::steady_clock::now();
    
    soup.piece_names = generatePieceNames(components);
    soup.piece_pivots = generatePivotPoints(components, soup);
    
    // Generate RBD metadata (could be stored in soup or returned separately)
    auto rbd_metadata = generateRBDMetadata(components, soup);
    (void)rbd_metadata;  // Would be stored or output
    
    auto meta_end = std::chrono::steady_clock::now();
    g_last_stats.metadata_time_ms = std::chrono::duration<double, std::milli>(
        meta_end - meta_start).count();
    
    // Step 4: Assign piece IDs
    EMBER_LOG_INFO("Stage 06: Assigning piece IDs...");
    assignPieceIds(soup, components);
    
    auto end_time = std::chrono::steady_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(
        end_time - start_time).count();
    
    EMBER_LOG_INFO("Stage 06 complete: %zu pieces, %.2f ms",
                   components.size(), total_time);
}

} // namespace stages
} // namespace ember
