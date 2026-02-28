/**
 * @file Stage05_Reconstruct.cpp
 * @brief Stage 05: Output Reconstruction implementation
 * 
 * This file implements the output reconstruction stage which builds
 * the final output mesh with attributes.
 * 
 * @author EMBER Boolean Engine Team
 * @version 1.0
 */

#include "Stage05_Reconstruct.h"
#include "ember/PolygonSoup.h"
#include "ember/Diagnostics.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <unordered_map>
#include <vector>

namespace ember {
namespace stages {

//=============================================================================
// COORDINATE DEQUANTIZATION
//=============================================================================

float dequantizeCoordinate(int32_t coord, int axis, const PolygonSoup& soup) {
    return static_cast<float>(soup.dequantize(coord, axis));
}

void materializeCoordinates(const PolygonSoup& soup,
                            std::vector<float>& output_vertices) {
    // Count total vertices needed
    size_t total_verts = 0;
    for (const auto& poly : soup.output_polygons) {
        total_verts += poly.vertex_indices.size();
    }
    
    output_vertices.resize(total_verts * 3);
    size_t idx = 0;
    
    for (const auto& poly : soup.output_polygons) {
        for (uint32_t vert_idx : poly.vertex_indices) {
            // Get the original triangle for coordinate lookup
            // In a full implementation, we'd have a vertex buffer
            // For now, use placeholder
            if (vert_idx < soup.triangles.size()) {
                const auto& tri = soup.triangles[vert_idx];
                // Use first vertex as placeholder
                output_vertices[idx * 3 + 0] = tri.v[0][0];
                output_vertices[idx * 3 + 1] = tri.v[0][1];
                output_vertices[idx * 3 + 2] = tri.v[0][2];
            } else {
                output_vertices[idx * 3 + 0] = 0.0f;
                output_vertices[idx * 3 + 1] = 0.0f;
                output_vertices[idx * 3 + 2] = 0.0f;
            }
            ++idx;
        }
    }
}

//=============================================================================
// NORMAL COMPUTATION
//=============================================================================

void computeFaceNormal(const Triangle& tri, float normal[3]) {
    // Edge vectors
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
    
    // Cross product
    normal[0] = e0[1] * e1[2] - e0[2] * e1[1];
    normal[1] = e0[2] * e1[0] - e0[0] * e1[2];
    normal[2] = e0[0] * e1[1] - e0[1] * e1[0];
}

void normalizeVector(float v[3]) {
    float len = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (len > 1e-10f) {
        v[0] /= len;
        v[1] /= len;
        v[2] /= len;
    }
}

void computeVertexNormals(PolygonSoup& soup) {
    // In a full implementation, this would:
    // 1. Compute face normals for all triangles
    // 2. Average normals at shared vertices
       // 3. Normalize result
    
    // For now, compute face normals only
    for (auto& tri : soup.triangles) {
        float normal[3];
        computeFaceNormal(tri, normal);
        normalizeVector(normal);
        // Store normal in triangle (would need additional storage)
    }
}

//=============================================================================
// BARYCENTRIC COORDINATES
//=============================================================================

bool computeBarycentricCoordinates(const float point[3],
                                   const Triangle& tri,
                                   float bary[3]) {
    // Vectors from v0 to v1 and v2
    float v0[3] = {tri.v[0][0], tri.v[0][1], tri.v[0][2]};
    float v1[3] = {tri.v[1][0], tri.v[1][1], tri.v[1][2]};
    float v2[3] = {tri.v[2][0], tri.v[2][1], tri.v[2][2]};
    
    float e0[3] = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
    float e1[3] = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};
    
    float pv[3] = {point[0] - v0[0], point[1] - v0[1], point[2] - v0[2]};
    
    // Compute dot products
    float d00 = e0[0]*e0[0] + e0[1]*e0[1] + e0[2]*e0[2];
    float d01 = e0[0]*e1[0] + e0[1]*e1[1] + e0[2]*e1[2];
    float d11 = e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2];
    float d20 = pv[0]*e0[0] + pv[1]*e0[1] + pv[2]*e0[2];
    float d21 = pv[0]*e1[0] + pv[1]*e1[1] + pv[2]*e1[2];
    
    float denom = d00 * d11 - d01 * d01;
    if (std::abs(denom) < 1e-10f) {
        bary[0] = 1.0f; bary[1] = 0.0f; bary[2] = 0.0f;
        return false;
    }
    
    bary[1] = (d11 * d20 - d01 * d21) / denom;
    bary[2] = (d00 * d21 - d01 * d20) / denom;
    bary[0] = 1.0f - bary[1] - bary[2];
    
    // Check if inside triangle
    return (bary[0] >= -1e-4f && bary[1] >= -1e-4f && bary[2] >= -1e-4f);
}

void interpolateAttributes(const Triangle& tri,
                           const float bary[3],
                           VertexAttributes& out) {
    (void)tri;  // Would use triangle attributes in full implementation
    
    // Simple interpolation (placeholder)
    out.normal[0] = 0.0f;
    out.normal[1] = 0.0f;
    out.normal[2] = 1.0f;
    out.uv[0] = bary[0];
    out.uv[1] = bary[1];
}

//=============================================================================
// VERTEX DEDUPLICATION
//=============================================================================

size_t deduplicateOutputVertices(PolygonSoup& soup, bool attribute_aware) {
    (void)attribute_aware;
    
    // Build vertex map
    std::unordered_map<VertexKey, uint32_t, VertexKeyHash> vertex_map;
    
    // In a full implementation, this would:
    // 1. Extract all vertices from output polygons
    // 2. Build hash map for deduplication
    // 3. Remap polygon indices
    // 4. Return count of duplicates removed
    
    // For now, return 0 (no deduplication)
    return 0;
}

//=============================================================================
// OUTPUT TRIANGLE CONSTRUCTION
//=============================================================================

size_t triangulateOutputPolygon(const OutputPolygon& poly,
                                const PolygonSoup& soup,
                                std::vector<Triangle>& triangles) {
    (void)soup;
    
    const auto& verts = poly.vertex_indices;
    size_t n = verts.size();
    
    if (n < 3) {
        return 0;  // Degenerate
    }
    
    if (n == 3) {
        // Already a triangle
        Triangle tri;
        tri.mesh_id = poly.mesh_id;
        tri.src_prim_idx = poly.original_tri;
        // Would set vertices here
        triangles.push_back(tri);
        return 1;
    }
    
    // Fan triangulation for n-gons
    size_t created = 0;
    for (size_t i = 1; i + 1 < n; ++i) {
        Triangle tri;
        tri.mesh_id = poly.mesh_id;
        tri.src_prim_idx = poly.original_tri;
        // Would set vertices (0, i, i+1)
        triangles.push_back(tri);
        ++created;
    }
    
    return created;
}

void buildOutputTriangles(PolygonSoup& soup) {
    std::vector<Triangle> output_tris;
    output_tris.reserve(soup.output_polygons.size());
    
    for (const auto& poly : soup.output_polygons) {
        triangulateOutputPolygon(poly, soup, output_tris);
    }
    
    // Replace triangles with output
    soup.triangles = std::move(output_tris);
    
    // Update mesh count
    soup.mesh_count = 1;  // Output is single mesh
    soup.mesh_tri_count = {static_cast<uint32_t>(soup.triangles.size())};
    
    // Clear output polygons (now in triangles)
    soup.output_polygons.clear();
}

//=============================================================================
// STATISTICS
//=============================================================================

static ReconstructionStats g_last_stats;

const ReconstructionStats& getLastReconstructionStats() {
    return g_last_stats;
}

//=============================================================================
// MAIN STAGE FUNCTION
//=============================================================================

void stage05_reconstruct(PolygonSoup& soup) {
    EMBER_SCOPED_TIMER("Stage 05: Reconstruction");
    
    g_last_stats = ReconstructionStats();
    auto start_time = std::chrono::steady_clock::now();
    
    g_last_stats.polygons_input = soup.output_polygons.size();
    
    if (soup.output_polygons.empty()) {
        EMBER_LOG_INFO("Stage 05: No polygons to reconstruct");
        return;
    }
    
    EMBER_LOG_INFO("Stage 05: Reconstructing %zu polygons",
                   soup.output_polygons.size());
    
    // Step 1: Deduplicate vertices
    EMBER_LOG_INFO("Stage 05: Deduplicating vertices...");
    auto dedup_start = std::chrono::steady_clock::now();
    
    g_last_stats.vertices_deduplicated = deduplicateOutputVertices(soup, false);
    
    auto dedup_end = std::chrono::steady_clock::now();
    g_last_stats.dedup_time_ms = std::chrono::duration<double, std::milli>(
        dedup_end - dedup_start).count();
    
    // Step 2: Compute normals
    EMBER_LOG_INFO("Stage 05: Computing normals...");
    computeVertexNormals(soup);
    
    // Step 3: Build output triangles
    EMBER_LOG_INFO("Stage 05: Building output triangles...");
    auto build_start = std::chrono::steady_clock::now();
    
    buildOutputTriangles(soup);
    
    auto build_end = std::chrono::steady_clock::now();
    g_last_stats.materialize_time_ms = std::chrono::duration<double, std::milli>(
        build_end - build_start).count();
    
    g_last_stats.triangles_output = soup.triangles.size();
    
    auto end_time = std::chrono::steady_clock::now();
    g_last_stats.total_time_ms = std::chrono::duration<double, std::milli>(
        end_time - start_time).count();
    
    EMBER_LOG_INFO("Stage 05 complete: %zu output triangles",
                   g_last_stats.triangles_output);
}

} // namespace stages
} // namespace ember
