/**
 * @file Stage05_Reconstruct.h
 * @brief Stage 05: Output Reconstruction for EMBER pipeline
 * 
 * This stage builds the final output mesh with:
 * - Output polygon construction
 * - Vertex attribute interpolation
 * - Coordinate materialization (exact to float conversion)
 * - Normal recomputation
 * - UV coordinate handling
 * 
 * @author EMBER Boolean Engine Team
 * @version 1.0
 */

#pragma once

#include "ember/PolygonSoup.h"
#include "ember/Diagnostics.h"
#include "ember/Plane.h"

#include <vector>
#include <cstdint>
#include <unordered_map>

namespace ember {
namespace stages {

//=============================================================================
// VERTEX DEDUPLICATION
//=============================================================================

/**
 * @brief Vertex key for deduplication
 * 
 * Combines position and attributes for exact matching.
 */
struct VertexKey {
    // Quantized position (for exact matching)
    int32_t px, py, pz;
    
    // Optional attributes (for attribute-aware deduplication)
    int32_t nx, ny, nz;  // Normal (quantized)
    int32_t u, v;        // UV coordinates (quantized)
    
    bool operator==(const VertexKey& other) const {
        return px == other.px && py == other.py && pz == other.pz;
    }
};

/**
 * @brief Hash function for VertexKey
 */
struct VertexKeyHash {
    size_t operator()(const VertexKey& v) const noexcept {
        // Combine hashes using bit mixing
        size_t h1 = static_cast<size_t>(v.px);
        size_t h2 = static_cast<size_t>(v.py);
        size_t h3 = static_cast<size_t>(v.pz);
        return h1 ^ (h2 << 10) ^ (h3 << 20);
    }
};

/**
 * @brief Deduplicate vertices in output
 * 
 * Merges vertices with identical positions (and optionally attributes).
 * Updates polygon vertex indices to reference deduplicated vertices.
 * 
 * @param[in,out] soup Polygon soup with output polygons
 * @param attribute_aware If true, also match attributes
 * @return Number of duplicate vertices removed
 */
size_t deduplicateOutputVertices(PolygonSoup& soup, bool attribute_aware = false);

//=============================================================================
// COORDINATE MATERIALIZATION
//=============================================================================

/**
 * @brief Materialize floating-point coordinates from exact representation
 * 
 * Converts integer/implicit coordinates to floating-point for output.
 * This is the only place where exact arithmetic converts to float.
 * 
 * @param soup Input polygon soup (with quantized vertices)
 * @param[out] output_vertices Output floating-point vertices
 */
void materializeCoordinates(const PolygonSoup& soup,
                            std::vector<float>& output_vertices);

/**
 * @brief Dequantize integer coordinate to float
 * 
 * @param coord Quantized integer coordinate
 * @param axis Axis index (0=X, 1=Y, 2=Z)
 * @param soup Polygon soup with quantization parameters
 * @return Dequantized floating-point coordinate
 */
float dequantizeCoordinate(int32_t coord, int axis, const PolygonSoup& soup);

//=============================================================================
// ATTRIBUTE INTERPOLATION
//=============================================================================

/**
 * @brief Vertex attributes structure
 */
struct VertexAttributes {
    float normal[3];     ///< Vertex normal
    float uv[2];         ///< Texture coordinates
    float color[4];      ///< Vertex color (RGBA)
    int material_id;     ///< Material identifier
    
    VertexAttributes() 
        : normal{0.0f, 0.0f, 1.0f},
          uv{0.0f, 0.0f},
          color{1.0f, 1.0f, 1.0f, 1.0f},
          material_id(0) {}
};

/**
 * @brief Interpolate attributes for a point inside a triangle
 * 
 * Uses barycentric coordinates for interpolation.
 * 
 * @param tri Source triangle with attributes
 * @param bary Barycentric coordinates (u, v, w)
 * @param[out] out Output attributes
 */
void interpolateAttributes(const Triangle& tri,
                           const float bary[3],
                           VertexAttributes& out);

/**
 * @brief Compute barycentric coordinates for a point in a triangle
 * 
 * @param point Query point
 * @param tri Triangle vertices
 * @param[out] bary Output barycentric coordinates
 * @return true if point is inside triangle
 */
bool computeBarycentricCoordinates(const float point[3],
                                   const Triangle& tri,
                                   float bary[3]);

//=============================================================================
// NORMAL COMPUTATION
//=============================================================================

/**
 * @brief Compute face normal for a triangle
 * 
 * @param tri Triangle
 * @param[out] normal Output normal (not normalized)
 */
void computeFaceNormal(const Triangle& tri, float normal[3]);

/**
 * @brief Compute vertex normals by averaging face normals
 * 
 * @param[in,out] soup Polygon soup
 */
void computeVertexNormals(PolygonSoup& soup);

/**
 * @brief Normalize a vector
 * 
 * @param[in,out] v Vector to normalize (in place)
 */
void normalizeVector(float v[3]);

//=============================================================================
// OUTPUT CONSTRUCTION
//=============================================================================

/**
 * @brief Build final output triangles from polygons
 * 
 * Converts output polygons to triangles (if needed) and builds
 * the final triangle list.
 * 
 * @param soup Input/output polygon soup
 */
void buildOutputTriangles(PolygonSoup& soup);

/**
 * @brief Triangulate an n-gon output polygon
 * 
 * Uses fan triangulation for polygons with > 3 vertices.
 * 
 * @param poly Output polygon
 * @param soup Polygon soup for vertex lookup
 * @param[out] triangles Output triangles
 * @return Number of triangles created
 */
size_t triangulateOutputPolygon(const OutputPolygon& poly,
                                const PolygonSoup& soup,
                                std::vector<Triangle>& triangles);

//=============================================================================
// MAIN STAGE FUNCTION
//=============================================================================

/**
 * @brief Stage 05: Output Reconstruction
 * 
 * Main entry point for output reconstruction.
 * 
 * Processing flow:
 *   1. Deduplicate output vertices
 *   2. Materialize floating-point coordinates
 *   3. Interpolate vertex attributes
 *   4. Compute/recompute normals
 *   5. Build final output triangles
 *   6. Validate output
 * 
 * @param soup Input/output polygon soup
 *             Input: output_polygons from classification
 *             Output: Final mesh with attributes
 */
void stage05_reconstruct(PolygonSoup& soup);

//=============================================================================
// STATISTICS
//=============================================================================

/**
 * @brief Reconstruction stage statistics
 */
struct ReconstructionStats {
    size_t polygons_input = 0;         ///< Input polygon count
    size_t triangles_output = 0;       ///< Output triangle count
    size_t vertices_deduplicated = 0;  ///< Duplicate vertices removed
    size_t attributes_interpolated = 0;///< Number of interpolated attributes
    double dedup_time_ms = 0.0;        ///< Deduplication time
    double materialize_time_ms = 0.0;  ///< Coordinate materialization time
    double total_time_ms = 0.0;        ///< Total reconstruction time
    
    void print() const {
        std::printf("[Reconstruct] Input: %zu polygons, Output: %zu triangles\n",
                    polygons_input, triangles_output);
        std::printf("[Reconstruct] Deduplicated: %zu vertices\n",
                    vertices_deduplicated);
        std::printf("[Reconstruct] Time: %.2f ms\n", total_time_ms);
    }
};

/**
 * @brief Get statistics from last reconstruction execution
 */
const ReconstructionStats& getLastReconstructionStats();

} // namespace stages
} // namespace ember
