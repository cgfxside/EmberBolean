/**
 * @file Stage06_ShatterPost.h
 * @brief Stage 06: Shatter Post-Processing for EMBER pipeline
 * 
 * This stage performs connected component extraction and RBD metadata
 * generation for shatter operations. It:
 * 
 * - Extracts connected components from the shattered mesh
 * - Generates piece metadata (names, pivot points)
 * - Assigns piece IDs to output polygons
 * - Computes bounding boxes and volumes for each piece
 * 
 * @author EMBER Boolean Engine Team
 * @version 1.0
 */

#pragma once

#include "ember/PolygonSoup.h"
#include "ember/Diagnostics.h"

#include <vector>
#include <cstdint>
#include <string>
#include <array>
#include <unordered_map>
#include <unordered_set>

namespace ember {
namespace stages {

//=============================================================================
// CONNECTED COMPONENT EXTRACTION
//=============================================================================

/**
 * @brief Connected component information
 */
struct ConnectedComponent {
    uint32_t id;                              ///< Component ID
    std::vector<uint32_t> triangle_indices;   ///< Triangles in this component
    std::array<float, 3> centroid;            ///< Centroid of component
    std::array<float, 3> bbox_min;            ///< Bounding box minimum
    std::array<float, 3> bbox_max;            ///< Bounding box maximum
    float volume;                             ///< Approximate volume
    float surface_area;                       ///< Surface area
    
    ConnectedComponent() 
        : id(0),
          centroid{0.0f, 0.0f, 0.0f},
          bbox_min{0.0f, 0.0f, 0.0f},
          bbox_max{0.0f, 0.0f, 0.0f},
          volume(0.0f),
          surface_area(0.0f) {}
};

/**
 * @brief Extract connected components from mesh
 * 
 * Uses triangle adjacency to find connected components.
 * Two triangles are adjacent if they share an edge.
 * 
 * @param soup Input polygon soup
 * @return Vector of connected components
 */
std::vector<ConnectedComponent> extractConnectedComponents(const PolygonSoup& soup);

/**
 * @brief Build triangle adjacency graph
 * 
 * @param soup Input polygon soup
 * @return Adjacency list (triangle index -> list of adjacent triangles)
 */
std::vector<std::vector<uint32_t>> buildTriangleAdjacency(const PolygonSoup& soup);

/**
 * @brief Edge key for adjacency building
 */
struct ShatterEdgeKey {
    uint32_t v0, v1;
    
    ShatterEdgeKey(uint32_t a, uint32_t b) {
        if (a < b) { v0 = a; v1 = b; }
        else { v0 = b; v1 = a; }
    }
    
    bool operator==(const ShatterEdgeKey& other) const {
        return v0 == other.v0 && v1 == other.v1;
    }
};

/**
 * @brief Hash function for ShatterEdgeKey
 */
struct ShatterEdgeKeyHash {
    size_t operator()(const ShatterEdgeKey& e) const noexcept {
        return (static_cast<size_t>(e.v0) << 16) ^ e.v1;
    }
};

//=============================================================================
// COMPONENT ANALYSIS
//=============================================================================

/**
 * @brief Compute centroid of a connected component
 * 
 * @param comp Connected component
 * @param soup Polygon soup
 */
void computeComponentCentroid(ConnectedComponent& comp, const PolygonSoup& soup);

/**
 * @brief Compute bounding box of a connected component
 * 
 * @param comp Connected component
 * @param soup Polygon soup
 */
void computeComponentBBox(ConnectedComponent& comp, const PolygonSoup& soup);

/**
 * @brief Compute volume of a connected component
 * 
 * Uses tetrahedralization or signed volume method.
 * 
 * @param comp Connected component
 * @param soup Polygon soup
 * @return Approximate volume
 */
float computeComponentVolume(const ConnectedComponent& comp, const PolygonSoup& soup);

/**
 * @brief Compute surface area of a connected component
 * 
 * @param comp Connected component
 * @param soup Polygon soup
 * @return Surface area
 */
float computeComponentSurfaceArea(const ConnectedComponent& comp, 
                                   const PolygonSoup& soup);

//=============================================================================
// PIECE METADATA GENERATION
//=============================================================================

/**
 * @brief Generate piece names
 * 
 * Creates human-readable names for each piece.
 * Format: "piece_001", "piece_002", etc.
 * 
 * @param components Connected components
 * @return Vector of piece names
 */
std::vector<std::string> generatePieceNames(
    const std::vector<ConnectedComponent>& components);

/**
 * @brief Generate pivot points for each piece
 * 
 * Pivot points are used for RBD simulation. Options:
 * - Centroid (default)
 * - Bounding box center
 * - User-specified point
 * 
 * @param components Connected components
 * @param soup Polygon soup
 * @return Vector of pivot points
 */
std::vector<std::array<float, 3>> generatePivotPoints(
    const std::vector<ConnectedComponent>& components,
    const PolygonSoup& soup);

/**
 * @brief RBD metadata for a piece
 */
struct RBDMetadata {
    float mass;                    ///< Piece mass (proportional to volume)
    std::array<float, 3> inertia;  ///< Diagonal inertia tensor
    float friction;                ///< Friction coefficient
    float restitution;             ///< Restitution (bounciness)
    
    RBDMetadata()
        : mass(1.0f),
          inertia{1.0f, 1.0f, 1.0f},
          friction(0.5f),
          restitution(0.3f) {}
};

/**
 * @brief Generate RBD metadata for pieces
 * 
 * @param components Connected components
 * @param soup Polygon soup
 * @return Vector of RBD metadata
 */
std::vector<RBDMetadata> generateRBDMetadata(
    const std::vector<ConnectedComponent>& components,
    const PolygonSoup& soup);

//=============================================================================
// PIECE ASSIGNMENT
//=============================================================================

/**
 * @brief Assign piece IDs to output polygons
 * 
 * Updates soup.piece_ids to map each output polygon to its piece.
 * 
 * @param[in,out] soup Polygon soup
 * @param components Connected components
 */
void assignPieceIds(PolygonSoup& soup, 
                    const std::vector<ConnectedComponent>& components);

//=============================================================================
// MAIN STAGE FUNCTION
//=============================================================================

/**
 * @brief Stage 06: Shatter Post-Processing
 * 
 * Main entry point for shatter post-processing.
 * 
 * Processing flow:
 *   1. Extract connected components
 *   2. Compute component properties (centroid, bbox, volume)
 *   3. Generate piece names
 *   4. Generate pivot points
 *   5. Generate RBD metadata
 *   6. Assign piece IDs to output
 * 
 * @param soup Input/output polygon soup
 *             Input: output_polygons from reconstruction
 *             Output: piece_ids, piece_names, piece_pivots populated
 */
void stage06_shatterPost(PolygonSoup& soup);

//=============================================================================
// STATISTICS
//=============================================================================

/**
 * @brief Shatter post-processing statistics
 */
struct ShatterStats {
    size_t pieces_total = 0;           ///< Total number of pieces
    size_t pieces_small = 0;           ///< Pieces below size threshold
    float total_volume = 0.0f;         ///< Total volume of all pieces
    float avg_piece_volume = 0.0f;     ///< Average piece volume
    double extraction_time_ms = 0.0;   ///< Component extraction time
    double metadata_time_ms = 0.0;     ///< Metadata generation time
    
    void print() const {
        std::printf("[Shatter] Pieces: %zu (small: %zu)\n",
                    pieces_total, pieces_small);
        std::printf("[Shatter] Total volume: %.4f, Avg: %.4f\n",
                    total_volume, avg_piece_volume);
        std::printf("[Shatter] Time: %.2f ms\n", 
                    extraction_time_ms + metadata_time_ms);
    }
};

/**
 * @brief Get statistics from last shatter execution
 */
const ShatterStats& getLastShatterStats();

} // namespace stages
} // namespace ember
