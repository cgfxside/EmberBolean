/**
 * @file ShatterPost.h
 * @brief Shatter post-processing for connected component extraction
 * 
 * This module handles post-processing of shatter operations:
 * - Connected component extraction using parallel DSU
 * - Piece pivot computation (center of mass)
 * - Watertight verification
 * - Name attribute generation
 * 
 * The key algorithm is parallel Disjoint Set Union (DSU) for finding
 * connected components in the shattered mesh.
 */

#pragma once

#include "PolygonSoup.h"
#include "MeshImport.h"
#include "Diagnostics.h"

namespace ember {

//=============================================================================
// CONFIGURATION
//=============================================================================

/**
 * @brief Configuration for shatter post-processing
 */
struct ShatterConfig {
    bool generate_name_attrib = true;   ///< Generate piece names
    bool generate_piece_pivots = true;  ///< Compute center of mass for each piece
    bool verify_watertight = true;      ///< Check if pieces are watertight
    std::string name_prefix = "piece_"; ///< Prefix for generated names
    
    /**
     * @brief Default configuration for standard shatter
     */
    static ShatterConfig standard() {
        ShatterConfig cfg;
        return cfg;
    }
    
    /**
     * @brief Configuration for minimal processing (faster)
     */
    static ShatterConfig minimal() {
        ShatterConfig cfg;
        cfg.verify_watertight = false;
        cfg.generate_piece_pivots = false;
        return cfg;
    }
    
    /**
     * @brief Configuration for full diagnostic output
     */
    static ShatterConfig diagnostic() {
        ShatterConfig cfg;
        cfg.verify_watertight = true;
        cfg.generate_piece_pivots = true;
        cfg.generate_name_attrib = true;
        return cfg;
    }
};

//=============================================================================
// CONNECTED COMPONENT EXTRACTION
//=============================================================================

/**
 * @brief Extract connected components from shattered output using parallel DSU
 * 
 * Uses a parallel Disjoint Set Union (Union-Find) algorithm to identify
 * connected components in the output mesh. Faces are connected if they
 * share an edge.
 * 
 * @param soup Polygon soup with output polygons (modified in-place)
 * @param config Shatter post-processing configuration
 */
void runShatterPost(PolygonSoup& soup, const ShatterConfig& config);

/**
 * @brief Extract connected components with diagnostic logging
 * 
 * @param soup Polygon soup with output polygons (modified in-place)
 * @param config Shatter post-processing configuration
 * @param log Diagnostic log for reporting
 */
void runShatterPost(PolygonSoup& soup, const ShatterConfig& config, DiagnosticLog& log);

//=============================================================================
// PIECE PIVOT COMPUTATION
//=============================================================================

/**
 * @brief Compute piece pivot (center of mass) for each connected component
 * 
 * The pivot is computed as the average of all vertex positions in the
 * piece. This is used by Houdini for transformation operations.
 * 
 * @param soup Polygon soup with piece_ids assigned (modified in-place)
 */
void computePiecePivots(PolygonSoup& soup);

/**
 * @brief Compute piece pivot with bounding box
 * 
 * @param soup Polygon soup with piece_ids assigned
 * @param piece_id Piece to compute pivot for
 * @param[out] pivot Output pivot point
 * @param[out] bbox_min Output bounding box minimum
 * @param[out] bbox_max Output bounding box maximum
 * @return true if piece exists and pivot was computed
 */
bool computePiecePivotAndBounds(const PolygonSoup& soup, int piece_id,
                                 std::array<float, 3>& pivot,
                                 std::array<float, 3>& bbox_min,
                                 std::array<float, 3>& bbox_max);

//=============================================================================
// WATERTIGHT VERIFICATION
//=============================================================================

/**
 * @brief Check if a piece is watertight (no boundary edges)
 * 
 * A mesh is watertight if every edge is shared by exactly two faces.
 * 
 * @param soup Polygon soup with output polygons
 * @param piece_id Piece to check
 * @param[out] boundary_edge_count Number of boundary edges found
 * @return true if piece is watertight
 */
bool isPieceWatertight(const PolygonSoup& soup, int piece_id, 
                       uint32_t& boundary_edge_count);

/**
 * @brief Verify watertightness for all pieces
 * 
 * @param soup Polygon soup with piece_ids assigned
 * @param log Diagnostic log for reporting
 * @return Number of non-watertight pieces
 */
uint32_t verifyAllPiecesWatertight(const PolygonSoup& soup, DiagnosticLog& log);

//=============================================================================
// STATISTICS
//=============================================================================

/**
 * @brief Statistics from shatter post-processing
 */
struct ShatterStats {
    int piece_count = 0;
    uint32_t total_polygons = 0;
    uint32_t watertight_pieces = 0;
    uint32_t non_watertight_pieces = 0;
    double avg_polygons_per_piece = 0.0;
    uint32_t min_polygons = 0;
    uint32_t max_polygons = 0;
    
    void print() const {
        std::printf("[EMBER Shatter] Pieces: %d, Total polygons: %u\n", 
                    piece_count, total_polygons);
        std::printf("[EMBER Shatter] Watertight: %u, Non-watertight: %u\n",
                    watertight_pieces, non_watertight_pieces);
        std::printf("[EMBER Shatter] Polygons/piece: avg=%.1f, min=%u, max=%u\n",
                    avg_polygons_per_piece, min_polygons, max_polygons);
    }
};

/**
 * @brief Compute shatter statistics
 * 
 * @param soup Polygon soup with piece_ids assigned
 * @param[out] stats Output statistics
 */
void computeShatterStats(const PolygonSoup& soup, ShatterStats& stats);

} // namespace ember
