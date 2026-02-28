/**
 * @file SeamExtract.h
 * @brief Seam extraction from intersection results
 * 
 * This module handles extraction of seam polylines from intersection
 * results. Seams are the curves where two meshes intersect, and they
 * are used for:
 * - Topology reconstruction
 * - Winding consistency checks
 * - Visual debugging
 * - Boolean operation boundary identification
 */

#pragma once

#include "PolygonSoup.h"
#include "Diagnostics.h"

namespace ember {

//=============================================================================
// CONFIGURATION
//=============================================================================

/**
 * @brief Configuration for seam extraction
 */
struct SeamConfig {
    bool extract_ab_seams = true;   ///< Extract A-B intersection seams
    bool extract_aa_seams = false;  ///< Extract A self-intersection seams
    bool extract_bb_seams = false;  ///< Extract B self-intersection seams
    
    double seam_tolerance = 1e-6;   ///< Tolerance for seam vertex merging
    bool chain_segments = true;     ///< Chain segments into continuous polylines
    bool mark_boundaries = true;    ///< Mark boundary vs interior seams
    
    /**
     * @brief Default configuration for standard Boolean operations
     */
    static SeamConfig standard() {
        SeamConfig cfg;
        cfg.extract_ab_seams = true;
        cfg.extract_aa_seams = false;
        cfg.extract_bb_seams = false;
        return cfg;
    }
    
    /**
     * @brief Configuration for self-intersection detection
     */
    static SeamConfig selfIntersection() {
        SeamConfig cfg;
        cfg.extract_ab_seams = true;
        cfg.extract_aa_seams = true;
        cfg.extract_bb_seams = true;
        return cfg;
    }
    
    /**
     * @brief Minimal configuration (faster, less detail)
     */
    static SeamConfig minimal() {
        SeamConfig cfg;
        cfg.chain_segments = false;
        cfg.mark_boundaries = false;
        return cfg;
    }
};

//=============================================================================
// SEAM SEGMENT STRUCTURES
//=============================================================================

/**
 * @brief Extended seam segment with additional metadata
 * 
 * This extends the base SeamSegment from PolygonSoup with additional
 * information for seam processing.
 */
struct SeamSegmentEx {
    uint32_t v0;              // First vertex index
    uint32_t v1;              // Second vertex index
    uint32_t tri0;            // Triangle from first mesh
    uint32_t tri1;            // Triangle from second mesh
    bool is_boundary;         // True if on mesh boundary
    uint8_t seam_type;        // 0=AB, 1=AA, 2=BB
    float length;             // Segment length (for sorting)
    
    SeamSegmentEx() : v0(0), v1(0), tri0(0), tri1(0), 
                      is_boundary(false), seam_type(0), length(0.0f) {}
    
    SeamSegmentEx(uint32_t a, uint32_t b, uint32_t t0, uint32_t t1)
        : v0(a), v1(b), tri0(t0), tri1(t1), 
          is_boundary(false), seam_type(0), length(0.0f) {}
};

/**
 * @brief Seam polyline - chained segments forming a continuous curve
 */
struct SeamPolyline {
    std::vector<uint32_t> vertices;  // Ordered vertex indices
    std::vector<uint32_t> tri0s;     // Triangles from mesh 0 per segment
    std::vector<uint32_t> tri1s;     // Triangles from mesh 1 per segment
    bool is_closed;                  // True if polyline forms a loop
    bool is_boundary;                // True if on mesh boundary
    uint8_t seam_type;               // 0=AB, 1=AA, 2=BB
    
    SeamPolyline() : is_closed(false), is_boundary(false), seam_type(0) {}
    
    /**
     * @brief Get number of segments in polyline
     */
    size_t segmentCount() const {
        return vertices.empty() ? 0 : vertices.size() - 1;
    }
    
    /**
     * @brief Get total vertex count
     */
    size_t vertexCount() const {
        return vertices.size();
    }
};

//=============================================================================
// SEAM EXTRACTION
//=============================================================================

/**
 * @brief Extract seam polylines from intersection results
 * 
 * Processes candidate pairs and intersection results to extract
 * seam segments, then chains them into continuous polylines.
 * 
 * @param soup Polygon soup with candidate pairs and triangles
 * @param config Seam extraction configuration
 */
void extractSeams(PolygonSoup& soup, const SeamConfig& config);

/**
 * @brief Extract seams with diagnostic logging
 * 
 * @param soup Polygon soup with candidate pairs
 * @param config Seam extraction configuration
 * @param log Diagnostic log for reporting
 */
void extractSeams(PolygonSoup& soup, const SeamConfig& config, DiagnosticLog& log);

/**
 * @brief Extract seam segments from a single candidate pair
 * 
 * @param soup Polygon soup with triangle data
 * @param pair Candidate pair to process
 * @param[out] segments Output seam segments
 * @return Number of segments extracted
 */
uint32_t extractSeamsFromPair(const PolygonSoup& soup, 
                               const CandidatePair& pair,
                               std::vector<SeamSegmentEx>& segments);

//=============================================================================
// SEGMENT CHAINING
//=============================================================================

/**
 * @brief Chain intersection segments into continuous polylines
 * 
 * Takes a collection of unordered seam segments and connects them
 * into continuous polylines. Handles:
 * - Open polylines (start and end at boundary)
 * - Closed loops (seams forming complete curves)
 * - Branching points (multiple seams meeting at a vertex)
 * 
 * @param segments Seam segments to chain (modified in-place, consumed)
 * @param[out] polylines Output chained polylines
 */
void chainSeamSegments(std::vector<SeamSegmentEx>& segments,
                       std::vector<SeamPolyline>& polylines);

/**
 * @brief Chain segments with explicit chaining configuration
 * 
 * @param segments Seam segments to chain
 * @param tolerance Vertex merging tolerance
 * @param allow_branching Whether to handle T-junctions
 * @param[out] polylines Output chained polylines
 */
void chainSeamSegmentsEx(const std::vector<SeamSegmentEx>& segments,
                         double tolerance,
                         bool allow_branching,
                         std::vector<SeamPolyline>& polylines);

//=============================================================================
// SEAM UTILITIES
//=============================================================================

/**
 * @brief Merge duplicate seam vertices within tolerance
 * 
 * @param soup Polygon soup with seam segments
 * @param tolerance Merge tolerance
 * @return Number of vertices merged
 */
uint32_t mergeSeamVertices(PolygonSoup& soup, double tolerance);

/**
 * @brief Compute seam segment length
 * 
 * @param seg Seam segment
 * @param soup Polygon soup with vertex positions
 * @return Segment length
 */
double computeSeamLength(const SeamSegment& seg, const PolygonSoup& soup);

/**
 * @brief Filter seams by length
 * 
 * @param segments Seam segments to filter
 * @param min_length Minimum length to keep
 * @param[out] filtered Output filtered segments
 */
void filterSeamsByLength(const std::vector<SeamSegmentEx>& segments,
                         double min_length,
                         std::vector<SeamSegmentEx>& filtered);

//=============================================================================
// SEAM STATISTICS
//=============================================================================

/**
 * @brief Statistics from seam extraction
 */
struct SeamStats {
    uint32_t total_segments = 0;
    uint32_t ab_segments = 0;
    uint32_t aa_segments = 0;
    uint32_t bb_segments = 0;
    uint32_t boundary_segments = 0;
    uint32_t total_polylines = 0;
    uint32_t closed_loops = 0;
    uint32_t open_chains = 0;
    double total_length = 0.0;
    
    void print() const {
        std::printf("[EMBER Seams] Segments: %u (AB:%u AA:%u BB:%u)\n",
                    total_segments, ab_segments, aa_segments, bb_segments);
        std::printf("[EMBER Seams] Polylines: %u (closed:%u open:%u)\n",
                    total_polylines, closed_loops, open_chains);
        std::printf("[EMBER Seams] Total length: %.6f\n", total_length);
    }
};

/**
 * @brief Compute seam statistics
 * 
 * @param segments Seam segments
 * @param polylines Chained polylines
 * @param[out] stats Output statistics
 */
void computeSeamStats(const std::vector<SeamSegmentEx>& segments,
                      const std::vector<SeamPolyline>& polylines,
                      SeamStats& stats);

} // namespace ember
