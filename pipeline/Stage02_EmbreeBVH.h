/**
 * @file Stage02_EmbreeBVH.h
 * @brief Stage 02: Embree Broad-Phase for EMBER pipeline
 * 
 * This stage builds a BVH (Bounding Volume Hierarchy) using Intel Embree
 * and performs broad-phase collision detection to find candidate triangle
 * pairs that may intersect.
 * 
 * The broad-phase significantly narrows down the number of triangle pairs
 * that need expensive exact intersection tests in Stage 03.
 * 
 * Key features:
 * - High-performance BVH construction using Embree
 * - Self-collision queries for target-cutter pairs
 * - Configurable intersection tolerance
 * - Candidate pair deduplication
 * 
 * @author EMBER Boolean Engine Team
 * @version 1.0
 */

#pragma once

#include "ember/PolygonSoup.h"
#include "ember/Diagnostics.h"

#include <cstdint>
#include <vector>

// Forward declaration for Embree device/scene
struct RTCDevice_t;
struct RTCScene_t;
typedef struct RTCDevice_t* RTCDevice;
typypedef struct RTCScene_t* RTCScene;

namespace ember {
namespace stages {

//=============================================================================
// EMBREE BVH BUILDER
//=============================================================================

/**
 * @brief BVH builder using Intel Embree
 * 
 * Wraps Embree's high-performance BVH construction for triangle meshes.
 * Provides efficient broad-phase collision detection.
 */
class EmbreeBVHBuilder {
public:
    /**
     * @brief Construct BVH builder
     */
    EmbreeBVHBuilder();
    
    /**
     * @brief Destructor cleans up Embree resources
     */
    ~EmbreeBVHBuilder();
    
    /**
     * @brief Build BVH from polygon soup
     * 
     * Creates an Embree scene containing all triangles from the soup.
     * Triangles are stored as separate geometries per mesh to enable
     * target-cutter collision queries.
     * 
     * @param soup Input polygon soup
     * @param log Diagnostic log for errors
     * @return true if BVH was built successfully
     */
    bool build(const PolygonSoup& soup, DiagnosticLog& log);
    
    /**
     * @brief Find candidate pairs between two meshes
     * 
     * Performs broad-phase collision detection between all triangles
     * in mesh0 and all triangles in mesh1.
     * 
     * @param mesh0 First mesh ID (typically target mesh = 0)
     * @param mesh1 Second mesh ID (typically cutter mesh)
     * @param[out] pairs Output candidate pairs
     * @return Number of candidate pairs found
     */
    size_t findCandidatePairs(uint32_t mesh0, uint32_t mesh1,
                              std::vector<CandidatePair>& pairs);
    
    /**
     * @brief Find all candidate pairs for multi-mesh input
     * 
     * Finds all target-cutter and cutter-cutter pairs.
     * 
     * @param soup Input polygon soup (for mesh organization)
     * @param[out] pairs Output candidate pairs
     * @return Number of candidate pairs found
     */
    size_t findAllCandidatePairs(const PolygonSoup& soup,
                                  std::vector<CandidatePair>& pairs);
    
    /**
     * @brief Get Embree device (for advanced usage)
     */
    RTCDevice device() const { return device_; }
    
    /**
     * @brief Get Embree scene (for advanced usage)
     */
    RTCScene scene() const { return scene_; }
    
    /**
     * @brief Check if BVH is valid
     */
    bool isValid() const { return scene_ != nullptr; }

private:
    RTCDevice device_;   ///< Embree device
    RTCScene scene_;     ///< Embree scene containing all geometries
    
    // Mesh geometry IDs in Embree scene
    std::vector<uint32_t> mesh_geom_ids_;
    
    // Triangle count per mesh (for index mapping)
    std::vector<uint32_t> mesh_tri_offsets_;
    
    // Disable copy/move
    EmbreeBVHBuilder(const EmbreeBVHBuilder&) = delete;
    EmbreeBVHBuilder& operator=(const EmbreeBVHBuilder&) = delete;
};

//=============================================================================
// CANDIDATE PAIR DEDUPLICATION
//=============================================================================

/**
 * @brief Remove duplicate candidate pairs
 * 
 * Ensures each unique triangle pair appears only once.
 * Uses sorting + unique for efficiency.
 * 
 * @param[in,out] pairs Input/output candidate pairs
 * @return Number of duplicates removed
 */
size_t deduplicateCandidatePairs(std::vector<CandidatePair>& pairs);

/**
 * @brief Sort candidate pairs for deterministic processing
 * 
 * Sorts by (tri0, tri1) for cache-friendly access patterns.
 * 
 * @param[in,out] pairs Input/output candidate pairs
 */
void sortCandidatePairs(std::vector<CandidatePair>& pairs);

//=============================================================================
// MAIN STAGE FUNCTION
//=============================================================================

/**
 * @brief Stage 02: Embree Broad-Phase
 * 
 * Main entry point for broad-phase collision detection.
 * Builds BVH and finds all candidate triangle pairs.
 * 
 * @param soup Input/output polygon soup
 *             Input: triangles from all meshes
 *             Output: candidate_pairs filled with potential intersections
 */
void stage02_broadPhase(PolygonSoup& soup);

//=============================================================================
// STATISTICS
//=============================================================================

/**
 * @brief Broad-phase statistics
 */
struct BroadPhaseStats {
    double build_time_ms = 0.0;       ///< BVH construction time
    double query_time_ms = 0.0;       ///< Collision query time
    size_t total_pairs = 0;           ///< Total candidate pairs found
    size_t duplicates_removed = 0;    ///< Duplicates eliminated
    size_t mesh_pairs_checked = 0;    ///< Number of mesh pairs queried
    
    void print() const {
        std::printf("[BroadPhase] Build: %.2f ms, Query: %.2f ms\n",
                    build_time_ms, query_time_ms);
        std::printf("[BroadPhase] Pairs: %zu (removed %zu duplicates)\n",
                    total_pairs, duplicates_removed);
    }
};

/**
 * @brief Get statistics from last broad-phase execution
 */
const BroadPhaseStats& getLastBroadPhaseStats();

} // namespace stages
} // namespace ember
