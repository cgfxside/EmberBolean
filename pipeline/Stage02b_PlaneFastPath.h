/**
 * @file Stage02b_PlaneFastPath.h
 * @brief Stage 02b: Plane Fast Path for Tier 1 planar cutters
 * 
 * This stage provides an analytical fast path for planar cutters (Tier 1).
 * Instead of building a BVH and performing expensive mesh-mesh intersection,
 * it classifies each triangle against the plane equation and splits
 * straddling triangles.
 * 
 * Key advantages:
 * - O(n) complexity vs O(n log n) for BVH
 * - No floating-point errors in classification (integer planes)
 * - Direct output without backend processing
 * - Supports grout planes (dual-plane pairs)
 * 
 * This is the optimal path when all cutters are planar and have no noise.
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

// Forward declaration for cutter dispatch
namespace ember {
    struct DispatchResult;
    struct PipelineConfig;
    struct CutterDescriptor;
}

namespace ember {
namespace stages {

//=============================================================================
// TRIANGLE CLASSIFICATION
//=============================================================================

/**
 * @brief Classification of a triangle relative to a plane
 * 
 * Each vertex is classified independently, then combined for the triangle.
 */
enum class PlaneClass : uint8_t {
    POSITIVE = 0,   ///< Above plane (ax + by + cz + d > 0)
    NEGATIVE = 1,   ///< Below plane (ax + by + cz + d < 0)
    ON_PLANE = 2    ///< On plane (ax + by + cz + d == 0)
};

/**
 * @brief Triangle classification result
 * 
 * Combines classifications for all three vertices.
 */
struct TriangleClass {
    PlaneClass v[3];  ///< Classification for each vertex
    
    /**
     * @brief Check if triangle is entirely on one side
     */
    bool isUniform() const {
        return (v[0] == v[1] && v[1] == v[2]);
    }
    
    /**
     * @brief Check if triangle straddles the plane
     */
    bool isStraddling() const {
        return !isUniform() && v[0] != PlaneClass::ON_PLANE &&
               v[1] != PlaneClass::ON_PLANE && v[2] != PlaneClass::ON_PLANE;
    }
    
    /**
     * @brief Check if triangle has vertices on the plane
     */
    bool hasOnPlane() const {
        return v[0] == PlaneClass::ON_PLANE ||
               v[1] == PlaneClass::ON_PLANE ||
               v[2] == PlaneClass::ON_PLANE;
    }
    
    /**
     * @brief Get dominant side (for non-straddling triangles)
     * @return POSITIVE or NEGATIVE for uniform triangles
     */
    PlaneClass dominantSide() const {
        if (v[0] != PlaneClass::ON_PLANE) return v[0];
        if (v[1] != PlaneClass::ON_PLANE) return v[1];
        return v[2];
    }
};

/**
 * @brief Classify a point against a plane using exact integer arithmetic
 * 
 * @param plane Integer plane equation
 * @param point 26-bit integer point coordinates
 * @return Plane classification
 */
PlaneClass classifyPoint(const Plane& plane, const int32_t point[3]);

/**
 * @brief Classify a triangle against a plane
 * 
 * @param tri Triangle to classify
 * @param plane Integer plane equation
 * @return Triangle classification
 */
TriangleClass classifyTriangle(const Triangle& tri, const Plane& plane);

//=============================================================================
// TRIANGLE SPLITTING
//=============================================================================

/**
 * @brief Result of splitting a triangle against a plane
 * 
 * A straddling triangle produces 3 triangles (two on one side, one on other).
 * A triangle with a vertex on the plane produces 2 triangles.
 */
struct SplitResult {
    static constexpr size_t MAX_TRIS = 4;
    
    Triangle triangles[MAX_TRIS];  ///< Output triangles
    size_t num_tris;               ///< Number of output triangles
    bool positive_side[MAX_TRIS];  ///< true if triangle is on positive side
    
    SplitResult() : num_tris(0) {
        for (size_t i = 0; i < MAX_TRIS; ++i) {
            positive_side[i] = false;
        }
    }
};

/**
 * @brief Split a triangle against a plane
 * 
 * For a straddling triangle, computes intersection points and produces
 * new triangles on each side of the plane.
 * 
 * @param tri Triangle to split
 * @param plane Integer plane equation
 * @return Split result with new triangles
 */
SplitResult splitTriangle(const Triangle& tri, const Plane& plane);

/**
 * @brief Compute intersection point of edge with plane
 * 
 * Uses linear interpolation for the intersection point.
 * The result is stored as floating-point coordinates.
 * 
 * @param v0 First edge vertex
 * @param v1 Second edge vertex
 * @param c0 Classification of v0
 * @param c1 Classification of v1
 * @param[out] out Intersection point (floating-point)
 */
void computeEdgePlaneIntersection(const float v0[3], const float v1[3],
                                   PlaneClass c0, PlaneClass c1,
                                   float out[3]);

//=============================================================================
// PLANE FAST PATH PROCESSING
//=============================================================================

/**
 * @brief Process all triangles against a single plane cutter
 * 
 * Classifies each triangle and splits straddling ones.
 * Produces output polygons for the Boolean result.
 * 
 * @param soup Input polygon soup
 * @param plane Plane equation for cutter
 * @param op Boolean operation
 * @param[out] output Output polygons
 */
void processPlaneCutter(const PolygonSoup& soup,
                        const Plane& plane,
                        BooleanOp op,
                        std::vector<OutputPolygon>& output);

/**
 * @brief Process all triangles against a grout cutter (dual planes)
 * 
 * For grout cutters, we have two planes (positive and negative offset).
 * The region between them is the "grout" region.
 * 
 * @param soup Input polygon soup
 * @param base_plane Base plane equation
 * @param grout_pos Positive offset plane
 * @param grout_neg Negative offset plane
 * @param op Boolean operation
 * @param[out] output Output polygons
 */
void processGroutCutter(const PolygonSoup& soup,
                        const Plane& base_plane,
                        const Plane& grout_pos,
                        const Plane& grout_neg,
                        BooleanOp op,
                        std::vector<OutputPolygon>& output);

//=============================================================================
// BOOLEAN FILTERING
//=============================================================================

/**
 * @brief Filter triangles based on Boolean operation
 * 
 * Given a triangle's classification relative to all cutters,
 * determine if it should be in the output.
 * 
 * @param side_per_cutter Side classification for each cutter
 * @param num_cutters Number of cutters
 * @param op Boolean operation
 * @return true if triangle should be in output
 */
bool shouldIncludeTriangle(const bool side_per_cutter[], 
                           size_t num_cutters,
                           BooleanOp op);

/**
 * @brief Compute Boolean result for union operation
 * 
 * Union: Keep triangles that are outside at least one cutter
 *        OR inside all cutters (intersection)
 * 
 * @param inside_bits Bitmask of cutters the triangle is inside
 * @param num_cutters Number of cutters
 * @return true if triangle should be kept
 */
bool filterUnion(uint32_t inside_bits, size_t num_cutters);

/**
 * @brief Compute Boolean result for intersection operation
 * 
 * Intersection: Keep triangles that are inside ALL cutters
 * 
 * @param inside_bits Bitmask of cutters the triangle is inside
 * @param num_cutters Number of cutters
 * @return true if triangle should be kept
 */
bool filterIntersection(uint32_t inside_bits, size_t num_cutters);

/**
 * @brief Compute Boolean result for difference operation
 * 
 * Difference: Keep triangles that are outside ALL cutters
 *             (target minus cutters)
 * 
 * @param inside_bits Bitmask of cutters the triangle is inside
 * @param num_cutters Number of cutters
 * @return true if triangle should be kept
 */
bool filterDifference(uint32_t inside_bits, size_t num_cutters);

//=============================================================================
// MAIN STAGE FUNCTION
//=============================================================================

/**
 * @brief Stage 02b: Plane Fast Path
 * 
 * Main entry point for Tier 1 planar cutter processing.
 * 
 * This stage:
 *   1. Processes each planar cutter
 *   2. Classifies all triangles against plane equations
 *   3. Splits straddling triangles
 *   4. Filters based on Boolean operation
 *   5. Produces output directly (no backend needed)
 * 
 * @param soup Input/output polygon soup
 * @param config Pipeline configuration
 * @param dispatch Cutter dispatch results with plane equations
 */
void stage02b_planeFastPath(PolygonSoup& soup,
                            const PipelineConfig& config,
                            const DispatchResult& dispatch);

//=============================================================================
// STATISTICS
//=============================================================================

/**
 * @brief Fast path statistics
 */
struct FastPathStats {
    size_t triangles_processed = 0;   ///< Total triangles classified
    size_t triangles_split = 0;       ///< Triangles that were split
    size_t triangles_kept = 0;        ///< Triangles in output
    size_t triangles_discarded = 0;   ///< Triangles filtered out
    double process_time_ms = 0.0;     ///< Total processing time
    
    void print() const {
        std::printf("[FastPath] Processed: %zu, Split: %zu\n",
                    triangles_processed, triangles_split);
        std::printf("[FastPath] Kept: %zu, Discarded: %zu\n",
                    triangles_kept, triangles_discarded);
        std::printf("[FastPath] Time: %.2f ms\n", process_time_ms);
    }
};

/**
 * @brief Get statistics from last fast path execution
 */
const FastPathStats& getLastFastPathStats();

} // namespace stages
} // namespace ember
