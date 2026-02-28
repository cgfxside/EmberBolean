/**
 * @file Stage04_Classify.h
 * @brief Stage 04: Classification for EMBER pipeline
 * 
 * This stage performs winding number evaluation and label propagation
 * to determine which output polygons should be kept based on the
 * Boolean operation.
 * 
 * Key operations:
 * - Winding number computation for each output polygon
 * - Label consistency propagation across shared edges
 * - Boolean operation filtering (union, intersection, difference, etc.)
 * - Handling of degenerate and boundary cases
 * 
 * @author EMBER Boolean Engine Team
 * @version 1.0
 */

#pragma once

#include "ember/PolygonSoup.h"
#include "ember/Diagnostics.h"
#include "backend/IBooleanBackend.h"

#include <vector>
#include <cstdint>
#include <unordered_map>

namespace ember {
namespace stages {

//=============================================================================
// WINDING NUMBER COMPUTATION
//=============================================================================

/**
 * @brief Winding number state for a polygon
 * 
 * The winding number indicates how many times the mesh wraps around
 * a point. Even = outside, Odd = inside.
 */
struct WindingState {
    int winding_number;      ///< Computed winding number
    bool is_inside;          ///< true if winding number is odd
    bool on_boundary;        ///< true if polygon is on mesh boundary
    bool computed;           ///< true if winding number has been computed
    
    WindingState() 
        : winding_number(0), 
          is_inside(false), 
          on_boundary(false),
          computed(false) {}
};

/**
 * @brief Compute winding number for a point relative to a mesh
 * 
 * Uses ray casting to count intersections with mesh triangles.
 * The winding number is the sum of signed intersections.
 * 
 * @param point Query point (floating-point coordinates)
 * @param mesh_triangles Triangles of the mesh to test against
 * @return Winding number (0 = outside, odd = inside)
 */
int computeWindingNumber(const float point[3], 
                         const std::vector<Triangle>& mesh_triangles);

/**
 * @brief Compute winding number for polygon centroid
 * 
 * @param poly Output polygon
 * @param soup Polygon soup with all triangles
 * @param mesh_id Mesh to test against
 * @return Winding state
 */
WindingState computePolygonWinding(const OutputPolygon& poly,
                                    const PolygonSoup& soup,
                                    uint32_t mesh_id);

/**
 * @brief Batch compute winding numbers for all output polygons
 * 
 * @param soup Input/output polygon soup
 * @return Vector of winding states (one per output polygon)
 */
std::vector<WindingState> computeAllWindingNumbers(const PolygonSoup& soup);

//=============================================================================
// LABEL PROPAGATION
//=============================================================================

/**
 * @brief Edge key for label propagation
 * 
 * Unordered edge used to find adjacent polygons.
 */
struct EdgeKey {
    uint32_t v0, v1;
    
    EdgeKey(uint32_t a, uint32_t b) {
        if (a < b) { v0 = a; v1 = b; }
        else { v0 = b; v1 = a; }
    }
    
    bool operator==(const EdgeKey& other) const {
        return v0 == other.v0 && v1 == other.v1;
    }
};

/**
 * @brief Hash function for EdgeKey
 */
struct EdgeKeyHash {
    size_t operator()(const EdgeKey& e) const noexcept {
        return (static_cast<size_t>(e.v0) << 16) ^ e.v1;
    }
};

/**
 * @brief Propagate labels across shared edges
 * 
 * Ensures consistent inside/outside classification for adjacent polygons.
 * This fixes issues where numerical errors cause inconsistent labels.
 * 
 * @param[in,out] winding_states Winding states to propagate
 * @param soup Polygon soup with output polygons
 */
void propagateLabels(std::vector<WindingState>& winding_states,
                     const PolygonSoup& soup);

/**
 * @brief Build edge-to-polygon adjacency map
 * 
 * @param soup Polygon soup with output polygons
 * @return Map from edge to list of polygon indices
 */
std::unordered_map<EdgeKey, std::vector<size_t>, EdgeKeyHash>
buildEdgeAdjacency(const PolygonSoup& soup);

//=============================================================================
// BOOLEAN FILTERING
//=============================================================================

/**
 * @brief Filter polygons based on Boolean operation
 * 
 * Applies the Boolean operation to determine which polygons to keep.
 * 
 * @param winding_states Winding states for each polygon
 * @param op Boolean operation
 * @param[out] keep_mask Output mask (true = keep polygon)
 * @return Number of polygons to keep
 */
size_t filterByBooleanOp(const std::vector<WindingState>& winding_states,
                         BooleanOp op,
                         std::vector<bool>& keep_mask);

/**
 * @brief Check if polygon should be kept for union operation
 * 
 * Union: Keep polygons that are outside at least one input mesh
 *        (equivalently: not inside all input meshes)
 * 
 * @param winding_per_mesh Winding state for each mesh
 * @param num_meshes Number of input meshes
 * @return true if polygon should be kept
 */
bool keepForUnion(const WindingState winding_per_mesh[], size_t num_meshes);

/**
 * @brief Check if polygon should be kept for intersection operation
 * 
 * Intersection: Keep polygons that are inside ALL input meshes
 * 
 * @param winding_per_mesh Winding state for each mesh
 * @param num_meshes Number of input meshes
 * @return true if polygon should be kept
 */
bool keepForIntersection(const WindingState winding_per_mesh[], 
                         size_t num_meshes);

/**
 * @brief Check if polygon should be kept for difference operation
 * 
 * Difference (A - B): Keep polygons from A that are outside B
 * 
 * @param winding_a Winding state relative to mesh A
 * @param winding_b Winding state relative to mesh B
 * @return true if polygon should be kept
 */
bool keepForDifference(const WindingState& winding_a, 
                       const WindingState& winding_b);

/**
 * @brief Check if polygon should be kept for symmetric difference (XOR)
 * 
 * XOR: Keep polygons that are inside an odd number of meshes
 * 
 * @param winding_per_mesh Winding state for each mesh
 * @param num_meshes Number of input meshes
 * @return true if polygon should be kept
 */
bool keepForXor(const WindingState winding_per_mesh[], size_t num_meshes);

/**
 * @brief Check if polygon should be kept for shatter operation
 * 
 * Shatter: Keep all pieces (no filtering)
 * 
 * @return true (always keep for shatter)
 */
inline bool keepForShatter() { return true; }

//=============================================================================
// BOUNDARY HANDLING
//=============================================================================

/**
 * @brief Handle boundary polygons
 * 
 * Polygons on the boundary of the result need special handling:
 * - Ensure consistent orientation
 * - Fix numerical issues
 * 
 * @param[in,out] winding_states Winding states to update
 * @param soup Polygon soup
 */
void handleBoundaryPolygons(std::vector<WindingState>& winding_states,
                            const PolygonSoup& soup);

//=============================================================================
// MAIN STAGE FUNCTION
//=============================================================================

/**
 * @brief Stage 04: Classification
 * 
 * Main entry point for classification stage.
 * 
 * Processing flow:
 *   1. Compute winding numbers for all output polygons
 *   2. Propagate labels for consistency
 *   3. Filter based on Boolean operation
 *   4. Handle boundary cases
 *   5. Update output polygon list
 * 
 * @param soup Input/output polygon soup
 *             Input: output_polygons from backend
 *             Output: filtered output_polygons
 * @param op Boolean operation for filtering
 */
void stage04_classify(PolygonSoup& soup, BooleanOp op);

//=============================================================================
// STATISTICS
//=============================================================================

/**
 * @brief Classification stage statistics
 */
struct ClassificationStats {
    size_t polygons_input = 0;         ///< Input polygon count
    size_t polygons_output = 0;        ///< Output polygon count
    size_t winding_computed = 0;       ///< Number of winding numbers computed
    size_t labels_propagated = 0;      ///< Number of label propagations
    size_t boundary_handled = 0;       ///< Number of boundary polygons handled
    double compute_time_ms = 0.0;      ///< Winding computation time
    double filter_time_ms = 0.0;       ///< Filtering time
    
    void print() const {
        std::printf("[Classify] Input: %zu, Output: %zu polygons\n",
                    polygons_input, polygons_output);
        std::printf("[Classify] Winding computed: %zu, Labels propagated: %zu\n",
                    winding_computed, labels_propagated);
        std::printf("[Classify] Time: %.2f ms\n", compute_time_ms + filter_time_ms);
    }
};

/**
 * @brief Get statistics from last classification execution
 */
const ClassificationStats& getLastClassificationStats();

} // namespace stages
} // namespace ember
