/**
 * @file LabelConsistency.h
 * @brief Label propagation for efficient Boolean classification
 * 
 * This module implements Sheng-style label propagation for efficient
 * classification of output faces. Instead of computing winding numbers
 * for every face independently, we propagate classifications through
 * face adjacency, exploiting spatial coherence.
 * 
 * The key insight is that adjacent faces typically have the same
 * classification (inside/outside), so we can propagate labels from
 * seed faces rather than computing from scratch for each face.
 */

#pragma once

#include "PolygonSoup.h"
#include "Diagnostics.h"

namespace ember {

//=============================================================================
// LABEL RESULT STRUCTURE
//=============================================================================

/**
 * @brief Result of label propagation
 * 
 * Contains per-face labels and metadata about the labeling process.
 */
struct LabelResult {
    std::vector<int> face_labels;   ///< Per output polygon label (-1 = unlabeled)
    int num_labels = 0;             ///< Number of distinct labels
    int num_propagated = 0;         ///< Number of faces labeled by propagation
    int num_seeds = 0;              ///< Number of seed faces used
    
    /**
     * @brief Check if a face has a valid label
     */
    bool hasLabel(size_t face_idx) const {
        return face_idx < face_labels.size() && face_labels[face_idx] >= 0;
    }
    
    /**
     * @brief Get label for a face
     */
    int getLabel(size_t face_idx) const {
        return (face_idx < face_labels.size()) ? face_labels[face_idx] : -1;
    }
};

//=============================================================================
// WINDING NUMBER CLASSIFICATION
//=============================================================================

/**
 * @brief Winding number classification for a point relative to a mesh
 * 
 * The winding number is a robust inside/outside test that handles
 * non-watertight meshes better than ray casting.
 * 
 * Winding number interpretation:
 * - 0: Point is outside the mesh
 * - ±1: Point is inside a watertight mesh
 * - Other values: Point is inside a non-watertight mesh
 */
struct WindingClassification {
    int winding_number = 0;     ///< Computed winding number
    bool is_inside = false;     ///< Simplified inside/outside
    double confidence = 1.0;    ///< Confidence in classification (0-1)
    
    /**
     * @brief Create from winding number
     */
    static WindingClassification fromWindingNumber(int wn) {
        WindingClassification wc;
        wc.winding_number = wn;
        wc.is_inside = (wn != 0);
        return wc;
    }
};

/**
 * @brief Compute winding number of a point relative to a mesh
 * 
 * Uses the solid angle method for robust computation.
 * 
 * @param soup Polygon soup with mesh triangles
 * @param mesh_id Mesh to test against
 * @param point Query point
 * @return Winding classification
 */
WindingClassification computeWindingNumber(const PolygonSoup& soup,
                                            int mesh_id,
                                            const double point[3]);

/**
 * @brief Compute winding number for a triangle centroid
 * 
 * @param soup Polygon soup with mesh triangles
 * @param mesh_id Mesh to test against
 * @param tri_idx Triangle index
 * @return Winding classification
 */
WindingClassification computeWindingNumberForTriangle(const PolygonSoup& soup,
                                                       int mesh_id,
                                                       uint32_t tri_idx);

//=============================================================================
// LABEL PROPAGATION
//=============================================================================

/**
 * @brief Sheng-style label propagation for efficient classification
 * 
 * Propagates winding number classifications through face adjacency.
 * This is much faster than computing winding numbers for every face
 * independently, especially for large meshes.
 * 
 * Algorithm:
 * 1. Build face adjacency graph
 * 2. Select seed faces (e.g., faces with known classification)
 * 3. BFS propagate labels to adjacent faces
 * 4. For unlabeled faces, compute winding number directly
 * 
 * @param soup Polygon soup with output polygons
 * @return LabelResult with per-face labels
 */
LabelResult propagateLabels(const PolygonSoup& soup);

/**
 * @brief Propagate labels with diagnostic logging
 * 
 * @param soup Polygon soup with output polygons
 * @param log Diagnostic log for reporting
 * @return LabelResult with per-face labels
 */
LabelResult propagateLabels(const PolygonSoup& soup, DiagnosticLog& log);

/**
 * @brief BFS propagation from seed faces
 * 
 * Propagates a label to all reachable faces through adjacency.
 * Stops at faces that are separated by seams or have different
 * winding classifications.
 * 
 * @param labels Face labels array (modified in-place)
 * @param start_face Starting face for propagation
 * @param label Label to propagate
 * @param soup Polygon soup with adjacency information
 * @return Number of faces labeled
 */
int propagateLabelBFS(std::vector<int>& labels, int start_face, 
                      int label, const PolygonSoup& soup);

/**
 * @brief BFS propagation with adjacency precomputed
 * 
 * @param labels Face labels array (modified in-place)
 * @param start_face Starting face for propagation
 * @param label Label to propagate
 * @param adjacency Precomputed face adjacency
 * @return Number of faces labeled
 */
int propagateLabelBFS(std::vector<int>& labels, int start_face, 
                      int label, const std::vector<std::vector<uint32_t>>& adjacency);

//=============================================================================
// SEED SELECTION
//=============================================================================

/**
 * @brief Strategy for selecting seed faces
 */
enum class SeedStrategy {
    Random,         ///< Random selection
    SpatialGrid,    ///< Regular grid sampling
    Extrema,        ///< Faces at bounding box extrema
    LargestComponent ///< From largest connected component
};

/**
 * @brief Select seed faces for label propagation
 * 
 * @param soup Polygon soup with output polygons
 * @param strategy Seed selection strategy
 * @param num_seeds Number of seeds to select
 * @param[out] seeds Output seed face indices
 */
void selectSeedFaces(const PolygonSoup& soup,
                     SeedStrategy strategy,
                     int num_seeds,
                     std::vector<uint32_t>& seeds);

/**
 * @brief Compute face adjacency for label propagation
 * 
 * Two faces are adjacent if they share an edge that is NOT a seam.
 * Seam edges separate different labels.
 * 
 * @param soup Polygon soup with output polygons
 * @param[out] adjacency Output adjacency lists
 */
void computeFaceAdjacencyForPropagation(const PolygonSoup& soup,
                                        std::vector<std::vector<uint32_t>>& adjacency);

//=============================================================================
// BOOLEAN CLASSIFICATION
//=============================================================================

/**
 * @brief Classification for Boolean operations
 * 
 * For mesh A and mesh B, each face can be classified as:
 * - A_in_B: Face of A inside B
 * - A_out_B: Face of A outside B
 * - B_in_A: Face of B inside A
 * - B_out_A: Face of B outside A
 */
enum class BooleanClass {
    Unknown = -1,
    A_in_B = 0,     // A inside B (for A - B, keep)
    A_out_B = 1,    // A outside B (for A ∩ B, discard)
    B_in_A = 2,     // B inside A (for A ∪ B, keep)
    B_out_A = 3     // B outside A (for A - B, discard)
};

/**
 * @brief Classify faces for Boolean operation
 * 
 * Uses label propagation to efficiently classify all faces
 * for a given Boolean operation.
 * 
 * @param soup Polygon soup with output polygons
 * @param op Boolean operation type
 * @param[out] classifications Output classifications
 * @return Number of faces classified
 */
int classifyForBoolean(const PolygonSoup& soup,
                       BooleanOp op,
                       std::vector<BooleanClass>& classifications);

/**
 * @brief Filter faces by Boolean classification
 * 
 * @param soup Polygon soup with output polygons
 * @param classifications Per-face classifications
 * @param keep_class Classes to keep
 * @param[out] filtered_indices Indices of faces to keep
 */
void filterByClassification(const PolygonSoup& soup,
                            const std::vector<BooleanClass>& classifications,
                            const std::vector<BooleanClass>& keep_class,
                            std::vector<uint32_t>& filtered_indices);

//=============================================================================
// STATISTICS
//=============================================================================

/**
 * @brief Statistics from label propagation
 */
struct LabelStats {
    int total_faces = 0;
    int labeled_faces = 0;
    int unlabeled_faces = 0;
    int num_labels = 0;
    int num_seeds = 0;
    double propagation_ratio = 0.0;
    double avg_component_size = 0.0;
    
    void print() const {
        std::printf("[EMBER Labels] Faces: %d total, %d labeled (%.1f%%)\n",
                    total_faces, labeled_faces, propagation_ratio * 100);
        std::printf("[EMBER Labels] Labels: %d, Seeds: %d, Avg component: %.1f\n",
                    num_labels, num_seeds, avg_component_size);
    }
};

/**
 * @brief Compute label propagation statistics
 * 
 * @param result Label propagation result
 * @param[out] stats Output statistics
 */
void computeLabelStats(const LabelResult& result, LabelStats& stats);

} // namespace ember
