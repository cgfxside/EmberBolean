/**
 * @file Stage01_Diagnostic.h
 * @brief Stage 01: Diagnostic Gatekeeper for EMBER pipeline
 * 
 * This stage validates and prepares input meshes for Boolean operations:
 * - Triangulates all polygons (n-gons → triangles)
 * - Removes degenerate triangles (zero area)
 * - Normalizes winding order (consistent orientation)
 * - Checks watertightness (closed mesh requirement)
 * - Performs cutter type detection (Tier 1 planar vs Tier 2 general)
 * 
 * The diagnostic gatekeeper ensures that downstream stages receive
 * clean, valid input and provides detailed diagnostics for debugging.
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
}

namespace ember {
namespace stages {

//=============================================================================
// TRIANGULATION
//=============================================================================

/**
 * @brief Triangulate all polygons in the input soup
 * 
 * Converts n-gons to triangles using ear-clipping or fan triangulation.
 * All output triangles are added to soup.triangles.
 * 
 * @param soup Input/output polygon soup
 * @param log Diagnostic log for tracking issues
 * @return Number of triangles created
 */
size_t triangulateAllPolygons(PolygonSoup& soup, DiagnosticLog& log);

/**
 * @brief Triangulate a single n-gon using fan triangulation
 * 
 * Creates triangles in fan pattern from vertex 0:
 *   (0,1,2), (0,2,3), (0,3,4), ...
 * 
 * @param vertices Array of n vertex indices
 * @param mesh_id Source mesh ID
 * @param src_prim_idx Original primitive index
 * @param soup Output soup to add triangles to
 * @return Number of triangles created
 */
size_t fanTriangulate(const std::vector<uint32_t>& vertices,
                      int mesh_id,
                      uint32_t src_prim_idx,
                      PolygonSoup& soup);

//=============================================================================
// DEGENERATE TRIANGLE REMOVAL
//=============================================================================

/**
 * @brief Remove degenerate triangles from the soup
 * 
 * A triangle is degenerate if it has zero area (collinear vertices
 * or duplicate vertices). These triangles are removed as they
 * cause numerical issues in predicates.
 * 
 * @param soup Input/output polygon soup
 * @param log Diagnostic log for tracking removed triangles
 * @return Number of triangles removed
 */
size_t removeDegenerateTriangles(PolygonSoup& soup, DiagnosticLog& log);

/**
 * @brief Check if a triangle is degenerate (zero area)
 * 
 * Uses exact integer arithmetic for robustness:
 * - Computes 2x signed area using cross product
 * - Returns true if area is exactly zero
 * 
 * @param tri Triangle to check
 * @return true if triangle is degenerate
 */
bool isTriangleDegenerate(const Triangle& tri);

/**
 * @brief Compute 2x signed area of triangle (for orientation)
 * 
 * @param v0 First vertex (integer coordinates)
 * @param v1 Second vertex
 * @param v2 Third vertex
 * @return 2x signed area (positive = CCW, negative = CW, zero = degenerate)
 */
int64_t triangleSignedArea2x(const int32_t v0[3], const int32_t v1[3], 
                              const int32_t v2[3]);

//=============================================================================
// WINDING NORMALIZATION
//=============================================================================

/**
 * @brief Normalize winding order for all triangles
 * 
 * Ensures consistent winding order (counter-clockwise when viewed
 * from outside). This is critical for correct plane equations
 * and winding number computation.
 * 
 * @param soup Input/output polygon soup
 * @param log Diagnostic log for tracking flipped triangles
 * @return Number of triangles flipped
 */
size_t normalizeWinding(PolygonSoup& soup, DiagnosticLog& log);

/**
 * @brief Compute expected winding order for a mesh
 * 
 * Analyzes mesh to determine if it should be CW or CCW.
 * Uses heuristics based on bounding box and face normals.
 * 
 * @param soup Input polygon soup
 * @param mesh_id Mesh to analyze
 * @return true if mesh should be CCW (positive normal pointing outward)
 */
bool computeExpectedWinding(const PolygonSoup& soup, uint32_t mesh_id);

//=============================================================================
// WATERTIGHTNESS CHECK
//=============================================================================

/**
 * @brief Check if a mesh is watertight (closed, no boundary edges)
 * 
 * A watertight mesh has every edge shared by exactly two triangles
 * with opposite orientation. This is required for correct winding
 * number computation.
 * 
 * @param soup Input polygon soup
 * @param mesh_id Mesh to check
 * @param log Diagnostic log for boundary edges found
 * @return true if mesh is watertight
 */
bool checkWatertightness(const PolygonSoup& soup, uint32_t mesh_id, 
                         DiagnosticLog& log);

/**
 * @brief Edge structure for watertightness checking
 */
struct Edge {
    uint32_t v0, v1;  // Vertex indices (ordered for hashing)
    
    Edge(uint32_t a, uint32_t b) {
        if (a < b) { v0 = a; v1 = b; }
        else { v0 = b; v1 = a; }
    }
    
    bool operator==(const Edge& other) const {
        return v0 == other.v0 && v1 == other.v1;
    }
};

/**
 * @brief Hash function for Edge
 */
struct EdgeHash {
    size_t operator()(const Edge& e) const noexcept {
        return (static_cast<size_t>(e.v0) << 16) ^ e.v1;
    }
};

//=============================================================================
// PLANE COMPUTATION
//=============================================================================

/**
 * @brief Compute plane equations for all triangles
 * 
 * Precomputes integer plane equations (ax + by + cz + d = 0) for
 * all triangles. These are used in predicates and classification.
 * 
 * @param soup Input/output polygon soup
 */
void computeTrianglePlanes(PolygonSoup& soup);

//=============================================================================
// CUTTER TYPE DETECTION
//=============================================================================

/**
 * @brief Detect cutter types and perform dispatch classification
 * 
 * Analyzes all connected components of cutter meshes and classifies:
 * - Tier 1: Planar cutters (all faces coplanar)
 * - Tier 1G: Planar with grout (dual plane pair)
 * - Tier 2: Noised meshes (tessellated noise)
 * - Tier 2: General meshes (arbitrary topology)
 * 
 * @param soup Input polygon soup
 * @param config Pipeline configuration
 * @param[out] dispatch Output dispatch result
 */
void detectCutterTypes(const PolygonSoup& soup, 
                       const PipelineConfig& config,
                       DispatchResult& dispatch);

/**
 * @brief Check if all triangles in a mesh are coplanar
 * 
 * @param soup Input polygon soup
 * @param mesh_id Mesh to check
 * @param[out] plane Output plane equation (if coplanar)
 * @return true if all triangles are coplanar
 */
bool isMeshCoplanar(const PolygonSoup& soup, uint32_t mesh_id, Plane& plane);

/**
 * @brief Compute connected components of a mesh
 * 
 * Groups triangles into connected components based on shared edges.
 * 
 * @param soup Input polygon soup
 * @param mesh_id Mesh to analyze
 * @return Vector of component assignments (one per triangle)
 */
std::vector<uint32_t> computeConnectedComponents(const PolygonSoup& soup, 
                                                  uint32_t mesh_id);

//=============================================================================
// MAIN STAGE FUNCTION
//=============================================================================

/**
 * @brief Stage 01: Diagnostic Gatekeeper
 * 
 * Main entry point for diagnostic stage. Performs all validation
 * and preparation steps, then performs cutter type detection.
 * 
 * Processing order:
 *   1. Triangulate all polygons
 *   2. Remove degenerate triangles
 *   3. Normalize winding order
 *   4. Compute plane equations
 *   5. Check watertightness
 *   6. Detect cutter types
 * 
 * @param soup Input/output polygon soup (modified in place)
 * @param config Pipeline configuration
 * @param[out] dispatch Output dispatch result for cutter classification
 * @return Diagnostic log with all issues found
 */
DiagnosticLog stage01_diagnostic(PolygonSoup& soup, 
                                  const PipelineConfig& config,
                                  DispatchResult& dispatch);

} // namespace stages
} // namespace ember
