/**
 * @file PolygonSoup.h
 * @brief Core mesh data structures for EMBER Boolean engine
 * 
 * EMBER - Exact Mesh Boolean Engine for Robust operations
 * Copyright (c) 2024-2025
 */

#pragma once

#include <cstdint>
#include <vector>
#include <array>
#include <string>
#include "IntegerTypes.h"

namespace ember {

//=============================================================================
// FORWARD DECLARATIONS
//=============================================================================

struct Plane;
struct IntPlane;
struct ImplicitPoint;
struct CandidatePair;
struct OutputPolygon;
struct SeamSegment;

//=============================================================================
// PLANE REPRESENTATIONS
//=============================================================================

/**
 * @brief Plane equation with 26-bit integer coefficients
 * 
 * Plane equation: a*x + b*y + c*z + d = 0
 * where (a, b, c) is the normal vector and d is the offset.
 * All coefficients are 26-bit signed integers.
 */
struct Plane {
    int32_t a, b, c, d;
    
    Plane() : a(0), b(0), c(0), d(0) {}
    Plane(int32_t a_, int32_t b_, int32_t c_, int32_t d_) 
        : a(a_), b(b_), c(c_), d(d_) {}
};

/**
 * @brief Plane equation with extended precision coefficients
 * 
 * Used for intermediate computations where 53-bit or 80-bit precision
 * is required before normalization to 26-bit.
 */
struct IntPlane {
    int64_t a, b, c;     // 53-bit normal (can be unnormalized)
    Int128 d;            // 80-bit offset (for exact dot products)
    
    IntPlane() : a(0), b(0), c(0), d(Int128::zero()) {}
    IntPlane(int64_t a_, int64_t b_, int64_t c_, const Int128& d_)
        : a(a_), b(b_), c(c_), d(d_) {}
};

//=============================================================================
// IMPLICIT POINT TYPES (Cherchi et al. 2020)
//=============================================================================

/**
 * @brief Type of implicit point for indirect predicates
 */
enum class ImplicitPointType {
    EXPLICIT,    // Original mesh vertex (quantized coordinates)
    LPI,         // Line-Plane Intersection
    TPI          // Triangle-Plane Intersection (3-plane intersection)
};

/**
 * @brief Implicit point descriptor for indirect predicates
 * 
 * Stores the defining geometric elements rather than materialized coordinates.
 * This allows exact evaluation without floating-point errors.
 */
struct ImplicitPoint {
    ImplicitPointType type;
    
    union {
        // EXPLICIT: index in the quantized coordinate array
        uint32_t explicit_idx;
        
        // LPI: line (line0, line1) intersected with plane (plane0, plane1, plane2)
        struct {
            uint32_t line0, line1;           // Line endpoints (vertex indices)
            uint32_t plane0, plane1, plane2; // Plane triangle vertices
        } lpi;
        
        // TPI: three planes intersect at a point
        struct {
            uint32_t plane0[3];  // First plane triangle
            uint32_t plane1[3];  // Second plane triangle
            uint32_t plane2[3];  // Third plane triangle
        } tpi;
    };
    
    ImplicitPoint() : type(ImplicitPointType::EXPLICIT), explicit_idx(0) {}
};

//=============================================================================
// TRIANGLE STRUCTURE
//=============================================================================

/**
 * @brief Single triangle with both float and quantized integer vertices
 */
struct Triangle {
    // Float vertices [vertex][axis] - axis 0=X, 1=Y, 2=Z
    std::array<std::array<float, 3>, 3> v;
    
    // Quantized 26-bit integer vertices [vertex][axis]
    std::array<std::array<int32_t, 3>, 3> iv;
    
    // Mesh ID: 0 = target mesh, 1..N = cutter meshes
    int mesh_id;
    
    // Original Houdini primitive index (for attribute tracking)
    uint32_t src_prim_idx;
    
    // Precomputed plane equation (from quantized vertices)
    Plane plane;
    
    Triangle() : mesh_id(0), src_prim_idx(0) {
        v = {{{{0,0,0}}, {{0,0,0}}, {{0,0,0}}}};
        iv = {{{{0,0,0}}, {{0,0,0}}, {{0,0,0}}}};
    }
};

//=============================================================================
// CANDIDATE PAIR (Broad-phase output)
//=============================================================================

/**
 * @brief Candidate triangle pair for intersection testing
 * 
 * Produced by Stage 02 (Embree BVH broad-phase) or Stage 02b (plane fast-path).
 * These pairs are passed to Stage 03 for exact intersection computation.
 */
struct CandidatePair {
    uint32_t tri_a;      // Triangle index in mesh A (target)
    uint32_t tri_b;      // Triangle index in mesh B (cutter)
    bool intersecting;   // True if triangles actually intersect (from narrow-phase)
    
    // Intersection segments (for CDT input)
    std::vector<std::array<float, 3>> intersection_points;
    
    CandidatePair() : tri_a(0), tri_b(0), intersecting(false) {}
    CandidatePair(uint32_t a, uint32_t b) : tri_a(a), tri_b(b), intersecting(false) {}
};

//=============================================================================
// OUTPUT POLYGON
//=============================================================================

/**
 * @brief Output polygon from Boolean operation
 * 
 * These are the final faces in the result mesh, with winding number
 * classification already applied.
 */
struct OutputPolygon {
    std::vector<uint32_t> vertex_indices;  // Indices into output vertex buffer
    int winding_number;                     // Winding number (for classification)
    int mesh_id;                            // Source mesh ID
    uint32_t src_prim_idx;                  // Original primitive index
    bool is_interior;                       // True if from intersection (cut face)
    
    OutputPolygon() : winding_number(0), mesh_id(0), src_prim_idx(0), is_interior(false) {}
};

//=============================================================================
// SEAM SEGMENT (for Seam operation)
//=============================================================================

/**
 * @brief Seam segment from intersection
 * 
 * Used for the "Seam" Boolean operation which extracts intersection curves.
 */
struct SeamSegment {
    std::array<float, 3> start;     // Start point
    std::array<float, 3> end;       // End point
    uint32_t tri_a;                 // Source triangle A
    uint32_t tri_b;                 // Source triangle B
    
    SeamSegment() : tri_a(0), tri_b(0) {
        start = {{0,0,0}};
        end = {{0,0,0}};
    }
};

//=============================================================================
// POLYGON SOUP (Main mesh container)
//=============================================================================

/**
 * @brief Container for triangle mesh data with EMBER-specific metadata
 * 
 * The PolygonSoup is the primary data structure throughout the EMBER pipeline.
 * It stores triangles with both floating-point and quantized integer vertices,
 * along with candidate pairs, output polygons, and shatter metadata.
 */
struct PolygonSoup {
    // Core triangle data
    std::vector<Triangle> triangles;
    
    // Integer vertices (deduplicated, indexed by triangles)
    // These are the 26-bit quantized coordinates used for exact predicates
    struct IntVertex {
        int32_t x, y, z;
        IntVertex() : x(0), y(0), z(0) {}
        IntVertex(int32_t x_, int32_t y_, int32_t z_) : x(x_), y(y_), z(z_) {}
    };
    std::vector<IntVertex> int_vertices;
    
    // Candidate pairs from broad-phase (Stage 02)
    std::vector<CandidatePair> candidate_pairs;
    
    // Output polygons from Boolean operation (Stage 05)
    std::vector<OutputPolygon> output_polygons;
    
    // Seam segments (for Seam operation)
    std::vector<SeamSegment> seam_segments;
    
    // Mesh metadata
    uint32_t mesh_count = 2;                    // Number of meshes (1 target + N cutters)
    std::vector<uint32_t> mesh_tri_count;       // Triangle count per mesh
    
    // Quantization parameters (set by MeshImport)
    double quantization_scale[3];               // Per-axis scales
    double quantization_inv_scale[3];           // Per-axis inverse scales
    double bbox_min[3], bbox_max[3];            // Bounding box
    bool per_axis_quantization;                 // True if using non-uniform scaling
    
    // Shatter results (Stage 06)
    int piece_count = 0;
    std::vector<int> piece_ids;                 // Piece ID per output polygon
    std::vector<std::string> piece_names;       // Named piece identifiers
    std::vector<std::array<float, 3>> piece_pivots;  // Center of mass per piece
    
    // Constructor
    PolygonSoup() : per_axis_quantization(false) {
        quantization_scale[0] = quantization_scale[1] = quantization_scale[2] = 1.0;
        quantization_inv_scale[0] = quantization_inv_scale[1] = quantization_inv_scale[2] = 1.0;
        bbox_min[0] = bbox_min[1] = bbox_min[2] = 0.0;
        bbox_max[0] = bbox_max[1] = bbox_max[2] = 0.0;
    }
    
    // Compute bounding box from triangle vertices
    void computeBoundingBox(double out_min[3], double out_max[3]) const {
        if (triangles.empty()) {
            for (int i = 0; i < 3; ++i) {
                out_min[i] = 0.0;
                out_max[i] = 0.0;
            }
            return;
        }
        
        // Initialize with first vertex
        for (int i = 0; i < 3; ++i) {
            out_min[i] = triangles[0].v[0][i];
            out_max[i] = triangles[0].v[0][i];
        }
        
        // Expand to contain all vertices
        for (const auto& tri : triangles) {
            for (int v = 0; v < 3; ++v) {
                for (int i = 0; i < 3; ++i) {
                    out_min[i] = std::min(out_min[i], static_cast<double>(tri.v[v][i]));
                    out_max[i] = std::max(out_max[i], static_cast<double>(tri.v[v][i]));
                }
            }
        }
    }
};

//=============================================================================
// BOOLEAN OPERATION ENUM
//=============================================================================

/**
 * @brief Boolean operation types
 */
enum class BooleanOp {
    Union,           // A ∪ B
    Intersection,    // A ∩ B
    DiffAB,          // A - B
    DiffBA,          // B - A
    Shatter,         // Split A into pieces using B
    Seam             // Extract intersection curves only
};

} // namespace ember
