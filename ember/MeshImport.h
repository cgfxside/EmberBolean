/**
 * @file MeshImport.h
 * @brief Mesh import and quantization for EMBER Boolean engine
 * 
 * EMBER - Exact Mesh Boolean Engine for Robust operations
 * Copyright (c) 2024-2025
 * 
 * This module handles:
 * - Per-axis quantization with aspect ratio detection
 * - Fixed-point coordinate conversion (26-bit precision)
 * - Plane equation computation from triangles
 * - Mesh validation and diagnostic reporting
 * 
 * The key innovation is per-axis quantization that preserves micro-details
 * on meshes with extreme aspect ratios (e.g., landscapes with small height
 * variations).
 */

#pragma once

#include <cstdint>
#include <cmath>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "PolygonSoup.h"

namespace ember {

//=============================================================================
// CONSTANTS
//=============================================================================

/**
 * @brief Fixed-point quantization parameters
 * 
 * We use 26-bit signed integers for coordinates, giving:
 * - Range: [-33,554,431, +33,554,431] (2^25 - 1)
 * - Precision: ~1 unit at scale 1.0
 * - Safety margin: 5 bits below 31-bit signed int limit
 */
constexpr int32_t EMBER_QUANT_BITS = 26;
constexpr int32_t EMBER_QUANT_MAX = (1 << (EMBER_QUANT_BITS - 1)) - 1;  // 33,554,431
constexpr int32_t EMBER_QUANT_MIN = -EMBER_QUANT_MAX;

/**
 * @brief Aspect ratio thresholds for quantization mode selection
 */
constexpr double EMBER_ASPECT_RATIO_PER_AXIS = 1000.0;   // Switch to per-axis
constexpr double EMBER_ASPECT_RATIO_WARNING = 10000.0;   // Log warning

/**
 * @brief Minimum dimension to avoid division by zero
 */
constexpr double EMBER_MIN_DIMENSION = 1e-12;

//=============================================================================
// FORWARD DECLARATIONS
//=============================================================================

struct QuantizationContext;
struct Plane;
struct PolygonSoup;
struct Triangle;

//=============================================================================
// DATA STRUCTURES
//=============================================================================

/**
 * @brief Quantization context for coordinate conversion
 * 
 * This structure captures the transformation between floating-point
 * and fixed-point coordinates. It supports both uniform and per-axis
 * quantization to handle extreme aspect ratios.
 */
struct QuantizationContext {
    double bbox_min[3];     ///< Bounding box minimum per axis
    double bbox_max[3];     ///< Bounding box maximum per axis
    double scale[3];        ///< Per-axis quantization scales
    double inv_scale[3];    ///< Per-axis inverse scales (for dequantization)
    bool per_axis;          ///< True if using non-uniform scaling
    double aspect_ratio_warning_threshold;  ///< Threshold for warnings
    double aspect_ratio_threshold;          ///< Threshold for per-axis mode
    
    /**
     * @brief Default constructor initializes to safe defaults
     */
    QuantizationContext()
        : per_axis(false)
        , aspect_ratio_warning_threshold(EMBER_ASPECT_RATIO_WARNING)
        , aspect_ratio_threshold(EMBER_ASPECT_RATIO_PER_AXIS)
    {
        for (int i = 0; i < 3; ++i) {
            bbox_min[i] = 0.0;
            bbox_max[i] = 0.0;
            scale[i] = 1.0;
            inv_scale[i] = 1.0;
        }
    }
    
    /**
     * @brief Quantize a single vertex
     */
    void quantize(double fx, double fy, double fz,
                  int32_t& ix, int32_t& iy, int32_t& iz) const;
    
    /**
     * @brief Dequantize an integer coordinate back to float
     */
    void dequantize(int32_t ix, int32_t iy, int32_t iz,
                    double& fx, double& fy, double& fz) const;
    
    /**
     * @brief Get the maximum aspect ratio among all axis pairs
     */
    double getMaxAspectRatio() const {
        double rx = bbox_max[0] - bbox_min[0];
        double ry = bbox_max[1] - bbox_min[1];
        double rz = bbox_max[2] - bbox_min[2];
        
        rx = std::max(rx, EMBER_MIN_DIMENSION);
        ry = std::max(ry, EMBER_MIN_DIMENSION);
        rz = std::max(rz, EMBER_MIN_DIMENSION);
        
        double max_ratio = 1.0;
        max_ratio = std::max(max_ratio, rx / ry);
        max_ratio = std::max(max_ratio, rx / rz);
        max_ratio = std::max(max_ratio, ry / rx);
        max_ratio = std::max(max_ratio, ry / rz);
        max_ratio = std::max(max_ratio, rz / rx);
        max_ratio = std::max(max_ratio, rz / ry);
        
        return max_ratio;
    }
};

//=============================================================================
// QUANTIZATION FUNCTIONS
//=============================================================================

/**
 * @brief Compute quantization parameters from polygon soup
 * 
 * Analyzes the bounding box to determine if per-axis quantization
 * is needed based on aspect ratios. Logs warnings for extreme ratios.
 * 
 * FIX 3.7: UNIFIED QUANTIZATION SCALE - SINGLE SOURCE OF TRUTH
 * This function computes THE quantization scale for the entire pipeline.
 * All quantization operations throughout EMBER must use this scale:
 *   - MeshImport::quantizeAll() for initial vertex quantization
 *   - ExactPredicates (via std::llround) for predicate slow-path
 *   - Any coordinate conversion in the pipeline
 * 
 * The scale is stored in PolygonSoup::quantization_scale[] and must be
 * used consistently by all components. Use verifyQuantizationScale() to
 * enforce this invariant.
 * 
 * @param soup Input polygon soup to analyze
 * @param aspect_ratio_threshold Threshold for switching to per-axis quantization
 * @return QuantizationContext with computed scales
 */
QuantizationContext computeQuantization(const PolygonSoup& soup, 
                                         double aspect_ratio_threshold = EMBER_ASPECT_RATIO_PER_AXIS);

/**
 * @brief Quantize all vertices in a polygon soup
 * 
 * Applies quantization to all triangles and stores results in
 * the iv[][] fields. Also updates soup quantization parameters.
 * 
 * FIX 3.7: This function stores the quantization scale in the soup,
 * which becomes the SINGLE SOURCE OF TRUTH for all pipeline operations.
 * 
 * @param soup Polygon soup to quantize (modified in-place)
 * @param ctx Quantization context with scales
 * @return Number of triangles successfully quantized
 */
uint32_t quantizeAll(PolygonSoup& soup, const QuantizationContext& ctx);

/**
 * @brief Verify that the quantization context matches the soup's stored scale
 * 
 * FIX 3.7: ENFORCES unified quantization scale across the pipeline.
 * Call this before any operation that uses quantization to ensure consistency.
 * 
 * @param ctx Quantization context to verify
 * @param soup Polygon soup with stored quantization scale
 * @return true if scales match (within tolerance), false otherwise
 */
bool verifyQuantizationScale(const QuantizationContext& ctx, const PolygonSoup& soup);

//=============================================================================
// PLANE COMPUTATION
//=============================================================================

/**
 * @brief Compute plane equation from three quantized vertices
 * 
 * Uses exact 53-bit integer arithmetic for cross product, then
 * normalizes to 26-bit coefficients.
 * 
 * @param v0 First vertex (quantized)
 * @param v1 Second vertex (quantized)
 * @param v2 Third vertex (quantized)
 * @return Plane structure with 26-bit coefficients
 */
Plane computePlaneFromTriangle(const int32_t v0[3],
                                const int32_t v1[3],
                                const int32_t v2[3]);

/**
 * @brief Compute plane from triangle using floating-point vertices
 * 
 * Convenience wrapper that first quantizes the vertices.
 * 
 * @param v0 First vertex (floating-point)
 * @param v1 Second vertex (floating-point)
 * @param v2 Third vertex (floating-point)
 * @return Plane structure with 26-bit coefficients
 */
Plane computePlaneFromTriangle(const double v0[3],
                                const double v1[3],
                                const double v2[3]);

//=============================================================================
// INTEGER SQRT HELPERS
//=============================================================================

/**
 * @brief Integer square root for int64_t (floor)
 * 
 * Uses Newton's method for fast convergence.
 */
int64_t isqrt(int64_t n);

/**
 * @brief Integer square root for Int128 (floor)
 * 
 * Uses Newton's method with guaranteed convergence for positive inputs.
 * Full 128-bit range, no truncation to int64_t.
 * PORTABLE: Uses ember::Int128 instead of GCC __int128.
 */
Int128 isqrt128(const Int128& n);

//=============================================================================
// NORMALIZATION
//=============================================================================

/**
 * @brief Normalize a 3D vector to 26-bit fixed-point
 * 
 * Takes 53-bit intermediate values and scales to 26-bit output
 * while preserving direction. Uses full Int128 arithmetic (portable).
 * 
 * @param nx, ny, nz Input vector components (53-bit)
 * @param out_nx, out_ny, out_nz Output 26-bit components
 */
void normalizeVector53To26(int64_t nx, int64_t ny, int64_t nz,
                           int32_t& out_nx, int32_t& out_ny, int32_t& out_nz);

//=============================================================================
// VALIDATION
//=============================================================================

/**
 * @brief Validate a polygon soup
 * 
 * Checks for:
 * - Degenerate triangles (zero area)
 * - Invalid coordinates (NaN, Inf)
 * - Quantization overflow
 * 
 * @param soup Polygon soup to validate
 * @param verbose If true, print diagnostic information
 * @return Number of invalid triangles found
 */
uint32_t validateSoup(const PolygonSoup& soup, bool verbose = false);

/**
 * @brief Check if a triangle is degenerate
 * 
 * A triangle is degenerate if its vertices are collinear or coincident.
 * 
 * @param v0 First vertex
 * @param v1 Second vertex
 * @param v2 Third vertex
 * @return True if degenerate
 */
bool isDegenerate(const int32_t v0[3], const int32_t v1[3], const int32_t v2[3]);

//=============================================================================
// DIAGNOSTICS
//=============================================================================

/**
 * @brief Print quantization context information
 */
void printQuantizationContext(const QuantizationContext& ctx, 
                               std::ostream& out = std::cout);

/**
 * @brief Print polygon soup statistics
 */
void printSoupStats(const PolygonSoup& soup,
                    std::ostream& out = std::cout);

/**
 * @brief Get EMBER version string
 */
const char* getEmberVersion();

} // namespace ember
