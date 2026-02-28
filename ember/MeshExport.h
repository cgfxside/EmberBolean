/**
 * @file MeshExport.h
 * @brief Mesh export and output coordinate materialization for EMBER Boolean engine
 * 
 * This module handles:
 * - Converting internal PolygonSoup representation back to output format
 * - Dequantizing 26-bit fixed-point coordinates to double-precision floats
 * - Validating output before returning to Houdini
 * - Materializing exact homogeneous coordinates into displayable vertices
 * 
 * The key operation is the conversion from exact arithmetic representations
 * (26-bit integers, homogeneous coordinates) back to floating-point for
 * downstream consumption.
 */

#pragma once

#include "PolygonSoup.h"
#include "MeshImport.h"
#include "Diagnostics.h"

namespace ember {

//=============================================================================
// OUTPUT COORDINATE MATERIALIZATION
//=============================================================================

/**
 * @brief Convert internal PolygonSoup representation back to output format
 * 
 * This is where exact homogeneous coordinates become double-precision floats.
 * All implicit points (LPI, TPI) are evaluated and converted to Cartesian
 * coordinates. This is the only place where exact arithmetic converts to
 * floating-point in the output pipeline.
 * 
 * @param soup Polygon soup with output polygons and implicit points
 * @param ctx Quantization context for dequantization
 */
void materializeOutputCoordinates(PolygonSoup& soup, const QuantizationContext& ctx);

/**
 * @brief Dequantize all vertices in output polygons
 * 
 * Converts 26-bit fixed-point integer coordinates back to floating-point
 * using the per-axis or uniform inverse scale factors from the quantization
 * context.
 * 
 * @param soup Polygon soup with quantized output polygons (modified in-place)
 * @param ctx Quantization context with inverse scales
 */
void dequantizeOutputPolygons(PolygonSoup& soup, const QuantizationContext& ctx);

/**
 * @brief Materialize a single implicit point to double-precision coordinates
 * 
 * Evaluates LPI (Line-Plane Intersection) or TPI (Triangle-Plane Intersection)
 * points and converts homogeneous coordinates to Cartesian.
 * 
 * @param point Implicit point to materialize
 * @param soup Source polygon soup with triangle data
 * @param ctx Quantization context
 * @param[out] out_x Output X coordinate
 * @param[out] out_y Output Y coordinate
 * @param[out] out_z Output Z coordinate
 * @return true if materialization succeeded
 */
bool materializeImplicitPoint(const ImplicitPoint& point,
                               const PolygonSoup& soup,
                               const QuantizationContext& ctx,
                               double& out_x, double& out_y, double& out_z);

//=============================================================================
// OUTPUT VALIDATION
//=============================================================================

/**
 * @brief Validate output before returning to Houdini
 * 
 * Checks for:
 * - Degenerate output polygons (zero area)
 * - Invalid vertex indices
 * - Inconsistent winding
 * - NaN/Inf coordinates
 * 
 * @param soup Output polygon soup to validate
 * @param log Diagnostic log for reporting issues
 * @return true if output passes all validation checks
 */
bool validateOutput(const PolygonSoup& soup, DiagnosticLog& log);

/**
 * @brief Validation statistics for output mesh
 */
struct OutputValidationStats {
    uint32_t total_polygons = 0;
    uint32_t degenerate_polygons = 0;
    uint32_t invalid_indices = 0;
    uint32_t nan_coordinates = 0;
    uint32_t flipped_polygons = 0;
    
    bool isValid() const {
        return degenerate_polygons == 0 && 
               invalid_indices == 0 && 
               nan_coordinates == 0;
    }
};

/**
 * @brief Validate output with detailed statistics
 * 
 * @param soup Output polygon soup to validate
 * @param log Diagnostic log for reporting issues
 * @param[out] stats Validation statistics
 * @return true if output passes all validation checks
 */
bool validateOutputDetailed(const PolygonSoup& soup, 
                            DiagnosticLog& log, 
                            OutputValidationStats& stats);

//=============================================================================
// OUTPUT ATTRIBUTES
//=============================================================================

/**
 * @brief Generate piece name attribute for shatter output
 * 
 * Creates names like "piece_0", "piece_1", etc. for each connected
 * component identified during shatter post-processing.
 * 
 * @param soup Polygon soup with piece_ids assigned
 * @param name_prefix Prefix for piece names (default: "piece_")
 */
void generatePieceNames(PolygonSoup& soup, const std::string& name_prefix = "piece_");

/**
 * @brief Generate pivot point attribute for each piece
 * 
 * Computes the center of mass (average vertex position) for each
 * connected component. Used by Houdini for transformation operations.
 * 
 * @param soup Polygon soup with piece_ids assigned
 * @param ctx Quantization context for coordinate conversion
 */
void generatePiecePivots(PolygonSoup& soup, const QuantizationContext& ctx);

//=============================================================================
// EXPORT FORMATS
//=============================================================================

/**
 * @brief Export polygon soup to raw vertex array
 * 
 * Flattens output polygons into a vertex array suitable for Houdini
 * or other downstream consumers.
 * 
 * @param soup Output polygon soup
 * @param[out] out_vertices Output vertex array [x0,y0,z0,x1,y1,z1,...]
 * @param[out] out_indices Output index array [p0i0,p0i1,p0i2,p1i0,...]
 * @param[out] out_vertex_counts Vertices per polygon [n0,n1,n2,...]
 */
void exportToRaw(const PolygonSoup& soup,
                 std::vector<double>& out_vertices,
                 std::vector<uint32_t>& out_indices,
                 std::vector<uint32_t>& out_vertex_counts);

/**
 * @brief Export polygon soup to indexed mesh format
 * 
 * Creates a deduplicated vertex list with indices, reducing memory
 * usage for meshes with shared vertices.
 * 
 * @param soup Output polygon soup
 * @param tolerance Vertex deduplication tolerance
 * @param[out] out_positions Unique vertex positions
 * @param[out] out_indices Polygon indices
 */
void exportToIndexed(const PolygonSoup& soup,
                     double tolerance,
                     std::vector<double>& out_positions,
                     std::vector<uint32_t>& out_indices);

} // namespace ember
