/**
 * @file MeshExport.cpp
 * @brief Implementation of mesh export and output coordinate materialization
 */

#include "MeshExport.h"
#include <cmath>
#include <unordered_map>
#include <algorithm>

namespace ember {

//=============================================================================
// IMPLICIT POINT MATERIALIZATION
//=============================================================================

/**
 * @brief Compute LPI (Line-Plane Intersection) point
 * 
 * Given edge (p0, p1) and plane (nx, ny, nz, d), computes:
 *   t = -(n·p0 + d) / (n·(p1 - p0))
 *   result = p0 + t * (p1 - p0)
 * 
 * Uses exact integer arithmetic for numerator and denominator,
 * then converts to double for the final division.
 */
static bool computeLPI(const int32_t p0[3], const int32_t p1[3],
                       const Plane& plane,
                       double& out_x, double& out_y, double& out_z) {
    // Compute n·p0 + d (numerator term)
    int64_t np0 = static_cast<int64_t>(plane.a) * p0[0] +
                  static_cast<int64_t>(plane.b) * p0[1] +
                  static_cast<int64_t>(plane.c) * p0[2] +
                  static_cast<int64_t>(plane.d);
    
    // Compute n·(p1 - p0) (denominator)
    int64_t dp[3] = {
        static_cast<int64_t>(p1[0]) - p0[0],
        static_cast<int64_t>(p1[1]) - p0[1],
        static_cast<int64_t>(p1[2]) - p0[2]
    };
    int64_t denom = static_cast<int64_t>(plane.a) * dp[0] +
                    static_cast<int64_t>(plane.b) * dp[1] +
                    static_cast<int64_t>(plane.c) * dp[2];
    
    // Check for parallel edge and plane
    if (denom == 0) {
        return false;  // Edge is parallel to plane
    }
    
    // Compute t = -np0 / denom
    double t = -static_cast<double>(np0) / static_cast<double>(denom);
    
    // Intersection point: p0 + t * (p1 - p0)
    out_x = static_cast<double>(p0[0]) + t * static_cast<double>(dp[0]);
    out_y = static_cast<double>(p0[1]) + t * static_cast<double>(dp[1]);
    out_z = static_cast<double>(p0[2]) + t * static_cast<double>(dp[2]);
    
    return true;
}

/**
 * @brief Compute TPI (Triangle-Plane Intersection) point
 * 
 * For coplanar triangle intersection, computes the intersection
 * of an edge from tri0 with the plane of tri1 (which is coplanar).
 * This is essentially the same as LPI for the coplanar case.
 */
static bool computeTPI(const Triangle& tri0, const Triangle& tri1,
                       uint8_t edge_idx,
                       double& out_x, double& out_y, double& out_z) {
    // Get edge vertices from tri0
    uint8_t v0_idx = edge_idx % 3;
    uint8_t v1_idx = (edge_idx + 1) % 3;
    
    const int32_t* p0 = tri0.iv[v0_idx];
    const int32_t* p1 = tri0.iv[v1_idx];
    
    // Use tri1's plane (coplanar case)
    return computeLPI(p0, p1, tri1.plane, out_x, out_y, out_z);
}

bool materializeImplicitPoint(const ImplicitPoint& point,
                               const PolygonSoup& soup,
                               const QuantizationContext& ctx,
                               double& out_x, double& out_y, double& out_z) {
    switch (point.type) {
        case ImplicitPointType::EXPLICIT: {
            // Direct lookup of original vertex
            uint32_t mesh_id = point.explicit_point.mesh_id;
            uint32_t tri_idx = point.explicit_point.tri_idx;
            uint8_t vert_idx = point.explicit_point.vert_idx;
            
            // Get triangle range for this mesh
            auto range = soup.getMeshTriangleRange(mesh_id);
            size_t global_tri_idx = range.first + tri_idx;
            
            if (global_tri_idx >= soup.triangles.size()) {
                return false;
            }
            
            const Triangle& tri = soup.triangles[global_tri_idx];
            if (vert_idx > 2) {
                return false;
            }
            
            // Dequantize vertex
            out_x = ctx.dequantize(tri.iv[vert_idx][0], 0);
            out_y = ctx.dequantize(tri.iv[vert_idx][1], 1);
            out_z = ctx.dequantize(tri.iv[vert_idx][2], 2);
            return true;
        }
        
        case ImplicitPointType::LPI: {
            // Line-Plane Intersection
            uint32_t tri0_idx = point.lpi.mesh0_tri;
            uint8_t edge_v0 = point.lpi.edge_v0;
            uint8_t edge_v1 = point.lpi.edge_v1;
            uint32_t tri1_idx = point.lpi.mesh1_tri;
            
            if (tri0_idx >= soup.triangles.size() || 
                tri1_idx >= soup.triangles.size()) {
                return false;
            }
            
            const Triangle& tri0 = soup.triangles[tri0_idx];
            const Triangle& tri1 = soup.triangles[tri1_idx];
            
            // Get edge vertices
            if (edge_v0 > 2 || edge_v1 > 2) {
                return false;
            }
            
            const int32_t* p0 = tri0.iv[edge_v0];
            const int32_t* p1 = tri0.iv[edge_v1];
            
            // Compute intersection with tri1's plane
            if (!computeLPI(p0, p1, tri1.plane, out_x, out_y, out_z)) {
                return false;
            }
            
            // Dequantize result
            out_x = ctx.dequantize(static_cast<int32_t>(out_x), 0);
            out_y = ctx.dequantize(static_cast<int32_t>(out_y), 1);
            out_z = ctx.dequantize(static_cast<int32_t>(out_z), 2);
            return true;
        }
        
        case ImplicitPointType::TPI: {
            // Triangle-Plane Intersection (coplanar case)
            uint32_t tri0_idx = point.tpi.mesh0_tri;
            uint32_t tri1_idx = point.tpi.mesh1_tri;
            uint8_t edge_idx = point.tpi.edge_idx;
            
            if (tri0_idx >= soup.triangles.size() || 
                tri1_idx >= soup.triangles.size()) {
                return false;
            }
            
            const Triangle& tri0 = soup.triangles[tri0_idx];
            const Triangle& tri1 = soup.triangles[tri1_idx];
            
            if (!computeTPI(tri0, tri1, edge_idx, out_x, out_y, out_z)) {
                return false;
            }
            
            // Dequantize result
            out_x = ctx.dequantize(static_cast<int32_t>(out_x), 0);
            out_y = ctx.dequantize(static_cast<int32_t>(out_y), 1);
            out_z = ctx.dequantize(static_cast<int32_t>(out_z), 2);
            return true;
        }
        
        default:
            return false;
    }
}

//=============================================================================
// OUTPUT COORDINATE MATERIALIZATION
//=============================================================================

void materializeOutputCoordinates(PolygonSoup& soup, const QuantizationContext& ctx) {
    EMBER_SCOPED_TIMER("Materialize output coordinates");
    
    // For now, dequantize the output polygons
    // In a full implementation, this would also evaluate implicit points
    dequantizeOutputPolygons(soup, ctx);
}

void dequantizeOutputPolygons(PolygonSoup& soup, const QuantizationContext& ctx) {
    // Update soup's quantization parameters from context
    for (int i = 0; i < 3; ++i) {
        soup.quantization_scale[i] = ctx.scale[i];
        soup.quantization_inv_scale[i] = ctx.inv_scale[i];
    }
    soup.per_axis_quantization = ctx.per_axis;
    
    // Note: Output polygons use vertex_indices into implicit point array
    // The actual dequantization happens when those points are materialized
    // This function updates the soup's quantization metadata
}

//=============================================================================
// OUTPUT VALIDATION
//=============================================================================

bool validateOutput(const PolygonSoup& soup, DiagnosticLog& log) {
    OutputValidationStats stats;
    return validateOutputDetailed(soup, log, stats);
}

bool validateOutputDetailed(const PolygonSoup& soup, 
                            DiagnosticLog& log, 
                            OutputValidationStats& stats) {
    EMBER_SCOPED_TIMER("Validate output");
    
    stats.total_polygons = static_cast<uint32_t>(soup.output_polygons.size());
    
    if (soup.output_polygons.empty()) {
        log.warn(DiagCategory::Warning, "Output contains no polygons");
        return true;  // Empty output is technically valid
    }
    
    // Track vertex indices to check for validity
    size_t max_vertex_idx = 0;
    for (const auto& poly : soup.output_polygons) {
        for (uint32_t idx : poly.vertex_indices) {
            max_vertex_idx = std::max(max_vertex_idx, static_cast<size_t>(idx));
        }
    }
    
    for (size_t i = 0; i < soup.output_polygons.size(); ++i) {
        const auto& poly = soup.output_polygons[i];
        
        // Check for degenerate polygons (less than 3 vertices)
        if (poly.vertex_indices.size() < 3) {
            ++stats.degenerate_polygons;
            log.warn(DiagCategory::Degenerate, 
                     "Polygon " + std::to_string(i) + " has fewer than 3 vertices");
            continue;
        }
        
        // Check for invalid vertex indices
        for (uint32_t idx : poly.vertex_indices) {
            if (idx > max_vertex_idx) {
                ++stats.invalid_indices;
                log.error(DiagCategory::Error,
                          "Polygon " + std::to_string(i) + " has invalid vertex index");
                break;
            }
        }
        
        // Check for flipped winding
        if (poly.flipped) {
            ++stats.flipped_polygons;
        }
    }
    
    // Validate triangle vertices for NaN/Inf
    for (size_t i = 0; i < soup.triangles.size(); ++i) {
        const auto& tri = soup.triangles[i];
        for (int v = 0; v < 3; ++v) {
            for (int c = 0; c < 3; ++c) {
                double val = tri.v[v][c];
                if (std::isnan(val) || std::isinf(val)) {
                    ++stats.nan_coordinates;
                    log.error(DiagCategory::Error,
                              "Triangle " + std::to_string(i) + " has invalid coordinate");
                }
            }
        }
    }
    
    // Log summary
    log.info("Output validation: " + std::to_string(stats.total_polygons) + 
             " polygons, " + std::to_string(stats.degenerate_polygons) + 
             " degenerate, " + std::to_string(stats.invalid_indices) + 
             " invalid indices, " + std::to_string(stats.nan_coordinates) + 
             " NaN/Inf");
    
    return stats.isValid();
}

//=============================================================================
// OUTPUT ATTRIBUTES
//=============================================================================

void generatePieceNames(PolygonSoup& soup, const std::string& name_prefix) {
    if (soup.piece_count <= 0) {
        return;
    }
    
    soup.piece_names.resize(soup.piece_count);
    for (int i = 0; i < soup.piece_count; ++i) {
        soup.piece_names[i] = name_prefix + std::to_string(i);
    }
}

void generatePiecePivots(PolygonSoup& soup, const QuantizationContext& ctx) {
    if (soup.piece_count <= 0 || soup.piece_ids.empty()) {
        return;
    }
    
    EMBER_SCOPED_TIMER("Generate piece pivots");
    
    // Accumulate vertex positions per piece
    std::vector<std::array<double, 3>> accum(soup.piece_count, {0.0, 0.0, 0.0});
    std::vector<uint32_t> counts(soup.piece_count, 0);
    
    // For each output polygon, accumulate its vertices to its piece
    for (size_t i = 0; i < soup.output_polygons.size(); ++i) {
        const auto& poly = soup.output_polygons[i];
        if (i >= soup.piece_ids.size()) {
            continue;
        }
        
        int piece_id = soup.piece_ids[i];
        if (piece_id < 0 || piece_id >= soup.piece_count) {
            continue;
        }
        
        // Get triangle for vertex positions
        if (poly.original_tri >= soup.triangles.size()) {
            continue;
        }
        
        const Triangle& tri = soup.triangles[poly.original_tri];
        
        // Accumulate triangle centroid
        for (int v = 0; v < 3; ++v) {
            accum[piece_id][0] += tri.v[v][0];
            accum[piece_id][1] += tri.v[v][1];
            accum[piece_id][2] += tri.v[v][2];
        }
        counts[piece_id] += 3;
    }
    
    // Compute averages
    soup.piece_pivots.resize(soup.piece_count);
    for (int i = 0; i < soup.piece_count; ++i) {
        if (counts[i] > 0) {
            soup.piece_pivots[i][0] = static_cast<float>(accum[i][0] / counts[i]);
            soup.piece_pivots[i][1] = static_cast<float>(accum[i][1] / counts[i]);
            soup.piece_pivots[i][2] = static_cast<float>(accum[i][2] / counts[i]);
        } else {
            soup.piece_pivots[i] = {0.0f, 0.0f, 0.0f};
        }
    }
}

//=============================================================================
// EXPORT FORMATS
//=============================================================================

void exportToRaw(const PolygonSoup& soup,
                 std::vector<double>& out_vertices,
                 std::vector<uint32_t>& out_indices,
                 std::vector<uint32_t>& out_vertex_counts) {
    EMBER_SCOPED_TIMER("Export to raw format");
    
    out_vertices.clear();
    out_indices.clear();
    out_vertex_counts.clear();
    
    // Reserve space
    size_t total_verts = 0;
    for (const auto& poly : soup.output_polygons) {
        total_verts += poly.vertex_indices.size();
    }
    
    out_vertices.reserve(total_verts * 3);
    out_indices.reserve(total_verts);
    out_vertex_counts.reserve(soup.output_polygons.size());
    
    // Export each polygon
    uint32_t current_idx = 0;
    for (const auto& poly : soup.output_polygons) {
        out_vertex_counts.push_back(static_cast<uint32_t>(poly.vertex_indices.size()));
        
        // Get source triangle for vertex positions
        if (poly.original_tri >= soup.triangles.size()) {
            continue;
        }
        
        const Triangle& tri = soup.triangles[poly.original_tri];
        
        // For triangles, use the original vertices
        // For n-gons, this is a simplified approach
        for (size_t v = 0; v < poly.vertex_indices.size() && v < 3; ++v) {
            out_vertices.push_back(tri.v[v][0]);
            out_vertices.push_back(tri.v[v][1]);
            out_vertices.push_back(tri.v[v][2]);
            out_indices.push_back(current_idx++);
        }
    }
}

void exportToIndexed(const PolygonSoup& soup,
                     double tolerance,
                     std::vector<double>& out_positions,
                     std::vector<uint32_t>& out_indices) {
    EMBER_SCOPED_TIMER("Export to indexed format");
    
    out_positions.clear();
    out_indices.clear();
    
    // Simple vertex deduplication using spatial hashing
    struct VertexKey {
        int32_t x, y, z;
        
        bool operator==(const VertexKey& other) const {
            return x == other.x && y == other.y && z == other.z;
        }
    };
    
    struct VertexKeyHash {
        size_t operator()(const VertexKey& k) const {
            // FNV-1a inspired hash
            uint64_t h = 14695981039346656037ULL;
            h ^= static_cast<uint64_t>(k.x);
            h *= 1099511628211ULL;
            h ^= static_cast<uint64_t>(k.y);
            h *= 1099511628211ULL;
            h ^= static_cast<uint64_t>(k.z);
            h *= 1099511628211ULL;
            return static_cast<size_t>(h);
        }
    };
    
    std::unordered_map<VertexKey, uint32_t, VertexKeyHash> vertex_map;
    
    // Quantization factor for deduplication
    double quant_factor = 1.0 / tolerance;
    
    for (const auto& poly : soup.output_polygons) {
        if (poly.original_tri >= soup.triangles.size()) {
            continue;
        }
        
        const Triangle& tri = soup.triangles[poly.original_tri];
        
        for (int v = 0; v < 3; ++v) {
            // Quantize position for deduplication
            VertexKey key;
            key.x = static_cast<int32_t>(tri.v[v][0] * quant_factor);
            key.y = static_cast<int32_t>(tri.v[v][1] * quant_factor);
            key.z = static_cast<int32_t>(tri.v[v][2] * quant_factor);
            
            auto it = vertex_map.find(key);
            if (it == vertex_map.end()) {
                // New vertex
                uint32_t new_idx = static_cast<uint32_t>(out_positions.size() / 3);
                vertex_map[key] = new_idx;
                out_positions.push_back(tri.v[v][0]);
                out_positions.push_back(tri.v[v][1]);
                out_positions.push_back(tri.v[v][2]);
                out_indices.push_back(new_idx);
            } else {
                // Existing vertex
                out_indices.push_back(it->second);
            }
        }
    }
}

} // namespace ember
