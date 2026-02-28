/**
 * @file SeamExtract.cpp
 * @brief Implementation of seam extraction and segment chaining
 */

#include "SeamExtract.h"
#include <cmath>
#include <unordered_map>
#include <queue>
#include <algorithm>

namespace ember {

//=============================================================================
// SEAM SEGMENT EXTRACTION FROM TRIANGLE PAIRS
//=============================================================================

/**
 * @brief Compute intersection of two triangles (simplified)
 * 
 * This is a placeholder for the full triangle-triangle intersection.
 * In a complete implementation, this would use exact predicates and
 * handle all intersection cases (coplanar, edge-plane, vertex-plane).
 * 
 * @param tri0 First triangle
 * @param tri1 Second triangle
 * @param[out] segments Intersection segments found
 * @return Number of segments (0, 1, or 2)
 */
static int intersectTriangles(const Triangle& tri0, const Triangle& tri1,
                               std::array<std::array<double, 3>, 2>* segments) {
    // Simplified intersection test using plane equations
    // A full implementation would use exact predicates
    
    // Check if triangles are coplanar
    int64_t dot = static_cast<int64_t>(tri0.plane.a) * tri1.plane.a +
                  static_cast<int64_t>(tri0.plane.b) * tri1.plane.b +
                  static_cast<int64_t>(tri0.plane.c) * tri1.plane.c;
    int64_t len0 = static_cast<int64_t>(tri0.plane.a) * tri0.plane.a +
                   static_cast<int64_t>(tri0.plane.b) * tri0.plane.b +
                   static_cast<int64_t>(tri0.plane.c) * tri0.plane.c;
    int64_t len1 = static_cast<int64_t>(tri1.plane.a) * tri1.plane.a +
                   static_cast<int64_t>(tri1.plane.b) * tri1.plane.b +
                   static_cast<int64_t>(tri1.plane.c) * tri1.plane.c;
    
    // Check if normals are parallel (coplanar case)
    int64_t cross_len_sq = len0 * len1 - dot * dot;
    
    if (cross_len_sq < 1000) {  // Threshold for coplanarity
        // Coplanar case - would need 2D intersection test
        // For now, return 0 segments (not implemented in simplified version)
        return 0;
    }
    
    // Non-coplanar: check for edge-plane intersections
    // Each edge of tri0 that crosses tri1's plane creates a potential intersection point
    int crossings0 = 0;
    std::array<std::array<double, 3>, 2> pts0;
    
    for (int e = 0; e < 3 && crossings0 < 2; ++e) {
        const int32_t* p0 = tri0.iv[e];
        const int32_t* p1 = tri0.iv[(e + 1) % 3];
        
        // Evaluate plane equation at edge endpoints
        int64_t d0 = tri1.plane.evaluate(p0[0], p0[1], p0[2]);
        int64_t d1 = tri1.plane.evaluate(p1[0], p1[1], p1[2]);
        
        // Check for sign change (edge crosses plane)
        if ((d0 < 0 && d1 > 0) || (d0 > 0 && d1 < 0)) {
            // Compute intersection point
            double t = static_cast<double>(-d0) / static_cast<double>(d1 - d0);
            pts0[crossings0][0] = p0[0] + t * (p1[0] - p0[0]);
            pts0[crossings0][1] = p0[1] + t * (p1[1] - p0[1]);
            pts0[crossings0][2] = p0[2] + t * (p1[2] - p0[2]);
            ++crossings0;
        }
    }
    
    if (crossings0 < 2) {
        return 0;  // No intersection or single point (vertex)
    }
    
    // Store the segment
    segments[0] = pts0;
    return 1;
}

/**
 * @brief Check if a point is inside a triangle (2D barycentric test)
 * 
 * Used for coplanar intersection testing.
 */
static bool pointInTriangle2D(const double p[2], 
                               const double a[2], 
                               const double b[2], 
                               const double c[2]) {
    double denom = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1]);
    if (std::abs(denom) < 1e-10) return false;
    
    double w1 = ((b[1] - c[1]) * (p[0] - c[0]) + (c[0] - b[0]) * (p[1] - c[1])) / denom;
    double w2 = ((c[1] - a[1]) * (p[0] - c[0]) + (a[0] - c[0]) * (p[1] - c[1])) / denom;
    double w3 = 1.0 - w1 - w2;
    
    return w1 >= -1e-10 && w2 >= -1e-10 && w3 >= -1e-10;
}

uint32_t extractSeamsFromPair(const PolygonSoup& soup, 
                               const CandidatePair& pair,
                               std::vector<SeamSegmentEx>& segments) {
    if (pair.tri0 >= soup.triangles.size() || pair.tri1 >= soup.triangles.size()) {
        return 0;
    }
    
    const Triangle& tri0 = soup.triangles[pair.tri0];
    const Triangle& tri1 = soup.triangles[pair.tri1];
    
    std::array<std::array<double, 3>, 2> intersection_pts;
    int num_segments = intersectTriangles(tri0, tri1, &intersection_pts);
    
    if (num_segments > 0) {
        // Create seam segment from intersection
        // In a full implementation, we would add vertices to the soup
        // and create proper vertex indices
        SeamSegmentEx seg;
        seg.tri0 = pair.tri0;
        seg.tri1 = pair.tri1;
        seg.seam_type = (tri0.mesh_id != tri1.mesh_id) ? 0 : 
                        (tri0.mesh_id == 0) ? 1 : 2;
        
        // Compute segment length
        double dx = intersection_pts[1][0] - intersection_pts[0][0];
        double dy = intersection_pts[1][1] - intersection_pts[0][1];
        double dz = intersection_pts[1][2] - intersection_pts[0][2];
        seg.length = static_cast<float>(std::sqrt(dx*dx + dy*dy + dz*dz));
        
        // Placeholder: would add vertices and set v0, v1 properly
        seg.v0 = static_cast<uint32_t>(segments.size() * 2);
        seg.v1 = static_cast<uint32_t>(segments.size() * 2 + 1);
        
        segments.push_back(seg);
    }
    
    return num_segments > 0 ? 1 : 0;
}

//=============================================================================
// MAIN SEAM EXTRACTION
//=============================================================================

void extractSeams(PolygonSoup& soup, const SeamConfig& config) {
    DiagnosticLog log;
    extractSeams(soup, config, log);
}

void extractSeams(PolygonSoup& soup, const SeamConfig& config, DiagnosticLog& log) {
    EMBER_SCOPED_TIMER("Extract seams");
    
    std::vector<SeamSegmentEx> segments;
    
    // Process all candidate pairs
    uint32_t processed = 0;
    for (const auto& pair : soup.candidate_pairs) {
        // Check if we should process this pair based on seam type
        bool is_ab = true;  // Simplified - would check mesh IDs
        bool is_aa = false;
        bool is_bb = false;
        
        if ((is_ab && config.extract_ab_seams) ||
            (is_aa && config.extract_aa_seams) ||
            (is_bb && config.extract_bb_seams)) {
            processed += extractSeamsFromPair(soup, pair, segments);
        }
    }
    
    log.info("Extracted " + std::to_string(segments.size()) + 
             " seam segments from " + std::to_string(soup.candidate_pairs.size()) + 
             " candidate pairs");
    
    // Merge duplicate vertices
    if (config.seam_tolerance > 0) {
        mergeSeamVertices(soup, config.seam_tolerance);
    }
    
    // Chain segments into polylines
    if (config.chain_segments) {
        std::vector<SeamPolyline> polylines;
        chainSeamSegments(segments, polylines);
        
        log.info("Chained into " + std::to_string(polylines.size()) + " polylines");
    }
    
    // Copy to soup's seam segments (simplified conversion)
    soup.seam_segments.clear();
    for (const auto& seg : segments) {
        SeamSegment s;
        s.v0 = seg.v0;
        s.v1 = seg.v1;
        s.tri0 = seg.tri0;
        s.tri1 = seg.tri1;
        s.is_boundary = seg.is_boundary;
        soup.seam_segments.push_back(s);
    }
}

//=============================================================================
// SEGMENT CHAINING
//=============================================================================

/**
 * @brief Build vertex-to-segments adjacency
 */
static void buildVertexSegmentAdjacency(
    const std::vector<SeamSegmentEx>& segments,
    std::unordered_map<uint32_t, std::vector<uint32_t>>& adjacency) {
    
    for (uint32_t i = 0; i < segments.size(); ++i) {
        adjacency[segments[i].v0].push_back(i);
        adjacency[segments[i].v1].push_back(i);
    }
}

void chainSeamSegments(std::vector<SeamSegmentEx>& segments,
                       std::vector<SeamPolyline>& polylines) {
    chainSeamSegmentsEx(segments, 1e-6, true, polylines);
}

void chainSeamSegmentsEx(const std::vector<SeamSegmentEx>& segments,
                         double tolerance,
                         bool allow_branching,
                         std::vector<SeamPolyline>& polylines) {
    EMBER_SCOPED_TIMER("Chain seam segments");
    
    polylines.clear();
    if (segments.empty()) {
        return;
    }
    
    // Build vertex-to-segments adjacency
    std::unordered_map<uint32_t, std::vector<uint32_t>> vertex_segments;
    buildVertexSegmentAdjacency(segments, vertex_segments);
    
    // Track used segments
    std::vector<bool> used(segments.size(), false);
    
    // Find endpoints (vertices with odd number of segments)
    std::vector<uint32_t> endpoints;
    for (const auto& pair : vertex_segments) {
        if (pair.second.size() % 2 == 1) {
            endpoints.push_back(pair.first);
        }
    }
    
    // Helper to extend a polyline from a starting vertex
    auto extendPolyline = [&](uint32_t start_vertex, bool follow_used) -> SeamPolyline {
        SeamPolyline polyline;
        polyline.vertices.push_back(start_vertex);
        
        uint32_t current = start_vertex;
        
        while (true) {
            // Find an unused segment connected to current vertex
            uint32_t next_segment = UINT32_MAX;
            uint32_t next_vertex = UINT32_MAX;
            
            auto it = vertex_segments.find(current);
            if (it == vertex_segments.end()) break;
            
            for (uint32_t seg_idx : it->second) {
                if (used[seg_idx] && !follow_used) continue;
                if (!used[seg_idx] && follow_used) continue;
                
                const auto& seg = segments[seg_idx];
                next_vertex = (seg.v0 == current) ? seg.v1 : seg.v0;
                next_segment = seg_idx;
                break;
            }
            
            if (next_segment == UINT32_MAX) break;
            
            // Mark segment as used and extend polyline
            used[next_segment] = true;
            
            // Check for closed loop
            if (next_vertex == start_vertex) {
                polyline.is_closed = true;
                break;
            }
            
            polyline.vertices.push_back(next_vertex);
            current = next_vertex;
        }
        
        return polyline;
    };
    
    // Chain from each endpoint (creates open polylines)
    for (uint32_t endpoint : endpoints) {
        // Check if this endpoint is already fully used
        bool has_unused = false;
        auto it = vertex_segments.find(endpoint);
        if (it != vertex_segments.end()) {
            for (uint32_t seg_idx : it->second) {
                if (!used[seg_idx]) {
                    has_unused = true;
                    break;
                }
            }
        }
        
        if (!has_unused) continue;
        
        SeamPolyline polyline = extendPolyline(endpoint, false);
        
        // Reverse and extend the other direction if possible
        if (polyline.vertices.size() > 1) {
            std::reverse(polyline.vertices.begin(), polyline.vertices.end());
            
            // Extend from the new start (original end)
            uint32_t new_start = polyline.vertices.back();
            SeamPolyline extension = extendPolyline(new_start, false);
            
            // Append extension (skip first vertex as it's already in polyline)
            for (size_t i = 1; i < extension.vertices.size(); ++i) {
                polyline.vertices.push_back(extension.vertices[i]);
            }
            polyline.is_closed = extension.is_closed;
        }
        
        if (polyline.vertices.size() >= 2) {
            polylines.push_back(std::move(polyline));
        }
    }
    
    // Find closed loops (all vertices have even degree)
    for (uint32_t i = 0; i < segments.size(); ++i) {
        if (used[i]) continue;
        
        // Start a new loop from this segment
        SeamPolyline loop;
        loop.is_closed = false;
        
        uint32_t start = segments[i].v0;
        uint32_t current = start;
        uint32_t prev = segments[i].v1;
        
        loop.vertices.push_back(start);
        used[i] = true;
        
        while (true) {
            // Find next segment in loop
            auto it = vertex_segments.find(current);
            if (it == vertex_segments.end()) break;
            
            uint32_t next_segment = UINT32_MAX;
            uint32_t next_vertex = UINT32_MAX;
            
            for (uint32_t seg_idx : it->second) {
                if (used[seg_idx]) continue;
                
                const auto& seg = segments[seg_idx];
                uint32_t other = (seg.v0 == current) ? seg.v1 : seg.v0;
                
                // Don't go back immediately
                if (other == prev && it->second.size() > 1) continue;
                
                next_vertex = other;
                next_segment = seg_idx;
                break;
            }
            
            if (next_segment == UINT32_MAX) break;
            
            used[next_segment] = true;
            prev = current;
            current = next_vertex;
            
            if (current == start) {
                loop.is_closed = true;
                break;
            }
            
            loop.vertices.push_back(current);
        }
        
        if (loop.vertices.size() >= 2) {
            polylines.push_back(std::move(loop));
        }
    }
}

//=============================================================================
// SEAM UTILITIES
//=============================================================================

uint32_t mergeSeamVertices(PolygonSoup& soup, double tolerance) {
    // Placeholder: would implement vertex merging using spatial hashing
    // For now, just return 0
    (void)soup;
    (void)tolerance;
    return 0;
}

double computeSeamLength(const SeamSegment& seg, const PolygonSoup& soup) {
    // Get vertex positions from triangles
    if (seg.tri0 >= soup.triangles.size()) {
        return 0.0;
    }
    
    const Triangle& tri = soup.triangles[seg.tri0];
    
    // Simplified: use first two vertices of triangle as proxy
    double dx = tri.v[1][0] - tri.v[0][0];
    double dy = tri.v[1][1] - tri.v[0][1];
    double dz = tri.v[1][2] - tri.v[0][2];
    
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

void filterSeamsByLength(const std::vector<SeamSegmentEx>& segments,
                         double min_length,
                         std::vector<SeamSegmentEx>& filtered) {
    filtered.clear();
    filtered.reserve(segments.size());
    
    for (const auto& seg : segments) {
        if (seg.length >= min_length) {
            filtered.push_back(seg);
        }
    }
}

//=============================================================================
// SEAM STATISTICS
//=============================================================================

void computeSeamStats(const std::vector<SeamSegmentEx>& segments,
                      const std::vector<SeamPolyline>& polylines,
                      SeamStats& stats) {
    stats.total_segments = static_cast<uint32_t>(segments.size());
    stats.total_polylines = static_cast<uint32_t>(polylines.size());
    
    // Count by type
    for (const auto& seg : segments) {
        switch (seg.seam_type) {
            case 0: ++stats.ab_segments; break;
            case 1: ++stats.aa_segments; break;
            case 2: ++stats.bb_segments; break;
        }
        if (seg.is_boundary) {
            ++stats.boundary_segments;
        }
        stats.total_length += seg.length;
    }
    
    // Count polylines
    for (const auto& poly : polylines) {
        if (poly.is_closed) {
            ++stats.closed_loops;
        } else {
            ++stats.open_chains;
        }
    }
}

} // namespace ember
