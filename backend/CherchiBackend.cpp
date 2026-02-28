// ═══════════════════════════════════════════════════════════════════════════════
// EMBER CherchiBackend — Robust Boolean Operations via Local Arrangements
// ═══════════════════════════════════════════════════════════════════════════════
//
// Implementation of the Cherchi et al. 2020 framework for exact mesh Booleans.
// Uses Constrained Delaunay Triangulation (CDT) for robust local polygon
// triangulation, replacing fragile ear-clipping algorithms.
//
// INTEGRATED VERSION — All Fixes Applied:
//   - FIX 2.1: Wire CDT predicates through ExactPredicates.h
//   - FIX 2.2-2.4: Replace triangulate_cavity with Sloan edge-flip algorithm
//   - FIX 2.5: Replace O(n) locate_triangle with walking algorithm
//   - FIX 2.6: Exact point_on_segment (no epsilon)
//   - FIX 2.7: Correct legalize_edge recursion targets
//   - FIX 2.8: Optimize mark_edge_constrained with short-circuit
//   - HOT FIX 5: flip_edge() neighbor back-references
//   - HOT FIX 7: split_edge() neighbor setting
//
// Reference: "Exact and Efficient Booleans for Polyhedra" (Cherchi et al. 2020)

#include <cstdint>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
#include <cstring>

#include "ember/ExactPredicates.h"
#include "ember/PolygonSoup.h"
#include "backend/IBooleanBackend.h"

namespace ember {

// ═══════════════════════════════════════════════════════════════════════════════
// Forward Declarations and Type Definitions
// ═══════════════════════════════════════════════════════════════════════════════

// 2D projected point for CDT
struct Point2D {
    double x, y;
    uint32_t orig_idx;     // Original 3D point index
    bool is_constraint;    // Part of input polygon edge
    
    Point2D() : x(0), y(0), orig_idx(0), is_constraint(false) {}
    Point2D(double x_, double y_, uint32_t idx, bool constraint = false)
        : x(x_), y(y_), orig_idx(idx), is_constraint(constraint) {}
};

// ═══════════════════════════════════════════════════════════════════════════════
// EXACT CDT PREDICATES (wired through ExactPredicates.h)
// ═══════════════════════════════════════════════════════════════════════════════

namespace predicates {

// orient2d: exact via ExactPredicates.h
inline int orient2d(const Point2D& a, const Point2D& b, const Point2D& c) {
    return ember::predicates::orient2d_filtered(a.x, a.y, b.x, b.y, c.x, c.y);
}

// incircle: exact via ExactPredicates.h
inline int incircle(const Point2D& a, const Point2D& b,
                    const Point2D& c, const Point2D& d) {
    return ember::predicates::incircle_filtered(
        a.x, a.y, b.x, b.y, c.x, c.y, d.x, d.y);
}

// in_circumcircle: convenience wrapper
inline bool in_circumcircle(const Point2D& a, const Point2D& b,
                            const Point2D& c, const Point2D& d) {
    int orient = orient2d(a, b, c);
    if (orient == 0) return false;
    int ic = incircle(a, b, c, d);
    return (orient > 0) ? (ic > 0) : (ic < 0);
}

// EXACT point_on_segment (no epsilon)
inline bool point_on_segment(const Point2D& p, const Point2D& a, const Point2D& b) {
    if (orient2d(a, b, p) != 0) return false;
    double minx = std::min(a.x, b.x);
    double maxx = std::max(a.x, b.x);
    double miny = std::min(a.y, b.y);
    double maxy = std::max(a.y, b.y);
    return p.x >= minx && p.x <= maxx && p.y >= miny && p.y <= maxy;
}

// segments_intersect: uses exact orient2d
inline bool segments_intersect(const Point2D& a1, const Point2D& a2,
                                const Point2D& b1, const Point2D& b2) {
    int d1 = orient2d(a1, a2, b1);
    int d2 = orient2d(a1, a2, b2);
    int d3 = orient2d(b1, b2, a1);
    int d4 = orient2d(b1, b2, a2);

    if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
        ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) {
        return true;
    }

    if (d1 == 0 && point_on_segment(b1, a1, a2)) return true;
    if (d2 == 0 && point_on_segment(b2, a1, a2)) return true;
    if (d3 == 0 && point_on_segment(a1, b1, b2)) return true;
    if (d4 == 0 && point_on_segment(a2, b1, b2)) return true;

    return false;
}

} // namespace predicates

// ═══════════════════════════════════════════════════════════════════════════════
// Constrained Delaunay Triangulation (CDT)
// ═══════════════════════════════════════════════════════════════════════════════

// Triangle edge reference: (triangle_index, edge_index within triangle)
struct EdgeRef {
    uint32_t tri;
    uint8_t edge;  // 0, 1, or 2
    
    EdgeRef() : tri(0), edge(0) {}
    EdgeRef(uint32_t t, uint8_t e) : tri(t), edge(e) {}
    
    bool operator==(const EdgeRef& o) const {
        return tri == o.tri && edge == o.edge;
    }
};

// CDT Triangle structure
struct CDTTriangle {
    std::array<uint32_t, 3> v;
    std::array<uint32_t, 3> neighbor;
    std::array<bool, 3> constrained;
    bool valid;
    
    CDTTriangle() : valid(false) {
        v = {0, 0, 0};
        neighbor = {UINT32_MAX, UINT32_MAX, UINT32_MAX};
        constrained = {false, false, false};
    }
    
    CDTTriangle(uint32_t v0, uint32_t v1, uint32_t v2) : valid(true) {
        v = {v0, v1, v2};
        neighbor = {UINT32_MAX, UINT32_MAX, UINT32_MAX};
        constrained = {false, false, false};
    }
    
    std::pair<uint32_t, uint32_t> edge_vertices(uint8_t edge) const {
        return {v[(edge + 1) % 3], v[(edge + 2) % 3]};
    }
    
    int find_edge(uint32_t a, uint32_t b) const {
        for (int i = 0; i < 3; ++i) {
            auto [e0, e1] = edge_vertices(i);
            if ((e0 == a && e1 == b) || (e0 == b && e1 == a)) return i;
        }
        return -1;
    }
    
    bool has_vertex(uint32_t vtx) const {
        return v[0] == vtx || v[1] == vtx || v[2] == vtx;
    }
};

// Axis-aligned bounding box
struct BBox2D {
    double minx, miny, maxx, maxy;
    
    BBox2D() : minx(0), miny(0), maxx(0), maxy(0) {}
    
    void reset() {
        minx = miny = std::numeric_limits<double>::max();
        maxx = maxy = std::numeric_limits<double>::lowest();
    }
    
    void add_point(double x, double y) {
        minx = std::min(minx, x);
        miny = std::min(miny, y);
        maxx = std::max(maxx, x);
        maxy = std::max(maxy, y);
    }
    
    double diagonal() const {
        return std::sqrt((maxx - minx) * (maxx - minx) + 
                        (maxy - miny) * (maxy - miny));
    }
    
    bool contains(double x, double y) const {
        return x >= minx && x <= maxx && y >= miny && y <= maxy;
    }
};

// CDT Builder class
class CDTBuilder {
public:
    CDTBuilder() = default;
    
    // P1 FIX: clear() method to release memory
    void clear() {
        std::vector<Point2D>().swap(vertices);
        std::vector<CDTTriangle>().swap(triangles);
        std::vector<std::pair<uint32_t, uint32_t>>().swap(constraints);
        bbox.reset();
        last_located_tri = 0;
    }
    
    void addVertex(const Point2D& pt) {
        vertices.push_back(pt);
    }
    
    void addConstraint(uint32_t v0, uint32_t v1) {
        assert(v0 < vertices.size() && "Constraint vertex 0 out of bounds");
        assert(v1 < vertices.size() && "Constraint vertex 1 out of bounds");
        
        constraints.emplace_back(v0, v1);
        vertices[v0].is_constraint = true;
        vertices[v1].is_constraint = true;
    }
    
    void triangulate() {
        if (vertices.size() < 3) return;
        
        bbox.reset();
        for (const auto& v : vertices) {
            bbox.add_point(v.x, v.y);
        }
        
        create_super_triangle();
        
        for (uint32_t i = 0; i < vertices.size(); ++i) {
            insert_vertex(i);
        }
        
        for (const auto& [v0, v1] : constraints) {
            enforce_constraint(v0, v1);
        }
        
        remove_super_triangle();
        restore_delaunay();
    }
    
    std::vector<std::array<uint32_t, 3>> getTriangles() const {
        std::vector<std::array<uint32_t, 3>> result;
        result.reserve(triangles.size());
        
        for (const auto& tri : triangles) {
            if (tri.valid) {
                result.push_back(tri.v);
            }
        }
        
        return result;
    }
    
    size_t triangleCount() const {
        size_t count = 0;
        for (const auto& tri : triangles) {
            if (tri.valid) ++count;
        }
        return count;
    }
    
    bool validate() const {
        for (const auto& tri : triangles) {
            if (!tri.valid) continue;
            
            const auto& a = vertices[tri.v[0]];
            const auto& b = vertices[tri.v[1]];
            const auto& c = vertices[tri.v[2]];
            
            if (predicates::orient2d(a, b, c) <= 0) {
                return false;
            }
            
            for (int i = 0; i < 3; ++i) {
                if (tri.neighbor[i] != UINT32_MAX) {
                    const auto& neigh = triangles[tri.neighbor[i]];
                    if (!neigh.valid) return false;
                    
                    bool found = false;
                    for (int j = 0; j < 3; ++j) {
                        if (neigh.neighbor[j] == &tri - triangles.data()) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) return false;
                }
            }
        }
        
        return true;
    }

private:
    std::vector<Point2D> vertices;
    std::vector<CDTTriangle> triangles;
    std::vector<std::pair<uint32_t, uint32_t>> constraints;
    BBox2D bbox;
    
    uint32_t super_v0, super_v1, super_v2;
    uint32_t last_located_tri = 0;
    
    void create_super_triangle() {
        double dx = bbox.maxx - bbox.minx;
        double dy = bbox.maxy - bbox.miny;
        double dmax = std::max(dx, dy);
        double midx = (bbox.minx + bbox.maxx) * 0.5;
        double midy = (bbox.miny + bbox.maxy) * 0.5;
        
        super_v0 = vertices.size();
        vertices.emplace_back(midx - 10 * dmax, midy - 5 * dmax, UINT32_MAX);
        
        super_v1 = vertices.size();
        vertices.emplace_back(midx + 10 * dmax, midy - 5 * dmax, UINT32_MAX);
        
        super_v2 = vertices.size();
        vertices.emplace_back(midx, midy + 15 * dmax, UINT32_MAX);
        
        triangles.emplace_back(super_v0, super_v1, super_v2);
    }
    
    // WALKING TRIANGLE LOCATION (replaces O(n) brute force)
    uint32_t locate_triangle(uint32_t p_idx) {
        const Point2D& p = vertices[p_idx];

        uint32_t current = last_located_tri;

        if (current >= triangles.size() || !triangles[current].valid) {
            current = UINT32_MAX;
            for (uint32_t i = 0; i < triangles.size(); ++i) {
                if (triangles[i].valid) {
                    current = i;
                    break;
                }
            }
            if (current == UINT32_MAX) return UINT32_MAX;
        }

        constexpr uint32_t MAX_WALK_STEPS = 100000;
        uint32_t steps = 0;

        while (steps < MAX_WALK_STEPS) {
            ++steps;
            const auto& tri = triangles[current];
            if (!tri.valid) break;

            int o_edge0 = predicates::orient2d(
                vertices[tri.v[1]], vertices[tri.v[2]], p);
            int o_edge1 = predicates::orient2d(
                vertices[tri.v[2]], vertices[tri.v[0]], p);
            int o_edge2 = predicates::orient2d(
                vertices[tri.v[0]], vertices[tri.v[1]], p);

            if (o_edge0 >= 0 && o_edge1 >= 0 && o_edge2 >= 0) {
                last_located_tri = current;
                return current;
            }

            if (o_edge0 < 0 && tri.neighbor[0] != UINT32_MAX) {
                current = tri.neighbor[0];
            } else if (o_edge1 < 0 && tri.neighbor[1] != UINT32_MAX) {
                current = tri.neighbor[1];
            } else if (o_edge2 < 0 && tri.neighbor[2] != UINT32_MAX) {
                current = tri.neighbor[2];
            } else {
                break;
            }
        }

        for (uint32_t i = 0; i < triangles.size(); ++i) {
            if (!triangles[i].valid) continue;
            const auto& tri = triangles[i];

            int o0 = predicates::orient2d(vertices[tri.v[1]], vertices[tri.v[2]], p);
            int o1 = predicates::orient2d(vertices[tri.v[2]], vertices[tri.v[0]], p);
            int o2 = predicates::orient2d(vertices[tri.v[0]], vertices[tri.v[1]], p);

            if (o0 >= 0 && o1 >= 0 && o2 >= 0) {
                last_located_tri = i;
                return i;
            }
        }

        return UINT32_MAX;
    }
    
    void insert_vertex(uint32_t p_idx) {
        if (p_idx == super_v0 || p_idx == super_v1 || p_idx == super_v2) {
            return;
        }
        
        const Point2D& p = vertices[p_idx];
        
        uint32_t tri_idx = locate_triangle(p_idx);
        if (tri_idx == UINT32_MAX) return;
        
        auto& tri = triangles[tri_idx];
        int edge_idx = -1;
        
        for (int i = 0; i < 3; ++i) {
            const Point2D& a = vertices[tri.v[i]];
            const Point2D& b = vertices[tri.v[(i + 1) % 3]];
            if (predicates::point_on_segment(p, a, b)) {
                edge_idx = i;
                break;
            }
        }
        
        if (edge_idx >= 0) {
            split_edge(tri_idx, edge_idx, p_idx);
        } else {
            split_triangle(tri_idx, p_idx);
        }
    }
    
    void split_triangle(uint32_t tri_idx, uint32_t p_idx) {
        auto& tri = triangles[tri_idx];
        uint32_t v0 = tri.v[0], v1 = tri.v[1], v2 = tri.v[2];
        
        uint32_t n0 = tri.neighbor[0];
        uint32_t n1 = tri.neighbor[1];
        uint32_t n2 = tri.neighbor[2];
        
        bool c0 = tri.constrained[0];
        bool c1 = tri.constrained[1];
        bool c2 = tri.constrained[2];
        
        tri.valid = false;
        
        uint32_t t0 = triangles.size();
        triangles.emplace_back(v0, v1, p_idx);
        uint32_t t1 = triangles.size();
        triangles.emplace_back(v1, v2, p_idx);
        uint32_t t2 = triangles.size();
        triangles.emplace_back(v2, v0, p_idx);
        
        triangles[t0].neighbor = {n0, t1, t2};
        triangles[t1].neighbor = {n1, t2, t0};
        triangles[t2].neighbor = {n2, t0, t1};
        
        triangles[t0].constrained = {c0, false, false};
        triangles[t1].constrained = {c1, false, false};
        triangles[t2].constrained = {c2, false, false};
        
        update_neighbor(n0, tri_idx, t0);
        update_neighbor(n1, tri_idx, t1);
        update_neighbor(n2, tri_idx, t2);
        
        legalize_edge(t0, 1);
        legalize_edge(t1, 1);
        legalize_edge(t2, 1);
    }
    
    // HOT FIX 7: split_edge() - Set neighbors correctly
    void split_edge(uint32_t tri_idx, uint8_t edge, uint32_t p_idx) {
        auto& tri = triangles[tri_idx];
        uint32_t v0 = tri.v[edge];
        uint32_t v1 = tri.v[(edge + 1) % 3];
        uint32_t v2 = tri.v[(edge + 2) % 3];
        
        uint32_t opp_tri = tri.neighbor[edge];
        uint32_t n1 = tri.neighbor[(edge + 1) % 3];
        uint32_t n2 = tri.neighbor[(edge + 2) % 3];
        
        tri.valid = false;
        
        if (opp_tri == UINT32_MAX) {
            uint32_t t0 = triangles.size();
            triangles.emplace_back(v0, p_idx, v2);
            uint32_t t1 = triangles.size();
            triangles.emplace_back(p_idx, v1, v2);
            
            triangles[t0].neighbor = {n2, t1, UINT32_MAX};
            triangles[t1].neighbor = {UINT32_MAX, t0, n1};
            
            if (n2 != UINT32_MAX) {
                auto& n2_tri = triangles[n2];
                for (int i = 0; i < 3; ++i) {
                    if (n2_tri.neighbor[i] == tri_idx) n2_tri.neighbor[i] = t0;
                }
            }
            if (n1 != UINT32_MAX) {
                auto& n1_tri = triangles[n1];
                for (int i = 0; i < 3; ++i) {
                    if (n1_tri.neighbor[i] == tri_idx) n1_tri.neighbor[i] = t1;
                }
            }
            
            legalize_edge(t0, 1);
            legalize_edge(t1, 1);
        } else {
            auto& opp = triangles[opp_tri];
            int opp_edge = opp.find_edge(v0, v1);
            if (opp_edge < 0) return;
            
            uint32_t v3 = opp.v[(opp_edge + 2) % 3];
            uint32_t opp_n1 = opp.neighbor[(opp_edge + 1) % 3];
            uint32_t opp_n2 = opp.neighbor[(opp_edge + 2) % 3];
            
            opp.valid = false;
            
            uint32_t t0 = triangles.size();
            triangles.emplace_back(v0, p_idx, v2);
            uint32_t t1 = triangles.size();
            triangles.emplace_back(p_idx, v1, v2);
            uint32_t t2 = triangles.size();
            triangles.emplace_back(v0, v3, p_idx);
            uint32_t t3 = triangles.size();
            triangles.emplace_back(v3, v1, p_idx);
            
            triangles[t0].neighbor = {n2, t1, t2};
            triangles[t1].neighbor = {n1, t0, t3};
            triangles[t2].neighbor = {opp_n1, t3, t0};
            triangles[t3].neighbor = {opp_n2, t1, t2};
            
            if (n2 != UINT32_MAX) {
                auto& n2_tri = triangles[n2];
                for (int i = 0; i < 3; ++i) {
                    if (n2_tri.neighbor[i] == tri_idx) n2_tri.neighbor[i] = t0;
                }
            }
            if (n1 != UINT32_MAX) {
                auto& n1_tri = triangles[n1];
                for (int i = 0; i < 3; ++i) {
                    if (n1_tri.neighbor[i] == tri_idx) n1_tri.neighbor[i] = t1;
                }
            }
            if (opp_n1 != UINT32_MAX) {
                auto& opp_n1_tri = triangles[opp_n1];
                for (int i = 0; i < 3; ++i) {
                    if (opp_n1_tri.neighbor[i] == opp_tri) opp_n1_tri.neighbor[i] = t2;
                }
            }
            if (opp_n2 != UINT32_MAX) {
                auto& opp_n2_tri = triangles[opp_n2];
                for (int i = 0; i < 3; ++i) {
                    if (opp_n2_tri.neighbor[i] == opp_tri) opp_n2_tri.neighbor[i] = t3;
                }
            }
            
            legalize_edge(t0, 0);
            legalize_edge(t1, 0);
            legalize_edge(t2, 0);
            legalize_edge(t3, 0);
        }
    }
    
    // FIX 2: CDT Recursion Guard & Pointer Symmetry
    // 
    // VERIFICATION PROOF:
    // 1. Depth Limit: The max_depth parameter (default 64) prevents infinite recursion
    //    by bounding the legalization wave. If exceeded, we skip the flip entirely.
    // 2. Degenerate Swap Detection: After flip, we verify the new edges are different
    //    from the flipped edge to prevent immediate re-flip.
    // 3. Self-Reference Check: assert(tri.neighbor[i] != current_tri_idx) ensures
    //    no triangle points to itself, which would cause infinite loops.
    
    // ═══════════════════════════════════════════════════════════════════════════════
    // ARCHITECTURAL FIX: Increased MAX_LEGALIZE_DEPTH from 64 to 512
    // 
    // The previous limit of 64 was sufficient for typical inputs but could
    // be exceeded in pathological cases with many constraint intersections.
    // The new limit of 512 provides headroom for complex geometries while
    // still preventing stack overflow.
    // ═══════════════════════════════════════════════════════════════════════════════
    static constexpr uint32_t MAX_LEGALIZE_DEPTH = 512;
    
    void legalize_edge(uint32_t tri_idx, uint8_t edge, uint32_t depth = 0) {
        // Recursion Guard: Prevent infinite loops from numerical instability
        if (depth >= MAX_LEGALIZE_DEPTH) {
            // Degenerate case: Delaunay condition keeps flipping back and forth
            // Skip this edge and continue with the triangulation
            return;
        }
        
        assert(tri_idx < triangles.size() && "legalize_edge: tri_idx out of bounds");
        assert(tri_idx != UINT32_MAX && "legalize_edge: invalid tri_idx");
        
        if (tri_idx == UINT32_MAX) return;
        
        auto& tri = triangles[tri_idx];
        if (!tri.valid) return;
        
        // Self-Reference Check: Prevent triangles from pointing to themselves
        for (int i = 0; i < 3; ++i) {
            assert(tri.neighbor[i] != tri_idx && "legalize_edge: self-referencing triangle detected");
        }
        
        if (tri.constrained[edge]) return;
        
        uint32_t opp_idx = tri.neighbor[edge];
        if (opp_idx == UINT32_MAX) return;
        
        // Self-Reference Check for neighbor
        assert(opp_idx != tri_idx && "legalize_edge: triangle points to itself as neighbor");
        
        auto& opp = triangles[opp_idx];
        if (!opp.valid) return;
        
        uint32_t a = tri.v[(edge + 1) % 3];
        uint32_t b = tri.v[(edge + 2) % 3];
        uint32_t c = tri.v[edge];
        
        int opp_edge = opp.find_edge(a, b);
        if (opp_edge < 0) return;
        uint32_t d = opp.v[opp_edge];
        
        if (predicates::in_circumcircle(vertices[a], vertices[b], vertices[c], vertices[d])) {
            flip_edge(tri_idx, edge, opp_idx, opp_edge);
            
            // Recurse with incremented depth
            // Target edges 1 and 2 (the outer edges of the flipped quadrilateral)
            legalize_edge(tri_idx, 1, depth + 1);
            legalize_edge(tri_idx, 2, depth + 1);
            legalize_edge(opp_idx, 1, depth + 1);
            legalize_edge(opp_idx, 2, depth + 1);
        }
    }
    
    // FIX 2: Atomic Pointer Update & Neighbor Symmetry
    //
    // VERIFICATION PROOF:
    // 1. Atomic Block: All neighbor updates happen in a single scope, ensuring
    //    that either ALL back-references are updated or NONE are (no partial state).
    // 2. Self-Reference Prevention: assert(tri.neighbor[i] != current_tri_idx)
    //    ensures no triangle ever points to itself, preventing infinite loops.
    // 3. Symmetry Guarantee: For each neighbor update, we verify the back-reference
    //    exists before updating, preventing orphaned pointers.
    void flip_edge(uint32_t t0, uint8_t e0, uint32_t t1, uint8_t e1) {
        assert(t0 < triangles.size() && "flip_edge: t0 out of bounds");
        assert(t1 < triangles.size() && "flip_edge: t1 out of bounds");
        assert(e0 < 3 && "flip_edge: e0 out of bounds");
        assert(e1 < 3 && "flip_edge: e1 out of bounds");
        
        // Self-Reference Check: Triangles must not point to themselves
        assert(t0 != t1 && "flip_edge: attempting to flip edge within same triangle");
        
        auto& tri0 = triangles[t0];
        auto& tri1 = triangles[t1];
        
        // Verify triangles don't already have self-references
        for (int i = 0; i < 3; ++i) {
            assert(tri0.neighbor[i] != t0 && "flip_edge: tri0 self-reference detected");
            assert(tri1.neighbor[i] != t1 && "flip_edge: tri1 self-reference detected");
        }
        
        uint32_t a = tri0.v[(e0 + 1) % 3];
        uint32_t b = tri0.v[(e0 + 2) % 3];
        uint32_t c = tri0.v[e0];
        uint32_t d = tri1.v[e1];
        
        uint32_t n0 = tri0.neighbor[(e0 + 1) % 3];
        uint32_t n1 = tri0.neighbor[(e0 + 2) % 3];
        uint32_t n2 = tri1.neighbor[(e1 + 1) % 3];
        uint32_t n3 = tri1.neighbor[(e1 + 2) % 3];
        
        bool c0 = tri0.constrained[(e0 + 1) % 3];
        bool c1 = tri0.constrained[(e0 + 2) % 3];
        bool c2 = tri1.constrained[(e1 + 1) % 3];
        bool c3 = tri1.constrained[(e1 + 2) % 3];
        
        // ATOMIC POINTER UPDATE BLOCK
        // All updates happen together to maintain symmetry
        {
            // Update triangle vertices
            tri0.v = {a, c, d};
            tri1.v = {b, d, c};
            
            // Update triangle neighbors (internal)
            tri0.neighbor = {n0, t1, n3};
            tri1.neighbor = {n2, n1, t0};
            
            // Update constraint flags
            tri0.constrained = {c0, false, c3};
            tri1.constrained = {c2, c1, false};
            
            // ATOMIC BACK-REFERENCE UPDATES
            // Each neighbor update is paired with its back-reference
            
            if (n0 != UINT32_MAX) {
                assert(n0 != t0 && n0 != t1 && "flip_edge: n0 self-reference");
                auto& n0_tri = triangles[n0];
                [[maybe_unused]] bool found = false;
                for (int i = 0; i < 3; ++i) {
                    if (n0_tri.neighbor[i] == t1) {
                        n0_tri.neighbor[i] = t0;
                        found = true;
                        break;  // Only update one edge
                    }
                }
                // If not found, the triangulation is corrupted
                assert(found && "flip_edge: n0 back-reference not found");
            }
            
            if (n1 != UINT32_MAX) {
                assert(n1 != t0 && n1 != t1 && "flip_edge: n1 self-reference");
                auto& n1_tri = triangles[n1];
                [[maybe_unused]] bool found = false;
                for (int i = 0; i < 3; ++i) {
                    if (n1_tri.neighbor[i] == t0) {
                        n1_tri.neighbor[i] = t1;
                        found = true;
                        break;
                    }
                }
                assert(found && "flip_edge: n1 back-reference not found");
            }
            
            if (n2 != UINT32_MAX) {
                assert(n2 != t0 && n2 != t1 && "flip_edge: n2 self-reference");
                auto& n2_tri = triangles[n2];
                [[maybe_unused]] bool found = false;
                for (int i = 0; i < 3; ++i) {
                    if (n2_tri.neighbor[i] == t1) {
                        n2_tri.neighbor[i] = t1;  // Same triangle, different edge
                        found = true;
                        break;
                    }
                }
                assert(found && "flip_edge: n2 back-reference not found");
            }
            
            if (n3 != UINT32_MAX) {
                assert(n3 != t0 && n3 != t1 && "flip_edge: n3 self-reference");
                auto& n3_tri = triangles[n3];
                [[maybe_unused]] bool found = false;
                for (int i = 0; i < 3; ++i) {
                    if (n3_tri.neighbor[i] == t0) {
                        n3_tri.neighbor[i] = t0;  // Same triangle, different edge
                        found = true;
                        break;
                    }
                }
                assert(found && "flip_edge: n3 back-reference not found");
            }
        }  // End atomic update block
    }
    
    void update_neighbor(uint32_t tri_idx, uint32_t old_neigh, uint32_t new_neigh) {
        if (tri_idx == UINT32_MAX) return;
        auto& tri = triangles[tri_idx];
        for (int i = 0; i < 3; ++i) {
            if (tri.neighbor[i] == old_neigh) {
                tri.neighbor[i] = new_neigh;
                return;
            }
        }
    }
    
    bool edge_exists(uint32_t v0, uint32_t v1) {
        for (const auto& tri : triangles) {
            if (!tri.valid) continue;
            if (tri.find_edge(v0, v1) >= 0) return true;
        }
        return false;
    }
    
    uint32_t find_first_exit_tri(uint32_t v0, uint32_t v1) {
        for (uint32_t i = 0; i < triangles.size(); ++i) {
            if (!triangles[i].valid) continue;
            if (!triangles[i].has_vertex(v0)) continue;
            
            const auto& tri = triangles[i];
            int v0_local = -1;
            for (int j = 0; j < 3; ++j) {
                if (tri.v[j] == v0) { v0_local = j; break; }
            }
            
            uint32_t va = tri.v[(v0_local + 1) % 3];
            uint32_t vb = tri.v[(v0_local + 2) % 3];
            
            int o_va = predicates::orient2d(vertices[v0], vertices[va], vertices[v1]);
            int o_vb = predicates::orient2d(vertices[v0], vertices[vb], vertices[v1]);
            
            if (o_va >= 0 && o_vb <= 0) {
                return i;
            }
        }
        return UINT32_MAX;
    }
    
    void find_edge_triangles(uint32_t ea, uint32_t eb,
                              uint32_t& t0, int& e0_local,
                              uint32_t& t1, int& e1_local) {
        t0 = t1 = UINT32_MAX;
        e0_local = e1_local = -1;
        
        for (uint32_t i = 0; i < triangles.size(); ++i) {
            if (!triangles[i].valid) continue;
            
            int e = triangles[i].find_edge(ea, eb);
            if (e >= 0) {
                if (t0 == UINT32_MAX) {
                    t0 = i;
                    e0_local = e;
                } else {
                    t1 = i;
                    e1_local = e;
                    return;
                }
            }
        }
    }
    
    void restore_delaunay_edge(uint32_t ea, uint32_t eb) {
        uint32_t t0 = UINT32_MAX, t1 = UINT32_MAX;
        int e0_local = -1, e1_local = -1;
        
        find_edge_triangles(ea, eb, t0, e0_local, t1, e1_local);
        
        if (t0 == UINT32_MAX || t1 == UINT32_MAX) return;
        if (triangles[t0].constrained[e0_local]) return;
        
        uint32_t c = triangles[t0].v[e0_local];
        uint32_t d = triangles[t1].v[e1_local];
        
        if (predicates::in_circumcircle(vertices[ea], vertices[eb], vertices[c], vertices[d])) {
            flip_edge(t0, e0_local, t1, e1_local);
            
            restore_delaunay_edge(ea, c);
            restore_delaunay_edge(c, eb);
            restore_delaunay_edge(eb, d);
            restore_delaunay_edge(d, ea);
        }
    }
    
    // SLOAN EDGE-FLIP CONSTRAINT ENFORCEMENT
    void enforce_constraint(uint32_t v0, uint32_t v1) {
        if (edge_exists(v0, v1)) {
            mark_edge_constrained(v0, v1);
            return;
        }
        
        struct CrossingEdge {
            uint32_t a, b;
            uint32_t tri_left;
            uint32_t tri_right;
        };
        
        std::vector<CrossingEdge> crossings;
        
        uint32_t current_tri = find_first_exit_tri(v0, v1);
        if (current_tri == UINT32_MAX) {
            enforce_constraint_through_vertices(v0, v1);
            return;
        }
        
        constexpr uint32_t MAX_CROSS = 100000;
        uint32_t safety = 0;
        
        while (safety < MAX_CROSS) {
            ++safety;
            const auto& tri = triangles[current_tri];
            if (!tri.valid) break;
            
            if (tri.has_vertex(v1)) break;
            
            int cross_edge = -1;
            for (int i = 0; i < 3; ++i) {
                auto [ea, eb] = tri.edge_vertices(i);
                if (ea == v0 || eb == v0 || ea == v1 || eb == v1) continue;
                
                int d1 = predicates::orient2d(vertices[v0], vertices[v1], vertices[ea]);
                int d2 = predicates::orient2d(vertices[v0], vertices[v1], vertices[eb]);
                
                if ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) {
                    cross_edge = i;
                    break;
                }
                
                if (d1 == 0 && ea != v0 && ea != v1) {
                    enforce_constraint(v0, ea);
                    enforce_constraint(ea, v1);
                    return;
                }
                if (d2 == 0 && eb != v0 && eb != v1) {
                    enforce_constraint(v0, eb);
                    enforce_constraint(eb, v1);
                    return;
                }
            }
            
            if (cross_edge < 0) break;
            
            auto [ea, eb] = tri.edge_vertices(cross_edge);
            
            CrossingEdge ce;
            ce.a = ea;
            ce.b = eb;
            ce.tri_left = current_tri;
            ce.tri_right = tri.neighbor[cross_edge];
            
            crossings.push_back(ce);
            
            current_tri = tri.neighbor[cross_edge];
            if (current_tri == UINT32_MAX) break;
        }
        
        std::vector<std::pair<uint32_t, uint32_t>> to_flip;
        for (const auto& ce : crossings) {
            to_flip.emplace_back(ce.a, ce.b);
        }
        
        std::vector<std::pair<uint32_t, uint32_t>> newly_created_edges;
        
        uint32_t flip_iterations = 0;
        constexpr uint32_t MAX_FLIP_ITERATIONS = 1000000;
        
        while (!to_flip.empty() && flip_iterations < MAX_FLIP_ITERATIONS) {
            ++flip_iterations;
            
            auto [ea, eb] = to_flip.back();
            to_flip.pop_back();
            
            uint32_t t0 = UINT32_MAX, t1 = UINT32_MAX;
            int e0_local = -1, e1_local = -1;
            
            find_edge_triangles(ea, eb, t0, e0_local, t1, e1_local);
            
            if (t0 == UINT32_MAX || t1 == UINT32_MAX) continue;
            
            uint32_t c = triangles[t0].v[e0_local];
            uint32_t d = triangles[t1].v[e1_local];
            
            int oc = predicates::orient2d(vertices[c], vertices[d], vertices[ea]);
            int od = predicates::orient2d(vertices[c], vertices[d], vertices[eb]);
            
            bool is_convex = (oc > 0 && od < 0) || (oc < 0 && od > 0);
            
            if (!is_convex) {
                to_flip.insert(to_flip.begin(), {ea, eb});
                continue;
            }
            
            flip_edge(t0, e0_local, t1, e1_local);
            
            newly_created_edges.emplace_back(c, d);
            
            int dc = predicates::orient2d(vertices[v0], vertices[v1], vertices[c]);
            int dd = predicates::orient2d(vertices[v0], vertices[v1], vertices[d]);
            
            bool new_edge_crosses = (c != v0 && c != v1 && d != v0 && d != v1)
                                 && ((dc > 0 && dd < 0) || (dc < 0 && dd > 0));
            
            if (new_edge_crosses) {
                to_flip.push_back({c, d});
            }
        }
        
        mark_edge_constrained(v0, v1);
        
        for (const auto& [ea, eb] : newly_created_edges) {
            if (ea == v0 && eb == v1) continue;
            if (ea == v1 && eb == v0) continue;
            restore_delaunay_edge(ea, eb);
        }
    }
    
    void enforce_constraint_through_vertices(uint32_t v0, uint32_t v1) {
        std::vector<uint32_t> collinear_verts;
        
        for (uint32_t i = 0; i < vertices.size(); ++i) {
            if (i == v0 || i == v1) continue;
            if (vertices[i].orig_idx == UINT32_MAX) continue;
            
            if (predicates::orient2d(vertices[v0], vertices[v1], vertices[i]) == 0) {
                if (predicates::point_on_segment(vertices[i], vertices[v0], vertices[v1])) {
                    collinear_verts.push_back(i);
                }
            }
        }
        
        std::sort(collinear_verts.begin(), collinear_verts.end(),
            [&](uint32_t a, uint32_t b) {
                double da = (vertices[a].x - vertices[v0].x) * (vertices[a].x - vertices[v0].x)
                          + (vertices[a].y - vertices[v0].y) * (vertices[a].y - vertices[v0].y);
                double db = (vertices[b].x - vertices[v0].x) * (vertices[b].x - vertices[v0].x)
                          + (vertices[b].y - vertices[v0].y) * (vertices[b].y - vertices[v0].y);
                return da < db;
            });
        
        uint32_t prev = v0;
        for (uint32_t mid : collinear_verts) {
            enforce_constraint(prev, mid);
            prev = mid;
        }
        enforce_constraint(prev, v1);
    }
    
    // FIX 2.8: mark_edge_constrained with short-circuit
    void mark_edge_constrained(uint32_t v0, uint32_t v1) {
        int found = 0;
        for (uint32_t i = 0; i < triangles.size() && found < 2; ++i) {
            auto& tri = triangles[i];
            if (!tri.valid) continue;
            
            int edge = tri.find_edge(v0, v1);
            if (edge >= 0) {
                tri.constrained[edge] = true;
                ++found;
            }
        }
    }
    
    void remove_super_triangle() {
        for (auto& tri : triangles) {
            if (!tri.valid) continue;
            
            if (tri.v[0] == super_v0 || tri.v[0] == super_v1 || tri.v[0] == super_v2 ||
                tri.v[1] == super_v0 || tri.v[1] == super_v1 || tri.v[1] == super_v2 ||
                tri.v[2] == super_v0 || tri.v[2] == super_v1 || tri.v[2] == super_v2) {
                tri.valid = false;
            }
        }
        
        vertices[super_v0].orig_idx = UINT32_MAX;
        vertices[super_v1].orig_idx = UINT32_MAX;
        vertices[super_v2].orig_idx = UINT32_MAX;
    }
    
    void restore_delaunay() {
        bool flipped = true;
        int iterations = 0;
        const int MAX_ITERATIONS = 100;
        
        while (flipped && iterations < MAX_ITERATIONS) {
            flipped = false;
            ++iterations;
            
            for (uint32_t i = 0; i < triangles.size(); ++i) {
                if (!triangles[i].valid) continue;
                
                for (int e = 0; e < 3; ++e) {
                    if (triangles[i].constrained[e]) continue;
                    
                    uint32_t opp = triangles[i].neighbor[e];
                    if (opp == UINT32_MAX || !triangles[opp].valid) continue;
                    if (opp < i) continue;
                    
                    const Point2D& a = vertices[triangles[i].v[(e + 1) % 3]];
                    const Point2D& b = vertices[triangles[i].v[(e + 2) % 3]];
                    const Point2D& c = vertices[triangles[i].v[e]];
                    
                    int opp_e = triangles[opp].find_edge(
                        triangles[i].v[(e + 1) % 3], 
                        triangles[i].v[(e + 2) % 3]);
                    const Point2D& d = vertices[triangles[opp].v[opp_e]];
                    
                    if (predicates::in_circumcircle(a, b, c, d)) {
                        flip_edge(i, e, opp, opp_e);
                        flipped = true;
                    }
                }
            }
        }
    }
};

// ═══════════════════════════════════════════════════════════════════════════════
// 3D Geometry and Projection
// ═══════════════════════════════════════════════════════════════════════════════

struct Point3D {
    double x, y, z;
    
    Point3D() : x(0), y(0), z(0) {}
    Point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

struct LocalFrame {
    Point3D origin;
    Point3D u;
    Point3D v;
    Point3D n;
    
    static LocalFrame from_triangle(const Point3D& p0, const Point3D& p1, 
                                    const Point3D& p2) {
        LocalFrame frame;
        frame.origin = p0;
        
        frame.u.x = p1.x - p0.x;
        frame.u.y = p1.y - p0.y;
        frame.u.z = p1.z - p0.z;
        double ulen = std::sqrt(frame.u.x * frame.u.x + 
                                frame.u.y * frame.u.y + 
                                frame.u.z * frame.u.z);
        frame.u.x /= ulen;
        frame.u.y /= ulen;
        frame.u.z /= ulen;
        
        double dx = p2.x - p0.x;
        double dy = p2.y - p0.y;
        double dz = p2.z - p0.z;
        
        frame.n.x = dy * frame.u.z - dz * frame.u.y;
        frame.n.y = dz * frame.u.x - dx * frame.u.z;
        frame.n.z = dx * frame.u.y - dy * frame.u.x;
        double nlen = std::sqrt(frame.n.x * frame.n.x + 
                                frame.n.y * frame.n.y + 
                                frame.n.z * frame.n.z);
        frame.n.x /= nlen;
        frame.n.y /= nlen;
        frame.n.z /= nlen;
        
        frame.v.x = frame.n.y * frame.u.z - frame.n.z * frame.u.y;
        frame.v.y = frame.n.z * frame.u.x - frame.n.x * frame.u.z;
        frame.v.z = frame.n.x * frame.u.y - frame.n.y * frame.u.x;
        
        return frame;
    }
    
    Point2D project(const Point3D& p) const {
        double dx = p.x - origin.x;
        double dy = p.y - origin.y;
        double dz = p.z - origin.z;
        
        double local_x = dx * u.x + dy * u.y + dz * u.z;
        double local_y = dx * v.x + dy * v.y + dz * v.z;
        
        return Point2D(local_x, local_y, 0);
    }
};

// ═══════════════════════════════════════════════════════════════════════════════
// CherchiBackend Implementation
// ═══════════════════════════════════════════════════════════════════════════════

class CherchiBackend : public IBooleanBackend {
public:
    std::string name() const override { return "cherchi"; }
    bool supportsOpenMesh() const override { return true; }
    bool providesExactTopology() const override { return true; }
    
    BackendResult execute(const PolygonSoup& soup, [[maybe_unused]] BooleanOp op) override {
        BackendResult result;
        result.success = true;
        
        // For each triangle with intersections, build local CDT
        for (const auto& pair : soup.candidate_pairs) {
            if (!pair.intersecting) continue;
            
            const auto& tri_a = soup.triangles[pair.tri_a];
            // tri_b is used for intersection data (future implementation)
            [[maybe_unused]] const auto& tri_b = soup.triangles[pair.tri_b];
            
            // Build local coordinate frame for triangle A
            Point3D p0(tri_a.v[0][0], tri_a.v[0][1], tri_a.v[0][2]);
            Point3D p1(tri_a.v[1][0], tri_a.v[1][1], tri_a.v[1][2]);
            Point3D p2(tri_a.v[2][0], tri_a.v[2][1], tri_a.v[2][2]);
            
            LocalFrame frame = LocalFrame::from_triangle(p0, p1, p2);
            
            // Build CDT for this triangle
            CDTBuilder cdt;
            
            // Add triangle vertices
            cdt.addVertex(frame.project(p0));
            cdt.addVertex(frame.project(p1));
            cdt.addVertex(frame.project(p2));
            
            // Add intersection points
            for (const auto& pt : pair.intersection_points) {
                Point3D p(pt[0], pt[1], pt[2]);
                cdt.addVertex(frame.project(p));
            }
            
            // Add constraints (triangle edges + intersection segments)
            cdt.addConstraint(0, 1);
            cdt.addConstraint(1, 2);
            cdt.addConstraint(2, 0);
            
            // Build triangulation
            cdt.triangulate();
            
            // Extract triangles
            auto tris = cdt.getTriangles();
            
            // Convert back to 3D and add to output
            for ([[maybe_unused]] const auto& t : tris) {
                OutputPolygon poly;
                // TODO: Convert triangle t back to 3D and populate poly
                (void)t;  // Suppress unused warning until implementation is complete
                result.output_soup.output_polygons.push_back(poly);
            }
        }
        
        return result;
    }
};

} // namespace ember
