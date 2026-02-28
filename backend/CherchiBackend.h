#pragma once
// ═══════════════════════════════════════════════════════════════════════════════
// EMBER CherchiBackend — Exact Boolean Operations via Local Arrangements
// ═══════════════════════════════════════════════════════════════════════════════
//
// Backend implementation using the Cherchi et al. 2020 framework for
// exact mesh Boolean operations.
//
// Key Features:
//   - Exact geometric predicates (orient2d, orient3d, insphere)
//   - Indirect predicates for implicit points (LPI, TPI)
//   - Constrained Delaunay Triangulation (CDT) for robust local triangulation
//   - No epsilon-based comparisons in topology decisions
//
// Reference: "Exact and Efficient Booleans for Polyhedra" (Cherchi et al. 2020)

#include "IBooleanBackend.h"

namespace ember {

class CherchiBackend : public IBooleanBackend {
public:
    // Backend identification
    std::string name() const override { return "cherchi"; }
    
    // Capabilities
    bool supportsOpenMesh() const override { return true; }
    bool providesExactTopology() const override { return true; }
    
    // Main execution method
    BackendResult execute(const PolygonSoup& soup, BooleanOp op) override;

private:
    // Internal types for the Cherchi algorithm
    struct ImplicitPoint;
    struct ArrangementPoint;
    struct TriangleIntersection;
    struct LocalArrangement;
    
    // Pipeline stages
    void detectIntersections(const PolygonSoup& soup,
                             std::vector<TriangleIntersection>& intersections);
    
    void buildLocalArrangements(const PolygonSoup& soup,
                                const std::vector<TriangleIntersection>& intersections,
                                std::vector<LocalArrangement>& arrangements);
    
    void classifyCells(const std::vector<LocalArrangement>& arrangements,
                       BooleanOp op,
                       std::vector<bool>& cell_inside);
    
    void extractResult(const PolygonSoup& soup,
                       const std::vector<LocalArrangement>& arrangements,
                       const std::vector<bool>& cell_inside,
                       BackendResult& result);
};

} // namespace ember
