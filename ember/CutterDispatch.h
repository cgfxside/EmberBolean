#pragma once
// ═══════════════════════════════════════════════════════════════
// EMBER CutterDispatch — Automatic Tier 1/Tier 2 Cutter Routing
// ═══════════════════════════════════════════════════════════════
//
// Sits inside Stage 01 (Diagnostic Gatekeeper).
// Analyzes each connected component of Input B and assigns
// a dispatch tier:
//
//   Tier 1:  Planar        → Stage 02b analytical fast-path
//   Tier 1G: PlanarGrout   → Stage 02b with dual-plane pairs
//   Tier 2:  NoisedMesh    → Tessellate + Stage 02/03 exact
//   Tier 2:  GeneralMesh   → Stage 02/03 exact (arbitrary mesh)
//
// The dispatch is fully automatic. The artist never selects a tier.
// Noise Amplitude > 0 promotes planar cutters to Tier 2.

#include <cstdint>
#include <vector>
#include "PolygonSoup.h"
#include "Plane.h"
#include "IntegerTypes.h"
#include "Parameters.h"

namespace ember {

// ─── Cutter Classification ───────────────────────────────────

enum class CutterType : uint8_t {
    Planar,          // All faces coplanar → single exact plane equation
    PlanarGrout,     // Planar + grout offset → dual-plane pair
    NoisedMesh,      // Tessellated noised cutter → triangle soup
    GeneralMesh,     // Arbitrary user mesh (Input 2) → triangle soup
};

// ─── Per-Cutter Descriptor ───────────────────────────────────

struct CutterDescriptor {
    CutterType  type;
    uint32_t    component_id;      // Connected component index in Input B
    uint32_t    mesh_id;           // Mesh ID for N-way WNV assignment

    // Tier 1 fields (valid when type == Planar or PlanarGrout)
    IntPlane    plane;             // 26-bit integer coefficients (a,b,c,d)
    IntPlane    grout_plane_pos;   // +offset plane (PlanarGrout only)
    IntPlane    grout_plane_neg;   // -offset plane (PlanarGrout only)

    // Tier 2 fields (valid when type == NoisedMesh or GeneralMesh)
    uint32_t    tri_start;         // Start index in PolygonSoup triangle buffer
    uint32_t    tri_count;         // Number of triangles for this cutter

    // Metadata for output tagging
    float       grout_width;       // 0 if no grout
    uint32_t    noise_seed;        // 0 if no noise
};

// ─── Aggregate Dispatch Result ───────────────────────────────

struct DispatchResult {
    std::vector<CutterDescriptor> cutters;

    // Aggregate flags for pipeline fast-path decisions
    bool all_planar;              // True if every cutter is Tier 1
    bool any_grout;               // True if any cutter has grout > 0
    bool any_noised;              // True if any cutter is Tier 2 (noised/general)
    bool force_exact;             // User override: force all to Tier 2

    uint32_t tier1_count;         // Number of analytical plane cutters
    uint32_t tier2_count;         // Number of mesh-mesh cutters

    // Tier 1 cutters: plane equations for Stage 02b
    std::vector<IntPlane> plane_equations;

    // Convenience: is this a pure fast-path cook?
    bool isFastPathOnly() const {
        return all_planar && !force_exact && tier2_count == 0;
    }

    // Convenience: does this cook need Stage 02 + 03?
    bool needsExactPath() const {
        return tier2_count > 0 || force_exact;
    }

    // Convenience: is this a mixed dispatch?
    bool isMixed() const {
        return tier1_count > 0 && tier2_count > 0 && !force_exact;
    }
};

// ─── Analysis Entry Point ────────────────────────────────────

// Analyze all connected components of Input B and classify each.
// Called from Stage 01 after triangulation + degenerate removal,
// BEFORE quantization.
DispatchResult analyzeCutters(
    const PolygonSoup& soup,
    const Parameters& params
);

// ─── Grout Plane Computation ─────────────────────────────────

// Compute dual grout planes from a base plane and grout width.
// Offset is applied along the plane normal in integer space.
void computeGroutPlanes(
    const IntPlane& base,
    float grout_width,
    double quantization_scale,
    IntPlane& out_pos,
    IntPlane& out_neg
);

} // namespace ember
