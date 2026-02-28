#pragma once
// ═══════════════════════════════════════════════════════════════
// EMBER CutterMeshGenerator — Tier 2 Noised Cutter Mesh Builder
// ═══════════════════════════════════════════════════════════════
//
// Generates tessellated, noise-displaced cutter meshes from
// plane equations. Used when CutterDispatch promotes a planar
// cutter to Tier 2 (noise_amplitude > 0).
//
// Pipeline position: Called AFTER CutterDispatch analysis,
// BEFORE 26-bit quantization. All work is in float space.
//
// Key invariant: Grout pairs use SYNCHRONIZED NOISE.
// Both sides share the same N(P), offset by ±grout_width/2.
// This guarantees the grout channel never inverts.

#include <cstdint>
#include <vector>
#include "CutterDispatch.h"
#include "PolygonSoup.h"
#include "Plane.h"

namespace ember {

// ─── Noise Configuration ─────────────────────────────────────

struct NoiseParams {
    float    amplitude_a;       // Layer A: large-scale displacement
    float    frequency_a;
    int      octaves_a;
    float    amplitude_b;       // Layer B: fine detail
    float    frequency_b;
    int      octaves_b;
    uint32_t seed;
};

// ─── Cutter Mesh Configuration ───────────────────────────────

struct CutterMeshConfig {
    float       segment_size;       // Max edge length for tessellation
    int         max_segments;       // Cap on subdivision per axis
    float       grout_width;        // 0 = no grout
    float       shell_noise_amp;    // Varies grout thickness spatially
    NoiseParams noise;
    bool        cull_external;      // Remove faces outside target hull
};

// ─── Generation Entry Point ──────────────────────────────────

// Generate a tessellated, noise-displaced cutter mesh for one
// cutter plane. The mesh is appended to out_soup and the
// descriptor is updated with tri_start/tri_count.
//
// Called in float space before quantization.
void generateNoisedCutterMesh(
    const IntPlane& base_plane,
    const CutterMeshConfig& config,
    const PolygonSoup& target_soup,     // For convex hull culling
    PolygonSoup& out_soup,              // Cutter tris appended here
    CutterDescriptor& out_desc          // Updated with tri_start/count
);

// ─── Batch Generation ────────────────────────────────────────

// Generate noised meshes for ALL Tier 2 NoisedMesh cutters
// in the dispatch result. Modifies soup and descriptors in-place.
void generateAllNoisedCutters(
    const CutterMeshConfig& config,
    const PolygonSoup& target_soup,
    PolygonSoup& soup,
    DispatchResult& dispatch
);

} // namespace ember
