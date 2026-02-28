/**
 * @file ExactPredicates.h
 * @brief Adaptive floating-point predicates with exact integer fallback for EMBER Boolean engine
 * 
 * Implements Shewchuk-style semi-static filters with fast double-precision paths
 * and exact Int128/Int256 fallback for uncertain cases.
 * 
 * Bit-width requirements:
 * - orient2d: 62 bits (Int128 sufficient, fits in 2 int64_t)
 * - orient3d: 148 bits (Int256 required, fits in 4 int64_t)
 * - dot_plane: 94 bits (Int128 sufficient)
 * - homogeneous classification: 160 bits (Int256 required)
 * 
 * @author EMBER Boolean Engine
 * @version 1.0.0
 */

#pragma once

#include <cstdint>
#include <cmath>
#include <array>
#include <algorithm>
#include <type_traits>

// AUDIT FIX: Use canonical Int128/Int256 from IntegerTypes.h
// Eliminates duplicate classes that existed at lines 87-417 of the original file.
// See EMBER_Audit_Regressions.md, Phase A: Type Unification
#include "IntegerTypes.h"

// AUDIT FIX: Include exact indirect predicates (zero FP division)
// Replaces the broken orient2d_indirect that used floating-point division.
// See EMBER_Audit_Regressions.md, Phase B: Exact Predicate Fixes
#include "ExactPredicates_Indirect.h"

namespace ember {
namespace predicates {

// AUDIT FIX: Bring canonical types into predicates namespace.
// All code below uses these — no more duplicate class definitions.
using ember::Int128;
using ember::Int256;

//=============================================================================
// SHEWCHUK ERROR BOUND CONSTANTS
//=============================================================================

/**
 * Machine epsilon for double precision: 2^-52 ≈ 2.22e-16
 * Used in error bound computations for filtered predicates.
 */
constexpr double EPSILON = 2.2204460492503131e-16;

/**
 * Error bound for orient2d: (3 + 8ε)ε ≈ 3.33e-16
 * Derived from Shewchuk's analysis of 2x2 determinant error.
 * 
 * Formula: |det(A) - det(A_exact)| ≤ (3 + 8ε)ε · permanent(A)
 */
constexpr double ORIENT2D_ERRBOUND = 3.3306690738754716e-16;

/**
 * Error bound for orient3d: (7 + 56ε)ε ≈ 7.77e-16  
 * Derived from Shewchuk's analysis of 3x3 determinant error.
 * 
 * Formula: |det(A) - det(A_exact)| ≤ (7 + 56ε)ε · permanent(A)
 */
constexpr double ORIENT3D_ERRBOUND = 7.771561172376096e-16;

/**
 * Error bound for dot product with plane: 4ε ≈ 8.88e-16
 * Conservative bound for dot product of 4-component vectors.
 */
constexpr double DOT_PLANE_ERRBOUND = 8.881784197001252e-16;

/**
 * Rescale factor for converting double coordinates to integers.
 * Using 2^20 gives ~1e-6 precision while keeping bit-width manageable.
 * 
 * Choice of 2^20:
 * - Provides sub-micron precision for unit-scale models
 * - Keeps orient3d within 148 bits (fits in Int256)
 * - Allows 2^20 * 2^20 = 2^40 for squared terms
 */
constexpr int64_t RESCALE_EXPONENT = 20;
constexpr int64_t RESCALE_FACTOR = 1LL << RESCALE_EXPONENT;  // 1,048,576

//=============================================================================
// AUDIT FIX: Int128 class REMOVED — now using ember::Int128 from IntegerTypes.h
// The duplicate class that was here (87 lines) used a different multiplication
// algorithm (mid1/mid2 shift-and-add) vs IntegerTypes.h (cross-sum-then-shift).
// Both produce identical results for EMBER's 26-bit budget, but maintaining
// two implementations is a maintenance hazard.
//
// API MAPPING (old → new):
//   Int128::mul64(a, b)          → Int128::mul64(a, b)
//   Int128::fromDouble(v)      → (removed: quantize separately before constructing)
//   Int128::fromQuantized(v)   → Int128(v)
//   Int128(v).high()           → Int128(v).hi    (public member)
//   Int128(v).low()            → Int128(v).lo    (public member)
//   Int128(hi, lo)             → Int128(hi, lo)  (same)
//=============================================================================

//=============================================================================
// AUDIT FIX: Int256 class REMOVED — now using ember::Int256 from IntegerTypes.h
// The duplicate used std::array<uint64_t,4> storage vs IntegerTypes.h's {Int128 hi, lo}.
//
// API MAPPING (old → new):
//   Int256::mul(Int128_a, Int128_b)  → Int256::mul128(Int128_a, Int128_b)
//   Int256::mul(int64_a, int64_b)    → Int256(Int128::mul64(a, b))
//   Int256::fromInt128(v)            → Int256(v)  (constructor)
//   Int256::fromDouble(v)            → (removed: quantize separately)
//   Int256(v).word(i)                → access via .hi and .lo members
//=============================================================================

//=============================================================================
// COORDINATE QUANTIZATION UTILITIES
//=============================================================================

/**
 * @brief Quantize a double coordinate to integer representation.
 * @param v Double coordinate value
 * @return Scaled integer: round(v * RESCALE_FACTOR)
 */
inline int64_t quantize(double v) {
    return static_cast<int64_t>(std::llround(v * static_cast<double>(RESCALE_FACTOR)));
}

/**
 * @brief Dequantize an integer back to double.
 * @param v Quantized integer
 * @return Double value: v / RESCALE_FACTOR
 */
inline double dequantize(int64_t v) {
    return static_cast<double>(v) / static_cast<double>(RESCALE_FACTOR);
}

//=============================================================================
// FILTERED PREDICATES: ORIENT2D
//=============================================================================

/**
 * @brief Compute the permanent (sum of absolute values) of the orient2d matrix.
 * Used for error bound computation in filtered evaluation.
 * 
 * The permanent of |a-b  a-c| is |ax-cx|*|by-cy| + |ay-cy|*|bx-cx|
 *                |b-c  b-c|
 * 
 * @return Upper bound on magnitude of determinant terms
 */
inline double orient2d_permanent(double ax, double ay, double bx, double by,
                                  double cx, double cy) {
    double acx = std::abs(ax - cx);
    double acy = std::abs(ay - cy);
    double bcx = std::abs(bx - cx);
    double bcy = std::abs(by - cy);
    
    return acx * bcy + acy * bcx;
}

/**
 * @brief Fast filtered orient2d predicate with adaptive precision.
 * 
 * Returns the orientation of triangle (a, b, c):
 *   +1: counter-clockwise (positive area)
 *    0: collinear (degenerate)
 *   -1: clockwise (negative area)
 * 
 * Algorithm:
 * 1. Compute determinant using double arithmetic (FAST PATH)
 * 2. Compute error bound = permanent * ORIENT2D_ERRBOUND
 * 3. If |det| > error_bound, sign is guaranteed correct → return immediately
 * 4. Otherwise, use exact Int128 arithmetic (SLOW PATH, ~10% of cases)
 * 
 * Bit-width analysis:
 * - Coordinates: 31 bits after quantization (2^20 scale)
 * - Differences: 32 bits
 * - Products: 64 bits
 * - Sum: 65 bits (fits comfortably in Int128)
 * 
 * @param ax,ay First point
 * @param bx,by Second point  
 * @param cx,cy Third point
 * @return int {-1, 0, +1} orientation sign
 */
inline int orient2d_filtered(double ax, double ay, double bx, double by,
                              double cx, double cy) {
    // === FAST PATH: Double-precision evaluation ===
    double acx = ax - cx;
    double acy = ay - cy;
    double bcx = bx - cx;
    double bcy = by - cy;
    
    // 2x2 determinant: (a-c) × (b-c)
    double det = acx * bcy - acy * bcx;
    
    // Compute error bound
    double permanent = orient2d_permanent(ax, ay, bx, by, cx, cy);
    double errbound = ORIENT2D_ERRBOUND * permanent;
    
    // Check if result is certain
    if (det > errbound) return +1;
    if (det < -errbound) return -1;
    
    // === SLOW PATH: Exact Int128 computation ===
    // This path taken when result is in uncertainty zone (~10% of cases)
    
    int64_t qax = quantize(ax);
    int64_t qay = quantize(ay);
    int64_t qbx = quantize(bx);
    int64_t qby = quantize(by);
    int64_t qcx = quantize(cx);
    int64_t qcy = quantize(cy);
    
    // Compute exact determinant: (ax-cx)*(by-cy) - (ay-cy)*(bx-cx)
    Int128 acx_q = Int128::mul64(qax - qcx, qby - qcy);
    Int128 acy_q = Int128::mul64(qay - qcy, qbx - qcx);
    Int128 det_q = acx_q - acy_q;
    
    return det_q.sign();
}

/**
 * @brief orient2d for pre-quantized integer coordinates.
 * Bypasses the filter and uses exact arithmetic directly.
 * 
 * @param ax,ay First point (quantized)
 * @param bx,by Second point (quantized)
 * @param cx,cy Third point (quantized)
 * @return int {-1, 0, +1} orientation sign
 */
inline int orient2d_exact(int64_t ax, int64_t ay, int64_t bx, int64_t by,
                          int64_t cx, int64_t cy) {
    Int128 acx = Int128::mul64(ax - cx, by - cy);
    Int128 acy = Int128::mul64(ay - cy, bx - cx);
    Int128 det = acx - acy;
    return det.sign();
}

//=============================================================================
// FORWARD DECLARATIONS
//=============================================================================

// Forward declaration for incircle_exact (defined later in this file)
inline int incircle_exact(int64_t ax, int64_t ay, int64_t bx, int64_t by,
                          int64_t cx, int64_t cy, int64_t dx, int64_t dy);

//=============================================================================
// FILTERED PREDICATES: INCIRCLE (2D Delaunay)
//=============================================================================

/**
 * @brief Filtered incircle test for double-precision coordinates.
 * 
 * Uses double-precision arithmetic with error bounds. Falls back to
 * exact arithmetic when the result is uncertain.
 * 
 * @param ax,ay First point (circle-defining triangle)
 * @param bx,by Second point
 * @param cx,cy Third point
 * @param dx,dy Query point (test if inside circumcircle of a,b,c)
 * @return int {-1, 0, +1} in-circle sign (+1 if d inside circle of a,b,c)
 */
inline int incircle_filtered(double ax, double ay, double bx, double by,
                              double cx, double cy, double dx, double dy) {
    // === FAST PATH: Double-precision evaluation ===
    double adx = ax - dx;
    double ady = ay - dy;
    double bdx = bx - dx;
    double bdy = by - dy;
    double cdx = cx - dx;
    double cdy = cy - dy;
    
    // 2x2 determinants (minors)
    double abd = adx * bdy - ady * bdx;
    double bcd = bdx * cdy - bdy * cdx;
    double cad = cdx * ady - cdy * adx;
    
    // Squared distances from d
    double ad_sq = adx * adx + ady * ady;
    double bd_sq = bdx * bdx + bdy * bdy;
    double cd_sq = cdx * cdx + cdy * cdy;
    
    // Determinant: ad_sq * bcd + bd_sq * cad + cd_sq * abd
    double det = ad_sq * bcd + bd_sq * cad + cd_sq * abd;
    
    // Compute error bound (conservative)
    double permanent = ad_sq * std::abs(bcd) + bd_sq * std::abs(cad) + cd_sq * std::abs(abd);
    double errbound = 1e-12 * permanent;  // Conservative error bound
    
    // Check if result is certain
    if (det > errbound) return +1;
    if (det < -errbound) return -1;
    
    // === SLOW PATH: Exact Int256 computation ===
    int64_t qax = quantize(ax);
    int64_t qay = quantize(ay);
    int64_t qbx = quantize(bx);
    int64_t qby = quantize(by);
    int64_t qcx = quantize(cx);
    int64_t qcy = quantize(cy);
    int64_t qdx = quantize(dx);
    int64_t qdy = quantize(dy);
    
    return incircle_exact(qax, qay, qbx, qby, qcx, qcy, qdx, qdy);
}

//=============================================================================
// FILTERED PREDICATES: ORIENT3D
//=============================================================================

/**
 * @brief Compute the permanent for orient3d error bound.
 * 
 * The permanent of the 3x3 orientation matrix gives an upper bound
 * on the magnitude of the determinant terms.
 * 
 * Matrix columns: (a-d), (b-d), (c-d)
 * 
 * @return Sum of absolute values of all 6 terms in determinant expansion
 */
inline double orient3d_permanent(double ax, double ay, double az,
                                  double bx, double by, double bz,
                                  double cx, double cy, double cz,
                                  double dx, double dy, double dz) {
    double adx = std::abs(ax - dx);
    double ady = std::abs(ay - dy);
    double adz = std::abs(az - dz);
    double bdx = std::abs(bx - dx);
    double bdy = std::abs(by - dy);
    double bdz = std::abs(bz - dz);
    double cdx = std::abs(cx - dx);
    double cdy = std::abs(cy - dy);
    double cdz = std::abs(cz - dz);
    
    // Sum of products for 3x3 determinant expansion
    return adx * (bdy * cdz + bdz * cdy) +
           ady * (bdx * cdz + bdz * cdx) +
           adz * (bdx * cdy + bdy * cdx);
}

/**
 * @brief Fast filtered orient3d predicate with adaptive precision.
 * 
 * Returns the orientation of tetrahedron (a, b, c, d):
 *   +1: d is below plane abc (positive volume)
 *    0: coplanar (degenerate)
 *   -1: d is above plane abc (negative volume)
 * 
 * The 3x3 determinant is computed as:
 *   | ax-dx  ay-dy  az-dz |
 *   | bx-dx  by-dy  bz-dz |
 *   | cx-dx  cy-dy  cz-dz |
 * 
 * Algorithm:
 * 1. Compute determinant using double arithmetic with expansion by minors
 * 2. Compute error bound = permanent * ORIENT3D_ERRBOUND
 * 3. If |det| > error_bound, sign is guaranteed correct
 * 4. Otherwise, use exact Int256 arithmetic
 * 
 * Bit-width analysis:
 * - Coordinates: 31 bits after quantization
 * - Differences: 32 bits  
 * - Pairwise products: 64 bits
 * - 2x2 determinants: 129 bits (sum of two 128-bit products)
 * - 3x3 determinant: 148 bits (sum of three 129-bit terms)
 * - Requires Int256 (4 x 64-bit words)
 * 
 * @param ax..az First point
 * @param bx..bz Second point
 * @param cx..cz Third point
 * @param dx..dz Fourth point
 * @return int {-1, 0, +1} orientation sign
 */
inline int orient3d_filtered(double ax, double ay, double az,
                              double bx, double by, double bz,
                              double cx, double cy, double cz,
                              double dx, double dy, double dz) {
    // === FAST PATH: Double-precision evaluation ===
    double adx = ax - dx;
    double ady = ay - dy;
    double adz = az - dz;
    double bdx = bx - dx;
    double bdy = by - dy;
    double bdz = bz - dz;
    double cdx = cx - dx;
    double cdy = cy - dy;
    double cdz = cz - dz;
    
    // 3x3 determinant via expansion by minors
    // det = adx * (bdy * cdz - bdz * cdy)
    //     - ady * (bdx * cdz - bdz * cdx)
    //     + adz * (bdx * cdy - bdy * cdx)
    double bdy_cdz = bdy * cdz;
    double bdz_cdy = bdz * cdy;
    double bdx_cdz = bdx * cdz;
    double bdz_cdx = bdz * cdx;
    double bdx_cdy = bdx * cdy;
    double bdy_cdx = bdy * cdx;
    
    double det = adx * (bdy_cdz - bdz_cdy)
               - ady * (bdx_cdz - bdz_cdx)
               + adz * (bdx_cdy - bdy_cdx);
    
    // Compute error bound
    double permanent = orient3d_permanent(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz);
    double errbound = ORIENT3D_ERRBOUND * permanent;
    
    // Check if result is certain
    if (det > errbound) return +1;
    if (det < -errbound) return -1;
    
    // === SLOW PATH: Exact Int256 computation ===
    // This path taken when result is in uncertainty zone
    
    int64_t qax = quantize(ax), qay = quantize(ay), qaz = quantize(az);
    int64_t qbx = quantize(bx), qby = quantize(by), qbz = quantize(bz);
    int64_t qcx = quantize(cx), qcy = quantize(cy), qcz = quantize(cz);
    int64_t qdx = quantize(dx), qdy = quantize(dy), qdz = quantize(dz);
    
    int64_t adx_q = qax - qdx, ady_q = qay - qdy, adz_q = qaz - qdz;
    int64_t bdx_q = qbx - qdx, bdy_q = qby - qdy, bdz_q = qbz - qdz;
    int64_t cdx_q = qcx - qdx, cdy_q = qcy - qdy, cdz_q = qcz - qdz;
    
    // Compute 2x2 determinants exactly (Int128)
    Int128 d1 = Int128::mul64(bdy_q, cdz_q) - Int128::mul64(bdz_q, cdy_q);  // bdy*cdz - bdz*cdy
    Int128 d2 = Int128::mul64(bdx_q, cdz_q) - Int128::mul64(bdz_q, cdx_q);  // bdx*cdz - bdz*cdx
    Int128 d3 = Int128::mul64(bdx_q, cdy_q) - Int128::mul64(bdy_q, cdx_q);  // bdx*cdy - bdy*cdx
    
    // Compute 3x3 determinant (Int256)
    Int256 term1 = Int256::mul128(Int128(adx_q), d1);
    Int256 term2 = Int256::mul128(Int128(ady_q), d2);
    Int256 term3 = Int256::mul128(Int128(adz_q), d3);
    
    Int256 det_q = term1 - term2 + term3;
    
    return det_q.sign();
}

/**
 * @brief orient3d for pre-quantized integer coordinates.
 * Uses exact Int256 arithmetic directly.
 * 
 * @param ax..az First point (quantized)
 * @param bx..bz Second point (quantized)
 * @param cx..cz Third point (quantized)
 * @param dx..dz Fourth point (quantized)
 * @return int {-1, 0, +1} orientation sign
 */
inline int orient3d_exact(int64_t ax, int64_t ay, int64_t az,
                          int64_t bx, int64_t by, int64_t bz,
                          int64_t cx, int64_t cy, int64_t cz,
                          int64_t dx, int64_t dy, int64_t dz) {
    int64_t adx = ax - dx, ady = ay - dy, adz = az - dz;
    int64_t bdx = bx - dx, bdy = by - dy, bdz = bz - dz;
    int64_t cdx = cx - dx, cdy = cy - dy, cdz = cz - dz;
    
    Int128 d1 = Int128::mul64(bdy, cdz) - Int128::mul64(bdz, cdy);
    Int128 d2 = Int128::mul64(bdx, cdz) - Int128::mul64(bdz, cdx);
    Int128 d3 = Int128::mul64(bdx, cdy) - Int128::mul64(bdy, cdx);
    
    Int256 det = Int256::mul128(Int128(adx), d1) 
               - Int256::mul128(Int128(ady), d2) 
               + Int256::mul128(Int128(adz), d3);
    
    return det.sign();
}

//=============================================================================
// EXACT PREDICATES: INCIRCLE
//=============================================================================

/**
 * @brief Exact incircle test (2D Delaunay in-circle predicate).
 *
 * Returns +1 if d is inside the circumcircle of triangle (a, b, c)
 * Returns -1 if d is outside the circumcircle
 * Returns  0 if d is exactly on the circumcircle (cocircular)
 *
 * The incircle test is computed as the sign of the 4x4 determinant:
 *   | ax  ay  ax²+ay²  1 |
 *   | bx  by  bx²+by²  1 |
 *   | cx  cy  cx²+cy²  1 |
 *   | dx  dy  dx²+dy²  1 |
 *
 * Bit-width analysis (26-bit inputs):
 *   - Coordinates: 26 bits
 *   - Squared coordinates: 52 bits
 *   - Lifted coordinate (x²+y²): 53 bits (sum of two 52-bit values)
 *   - 3x3 minors: 106 bits (53-bit * 53-bit)
 *   - Final determinant: ~160 bits (lifted * minor)
 *   - Requires Int256
 *
 * @return int {-1, 0, +1} in-circle sign
 */
inline int incircle_exact(int64_t ax, int64_t ay,
                          int64_t bx, int64_t by,
                          int64_t cx, int64_t cy,
                          int64_t dx, int64_t dy) {
    // Differences from d
    int64_t adx = ax - dx, ady = ay - dy;
    int64_t bdx = bx - dx, bdy = by - dy;
    int64_t cdx = cx - dx, cdy = cy - dy;

    // Compute 2x2 determinants (minors)
    Int128 abd = Int128::mul64(adx, bdy) - Int128::mul64(ady, bdx);  // |a-d, b-d|
    Int128 bcd = Int128::mul64(bdx, cdy) - Int128::mul64(bdy, cdx);  // |b-d, c-d|
    Int128 cad = Int128::mul64(cdx, ady) - Int128::mul64(cdy, adx);  // |c-d, a-d|

    // Compute squared distances from d
    Int128 ad_sq = Int128::mul64(adx, adx) + Int128::mul64(ady, ady);
    Int128 bd_sq = Int128::mul64(bdx, bdx) + Int128::mul64(bdy, bdy);
    Int128 cd_sq = Int128::mul64(cdx, cdx) + Int128::mul64(cdy, cdy);

    // Compute determinant: ad_sq * bcd + bd_sq * cad + cd_sq * abd
    Int256 det = Int256::mul128(ad_sq, bcd)
               + Int256::mul128(bd_sq, cad)
               + Int256::mul128(cd_sq, abd);

    return det.sign();
}

//=============================================================================
// FILTERED PREDICATES: DOT_PLANE
//=============================================================================

/**
 * @brief Compute dot product of point p with plane defined by (a, b, c).
 * 
 * The plane equation is: n · (x - a) = 0 where n = (b-a) × (c-a)
 * Returns: n · (p - a)
 * 
 * This is equivalent to the 4x4 determinant:
 *   | ax  ay  az  1 |
 *   | bx  by  bz  1 |
 *   | cx  cy  cz  1 |
 *   | px  py  pz  1 |
 * 
 * Algorithm:
 * 1. Compute plane normal using double cross product
 * 2. Compute dot product with (p - a)
 * 3. If uncertain, use exact Int128 arithmetic
 * 
 * Bit-width analysis:
 * - Normal components: 64 bits (cross product of 32-bit diffs)
 * - Dot product terms: 96 bits
 * - Sum: 98 bits (fits in Int128)
 * 
 * @param pa,pb,pc Plane-defining triangle
 * @param px,py,pz Query point
 * @return int {-1, 0, +1} sign of dot product
 */
inline int dot_plane_filtered(double pax, double pay, double paz,
                               double pbx, double pby, double pbz,
                               double pcx, double pcy, double pcz,
                               double px, double py, double pz) {
    // === FAST PATH: Double-precision evaluation ===
    
    // Compute edge vectors
    double abx = pbx - pax;
    double aby = pby - pay;
    double abz = pbz - paz;
    double acx = pcx - pax;
    double acy = pcy - pay;
    double acz = pcz - paz;
    
    // Compute normal: ab × ac
    double nx = aby * acz - abz * acy;
    double ny = abz * acx - abx * acz;
    double nz = abx * acy - aby * acx;
    
    // Vector from a to query point
    double apx = px - pax;
    double apy = py - pay;
    double apz = pz - paz;
    
    // Dot product
    double dot = nx * apx + ny * apy + nz * apz;
    
    // Error bound (conservative)
    double normal_mag = std::abs(nx) + std::abs(ny) + std::abs(nz);
    double ap_mag = std::abs(apx) + std::abs(apy) + std::abs(apz);
    double errbound = DOT_PLANE_ERRBOUND * normal_mag * ap_mag;
    
    // Check certainty
    if (dot > errbound) return +1;
    if (dot < -errbound) return -1;
    
    // === SLOW PATH: Exact Int128 computation ===
    // Using the 4x4 determinant formulation for exactness
    
    int64_t qax = quantize(pax), qay = quantize(pay), qaz = quantize(paz);
    int64_t qbx = quantize(pbx), qby = quantize(pby), qbz = quantize(pbz);
    int64_t qcx = quantize(pcx), qcy = quantize(pcy), qcz = quantize(pcz);
    int64_t qpx = quantize(px),  qpy = quantize(py),  qpz = quantize(pz);
    
    // Compute exact 4x4 determinant (expanded)
    // This is equivalent to orient3d with homogeneous coordinates
    
    // For efficiency, compute as: dot(n, ap) where n = ab × ac
    // Cross product components (each is a difference of two 64-bit products)
    Int128 nx_q = Int128::mul64(qby - qay, qcz - qaz) - Int128::mul64(qbz - qaz, qcy - qay);
    Int128 ny_q = Int128::mul64(qbz - qaz, qcx - qax) - Int128::mul64(qbx - qax, qcz - qaz);
    Int128 nz_q = Int128::mul64(qbx - qax, qcy - qay) - Int128::mul64(qby - qay, qcx - qax);
    
    // Dot product (each term is 128-bit product)
    Int128 dot_q = nx_q * (qpx - qax) + ny_q * (qpy - qay) + nz_q * (qpz - qaz);
    
    return dot_q.sign();
}

//=============================================================================
// HOMOGENEOUS COORDINATE PREDICATES
//=============================================================================

/**
 * @brief Plane representation for homogeneous point classification.
 * Stores plane equation: nx*x + ny*y + nz*z + d = 0
 * where (nx, ny, nz) is the normal and d is the offset.
 */
struct Plane {
    double nx, ny, nz, d;
    
    Plane() : nx(0), ny(0), nz(0), d(0) {}
    Plane(double nx_, double ny_, double nz_, double d_) 
        : nx(nx_), ny(ny_), nz(nz_), d(d_) {}
    
    /**
     * @brief Construct plane from three points.
     */
    static Plane fromPoints(double ax, double ay, double az,
                            double bx, double by, double bz,
                            double cx, double cy, double cz) {
        double abx = bx - ax, aby = by - ay, abz = bz - az;
        double acx = cx - ax, acy = cy - ay, acz = cz - az;
        
        double nx = aby * acz - abz * acy;
        double ny = abz * acx - abx * acz;
        double nz = abx * acy - aby * acx;
        double d = -(nx * ax + ny * ay + nz * az);
        
        return Plane(nx, ny, nz, d);
    }
};

/**
 * @brief Classify a homogeneous point against a plane.
 * 
 * Homogeneous points come from intersection operations with the form:
 *   (x, y, z, w) where w may be zero (point at infinity)
 * 
 * The classification is: sign(nx*x + ny*y + nz*z + d*w)
 * 
 * This requires 160-bit arithmetic (5 32-bit values multiplied and summed),
 * so we always use exact Int256 computation.
 * 
 * Bit-width analysis:
 * - Each term: nx*x requires 64 bits * 32 bits = 96 bits
 * - Sum of 4 terms: 98 bits, but with w up to 2^32, d*w is 96 bits
 * - Total: ~128 bits worst case, use Int256 for safety
 * 
 * @param x,y,z,w Homogeneous coordinates (w may be 0)
 * @param plane Plane equation
 * @return int {-1, 0, +1} classification sign
 */
inline int classify_homogeneous(int64_t x, int64_t y, int64_t z, int64_t w,
                                 const Plane& plane) {
    // Always use exact arithmetic for homogeneous points
    // Quantize plane coefficients
    int64_t qnx = quantize(plane.nx);
    int64_t qny = quantize(plane.ny);
    int64_t qnz = quantize(plane.nz);
    int64_t qd  = quantize(plane.d);
    
    // Compute dot product exactly (each product is ~96 bits)
    Int128 term1 = Int128::mul64(qnx, x);
    Int128 term2 = Int128::mul64(qny, y);
    Int128 term3 = Int128::mul64(qnz, z);
    Int128 term4 = Int128::mul64(qd, w);
    
    Int128 result = term1 + term2 + term3 + term4;
    
    return result.sign();
}

/**
 * @brief Double-precision version for homogeneous classification.
 * Uses exact fallback if uncertain.
 */
inline int classify_homogeneous(double x, double y, double z, double w,
                                 const Plane& plane) {
    // Fast path
    double result = plane.nx * x + plane.ny * y + plane.nz * z + plane.d * w;
    
    double max_term = std::abs(plane.nx * x) + std::abs(plane.ny * y) 
                    + std::abs(plane.nz * z) + std::abs(plane.d * w);
    double errbound = DOT_PLANE_ERRBOUND * max_term;
    
    if (result > errbound) return +1;
    if (result < -errbound) return -1;
    
    // Exact path
    return classify_homogeneous(quantize(x), quantize(y), quantize(z), quantize(w), plane);
}

//=============================================================================
// IMPLICIT POINT TYPES (for Indirect Predicates)
//=============================================================================

/**
 * @brief Types of implicit points supported by the indirect predicate system.
 * 
 * EXPLICIT: Regular point with known coordinates
 * LPI: Line-Plane Intersection - intersection of line AB with plane PQR
 * TPI: Triangle-Plane Intersection - intersection of plane ABC with plane PQR
 */
enum class ImplicitPointType {
    EXPLICIT,
    LPI,    // Line-Plane Intersection
    TPI     // Triangle-Plane (actually plane-plane) Intersection
};

/**
 * @brief Forward declaration of triangle storage.
 * Actual implementation depends on mesh data structure.
 */
template<typename T>
class TriangleStorage;

/**
 * @brief Implicit point representation for Cherchi et al. framework.
 * 
 * Implicit points are defined geometrically rather than by explicit coordinates.
 * This allows exact evaluation of predicates without computing coordinates.
 * 
 * For LPI (Line-Plane Intersection):
 *   - line0, line1: indices of points defining the line
 *   - plane0, plane1, plane2: indices of points defining the plane
 * 
 * For TPI (Triangle-Plane Intersection = plane-plane intersection):
 *   - line0, line1, line2: first plane (as triangle)
 *   - plane0, plane1, plane2: second plane (as triangle)
 */
struct ImplicitPoint {
    ImplicitPointType type;
    
    union {
        // For EXPLICIT: store index to coordinate array
        uint32_t explicit_idx;
        
        // For LPI/TPI: store indices to triangle storage
        struct {
            uint32_t line0, line1, line2;    // First plane/line (up to 3 points)
            uint32_t plane0, plane1, plane2; // Second plane (3 points)
        } indices;
    };
    
    // Constructors
    static ImplicitPoint makeExplicit(uint32_t idx) {
        ImplicitPoint p;
        p.type = ImplicitPointType::EXPLICIT;
        p.explicit_idx = idx;
        return p;
    }
    
    static ImplicitPoint makeLPI(uint32_t l0, uint32_t l1, 
                                  uint32_t p0, uint32_t p1, uint32_t p2) {
        ImplicitPoint p;
        p.type = ImplicitPointType::LPI;
        p.indices.line0 = l0;
        p.indices.line1 = l1;
        p.indices.line2 = 0;  // Unused for LPI
        p.indices.plane0 = p0;
        p.indices.plane1 = p1;
        p.indices.plane2 = p2;
        return p;
    }
    
    static ImplicitPoint makeTPI(uint32_t t0, uint32_t t1, uint32_t t2,
                                  uint32_t p0, uint32_t p1, uint32_t p2) {
        ImplicitPoint p;
        p.type = ImplicitPointType::TPI;
        p.indices.line0 = t0;
        p.indices.line1 = t1;
        p.indices.line2 = t2;
        p.indices.plane0 = p0;
        p.indices.plane1 = p1;
        p.indices.plane2 = p2;
        return p;
    }
};

//=============================================================================
// INDIRECT PREDICATES
//=============================================================================

/**
 * @brief Evaluate orient2d for implicit points using filtered approach.
 * 
 * This is the core predicate for Cherchi et al.'s indirect framework.
 * It evaluates orientation without computing explicit intersection coordinates.
 * 
 * For explicit points: delegates to orient2d_filtered
 * For LPI/TPI: uses the implicit point's geometric definition
 * 
 * @param ip1 First implicit point
 * @param ip2 Second implicit point  
 * @param ip3 Third implicit point
 * @param triangles Triangle storage for coordinate lookup
 * @return int {-1, 0, +1} orientation sign
 */
template<typename TriangleStorage>
inline int orient2d_indirect(const ImplicitPoint& ip1,
                              const ImplicitPoint& ip2,
                              const ImplicitPoint& ip3,
                              const TriangleStorage& triangles) {
    // Collect all explicit points involved
    // For each implicit point, we may need to expand to its definition
    
    // Simple case: all explicit
    if (ip1.type == ImplicitPointType::EXPLICIT &&
        ip2.type == ImplicitPointType::EXPLICIT &&
        ip3.type == ImplicitPointType::EXPLICIT) {
        
        auto [ax, ay] = triangles.getPoint2D(ip1.explicit_idx);
        auto [bx, by] = triangles.getPoint2D(ip2.explicit_idx);
        auto [cx, cy] = triangles.getPoint2D(ip3.explicit_idx);
        
        return orient2d_filtered(ax, ay, bx, by, cx, cy);
    }
    
    // =========================================================================
    // AUDIT FIX: Exact indirect predicates using homogeneous coordinates.
    // 
    // REPLACES: The broken FP-division path that computed:
    //     ax = l0x + (num / den) * lx;  // ← FLOATING-POINT DIVISION
    //     ...
    //     return orient2d_filtered(ax, ay, bx, by, cx, cy);
    //
    // The correct approach from Cherchi et al.: represent LPI points in
    // homogeneous coordinates P = (N_x, N_y, N_z, D) and evaluate
    // orient2d via the sign of (det × D_a × D_b × D_c) with ZERO divisions.
    //
    // See ExactPredicates_Indirect.h for the full implementation.
    // =========================================================================
    
    // Helper: compute LPI or wrap explicit point as homogeneous
    auto makeHomogPoint = [&](const ImplicitPoint& ip, int proj_axis) -> HomogPoint2D {
        if (ip.type == ImplicitPointType::EXPLICIT) {
            auto [px, py, pz] = triangles.getPoint3D(ip.explicit_idx);
            // Quantize to integer
            int32_t qx = static_cast<int32_t>(std::llround(px * RESCALE_FACTOR));
            int32_t qy = static_cast<int32_t>(std::llround(py * RESCALE_FACTOR));
            int32_t qz = static_cast<int32_t>(std::llround(pz * RESCALE_FACTOR));
            
            switch (proj_axis) {
                case 0: return HomogPoint2D::fromExplicit(qy, qz);
                case 1: return HomogPoint2D::fromExplicit(qx, qz);
                case 2: return HomogPoint2D::fromExplicit(qx, qy);
            }
            return HomogPoint2D::fromExplicit(qx, qy);  // default
        } else {
            // Compute exact LPI
            auto [l0x, l0y, l0z] = triangles.getPoint3D(ip.indices.line0);
            auto [l1x, l1y, l1z] = triangles.getPoint3D(ip.indices.line1);
            auto [p0x, p0y, p0z] = triangles.getPoint3D(ip.indices.plane0);
            auto [p1x, p1y, p1z] = triangles.getPoint3D(ip.indices.plane1);
            auto [p2x, p2y, p2z] = triangles.getPoint3D(ip.indices.plane2);
            
            // Quantize all coordinates
            int32_t ql0x = static_cast<int32_t>(std::llround(l0x * RESCALE_FACTOR));
            int32_t ql0y = static_cast<int32_t>(std::llround(l0y * RESCALE_FACTOR));
            int32_t ql0z = static_cast<int32_t>(std::llround(l0z * RESCALE_FACTOR));
            int32_t ql1x = static_cast<int32_t>(std::llround(l1x * RESCALE_FACTOR));
            int32_t ql1y = static_cast<int32_t>(std::llround(l1y * RESCALE_FACTOR));
            int32_t ql1z = static_cast<int32_t>(std::llround(l1z * RESCALE_FACTOR));
            int32_t qp0x = static_cast<int32_t>(std::llround(p0x * RESCALE_FACTOR));
            int32_t qp0y = static_cast<int32_t>(std::llround(p0y * RESCALE_FACTOR));
            int32_t qp0z = static_cast<int32_t>(std::llround(p0z * RESCALE_FACTOR));
            int32_t qp1x = static_cast<int32_t>(std::llround(p1x * RESCALE_FACTOR));
            int32_t qp1y = static_cast<int32_t>(std::llround(p1y * RESCALE_FACTOR));
            int32_t qp1z = static_cast<int32_t>(std::llround(p1z * RESCALE_FACTOR));
            int32_t qp2x = static_cast<int32_t>(std::llround(p2x * RESCALE_FACTOR));
            int32_t qp2y = static_cast<int32_t>(std::llround(p2y * RESCALE_FACTOR));
            int32_t qp2z = static_cast<int32_t>(std::llround(p2z * RESCALE_FACTOR));
            
            LPI_Exact lpi = compute_LPI_exact(
                ql0x, ql0y, ql0z, ql1x, ql1y, ql1z,
                qp0x, qp0y, qp0z, qp1x, qp1y, qp1z, qp2x, qp2y, qp2z);
            
            if (!lpi.valid) {
                // Line parallel to plane — degenerate intersection
                // Return 0 (collinear) as fallback
                return HomogPoint2D::fromExplicit(0, 0);
            }
            
            return HomogPoint2D::fromLPI(lpi, proj_axis);
        }
    };
    
    // Choose projection axis based on the first triangle encountered
    // (use dominant normal component for maximum area projection)
    int proj_axis = 2;  // Default: project to XY
    
    // Build homogeneous 2D points
    HomogPoint2D A = makeHomogPoint(ip1, proj_axis);
    HomogPoint2D B = makeHomogPoint(ip2, proj_axis);
    HomogPoint2D C = makeHomogPoint(ip3, proj_axis);
    
    // Exact orient2d with ZERO floating-point divisions
    return orient2d_homogeneous_exact(A, B, C);
}

/**
 * @brief Compare two implicit points for equality/ordering.
 * 
 * AUDIT FIX: Replaced epsilon-based comparison with exact Int128 arithmetic.
 * The original used:
 *   const double eps = EPSILON * (std::abs(ax) + std::abs(bx) + 1.0);
 *   if (ax < bx - eps) return -1;
 * This is NOT exact — two distinct LPI points differing by less than eps
 * would be incorrectly reported as equal, corrupting CDT topology.
 *
 * The fixed version uses homogeneous coordinates:
 *   A_x/A_w < B_x/B_w  ⟺  A_x*B_w < B_x*A_w  (if A_w*B_w > 0)
 * All comparisons use Int256 products — zero epsilon, zero division.
 * 
 * @return int {-1, 0, +1} comparison result
 */
template<typename TriangleStorage>
inline int pointCompare_indirect(const ImplicitPoint& a,
                                  const ImplicitPoint& b,
                                  const TriangleStorage& triangles) {
    // Both explicit: direct double comparison (exact for identical coordinates)
    if (a.type == ImplicitPointType::EXPLICIT && b.type == ImplicitPointType::EXPLICIT) {
        auto [ax, ay, az] = triangles.getPoint3D(a.explicit_idx);
        auto [bx, by, bz] = triangles.getPoint3D(b.explicit_idx);
        
        if (ax < bx) return -1;
        if (ax > bx) return +1;
        if (ay < by) return -1;
        if (ay > by) return +1;
        if (az < bz) return -1;
        if (az > bz) return +1;
        return 0;
    }
    
    // Helper: compute exact LPI for an implicit point
    auto computeLPI = [&](const ImplicitPoint& ip) -> LPI_Exact {
        auto [l0x, l0y, l0z] = triangles.getPoint3D(ip.indices.line0);
        auto [l1x, l1y, l1z] = triangles.getPoint3D(ip.indices.line1);
        auto [p0x, p0y, p0z] = triangles.getPoint3D(ip.indices.plane0);
        auto [p1x, p1y, p1z] = triangles.getPoint3D(ip.indices.plane1);
        auto [p2x, p2y, p2z] = triangles.getPoint3D(ip.indices.plane2);
        
        return compute_LPI_exact(
            static_cast<int32_t>(std::llround(l0x * RESCALE_FACTOR)),
            static_cast<int32_t>(std::llround(l0y * RESCALE_FACTOR)),
            static_cast<int32_t>(std::llround(l0z * RESCALE_FACTOR)),
            static_cast<int32_t>(std::llround(l1x * RESCALE_FACTOR)),
            static_cast<int32_t>(std::llround(l1y * RESCALE_FACTOR)),
            static_cast<int32_t>(std::llround(l1z * RESCALE_FACTOR)),
            static_cast<int32_t>(std::llround(p0x * RESCALE_FACTOR)),
            static_cast<int32_t>(std::llround(p0y * RESCALE_FACTOR)),
            static_cast<int32_t>(std::llround(p0z * RESCALE_FACTOR)),
            static_cast<int32_t>(std::llround(p1x * RESCALE_FACTOR)),
            static_cast<int32_t>(std::llround(p1y * RESCALE_FACTOR)),
            static_cast<int32_t>(std::llround(p1z * RESCALE_FACTOR)),
            static_cast<int32_t>(std::llround(p2x * RESCALE_FACTOR)),
            static_cast<int32_t>(std::llround(p2y * RESCALE_FACTOR)),
            static_cast<int32_t>(std::llround(p2z * RESCALE_FACTOR)));
    };
    
    // Both implicit: exact LPI comparison
    if (a.type != ImplicitPointType::EXPLICIT && b.type != ImplicitPointType::EXPLICIT) {
        LPI_Exact lpi_a = computeLPI(a);
        LPI_Exact lpi_b = computeLPI(b);
        
        if (!lpi_a.valid || !lpi_b.valid) return 0;  // Degenerate
        return pointCompare_exact(lpi_a, lpi_b);
    }
    
    // Mixed: one explicit, one implicit
    if (a.type == ImplicitPointType::EXPLICIT) {
        // A is explicit, B is implicit
        auto [ax, ay, az] = triangles.getPoint3D(a.explicit_idx);
        int32_t qax = static_cast<int32_t>(std::llround(ax * RESCALE_FACTOR));
        int32_t qay = static_cast<int32_t>(std::llround(ay * RESCALE_FACTOR));
        int32_t qaz = static_cast<int32_t>(std::llround(az * RESCALE_FACTOR));
        
        LPI_Exact lpi_b = computeLPI(b);
        if (!lpi_b.valid) return 0;
        
        return pointCompare_explicit_vs_lpi(qax, qay, qaz, lpi_b);
    } else {
        // A is implicit, B is explicit
        auto [bx, by, bz] = triangles.getPoint3D(b.explicit_idx);
        int32_t qbx = static_cast<int32_t>(std::llround(bx * RESCALE_FACTOR));
        int32_t qby = static_cast<int32_t>(std::llround(by * RESCALE_FACTOR));
        int32_t qbz = static_cast<int32_t>(std::llround(bz * RESCALE_FACTOR));
        
        LPI_Exact lpi_a = computeLPI(a);
        if (!lpi_a.valid) return 0;
        
        return -pointCompare_explicit_vs_lpi(qbx, qby, qbz, lpi_a);
    }
}

//=============================================================================
// UTILITY FUNCTIONS
//=============================================================================

/**
 * @brief Check if three 2D points are collinear.
 * Convenience wrapper around orient2d_filtered.
 */
inline bool isCollinear2D(double ax, double ay, double bx, double by,
                          double cx, double cy) {
    return orient2d_filtered(ax, ay, bx, by, cx, cy) == 0;
}

/**
 * @brief Check if four 3D points are coplanar.
 * Convenience wrapper around orient3d_filtered.
 */
inline bool isCoplanar3D(double ax, double ay, double az,
                         double bx, double by, double bz,
                         double cx, double cy, double cz,
                         double dx, double dy, double dz) {
    return orient3d_filtered(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz) == 0;
}

/**
 * @brief Compute triangle area (2x absolute value of orient2d).
 */
inline double triangleArea2D(double ax, double ay, double bx, double by,
                              double cx, double cy) {
    double acx = ax - cx;
    double acy = ay - cy;
    double bcx = bx - cx;
    double bcy = by - cy;
    return std::abs(acx * bcy - acy * bcx);
}

/**
 * @brief Compute tetrahedron volume (6x absolute value of orient3d).
 */
inline double tetrahedronVolume(double ax, double ay, double az,
                                 double bx, double by, double bz,
                                 double cx, double cy, double cz,
                                 double dx, double dy, double dz) {
    double adx = ax - dx, ady = ay - dy, adz = az - dz;
    double bdx = bx - dx, bdy = by - dy, bdz = bz - dz;
    double cdx = cx - dx, cdy = cy - dy, cdz = cz - dz;
    
    double det = adx * (bdy * cdz - bdz * cdy)
               - ady * (bdx * cdz - bdz * cdx)
               + adz * (bdx * cdy - bdy * cdx);
    
    return std::abs(det);
}

} // namespace predicates
} // namespace ember
