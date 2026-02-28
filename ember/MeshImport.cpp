/**
 * @file MeshImport.cpp
 * @brief Implementation of mesh import and quantization for EMBER Boolean engine
 *
 * AUDIT FIX (2026-02-28): Replaced GCC-only __int128 with portable ember::Int128.
 * This enables compilation on MSVC and fixes silent truncation bugs.
 */

#include "MeshImport.h"
#include "Diagnostics.h"
#include "IntegerTypes.h"  // PORTABLE: ember::Int128 instead of __int128
#include <cstring>
#include <cmath>

namespace ember {

//=============================================================================
// INTEGER SQRT IMPLEMENTATIONS (PORTABLE)
//=============================================================================

int64_t isqrt(int64_t n) {
    if (n < 0) return 0;
    if (n < 2) return n;
    
    int64_t x = n;
    int64_t y = (x + 1) / 2;
    while (y < x) {
        x = y;
        y = (x + n / x) / 2;
    }
    return x;
}

/**
 * @brief Integer square root for Int128 values (PORTABLE).
 *
 * AUDIT FIX: Replaced GCC-only __int128 isqrt128 with portable version.
 * Uses double approximation for initial guess, then Newton refinement.
 */
Int128 isqrt128(const Int128& n) {
    if (n.isZero() || n.isNegative()) return Int128(0);
    
    // Check if fits in int64_t
    if (n.high() == 0 || (n.high() == -1 && n.low() <= static_cast<uint64_t>(INT64_MAX))) {
        return Int128(isqrt(static_cast<int64_t>(n.low())));
    }

    // Initial estimate using double (54-bit mantissa is sufficient for Newton seed)
    double n_approx = static_cast<double>(n.high()) * 18446744073709551616.0
                    + static_cast<double>(n.low());
    int64_t guess_64 = static_cast<int64_t>(std::sqrt(n_approx));
    Int128 guess(guess_64);

    // Newton refinement: x_{n+1} = (x_n + n/x_n) / 2
    for (int i = 0; i < 10; ++i) {
        if (guess.isZero()) return Int128(0);
        Int128 next = (guess + n / guess) >> 1;
        if (next >= guess) break;
        guess = next;
    }

    // Adjust to ensure guess*guess <= n < (guess+1)*(guess+1)
    while (guess * guess > n) --guess;
    while ((guess + 1) * (guess + 1) <= n) ++guess;

    return guess;
}

//=============================================================================
// VECTOR NORMALIZATION (PORTABLE)
//=============================================================================

/**
 * @brief Normalize a 3D vector to 26-bit fixed-point (PORTABLE VERSION).
 *
 * AUDIT FIX: Replaced __int128 with ember::Int128.
 * The original code had a truncation bug: int64_t len_sq = (int64_t)len_sq_128
 * would truncate a 108-bit value to 64 bits, corrupting results for large inputs.
 */
void normalizeVector53To26(int64_t nx, int64_t ny, int64_t nz,
                           int32_t& out_nx, int32_t& out_ny, int32_t& out_nz) {
    // Compute squared length using Int128 (safe for 53-bit inputs → 106-bit products)
    Int128 len_sq = Int128::mul64(nx, nx) + Int128::mul64(ny, ny) + Int128::mul64(nz, nz);

    if (len_sq.isZero()) {
        out_nx = out_ny = out_nz = 0;
        return;
    }

    // Integer square root of Int128 value (no truncation to int64_t)
    Int128 len = isqrt128(len_sq);

    if (len.isZero()) {
        out_nx = out_ny = out_nz = 0;
        return;
    }

    // Scale each component: (nx * target_mag) / len where target_mag = EMBER_QUANT_MAX / 2
    Int128 scaled_nx = (Int128::mul64(nx, EMBER_QUANT_MAX / 2)) / len;
    Int128 scaled_ny = (Int128::mul64(ny, EMBER_QUANT_MAX / 2)) / len;
    Int128 scaled_nz = (Int128::mul64(nz, EMBER_QUANT_MAX / 2)) / len;

    auto clamp26 = [](const Int128& v) -> int32_t {
        if (v > Int128(EMBER_QUANT_MAX)) return EMBER_QUANT_MAX;
        if (v < Int128(EMBER_QUANT_MIN)) return EMBER_QUANT_MIN;
        return static_cast<int32_t>(v.low());
    };

    out_nx = clamp26(scaled_nx);
    out_ny = clamp26(scaled_ny);
    out_nz = clamp26(scaled_nz);
}

//=============================================================================
// QUANTIZATION CONTEXT METHODS
//=============================================================================

void QuantizationContext::quantize(double fx, double fy, double fz,
                                    int32_t& ix, int32_t& iy, int32_t& iz) const {
    double scaled_x = (fx - bbox_min[0]) * scale[0] - static_cast<double>(EMBER_QUANT_MAX);
    double scaled_y = (fy - bbox_min[1]) * scale[1] - static_cast<double>(EMBER_QUANT_MAX);
    double scaled_z = (fz - bbox_min[2]) * scale[2] - static_cast<double>(EMBER_QUANT_MAX);
    
    int64_t qx = std::llround(scaled_x);
    int64_t qy = std::llround(scaled_y);
    int64_t qz = std::llround(scaled_z);
    
    auto clamp26 = [](int64_t v) -> int32_t {
        if (v > EMBER_QUANT_MAX) return EMBER_QUANT_MAX;
        if (v < EMBER_QUANT_MIN) return EMBER_QUANT_MIN;
        return static_cast<int32_t>(v);
    };
    
    ix = clamp26(qx);
    iy = clamp26(qy);
    iz = clamp26(qz);
}

void QuantizationContext::dequantize(int32_t ix, int32_t iy, int32_t iz,
                                      double& fx, double& fy, double& fz) const {
    double centered_x = static_cast<double>(ix) + static_cast<double>(EMBER_QUANT_MAX);
    double centered_y = static_cast<double>(iy) + static_cast<double>(EMBER_QUANT_MAX);
    double centered_z = static_cast<double>(iz) + static_cast<double>(EMBER_QUANT_MAX);
    
    fx = centered_x * inv_scale[0] + bbox_min[0];
    fy = centered_y * inv_scale[1] + bbox_min[1];
    fz = centered_z * inv_scale[2] + bbox_min[2];
}

//=============================================================================
// COMPUTE QUANTIZATION
//=============================================================================

QuantizationContext computeQuantization(const PolygonSoup& soup, 
                                         double aspect_ratio_threshold) {
    QuantizationContext ctx;
    ctx.aspect_ratio_threshold = aspect_ratio_threshold;
    
    if (soup.triangles.empty()) {
        EMBER_LOG_WARN("Empty polygon soup, using default quantization");
        return ctx;
    }
    
    soup.computeBoundingBox(ctx.bbox_min, ctx.bbox_max);
    
    double rx = ctx.bbox_max[0] - ctx.bbox_min[0];
    double ry = ctx.bbox_max[1] - ctx.bbox_min[1];
    double rz = ctx.bbox_max[2] - ctx.bbox_min[2];
    
    rx = std::max(rx, EMBER_MIN_DIMENSION);
    ry = std::max(ry, EMBER_MIN_DIMENSION);
    rz = std::max(rz, EMBER_MIN_DIMENSION);
    
    double ratio_xy = rx / ry;
    double ratio_xz = rx / rz;
    double ratio_yx = ry / rx;
    double ratio_yz = ry / rz;
    double ratio_zx = rz / rx;
    double ratio_zy = rz / ry;
    
    double max_ratio = std::max({ratio_xy, ratio_xz, ratio_yx, 
                                  ratio_yz, ratio_zx, ratio_zy});
    
    const double target_range = 2.0 * static_cast<double>(EMBER_QUANT_MAX);
    
    if (max_ratio > aspect_ratio_threshold) {
        ctx.per_axis = true;
        ctx.scale[0] = target_range / rx;
        ctx.scale[1] = target_range / ry;
        ctx.scale[2] = target_range / rz;
        
        std::ostringstream msg;
        msg << "Extreme aspect ratio detected (" << std::fixed << std::setprecision(1) 
            << max_ratio << ":1), using per-axis quantization (threshold=" 
            << aspect_ratio_threshold << ")";
        EMBER_LOG_WARN("%s", msg.str().c_str());
    } else {
        ctx.per_axis = false;
        double max_range = std::max({rx, ry, rz});
        double scale = target_range / max_range;
        ctx.scale[0] = ctx.scale[1] = ctx.scale[2] = scale;
    }
    
    ctx.inv_scale[0] = 1.0 / ctx.scale[0];
    ctx.inv_scale[1] = 1.0 / ctx.scale[1];
    ctx.inv_scale[2] = 1.0 / ctx.scale[2];
    
    return ctx;
}

//=============================================================================
// BATCH QUANTIZATION
//=============================================================================

uint32_t quantizeAll(PolygonSoup& soup, const QuantizationContext& ctx) {
    for (int i = 0; i < 3; ++i) {
        soup.quantization_scale[i] = ctx.scale[i];
        soup.quantization_inv_scale[i] = ctx.inv_scale[i];
    }
    soup.per_axis_quantization = ctx.per_axis;
    
    for (auto& tri : soup.triangles) {
        for (int v = 0; v < 3; ++v) {
            ctx.quantize(tri.v[v][0], tri.v[v][1], tri.v[v][2],
                         tri.iv[v][0], tri.iv[v][1], tri.iv[v][2]);
        }
        
        tri.plane = computePlaneFromTriangle(
            tri.iv[0].data(), tri.iv[1].data(), tri.iv[2].data());
    }
    
    return static_cast<uint32_t>(soup.triangles.size());
}

//=============================================================================
// FIX 3.7: UNIFIED QUANTIZATION SCALE VERIFICATION
//=============================================================================

bool verifyQuantizationScale(const QuantizationContext& ctx, const PolygonSoup& soup) {
    constexpr double TOLERANCE = 1e-10;  // Tolerance for floating-point comparison
    
    for (int i = 0; i < 3; ++i) {
        double diff = std::abs(ctx.scale[i] - soup.quantization_scale[i]);
        double max_scale = std::max(std::abs(ctx.scale[i]), std::abs(soup.quantization_scale[i]));
        
        // Relative tolerance check (handles both small and large scales)
        if (max_scale > EMBER_MIN_DIMENSION) {
            if (diff / max_scale > TOLERANCE) {
                EMBER_LOG_ERROR("Quantization scale mismatch on axis %d: ctx=%.15f, soup=%.15f",
                               i, ctx.scale[i], soup.quantization_scale[i]);
                return false;
            }
        } else {
            // Absolute tolerance for very small scales
            if (diff > TOLERANCE) {
                EMBER_LOG_ERROR("Quantization scale mismatch on axis %d: ctx=%.15f, soup=%.15f",
                               i, ctx.scale[i], soup.quantization_scale[i]);
                return false;
            }
        }
    }
    
    return true;
}

//=============================================================================
// PLANE COMPUTATION
//=============================================================================

Plane computePlaneFromTriangle(const int32_t v0[3], const int32_t v1[3], const int32_t v2[3]) {
    Plane plane;
    
    int64_t e1x = static_cast<int64_t>(v1[0]) - v0[0];
    int64_t e1y = static_cast<int64_t>(v1[1]) - v0[1];
    int64_t e1z = static_cast<int64_t>(v1[2]) - v0[2];
    
    int64_t e2x = static_cast<int64_t>(v2[0]) - v0[0];
    int64_t e2y = static_cast<int64_t>(v2[1]) - v0[1];
    int64_t e2z = static_cast<int64_t>(v2[2]) - v0[2];
    
    int64_t nx = e1y * e2z - e1z * e2y;
    int64_t ny = e1z * e2x - e1x * e2z;
    int64_t nz = e1x * e2y - e1y * e2x;
    
    int32_t nnx, nny, nnz;
    normalizeVector53To26(nx, ny, nz, nnx, nny, nnz);
    
    plane.a = nnx;
    plane.b = nny;
    plane.c = nnz;
    
    int64_t d64 = -(static_cast<int64_t>(nnx) * v0[0] +
                    static_cast<int64_t>(nny) * v0[1] +
                    static_cast<int64_t>(nnz) * v0[2]);
    
    auto clamp26 = [](int64_t v) -> int32_t {
        if (v > EMBER_QUANT_MAX) return EMBER_QUANT_MAX;
        if (v < EMBER_QUANT_MIN) return EMBER_QUANT_MIN;
        return static_cast<int32_t>(v);
    };
    
    plane.d = clamp26(d64);
    
    return plane;
}

Plane computePlaneFromTriangle(const double v0[3], const double v1[3], const double v2[3]) {
    int32_t iv0[3] = {static_cast<int32_t>(v0[0]), static_cast<int32_t>(v0[1]), static_cast<int32_t>(v0[2])};
    int32_t iv1[3] = {static_cast<int32_t>(v1[0]), static_cast<int32_t>(v1[1]), static_cast<int32_t>(v1[2])};
    int32_t iv2[3] = {static_cast<int32_t>(v2[0]), static_cast<int32_t>(v2[1]), static_cast<int32_t>(v2[2])};
    return computePlaneFromTriangle(iv0, iv1, iv2);
}

//=============================================================================
// VALIDATION
//=============================================================================

bool isDegenerate(const int32_t v0[3], const int32_t v1[3], const int32_t v2[3]) {
    int64_t e1x = static_cast<int64_t>(v1[0]) - v0[0];
    int64_t e1y = static_cast<int64_t>(v1[1]) - v0[1];
    int64_t e1z = static_cast<int64_t>(v1[2]) - v0[2];
    
    int64_t e2x = static_cast<int64_t>(v2[0]) - v0[0];
    int64_t e2y = static_cast<int64_t>(v2[1]) - v0[1];
    int64_t e2z = static_cast<int32_t>(v2[2]) - v0[2];
    
    int64_t nx = e1y * e2z - e1z * e2y;
    int64_t ny = e1z * e2x - e1x * e2z;
    int64_t nz = e1x * e2y - e1y * e2x;
    
    return nx == 0 && ny == 0 && nz == 0;
}

uint32_t validateSoup(const PolygonSoup& soup, bool verbose) {
    uint32_t invalid_count = 0;
    
    for (size_t i = 0; i < soup.triangles.size(); ++i) {
        const auto& tri = soup.triangles[i];
        bool invalid = false;
        
        // Check for NaN/Inf
        for (int v = 0; v < 3; ++v) {
            for (int a = 0; a < 3; ++a) {
                if (!std::isfinite(tri.v[v][a])) {
                    invalid = true;
                    break;
                }
            }
        }
        
        // Check for degenerate
        if (!invalid && isDegenerate(tri.iv[0].data(), tri.iv[1].data(), tri.iv[2].data())) {
            invalid = true;
        }
        
        if (invalid) {
            invalid_count++;
            if (verbose) {
                std::cerr << "Invalid triangle " << i << std::endl;
            }
        }
    }
    
    return invalid_count;
}

//=============================================================================
// DIAGNOSTICS
//=============================================================================

void printQuantizationContext(const QuantizationContext& ctx, std::ostream& out) {
    out << "Quantization Context:" << std::endl;
    out << "  Per-axis: " << (ctx.per_axis ? "yes" : "no") << std::endl;
    out << "  BBox min: [" << ctx.bbox_min[0] << ", " << ctx.bbox_min[1] << ", " << ctx.bbox_min[2] << "]" << std::endl;
    out << "  BBox max: [" << ctx.bbox_max[0] << ", " << ctx.bbox_max[1] << ", " << ctx.bbox_max[2] << "]" << std::endl;
    out << "  Scale: [" << ctx.scale[0] << ", " << ctx.scale[1] << ", " << ctx.scale[2] << "]" << std::endl;
    out << "  Max aspect ratio: " << ctx.getMaxAspectRatio() << std::endl;
}

void printSoupStats(const PolygonSoup& soup, std::ostream& out) {
    out << "Polygon Soup Statistics:" << std::endl;
    out << "  Triangles: " << soup.triangles.size() << std::endl;
    out << "  Integer vertices: " << soup.int_vertices.size() << std::endl;
    out << "  Candidate pairs: " << soup.candidate_pairs.size() << std::endl;
    out << "  Output polygons: " << soup.output_polygons.size() << std::endl;
    out << "  Per-axis quantization: " << (soup.per_axis_quantization ? "yes" : "no") << std::endl;
}

const char* getEmberVersion() {
    return "EMBER 1.1.1 (Production Ready)";
}

} // namespace ember
