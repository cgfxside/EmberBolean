/**
 * @file CutterMeshGenerator.cpp
 * @brief Implementation of noise mesh generator for EMBER Boolean engine
 */

#include "CutterMeshGenerator.h"
#include "Diagnostics.h"
#include <cmath>
#include <random>
#include <algorithm>

namespace ember {

//=============================================================================
// SIMPLEX NOISE IMPLEMENTATION
//=============================================================================

// Permutation table for simplex noise
static const uint8_t SIMPLEX_PERM[256] = {
    151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225,
    140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148,
    247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32,
    57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175,
    74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122,
    60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54,
    65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169,
    200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64,
    52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212,
    207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213,
    119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
    129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104,
    218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241,
    81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157,
    184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93,
    222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
};

// Gradient table for 3D
static const int8_t GRAD3[16][3] = {
    {1, 1, 0}, {-1, 1, 0}, {1, -1, 0}, {-1, -1, 0},
    {1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1},
    {0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1},
    {1, 1, 0}, {-1, 1, 0}, {0, -1, 1}, {0, -1, -1}
};

// Fast floor function
static inline int fastFloor(float x) {
    int xi = static_cast<int>(x);
    return x < xi ? xi - 1 : xi;
}

// Dot product with gradient
static inline float dotGrad(int gi, float x, float y, float z) {
    return GRAD3[gi & 15][0] * x + GRAD3[gi & 15][1] * y + GRAD3[gi & 15][2] * z;
}

float simplexNoise3D(float x, float y, float z, uint32_t seed) {
    // Skewing and unskewing factors for 3D
    static const float F3 = 1.0f / 3.0f;
    static const float G3 = 1.0f / 6.0f;
    
    // Skew the input space to determine which simplex cell we're in
    float s = (x + y + z) * F3;
    int i = fastFloor(x + s);
    int j = fastFloor(y + s);
    int k = fastFloor(z + s);
    
    float t = (i + j + k) * G3;
    float X0 = i - t;
    float Y0 = j - t;
    float Z0 = k - t;
    float x0 = x - X0;
    float y0 = y - Y0;
    float z0 = z - Z0;
    
    // Determine which simplex we are in
    int i1, j1, k1;
    int i2, j2, k2;
    
    if (x0 >= y0) {
        if (y0 >= z0) {
            i1 = 1; j1 = 0; k1 = 0;
            i2 = 1; j2 = 1; k2 = 0;
        } else if (x0 >= z0) {
            i1 = 1; j1 = 0; k1 = 0;
            i2 = 1; j2 = 0; k2 = 1;
        } else {
            i1 = 0; j1 = 0; k1 = 1;
            i2 = 1; j2 = 0; k2 = 1;
        }
    } else {
        if (y0 < z0) {
            i1 = 0; j1 = 0; k1 = 1;
            i2 = 0; j2 = 1; k2 = 1;
        } else if (x0 < z0) {
            i1 = 0; j1 = 1; k1 = 0;
            i2 = 0; j2 = 1; k2 = 1;
        } else {
            i1 = 0; j1 = 1; k1 = 0;
            i2 = 1; j2 = 1; k2 = 0;
        }
    }
    
    float x1 = x0 - i1 + G3;
    float y1 = y0 - j1 + G3;
    float z1 = z0 - k1 + G3;
    float x2 = x0 - i2 + 2.0f * G3;
    float y2 = y0 - j2 + 2.0f * G3;
    float z2 = z0 - k2 + 2.0f * G3;
    float x3 = x0 - 1.0f + 3.0f * G3;
    float y3 = y0 - 1.0f + 3.0f * G3;
    float z3 = z0 - 1.0f + 3.0f * G3;
    
    // Hash coordinates of the corners
    uint32_t ii = static_cast<uint32_t>(i) & 255;
    uint32_t jj = static_cast<uint32_t>(j) & 255;
    uint32_t kk = static_cast<uint32_t>(k) & 255;
    
    // Mix in seed
    uint32_t s0 = (seed + 0x9e3779b9) & 255;
    
    int gi0 = SIMPLEX_PERM[(ii + SIMPLEX_PERM[(jj + SIMPLEX_PERM[kk & 255]) & 255]) & 255] ^ s0;
    int gi1 = SIMPLEX_PERM[(ii + i1 + SIMPLEX_PERM[(jj + j1 + SIMPLEX_PERM[(kk + k1) & 255]) & 255]) & 255] ^ s0;
    int gi2 = SIMPLEX_PERM[(ii + i2 + SIMPLEX_PERM[(jj + j2 + SIMPLEX_PERM[(kk + k2) & 255]) & 255]) & 255] ^ s0;
    int gi3 = SIMPLEX_PERM[(ii + 1 + SIMPLEX_PERM[(jj + 1 + SIMPLEX_PERM[(kk + 1) & 255]) & 255]) & 255] ^ s0;
    
    // Calculate contributions from the four corners
    float n0 = 0.0f, n1 = 0.0f, n2 = 0.0f, n3 = 0.0f;
    
    float t0 = 0.6f - x0 * x0 - y0 * y0 - z0 * z0;
    if (t0 >= 0.0f) {
        t0 *= t0;
        n0 = t0 * t0 * dotGrad(gi0, x0, y0, z0);
    }
    
    float t1 = 0.6f - x1 * x1 - y1 * y1 - z1 * z1;
    if (t1 >= 0.0f) {
        t1 *= t1;
        n1 = t1 * t1 * dotGrad(gi1, x1, y1, z1);
    }
    
    float t2 = 0.6f - x2 * x2 - y2 * y2 - z2 * z2;
    if (t2 >= 0.0f) {
        t2 *= t2;
        n2 = t2 * t2 * dotGrad(gi2, x2, y2, z2);
    }
    
    float t3 = 0.6f - x3 * x3 - y3 * y3 - z3 * z3;
    if (t3 >= 0.0f) {
        t3 *= t3;
        n3 = t3 * t3 * dotGrad(gi3, x3, y3, z3);
    }
    
    // Add contributions from each corner and scale to [-1, 1] range
    return 32.0f * (n0 + n1 + n2 + n3);
}

//=============================================================================
// FRACTAL NOISE
//=============================================================================

float fractalNoise(float x, float y, float z, const NoiseParams& params) {
    float result = 0.0f;
    
    // Primary noise layer (A)
    if (params.amplitude_a > 0.0f && params.octaves_a > 0) {
        float amp = params.amplitude_a;
        float freq = params.frequency_a;
        for (int i = 0; i < params.octaves_a; ++i) {
            result += amp * simplexNoise3D(x * freq, y * freq, z * freq, params.seed + i);
            amp *= 0.5f;
            freq *= 2.0f;
        }
    }
    
    // Secondary noise layer (B)
    if (params.amplitude_b > 0.0f && params.octaves_b > 0) {
        float amp = params.amplitude_b;
        float freq = params.frequency_b;
        for (int i = 0; i < params.octaves_b; ++i) {
            result += amp * simplexNoise3D(x * freq, y * freq, z * freq, params.seed + 100 + i);
            amp *= 0.5f;
            freq *= 2.0f;
        }
    }
    
    return result;
}

//=============================================================================
// TARGET BOUNDS COMPUTATION
//=============================================================================

AABB computeTargetBounds(const PolygonSoup& target_soup) {
    AABB bounds;
    
    // Get target mesh range (mesh_id = 0)
    auto range = target_soup.getMeshTriangleRange(0);
    
    for (size_t i = range.first; i < range.second; ++i) {
        const auto& tri = target_soup.triangles[i];
        for (int v = 0; v < 3; ++v) {
            double p[3] = {tri.v[v][0], tri.v[v][1], tri.v[v][2]};
            bounds.expand(p);
        }
    }
    
    return bounds;
}

//=============================================================================
// PLANE BASIS COMPUTATION
//=============================================================================

/**
 * @brief Compute orthonormal basis for a plane
 * 
 * Given a plane normal, computes two tangent vectors that form
 * an orthonormal basis with the normal.
 */
static void computePlaneBasis(const ExtPlane& plane,
                               float out_tangent[3],
                               float out_bitangent[3],
                               float out_normal[3]) {
    // Normalize normal
    float nx = static_cast<float>(plane.a);
    float ny = static_cast<float>(plane.b);
    float nz = static_cast<float>(plane.c);
    float nlen = std::sqrt(nx * nx + ny * ny + nz * nz);
    
    if (nlen < 1e-10f) {
        // Degenerate - use default basis
        out_normal[0] = 0.0f; out_normal[1] = 0.0f; out_normal[2] = 1.0f;
        out_tangent[0] = 1.0f; out_tangent[1] = 0.0f; out_tangent[2] = 0.0f;
        out_bitangent[0] = 0.0f; out_bitangent[1] = 1.0f; out_bitangent[2] = 0.0f;
        return;
    }
    
    nx /= nlen;
    ny /= nlen;
    nz /= nlen;
    
    out_normal[0] = nx;
    out_normal[1] = ny;
    out_normal[2] = nz;
    
    // Choose tangent based on normal direction
    // Use the axis with smallest normal component as reference
    if (std::abs(nx) < std::abs(ny) && std::abs(nx) < std::abs(nz)) {
        // X is smallest - use X axis as reference
        out_tangent[0] = 0.0f;
        out_tangent[1] = -nz;
        out_tangent[2] = ny;
    } else if (std::abs(ny) < std::abs(nz)) {
        // Y is smallest - use Y axis as reference
        out_tangent[0] = nz;
        out_tangent[1] = 0.0f;
        out_tangent[2] = -nx;
    } else {
        // Z is smallest - use Z axis as reference
        out_tangent[0] = -ny;
        out_tangent[1] = nx;
        out_tangent[2] = 0.0f;
    }
    
    // Normalize tangent
    float tlen = std::sqrt(out_tangent[0] * out_tangent[0] +
                           out_tangent[1] * out_tangent[1] +
                           out_tangent[2] * out_tangent[2]);
    if (tlen > 1e-10f) {
        out_tangent[0] /= tlen;
        out_tangent[1] /= tlen;
        out_tangent[2] /= tlen;
    }
    
    // Bitangent = normal × tangent
    out_bitangent[0] = ny * out_tangent[2] - nz * out_tangent[1];
    out_bitangent[1] = nz * out_tangent[0] - nx * out_tangent[2];
    out_bitangent[2] = nx * out_tangent[1] - ny * out_tangent[0];
}

//=============================================================================
// PLANE MESH CREATION
//=============================================================================

static void createPlaneMesh(const ExtPlane& plane, const AABB& bounds,
                            float margin, float segment_size, int max_segments,
                            std::vector<std::array<float, 3>>& vertices,
                            std::vector<std::array<uint32_t, 3>>& triangles) {
    
    vertices.clear();
    triangles.clear();
    
    // Compute expanded bounds with margin
    AABB expanded = bounds;
    expanded.expand(margin);
    
    // Compute plane basis
    float tangent[3], bitangent[3], normal[3];
    computePlaneBasis(plane, tangent, bitangent, normal);
    
    // Find a point on the plane
    // From plane equation ax + by + cz + d = 0
    // We can use: if |c| > epsilon, z = (-d - ax - by) / c
    double plane_point[3];
    double a = plane.a, b = plane.b, c = plane.c, d = plane.d;
    double nlen_sq = a * a + b * b + c * c;
    
    if (std::abs(c) > 1e-10) {
        plane_point[0] = 0.0;
        plane_point[1] = 0.0;
        plane_point[2] = -d / c;
    } else if (std::abs(b) > 1e-10) {
        plane_point[0] = 0.0;
        plane_point[1] = -d / b;
        plane_point[2] = 0.0;
    } else if (std::abs(a) > 1e-10) {
        plane_point[0] = -d / a;
        plane_point[1] = 0.0;
        plane_point[2] = 0.0;
    } else {
        // Degenerate plane
        return;
    }
    
    // Project expanded bounds onto plane tangent/bitangent space
    double corners[8][3];
    int idx = 0;
    for (int ix = 0; ix < 2; ++ix) {
        for (int iy = 0; iy < 2; ++iy) {
            for (int iz = 0; iz < 2; ++iz) {
                corners[idx][0] = (ix == 0) ? expanded.min[0] : expanded.max[0];
                corners[idx][1] = (iy == 0) ? expanded.min[1] : expanded.max[1];
                corners[idx][2] = (iz == 0) ? expanded.min[2] : expanded.max[2];
                ++idx;
            }
        }
    }
    
    // Project corners to plane coordinates
    float min_u = std::numeric_limits<float>::max();
    float max_u = std::numeric_limits<float>::lowest();
    float min_v = std::numeric_limits<float>::max();
    float max_v = std::numeric_limits<float>::lowest();
    
    for (int i = 0; i < 8; ++i) {
        double dx = corners[i][0] - plane_point[0];
        double dy = corners[i][1] - plane_point[1];
        double dz = corners[i][2] - plane_point[2];
        
        float u = static_cast<float>(dx * tangent[0] + dy * tangent[1] + dz * tangent[2]);
        float v = static_cast<float>(dx * bitangent[0] + dy * bitangent[1] + dz * bitangent[2]);
        
        min_u = std::min(min_u, u);
        max_u = std::max(max_u, u);
        min_v = std::min(min_v, v);
        max_v = std::max(max_v, v);
    }
    
    // Add extra margin in UV space
    float uv_margin = margin;
    min_u -= uv_margin;
    max_u += uv_margin;
    min_v -= uv_margin;
    max_v += uv_margin;
    
    // Compute segment count
    float extent_u = max_u - min_u;
    float extent_v = max_v - min_v;
    
    int segs_u = 1;
    int segs_v = 1;
    
    if (segment_size > 0.0f) {
        segs_u = std::max(1, std::min(max_segments, static_cast<int>(std::ceil(extent_u / segment_size))));
        segs_v = std::max(1, std::min(max_segments, static_cast<int>(std::ceil(extent_v / segment_size))));
    }
    
    // Create grid vertices
    vertices.reserve((segs_u + 1) * (segs_v + 1));
    
    for (int iv = 0; iv <= segs_v; ++iv) {
        float v = min_v + (max_v - min_v) * (static_cast<float>(iv) / segs_v);
        for (int iu = 0; iu <= segs_u; ++iu) {
            float u = min_u + (max_u - min_u) * (static_cast<float>(iu) / segs_u);
            
            // Compute 3D position: p = plane_point + u * tangent + v * bitangent
            std::array<float, 3> pos;
            pos[0] = static_cast<float>(plane_point[0] + u * tangent[0] + v * bitangent[0]);
            pos[1] = static_cast<float>(plane_point[1] + u * tangent[1] + v * bitangent[1]);
            pos[2] = static_cast<float>(plane_point[2] + u * tangent[2] + v * bitangent[2]);
            
            vertices.push_back(pos);
        }
    }
    
    // Create triangles (as quads split into two triangles each)
    triangles.reserve(segs_u * segs_v * 2);
    
    for (int iv = 0; iv < segs_v; ++iv) {
        for (int iu = 0; iu < segs_u; ++iu) {
            uint32_t i00 = iv * (segs_u + 1) + iu;
            uint32_t i10 = iv * (segs_u + 1) + (iu + 1);
            uint32_t i01 = (iv + 1) * (segs_u + 1) + iu;
            uint32_t i11 = (iv + 1) * (segs_u + 1) + (iu + 1);
            
            // First triangle: (0,0), (1,0), (0,1)
            triangles.push_back({i00, i10, i01});
            
            // Second triangle: (1,0), (1,1), (0,1)
            triangles.push_back({i10, i11, i01});
        }
    }
    
    EMBER_LOG_DEBUG("Created plane mesh: %zu vertices, %zu triangles (%dx%d grid)",
                    vertices.size(), triangles.size(), segs_u, segs_v);
}

//=============================================================================
// NOISE DISPLACEMENT
//=============================================================================

static void applyNoiseDisplacement(std::vector<std::array<float, 3>>& vertices,
                                    const ExtPlane& plane, const NoiseParams& params,
                                    float grout_offset) {
    
    if (!params.hasNoise() && grout_offset == 0.0f) {
        return;
    }
    
    // Get plane normal
    float nx = static_cast<float>(plane.a);
    float ny = static_cast<float>(plane.b);
    float nz = static_cast<float>(plane.c);
    float nlen = std::sqrt(nx * nx + ny * ny + nz * nz);
    
    if (nlen < 1e-10f) {
        return;
    }
    
    nx /= nlen;
    ny /= nlen;
    nz /= nlen;
    
    // Apply displacement to each vertex
    for (auto& v : vertices) {
        // Evaluate noise at vertex position
        float displacement = fractalNoise(v[0], v[1], v[2], params);
        
        // Add grout offset
        displacement += grout_offset;
        
        // Displace along normal
        v[0] += nx * displacement;
        v[1] += ny * displacement;
        v[2] += nz * displacement;
    }
}

//=============================================================================
// SHELL EXTRUSION
//=============================================================================

static void extrudeShell(const std::vector<std::array<float, 3>>& vertices,
                         const std::vector<std::array<uint32_t, 3>>& triangles,
                         const ExtPlane& plane, float grout_width,
                         std::vector<std::array<float, 3>>& out_vertices,
                         std::vector<std::array<uint32_t, 3>>& out_triangles) {
    
    if (grout_width <= 0.0f) {
        out_vertices = vertices;
        out_triangles = triangles;
        return;
    }
    
    // Get plane normal
    float nx = static_cast<float>(plane.a);
    float ny = static_cast<float>(plane.b);
    float nz = static_cast<float>(plane.c);
    float nlen = std::sqrt(nx * nx + ny * ny + nz * nz);
    
    if (nlen < 1e-10f) {
        out_vertices = vertices;
        out_triangles = triangles;
        return;
    }
    
    nx /= nlen;
    ny /= nlen;
    nz /= nlen;
    
    float half_width = grout_width * 0.5f;
    
    // Count output size
    size_t num_verts = vertices.size();
    size_t num_tris = triangles.size();
    
    // Output needs: top vertices + bottom vertices + side vertices
    // For simplicity, we create top and bottom surfaces
    out_vertices.reserve(num_verts * 2);
    out_triangles.reserve(num_tris * 2);
    
    // Add top surface (displaced in +normal direction)
    for (const auto& v : vertices) {
        std::array<float, 3> top_v = {
            v[0] + nx * half_width,
            v[1] + ny * half_width,
            v[2] + nz * half_width
        };
        out_vertices.push_back(top_v);
    }
    
    // Add bottom surface (displaced in -normal direction)
    for (const auto& v : vertices) {
        std::array<float, 3> bottom_v = {
            v[0] - nx * half_width,
            v[1] - ny * half_width,
            v[2] - nz * half_width
        };
        out_vertices.push_back(bottom_v);
    }
    
    // Add top triangles (same winding)
    for (const auto& tri : triangles) {
        out_triangles.push_back(tri);
    }
    
    // Add bottom triangles (reversed winding)
    for (const auto& tri : triangles) {
        out_triangles.push_back({tri[0] + static_cast<uint32_t>(num_verts),
                                  tri[2] + static_cast<uint32_t>(num_verts),
                                  tri[1] + static_cast<uint32_t>(num_verts)});
    }
    
    // Note: Side walls would require edge information, which we skip for simplicity
    // The Boolean engine will handle the open mesh correctly
    
    EMBER_LOG_DEBUG("Extruded shell: %zu vertices, %zu triangles",
                    out_vertices.size(), out_triangles.size());
}

//=============================================================================
// EXTERNAL FACE CULLING
//=============================================================================

static void cullExternalFaces(std::vector<std::array<float, 3>>& vertices,
                               std::vector<std::array<uint32_t, 3>>& triangles,
                               const AABB& target_bounds) {
    
    if (triangles.empty()) {
        return;
    }
    
    // Expand target bounds by small margin for numerical safety
    AABB expanded = target_bounds;
    expanded.expand(0.001f);
    
    // Filter triangles
    std::vector<std::array<uint32_t, 3>> kept_triangles;
    kept_triangles.reserve(triangles.size());
    
    for (const auto& tri : triangles) {
        // Compute triangle centroid
        double cx = (vertices[tri[0]][0] + vertices[tri[1]][0] + vertices[tri[2]][0]) / 3.0;
        double cy = (vertices[tri[0]][1] + vertices[tri[1]][1] + vertices[tri[2]][1]) / 3.0;
        double cz = (vertices[tri[0]][2] + vertices[tri[1]][2] + vertices[tri[2]][2]) / 3.0;
        
        double centroid[3] = {cx, cy, cz};
        
        // Keep triangle if centroid is inside expanded bounds
        if (expanded.contains(centroid)) {
            kept_triangles.push_back(tri);
        }
    }
    
    size_t culled = triangles.size() - kept_triangles.size();
    if (culled > 0) {
        EMBER_LOG_DEBUG("Culled %zu external triangles (kept %zu)",
                        culled, kept_triangles.size());
    }
    
    triangles = std::move(kept_triangles);
}

//=============================================================================
// PLANE FROM TRIANGLE FOR MESH GENERATION
//=============================================================================

/**
 * @brief Compute plane from triangle vertices for mesh generation
 */
static ExtPlane planeFromTriangleMesh(const float v0[3], const float v1[3], const float v2[3]) {
    // Convert to integers with reasonable scale
    const float scale = 1000.0f;
    int32_t iv0[3] = {
        static_cast<int32_t>(v0[0] * scale),
        static_cast<int32_t>(v0[1] * scale),
        static_cast<int32_t>(v0[2] * scale)
    };
    int32_t iv1[3] = {
        static_cast<int32_t>(v1[0] * scale),
        static_cast<int32_t>(v1[1] * scale),
        static_cast<int32_t>(v1[2] * scale)
    };
    int32_t iv2[3] = {
        static_cast<int32_t>(v2[0] * scale),
        static_cast<int32_t>(v2[1] * scale),
        static_cast<int32_t>(v2[2] * scale)
    };
    
    return planeFromTriangleExt(iv0, iv1, iv2);
}

//=============================================================================
// MAIN MESH GENERATION
//=============================================================================

void generateNoisedCutterMesh(const ExtPlane& base_plane,
                               const CutterMeshConfig& config,
                               const PolygonSoup& target_soup,
                               PolygonSoup& out_soup,
                               CutterDescriptor& out_desc) {
    
    EMBER_SCOPED_TIMER("generateNoisedCutterMesh");
    
    EMBER_LOG_INFO("Generating noised cutter mesh for mesh %u, component %u",
                   out_desc.mesh_id, out_desc.component_id);
    
    // 1. Compute target bounds
    AABB target_bounds = computeTargetBounds(target_soup);
    
    if (!target_bounds.isValid()) {
        EMBER_LOG_WARN("Invalid target bounds - skipping mesh generation");
        return;
    }
    
    // 2. Create base plane mesh
    std::vector<std::array<float, 3>> vertices;
    std::vector<std::array<uint32_t, 3>> triangles;
    
    float margin = config.noise.maxDisplacement() + config.grout_width * 0.5f + 0.1f;
    
    createPlaneMesh(base_plane, target_bounds, margin,
                    config.segment_size, config.max_segments,
                    vertices, triangles);
    
    if (vertices.empty() || triangles.empty()) {
        EMBER_LOG_WARN("Empty plane mesh generated");
        return;
    }
    
    // 3. Apply noise displacement
    float grout_offset = config.hasGrout() ? 0.0f : 0.0f;
    applyNoiseDisplacement(vertices, base_plane, config.noise, grout_offset);
    
    // 4. Extrude shell if grout enabled
    std::vector<std::array<float, 3>> final_vertices;
    std::vector<std::array<uint32_t, 3>> final_triangles;
    
    if (config.hasGrout()) {
        extrudeShell(vertices, triangles, base_plane, config.grout_width,
                     final_vertices, final_triangles);
    } else {
        final_vertices = std::move(vertices);
        final_triangles = std::move(triangles);
    }
    
    // 5. Cull external faces
    if (config.cull_external) {
        cullExternalFaces(final_vertices, final_triangles, target_bounds);
    }
    
    if (final_triangles.empty()) {
        EMBER_LOG_WARN("All triangles culled - no mesh generated");
        return;
    }
    
    // 6. Append to output soup
    uint32_t tri_start = static_cast<uint32_t>(out_soup.triangles.size());
    
    for (const auto& tri : final_triangles) {
        Triangle new_tri;
        
        // Set vertices
        for (int v = 0; v < 3; ++v) {
            new_tri.v[v][0] = final_vertices[tri[v]][0];
            new_tri.v[v][1] = final_vertices[tri[v]][1];
            new_tri.v[v][2] = final_vertices[tri[v]][2];
            
            // Quantize vertices
            for (int i = 0; i < 3; ++i) {
                new_tri.iv[v][i] = out_soup.quantize(new_tri.v[v][i], i);
            }
        }
        
        // Set mesh ID and compute plane
        new_tri.mesh_id = static_cast<int>(out_desc.mesh_id);
        new_tri.src_prim_idx = 0;
        new_tri.plane = planeFromTriangleMesh(new_tri.v[0].data(), new_tri.v[1].data(), new_tri.v[2].data());
        
        out_soup.triangles.push_back(new_tri);
    }
    
    uint32_t tri_count = static_cast<uint32_t>(out_soup.triangles.size()) - tri_start;
    
    // 7. Update descriptor
    out_desc.tri_start = tri_start;
    out_desc.tri_count = tri_count;
    
    EMBER_LOG_INFO("Generated noised cutter mesh: %u triangles (%u -> %u)",
                   tri_count, tri_start, tri_start + tri_count);
}

//=============================================================================
// BATCH GENERATION
//=============================================================================

void generateAllNoisedCutters(const CutterMeshConfig& config,
                               const PolygonSoup& target_soup,
                               PolygonSoup& soup,
                               DispatchResult& dispatch) {
    
    EMBER_SCOPED_TIMER("generateAllNoisedCutters");
    
    EMBER_LOG_INFO("Generating all Tier 2 NoisedMesh cutters");
    
    int generated_count = 0;
    
    for (auto& cutter : dispatch.cutters) {
        if (cutter.tier == CutterTier::Tier2_NoisedMesh) {
            generateNoisedCutterMesh(cutter.base_plane, config, target_soup, soup, cutter);
            ++generated_count;
        }
    }
    
    EMBER_LOG_INFO("Generated %d NoisedMesh cutters", generated_count);
}

} // namespace ember
