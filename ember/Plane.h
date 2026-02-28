/**
 * @file Plane.h
 * @brief Extended plane operations and homogeneous points for EMBER
 * 
 * This file provides:
 * - Plane utility functions (working with Plane from PolygonSoup.h)
 * - IntPlane utility functions (working with IntPlane from PolygonSoup.h)
 * - HomogeneousPoint for exact intersection representation
 * 
 * The Plane and IntPlane structs are defined in PolygonSoup.h.
 * This header adds operations and additional types.
 * 
 * @author EMBER Boolean Engine Team
 * @version 1.0
 */

#pragma once

#include <array>
#include <cstdint>
#include <cmath>
#include <limits>
#include "IntegerTypes.h"
#include "PolygonSoup.h"

namespace ember {

//=============================================================================
// PLANE UTILITIES
//=============================================================================

/**
 * @brief Utility functions for Plane (from PolygonSoup.h)
 */
namespace PlaneUtils {

/**
 * @brief Evaluate plane equation at integer point
 * @param p The plane
 * @param x X coordinate (26-bit)
 * @param y Y coordinate (26-bit)
 * @param z Z coordinate (26-bit)
 * @return ax + by + cz + d (fits in 80-bit signed)
 */
inline int64_t evaluate(const Plane& p, int32_t x, int32_t y, int32_t z) {
    return static_cast<int64_t>(p.a) * x + 
           static_cast<int64_t>(p.b) * y + 
           static_cast<int64_t>(p.c) * z + p.d;
}

/**
 * @brief Evaluate plane equation at floating-point point
 * @return Double-precision result (may have rounding errors)
 */
inline double evaluate(const Plane& p, double x, double y, double z) {
    return static_cast<double>(p.a) * x + 
           static_cast<double>(p.b) * y + 
           static_cast<double>(p.c) * z + 
           static_cast<double>(p.d);
}

/**
 * @brief Check if plane is valid (non-zero normal)
 */
inline bool isValid(const Plane& p) {
    return p.a != 0 || p.b != 0 || p.c != 0;
}

/**
 * @brief Get squared length of normal vector
 */
inline int64_t normalLengthSq(const Plane& p) {
    return static_cast<int64_t>(p.a) * p.a + 
           static_cast<int64_t>(p.b) * p.b + 
           static_cast<int64_t>(p.c) * p.c;
}

/**
 * @brief Check if point is on positive side of plane
 * @return true if ax + by + cz + d > 0
 */
inline bool isPositive(const Plane& p, int32_t x, int32_t y, int32_t z) {
    return evaluate(p, x, y, z) > 0;
}

/**
 * @brief Check if point is on negative side of plane
 * @return true if ax + by + cz + d < 0
 */
inline bool isNegative(const Plane& p, int32_t x, int32_t y, int32_t z) {
    return evaluate(p, x, y, z) < 0;
}

/**
 * @brief Check if point is on plane
 * @return true if ax + by + cz + d == 0
 */
inline bool contains(const Plane& p, int32_t x, int32_t y, int32_t z) {
    return evaluate(p, x, y, z) == 0;
}

/**
 * @brief Compute orientation of point relative to plane
 * @return -1 (negative), 0 (on plane), or +1 (positive)
 */
inline int orientation(const Plane& p, int32_t x, int32_t y, int32_t z) {
    auto val = evaluate(p, x, y, z);
    return (val > 0) - (val < 0);
}

/**
 * @brief Flip plane normal (multiply all coefficients by -1)
 */
inline void flip(Plane& p) {
    p.a = -p.a;
    p.b = -p.b;
    p.c = -p.c;
    p.d = -p.d;
}

} // namespace PlaneUtils

//=============================================================================
// INTPLANE UTILITIES
//=============================================================================

/**
 * @brief Utility functions for IntPlane (from PolygonSoup.h)
 */
namespace IntPlaneUtils {

/**
 * @brief Check if IntPlane is valid (non-zero normal)
 */
inline bool isValid(const IntPlane& p) {
    return p.a != 0 || p.b != 0 || p.c != 0;
}

/**
 * @brief Evaluate IntPlane at integer point (returns Int256 for safety)
 * 
 * Bit-width: 53 + 26 + log2(4) = 81 bits per product
 * Sum of 3 products + d: 81 + log2(3) + 128 ≈ 211 bits (fits in Int256)
 */
inline Int256 evaluate(const IntPlane& p, int64_t x, int64_t y, int64_t z) {
    Int256 ax = Int256::mul128(Int128(p.a), Int128(x));
    Int256 by = Int256::mul128(Int128(p.b), Int128(y));
    Int256 cz = Int256::mul128(Int128(p.c), Int128(z));
    return ax + by + cz + Int256(p.d);
}

} // namespace IntPlaneUtils

//=============================================================================
// HOMOGENEOUS POINT (80-BIT COORDINATES)
//=============================================================================

/**
 * @brief 3D point in homogeneous coordinates (x, y, z, w)
 * 
 * Homogeneous coordinates represent a 3D point (X, Y, Z) as (x, y, z, w)
 * where X = x/w, Y = y/w, Z = z/w.
 * 
 * EMBER uses homogeneous points for exact representation of:
 * - LPI (Line-Plane Intersection) points: up to 80-bit per component
 * - TPI (Triangle-Plane Intersection) points: up to 160-bit per component
 * 
 * Using Int128 provides headroom for both cases with a single type.
 * 
 * @note When w = 0, the point represents a direction (point at infinity)
 */
struct HomogeneousPoint {
    Int128 x;  ///< X numerator
    Int128 y;  ///< Y numerator
    Int128 z;  ///< Z numerator
    Int128 w;  ///< Denominator
    
    /**
     * @brief Default constructor creates point at origin
     */
    HomogeneousPoint() : x(0), y(0), z(0), w(1) {}
    
    /**
     * @brief Construct from explicit integer coordinates (w = 1)
     */
    HomogeneousPoint(int64_t x_, int64_t y_, int64_t z_) 
        : x(x_), y(y_), z(z_), w(1) {}
    
    /**
     * @brief Construct from Int128 components
     */
    HomogeneousPoint(const Int128& x_, const Int128& y_, 
                     const Int128& z_, const Int128& w_)
        : x(x_), y(y_), z(z_), w(w_) {}
    
    /**
     * @brief Check if point is at infinity (w = 0)
     */
    bool isAtInfinity() const {
        return w.isZero();
    }
    
    /**
     * @brief Check if point is valid (finite and w != 0)
     */
    bool isValid() const {
        return !w.isZero();
    }
    
    /**
     * @brief Check if point is the origin
     */
    bool isOrigin() const {
        return x.isZero() && y.isZero() && z.isZero() && !w.isZero();
    }
    
    /**
     * @brief Convert to double-precision Cartesian coordinates
     * 
     * Performs division x/w, y/w, z/w with appropriate rounding.
     * 
     * @param[out] ox Output X coordinate
     * @param[out] oy Output Y coordinate
     * @param[out] oz Output Z coordinate
     * 
     * @warning May produce inf/nan if w = 0 or values overflow double range
     */
    void toDouble(double& ox, double& oy, double& oz) const {
        if (w.isZero()) {
            ox = oy = oz = std::numeric_limits<double>::infinity();
            return;
        }
        double wd = static_cast<double>(w);
        ox = static_cast<double>(x) / wd;
        oy = static_cast<double>(y) / wd;
        oz = static_cast<double>(z) / wd;
    }
    
    /**
     * @brief Get Cartesian X coordinate as double
     */
    double xDouble() const {
        return w.isZero() ? std::numeric_limits<double>::infinity() 
                          : static_cast<double>(x) / static_cast<double>(w);
    }
    
    /**
     * @brief Get Cartesian Y coordinate as double
     */
    double yDouble() const {
        return w.isZero() ? std::numeric_limits<double>::infinity() 
                          : static_cast<double>(y) / static_cast<double>(w);
    }
    
    /**
     * @brief Get Cartesian Z coordinate as double
     */
    double zDouble() const {
        return w.isZero() ? std::numeric_limits<double>::infinity() 
                          : static_cast<double>(z) / static_cast<double>(w);
    }
};

//=============================================================================
// PLANE CONSTRUCTION FROM POINTS
//=============================================================================

/**
 * @brief Construct a plane from three points
 * 
 * Computes the plane equation ax + by + cz + d = 0 where:
 * - (a, b, c) is the cross product of (p1-p0) and (p2-p0)
 * - d = -(a*p0.x + b*p0.y + c*p0.z)
 * 
 * @param p0 First point (vertex indices)
 * @param p1 Second point (vertex indices)
 * @param p2 Third point (vertex indices)
 * @param vertices Array of vertex coordinates
 * @return Plane with 26-bit coefficients
 */
Plane planeFromPoints(uint32_t p0, uint32_t p1, uint32_t p2, 
                      const std::vector<PolygonSoup::IntVertex>& vertices);

/**
 * @brief Construct a plane from three explicit points
 */
Plane planeFromPoints(const PolygonSoup::IntVertex& p0, 
                      const PolygonSoup::IntVertex& p1, 
                      const PolygonSoup::IntVertex& p2);

} // namespace ember
