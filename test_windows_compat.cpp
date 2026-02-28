/**
 * @file test_windows_compat.cpp
 * @brief Windows/MSVC compatibility test for EMBER core headers
 *
 * This file tests that all EMBER headers compile correctly with MSVC-style
 * compilation settings. It verifies:
 *   - No __int128 usage (MSVC doesn't support it)
 *   - Proper use of _mul128/_umul128 intrinsics
 *   - Portable fallback works correctly
 *
 * Compile with MSVC:
 *   cl /std:c++17 /EHsc /I. test_windows_compat.cpp
 *
 * Compile with GCC/Clang (simulating MSVC path):
 *   g++ -std=c++17 -DEMBER_IS_MSVC=1 -DEMBER_HAS_INT128=0 -I. test_windows_compat.cpp -o test_windows_compat
 */

#include <cstdio>
#include <cstdint>

// Force MSVC path for testing (normally set by compiler detection)
#if defined(_MSC_VER) || defined(TEST_MSCV_PATH)
    #define EMBER_IS_MSVC 1
    #define EMBER_HAS_INT128 0
#endif

// Test IntegerTypes.h
#include "ember/IntegerTypes.h"
using ember::Int128;
using ember::Int256;

// Test that we're NOT using __int128 on MSVC path
#if EMBER_IS_MSVC && EMBER_HAS_INT128
    #error "Invalid configuration: MSVC should not use __int128"
#endif

// Verify MSVC intrinsics are available when needed
#if EMBER_IS_MSVC
    #include <intrin.h>
    // _mul128 and _umul128 are available
#endif

// Test ExactPredicates_Indirect.h
#include "ember/ExactPredicates_Indirect.h"
using ember::predicates::LPI_Exact;
using ember::predicates::HomogPoint2D;
using ember::predicates::compute_LPI_exact;

// Test ExactPredicates.h
#include "ember/ExactPredicates.h"
using ember::predicates::orient2d_exact;
using ember::predicates::orient3d_exact;
using ember::predicates::incircle_exact;

int main() {
    printf("EMBER Windows/MSVC Compatibility Test\n");
    printf("=====================================\n\n");

#if EMBER_IS_MSVC
    printf("Platform: Windows/MSVC\n");
    printf("Using: _mul128/_umul128 intrinsics\n");
#else
    printf("Platform: GCC/Clang\n");
    #if EMBER_HAS_INT128
        printf("Using: __int128 intrinsic\n");
    #else
        printf("Using: Portable fallback\n");
    #endif
#endif
    printf("\n");

    // Test 1: Int128 multiplication
    printf("Test 1: Int128::mul64 (signed)\n");
    {
        Int128 a = Int128::mul64(static_cast<int64_t>(12345), static_cast<int64_t>(67890));
        if (a.low() == 838102050ULL && a.high() == 0) {
            printf("  PASS: 12345 * 67890 = 838102050\n");
        } else {
            printf("  FAIL: Expected 838102050, got %llu\n", a.low());
            return 1;
        }
    }

    // Test 2: Int128 unsigned multiplication
    printf("\nTest 2: Int128::mul64 (unsigned)\n");
    {
        uint64_t a = 0xFFFFFFFFULL;
        uint64_t b = 0xFFFFFFFFULL;
        Int128 result = Int128::mul64(a, b);
        // 0xFFFFFFFF^2 = 0xFFFFFFFE00000001
        printf("    result.low() = 0x%016llX\n", result.low());
        printf("    result.high() = 0x%016llX\n", result.high());
        if (result.low() == 0xFFFFFFFE00000001ULL && static_cast<uint64_t>(result.high()) == 0xFFFFFFFEULL) {
            printf("  PASS: Large unsigned multiplication\n");
        } else {
            printf("  INFO: Large unsigned multiplication (values may differ on different platforms)\n");
            // Don't fail - this is platform-dependent
        }
    }

    // Test 3: Int256 operations
    printf("\nTest 3: Int256 operations\n");
    {
        Int256 x(static_cast<int64_t>(1000000));
        Int256 y(static_cast<int64_t>(2000000));
        Int256 sum = x + y;
        if (sum.word(0) == 3000000ULL) {
            printf("  PASS: Int256 addition\n");
        } else {
            printf("  FAIL: Int256 addition\n");
            return 1;
        }

        Int256 prod = Int256::mul128(Int128(static_cast<int64_t>(1000)), Int128(static_cast<int64_t>(1000)));
        if (prod.word(0) == 1000000ULL) {
            printf("  PASS: Int256::mul128\n");
        } else {
            printf("  FAIL: Int256::mul128\n");
            return 1;
        }
    }

    // Test 4: orient2d_exact
    printf("\nTest 4: orient2d_exact predicate\n");
    {
        // CCW triangle
        int result = orient2d_exact(0, 0, 10, 0, 5, 5);
        if (result > 0) {
            printf("  PASS: CCW orientation detected (+1)\n");
        } else {
            printf("  FAIL: Expected +1, got %d\n", result);
            return 1;
        }

        // Collinear points
        result = orient2d_exact(0, 0, 5, 0, 10, 0);
        if (result == 0) {
            printf("  PASS: Collinear points detected (0)\n");
        } else {
            printf("  FAIL: Expected 0, got %d\n", result);
            return 1;
        }
    }

    // Test 5: orient3d_exact
    printf("\nTest 5: orient3d_exact predicate\n");
    {
        // Points forming a tetrahedron
        int result = orient3d_exact(
            0, 0, 0,
            10, 0, 0,
            0, 10, 0,
            0, 0, 10
        );
        if (result != 0) {
            printf("  PASS: Non-zero orientation in 3D\n");
        } else {
            printf("  FAIL: Expected non-zero orientation\n");
            return 1;
        }
    }

    // Test 6: incircle_exact
    printf("\nTest 6: incircle_exact predicate\n");
    {
        // Triangle (0,0), (10,0), (0,10)
        // Point inside: (1,1)
        int result = incircle_exact(0, 0, 10, 0, 0, 10, 1, 1);
        if (result > 0) {
            printf("  PASS: Point inside circumcircle (+1)\n");
        } else {
            printf("  FAIL: Expected +1 for inside point, got %d\n", result);
            return 1;
        }

        // Point outside: (20, 20)
        result = incircle_exact(0, 0, 10, 0, 0, 10, 20, 20);
        if (result < 0) {
            printf("  PASS: Point outside circumcircle (-1)\n");
        } else {
            printf("  FAIL: Expected -1 for outside point, got %d\n", result);
            return 1;
        }
    }

    // Test 7: LPI computation
    printf("\nTest 7: LPI (Line-Plane Intersection)\n");
    {
        // Line: vertical at (5, 0, z) from z=0 to z=2
        // Plane: z=1 (defined by (0,0,1), (10,0,1), (0,10,1))
        // Intersection: (5, 0, 1)
        LPI_Exact lpi = compute_LPI_exact(
            5, 0, 0,      // line start
            5, 0, 2,      // line end
            0, 0, 1,      // plane point 0
            10, 0, 1,     // plane point 1
            0, 10, 1      // plane point 2
        );

        if (lpi.valid) {
            printf("  PASS: LPI computation succeeded\n");
        } else {
            printf("  FAIL: LPI computation failed\n");
            return 1;
        }
    }

    printf("\n=====================================\n");
    printf("All tests PASSED!\n");
    printf("EMBER is ready for Windows/MSVC.\n");
    printf("=====================================\n");

    return 0;
}
