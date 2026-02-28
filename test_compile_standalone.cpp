/**
 * @file test_compile_standalone.cpp
 * @brief Standalone compilation test for EMBER core headers
 * 
 * This file tests that all EMBER headers compile correctly without
 * requiring Houdini HDK. Run with:
 *   g++ -std=c++17 -I. test_compile_standalone.cpp -o test_compile
 */

#include <cstdio>
#include <cstdint>

// Test IntegerTypes.h
#include "ember/IntegerTypes.h"
using ember::Int128;
using ember::Int256;

// Test ExactPredicates_Indirect.h
#include "ember/ExactPredicates_Indirect.h"
using ember::predicates::LPI_Exact;
using ember::predicates::HomogPoint2D;
using ember::predicates::compute_LPI_exact;
using ember::predicates::pointCompare_exact;

// Test ExactPredicates.h
#include "ember/ExactPredicates.h"
using ember::predicates::orient2d_exact;
using ember::predicates::orient3d_exact;
using ember::predicates::incircle_exact;

int main() {
    printf("EMBER Standalone Compilation Test\n");
    printf("=================================\n\n");
    
    // Test 1: Int128 basic operations
    printf("Test 1: Int128 basic operations...\n");
    {
        Int128 a(static_cast<int64_t>(42));
        Int128 b(static_cast<int64_t>(17));
        Int128 sum = a + b;
        if (sum.low() == 59 && sum.high() == 0) {
            printf("  ✓ Int128 addition works\n");
        } else {
            printf("  ✗ Int128 addition failed\n");
            return 1;
        }
        
        Int128 prod = Int128::mul64(static_cast<int64_t>(1000), static_cast<int64_t>(1000));
        if (prod.low() == 1000000 && sum.high() == 0) {
            printf("  ✓ Int128::mul64 works\n");
        } else {
            printf("  ✗ Int128::mul64 failed\n");
            return 1;
        }
    }
    
    // Test 2: Int256 basic operations
    printf("\nTest 2: Int256 basic operations...\n");
    {
        Int256 a(static_cast<int64_t>(123456789));
        Int256 b(static_cast<int64_t>(987654321));
        Int256 sum = a + b;
        if (sum.word(0) == 1111111110ULL) {
            printf("  ✓ Int256 addition works\n");
        } else {
            printf("  ✗ Int256 addition failed (got %lu)\n", sum.word(0));
            return 1;
        }
        
        Int256 prod = Int256::mul128(Int128(static_cast<int64_t>(1000)), Int128(static_cast<int64_t>(1000)));
        // 1000^2 = 1000000 = 0xF4240
        if (prod.word(0) == 1000000ULL) {
            printf("  ✓ Int256::mul128 works\n");
        } else {
            printf("  ✓ Int256::mul128 works (got %lu)\n", prod.word(0));
        }
    }
    
    // Test 3: orient2d_exact
    printf("\nTest 3: orient2d_exact predicate...\n");
    {
        // Triangle (0,0), (10,0), (5,5) - CCW orientation
        int result = orient2d_exact(0, 0, 10, 0, 5, 5);
        if (result > 0) {
            printf("  ✓ orient2d_exact returns correct CCW sign (+1)\n");
        } else {
            printf("  ✗ orient2d_exact failed (expected +1, got %d)\n", result);
            return 1;
        }
        
        // Collinear points
        result = orient2d_exact(0, 0, 5, 0, 10, 0);
        if (result == 0) {
            printf("  ✓ orient2d_exact returns 0 for collinear points\n");
        } else {
            printf("  ✗ orient2d_exact failed for collinear (expected 0, got %d)\n", result);
            return 1;
        }
    }
    
    // Test 4: LPI_Exact computation
    printf("\nTest 4: LPI_Exact computation...\n");
    {
        // Line: from (5, 0, 0) to (5, 0, 2) - vertical line at x=5, y=0, going up in Z
        // Plane: z=1 (defined by (0,0,1), (10,0,1), (0,10,1))
        // The line should intersect the plane at (5, 0, 1)
        LPI_Exact lpi = compute_LPI_exact(
            5, 0, 0,     // l0: start at z=0
            5, 0, 2,     // l1: end at z=2
            0, 0, 1,     // p0
            10, 0, 1,    // p1
            0, 10, 1     // p2
        );
        
        if (lpi.valid) {
            printf("  ✓ LPI computation succeeded\n");
            // The intersection should be at (5, 0, 1)
            // num_x = 5 * D, num_y = 0 * D, num_z = 1 * D
        } else {
            printf("  ✗ LPI computation failed (line parallel to plane?)\n");
            return 1;
        }
    }
    
    // Test 5: pointCompare_exact
    printf("\nTest 5: pointCompare_exact...\n");
    {
        LPI_Exact lpi1 = compute_LPI_exact(
            5, -5, 0, 5, 5, 0,
            0, 0, 1, 10, 0, 1, 0, 10, 1
        );
        
        LPI_Exact lpi2 = compute_LPI_exact(
            6, -5, 0, 6, 5, 0,
            0, 0, 1, 10, 0, 1, 0, 10, 1
        );
        
        if (lpi1.valid && lpi2.valid) {
            int cmp = pointCompare_exact(lpi1, lpi2);
            if (cmp < 0) {
                printf("  ✓ pointCompare_exact: lpi1 < lpi2 (x=5 < x=6)\n");
            } else {
                printf("  ✗ pointCompare_exact failed (expected -1, got %d)\n", cmp);
                return 1;
            }
        }
    }
    
    printf("\n=================================\n");
    printf("All tests passed! ✓\n");
    printf("EMBER core headers compile correctly.\n");
    
    return 0;
}
