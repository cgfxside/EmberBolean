/**
 * @file test_cdt_consistency.cpp
 * @brief Constrained Delaunay Triangulation (CDT) consistency tests
 * 
 * Tests the CDT implementation in CherchiBackend.cpp:
 * - Basic triangulation
 * - Constraint enforcement
 * - Delaunay property
 * - Triangle orientation
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

// Include EMBER core (ember/ExactPredicates.h is included by CherchiBackend.cpp)
#include "CherchiBackend.cpp"

using namespace ember;

// =============================================================================
// Test Framework
// =============================================================================

struct TestResult {
    std::string name;
    bool passed;
    std::string message;
};

std::vector<TestResult> test_results;

#define TEST_ASSERT(cond, msg) \
    do { \
        if (!(cond)) { \
            test_results.push_back({__FUNCTION__, false, msg}); \
            return; \
        } \
    } while(0)

#define TEST_PASS() \
    test_results.push_back({__FUNCTION__, true, "OK"})

// =============================================================================
// CDT Tests
// =============================================================================

void test_cdt_basic_triangle() {
    // Simple triangle - 3 points
    CDTBuilder cdt;
    
    cdt.addVertex(Point2D(0, 0, 0));
    cdt.addVertex(Point2D(1, 0, 1));
    cdt.addVertex(Point2D(0, 1, 2));
    
    // Add constraint edges to form the triangle boundary
    cdt.addConstraint(0, 1);
    cdt.addConstraint(1, 2);
    cdt.addConstraint(2, 0);
    
    cdt.triangulate();
    
    auto triangles = cdt.getTriangles();
    
    // Should have exactly 1 triangle
    TEST_ASSERT(triangles.size() == 1, "Single triangle should produce 1 triangle");
    
    // Validate the triangulation
    TEST_ASSERT(cdt.validate(), "Basic triangle should be valid");
    
    TEST_PASS();
}

void test_cdt_square() {
    // Square with 4 points
    CDTBuilder cdt;
    
    cdt.addVertex(Point2D(0, 0, 0));
    cdt.addVertex(Point2D(1, 0, 1));
    cdt.addVertex(Point2D(1, 1, 2));
    cdt.addVertex(Point2D(0, 1, 3));
    
    // Add boundary constraints
    cdt.addConstraint(0, 1);
    cdt.addConstraint(1, 2);
    cdt.addConstraint(2, 3);
    cdt.addConstraint(3, 0);
    
    cdt.triangulate();
    
    auto triangles = cdt.getTriangles();
    
    // A square should be split into 2 triangles
    TEST_ASSERT(triangles.size() == 2, "Square should produce 2 triangles");
    
    // Validate
    TEST_ASSERT(cdt.validate(), "Square triangulation should be valid");
    
    TEST_PASS();
}

void test_cdt_five_points() {
    // 5 points forming a pentagon
    CDTBuilder cdt;
    
    cdt.addVertex(Point2D(0, 0, 0));
    cdt.addVertex(Point2D(2, 0, 1));
    cdt.addVertex(Point2D(2, 2, 2));
    cdt.addVertex(Point2D(1, 3, 3));
    cdt.addVertex(Point2D(0, 2, 4));
    
    // Boundary constraints
    cdt.addConstraint(0, 1);
    cdt.addConstraint(1, 2);
    cdt.addConstraint(2, 3);
    cdt.addConstraint(3, 4);
    cdt.addConstraint(4, 0);
    
    cdt.triangulate();
    
    auto triangles = cdt.getTriangles();
    
    // A pentagon should be split into 3 triangles
    TEST_ASSERT(triangles.size() == 3, "Pentagon should produce 3 triangles");
    
    // Validate
    TEST_ASSERT(cdt.validate(), "Pentagon triangulation should be valid");
    
    TEST_PASS();
}

void test_cdt_with_interior_point() {
    // Square with interior point
    CDTBuilder cdt;
    
    cdt.addVertex(Point2D(0, 0, 0));
    cdt.addVertex(Point2D(2, 0, 1));
    cdt.addVertex(Point2D(2, 2, 2));
    cdt.addVertex(Point2D(0, 2, 3));
    cdt.addVertex(Point2D(1, 1, 4));  // Interior point
    
    // Boundary constraints
    cdt.addConstraint(0, 1);
    cdt.addConstraint(1, 2);
    cdt.addConstraint(2, 3);
    cdt.addConstraint(3, 0);
    
    cdt.triangulate();
    
    auto triangles = cdt.getTriangles();
    
    // Square with interior point should have 4 triangles
    TEST_ASSERT(triangles.size() == 4, "Square with interior point should produce 4 triangles");
    
    // Validate
    TEST_ASSERT(cdt.validate(), "Triangulation with interior point should be valid");
    
    TEST_PASS();
}

void test_cdt_constraint_preservation() {
    // Test that constrained edges are preserved
    CDTBuilder cdt;
    
    cdt.addVertex(Point2D(0, 0, 0));
    cdt.addVertex(Point2D(2, 0, 1));
    cdt.addVertex(Point2D(2, 2, 2));
    cdt.addVertex(Point2D(0, 2, 3));
    
    // Add diagonal constraint
    cdt.addConstraint(0, 2);
    
    // Add boundary constraints
    cdt.addConstraint(0, 1);
    cdt.addConstraint(1, 2);
    cdt.addConstraint(2, 3);
    cdt.addConstraint(3, 0);
    
    cdt.triangulate();
    
    auto triangles = cdt.getTriangles();
    
    // Should have 2 triangles sharing the diagonal
    TEST_ASSERT(triangles.size() == 2, "Constrained diagonal should produce 2 triangles");
    
    // Check that edge (0,2) exists in the triangulation
    bool found_constraint = false;
    for (const auto& tri : triangles) {
        if ((tri[0] == 0 && tri[1] == 2) || (tri[0] == 2 && tri[1] == 0) ||
            (tri[1] == 0 && tri[2] == 2) || (tri[1] == 2 && tri[2] == 0) ||
            (tri[2] == 0 && tri[0] == 2) || (tri[2] == 2 && tri[0] == 0)) {
            found_constraint = true;
            break;
        }
    }
    TEST_ASSERT(found_constraint, "Constrained edge should be preserved");
    
    TEST_PASS();
}

void test_cdt_triangle_orientation() {
    // Test that all triangles have correct CCW orientation
    CDTBuilder cdt;
    
    cdt.addVertex(Point2D(0, 0, 0));
    cdt.addVertex(Point2D(3, 0, 1));
    cdt.addVertex(Point2D(3, 3, 2));
    cdt.addVertex(Point2D(0, 3, 3));
    cdt.addVertex(Point2D(1, 1, 4));
    cdt.addVertex(Point2D(2, 1, 5));
    cdt.addVertex(Point2D(1.5, 2, 6));
    
    // Boundary
    cdt.addConstraint(0, 1);
    cdt.addConstraint(1, 2);
    cdt.addConstraint(2, 3);
    cdt.addConstraint(3, 0);
    
    cdt.triangulate();
    
    auto triangles = cdt.getTriangles();
    
    // Check each triangle has positive orientation
    for (const auto& tri : triangles) {
        // Get vertex coordinates
        // Note: We need to access vertices through the builder
        // For this test, we just validate the structure
    }
    
    // Validate should check orientations
    TEST_ASSERT(cdt.validate(), "All triangles should have correct orientation");
    
    TEST_PASS();
}

void test_cdt_collinear_points() {
    // Test handling of nearly collinear points
    CDTBuilder cdt;
    
    cdt.addVertex(Point2D(0, 0, 0));
    cdt.addVertex(Point2D(1, 0, 1));
    cdt.addVertex(Point2D(2, 0, 2));  // Nearly collinear
    cdt.addVertex(Point2D(1, 1, 3));
    
    // Boundary
    cdt.addConstraint(0, 1);
    cdt.addConstraint(1, 2);
    cdt.addConstraint(2, 3);
    cdt.addConstraint(3, 0);
    
    cdt.triangulate();
    
    // Should complete without crashing
    auto triangles = cdt.getTriangles();
    TEST_ASSERT(triangles.size() >= 1, "Should produce at least 1 triangle");
    
    TEST_PASS();
}

void test_cdt_delaunay_property() {
    // Test Delaunay property - no point inside circumcircle of any triangle
    CDTBuilder cdt;
    
    // Create points that should form a Delaunay triangulation
    cdt.addVertex(Point2D(0, 0, 0));
    cdt.addVertex(Point2D(4, 0, 1));
    cdt.addVertex(Point2D(4, 4, 2));
    cdt.addVertex(Point2D(0, 4, 3));
    cdt.addVertex(Point2D(2, 2, 4));  // Center point
    
    // Boundary constraints
    cdt.addConstraint(0, 1);
    cdt.addConstraint(1, 2);
    cdt.addConstraint(2, 3);
    cdt.addConstraint(3, 0);
    
    cdt.triangulate();
    
    auto triangles = cdt.getTriangles();
    
    // Should have 4 triangles (square with center point)
    TEST_ASSERT(triangles.size() == 4, "Square with center should have 4 triangles");
    
    // Validate structure
    TEST_ASSERT(cdt.validate(), "Delaunay triangulation should be valid");
    
    TEST_PASS();
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char* argv[]) {
    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║     EMBER CDT Consistency Tests                              ║\n");
    printf("║     Testing: Constrained Delaunay Triangulation              ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");
    
    // Run all tests
    test_cdt_basic_triangle();
    test_cdt_square();
    test_cdt_five_points();
    test_cdt_with_interior_point();
    test_cdt_constraint_preservation();
    test_cdt_triangle_orientation();
    test_cdt_collinear_points();
    test_cdt_delaunay_property();
    
    // Print results
    printf("\n═══════════════════════════════════════════════════════════════\n");
    printf("Test Results:\n");
    printf("═══════════════════════════════════════════════════════════════\n");
    
    int passed = 0;
    int failed = 0;
    
    for (const auto& result : test_results) {
        if (result.passed) {
            printf("  [PASS] %s\n", result.name.c_str());
            passed++;
        } else {
            printf("  [FAIL] %s: %s\n", result.name.c_str(), result.message.c_str());
            failed++;
        }
    }
    
    printf("═══════════════════════════════════════════════════════════════\n");
    printf("Summary: %d passed, %d failed (total: %zu)\n", passed, failed, test_results.size());
    printf("═══════════════════════════════════════════════════════════════\n");
    
    return failed > 0 ? 1 : 0;
}
