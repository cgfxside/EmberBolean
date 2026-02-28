/**
 * @file test_integers.cpp
 * @brief Unit tests for Int128 and Int256
 */

#include "IntegerTypes.h"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace ember;

void test_int128_basic() {
    std::cout << "Testing Int128 basic operations..." << std::endl;
    
    // Construction
    Int128 a;
    assert(a.isZero());
    assert(a.sign() == 0);
    
    Int128 b(int64_t(42));
    assert(b.lo == uint64_t(42));
    assert(b.hi == int64_t(0));
    
    Int128 c(int64_t(-42));
    assert(c.isNegative());
    assert(c.sign() == -1);
    
    // Comparison
    assert(b > c);
    assert(c < b);
    assert(b == Int128(int64_t(42)));
    assert(b != c);
    
    // Addition
    Int128 d = b + Int128(int64_t(8));
    assert(d == Int128(int64_t(50)));
    
    // Subtraction
    Int128 e = b - Int128(int64_t(10));
    assert(e == Int128(int64_t(32)));
    
    // Negation
    Int128 f = -b;
    assert(f == c);
    assert(-f == b);
    
    // Increment/Decrement
    Int128 g(int64_t(5));
    assert(++g == Int128(int64_t(6)));
    assert(g++ == Int128(int64_t(6)));
    assert(g == Int128(int64_t(7)));
    assert(--g == Int128(int64_t(6)));
    assert(g-- == Int128(int64_t(6)));
    assert(g == Int128(int64_t(5)));
    
    std::cout << "  Int128 basic operations: PASSED" << std::endl;
}

void test_int128_mul64() {
    std::cout << "Testing Int128::mul64..." << std::endl;
    
    // Basic multiplication
    Int128 p1 = Int128::mul64(int64_t(5), int64_t(7));
    assert(p1 == Int128(int64_t(35)));
    
    // Negative multiplication
    Int128 p2 = Int128::mul64(int64_t(-5), int64_t(7));
    assert(p2 == Int128(int64_t(-35)));
    
    Int128 p3 = Int128::mul64(int64_t(5), int64_t(-7));
    assert(p3 == Int128(int64_t(-35)));
    
    Int128 p4 = Int128::mul64(int64_t(-5), int64_t(-7));
    assert(p4 == Int128(int64_t(35)));
    
    // Large values (26-bit * 26-bit = 52-bit, fits in 80-bit budget)
    int64_t large1 = (1LL << 25) - 1;  // 2^25 - 1 = 33,554,431
    int64_t large2 = (1LL << 25) - 1;
    Int128 p5 = Int128::mul64(large1, large2);
    // Verify by computing expected result
    __int128 expected = static_cast<__int128>(large1) * large2;
    assert(p5.lo == static_cast<uint64_t>(expected));
    assert(p5.hi == static_cast<int64_t>(expected >> 64));
    
    // Edge case: INT64_MIN
    Int128 p6 = Int128::mul64(INT64_MIN, int64_t(0));
    assert(p6.isZero());
    
    Int128 p7 = Int128::mul64(INT64_MIN, int64_t(1));
    assert(p7 == Int128(INT64_MIN));
    
    std::cout << "  Int128::mul64: PASSED" << std::endl;
}

void test_int128_bitwise() {
    std::cout << "Testing Int128 bitwise operations..." << std::endl;
    
    Int128 a(int64_t(0x0F0F0F0F0F0F0F0FULL), uint64_t(0xF0F0F0F0F0F0F0F0ULL));
    Int128 b(int64_t(0xAAAAAAAAAAAAAAAAULL), uint64_t(0x5555555555555555ULL));
    
    // AND
    Int128 c = a & b;
    assert(c.lo == uint64_t(0x5050505050505050ULL));
    assert(static_cast<uint64_t>(c.hi) == uint64_t(0x0A0A0A0A0A0A0A0AULL));
    
    // OR
    Int128 d = a | b;
    assert(d.lo == uint64_t(0xF5F5F5F5F5F5F5F5ULL));
    assert(static_cast<uint64_t>(d.hi) == uint64_t(0xAFAFAFAFAFAFAFAFULL));
    
    // XOR
    Int128 e = a ^ b;
    assert(e.lo == uint64_t(0xA5A5A5A5A5A5A5A5ULL));
    assert(static_cast<uint64_t>(e.hi) == uint64_t(0xA5A5A5A5A5A5A5A5ULL));
    
    // NOT
    Int128 f = ~a;
    assert(f.lo == ~uint64_t(0xF0F0F0F0F0F0F0F0ULL));
    assert(f.hi == ~int64_t(0x0F0F0F0F0F0F0F0FULL));
    
    // Shift left
    Int128 g(int64_t(1));
    Int128 h = g << 64;
    assert(h.lo == uint64_t(0));
    assert(h.hi == int64_t(1));
    
    Int128 i = g << 65;
    assert(i.lo == uint64_t(0));
    assert(i.hi == int64_t(2));
    
    // Shift right (arithmetic)
    Int128 j(int64_t(-4));
    Int128 k = j >> 1;
    assert(k == Int128(int64_t(-2)));
    
    Int128 l(int64_t(-8));
    Int128 m = l >> 2;
    assert(m == Int128(int64_t(-2)));
    
    std::cout << "  Int128 bitwise operations: PASSED" << std::endl;
}

void test_int256_basic() {
    std::cout << "Testing Int256 basic operations..." << std::endl;
    
    // Construction
    Int256 a;
    assert(a.isZero());
    assert(a.sign() == 0);
    
    Int256 b(Int128(int64_t(42)));
    assert(b.lo == Int128(int64_t(42)));
    assert(b.hi == Int128::zero());
    
    Int256 c(Int128(int64_t(-42)));
    assert(c.isNegative());
    assert(c.sign() == -1);
    
    // Comparison
    assert(b > c);
    assert(c < b);
    assert(b == Int256(Int128(int64_t(42))));
    assert(b != c);
    
    // Addition
    Int256 d = b + Int256(Int128(int64_t(8)));
    assert(d == Int256(Int128(int64_t(50))));
    
    // Subtraction
    Int256 e = b - Int256(Int128(int64_t(10)));
    assert(e == Int256(Int128(int64_t(32))));
    
    // Negation
    Int256 f = -b;
    assert(f == c);
    assert(-f == b);
    
    std::cout << "  Int256 basic operations: PASSED" << std::endl;
}

void test_int256_mul128() {
    std::cout << "Testing Int256::mul128..." << std::endl;
    
    // Basic multiplication
    Int256 p1 = Int256::mul128(Int128(int64_t(5)), Int128(int64_t(7)));
    assert(p1.lo == Int128(int64_t(35)));
    assert(p1.hi == Int128::zero());
    
    // Negative multiplication
    Int256 p2 = Int256::mul128(Int128(int64_t(-5)), Int128(int64_t(7)));
    assert(p2.lo == Int128(int64_t(-35)));
    assert(p2.hi == Int128(int64_t(-1)));  // Sign extension
    
    Int256 p3 = Int256::mul128(Int128(int64_t(5)), Int128(int64_t(-7)));
    assert(p3.lo == Int128(int64_t(-35)));
    assert(p3.hi == Int128(int64_t(-1)));
    
    Int256 p4 = Int256::mul128(Int128(int64_t(-5)), Int128(int64_t(-7)));
    assert(p4.lo == Int128(int64_t(35)));
    assert(p4.hi == Int128::zero());
    
    // Large values (128-bit * 128-bit)
    Int128 large1(int64_t(1) << 62, uint64_t(0));
    Int128 large2(int64_t(1) << 62, uint64_t(0));
    Int256 p5 = Int256::mul128(large1, large2);
    // (2^126)^2 = 2^252
    // p5.hi is Int128, so p5.hi.lo is uint64_t, p5.hi.hi is int64_t
    assert(p5.hi.lo == uint64_t(0));
    assert(p5.hi.hi == int64_t(1) << 60);  // 2^124 in high 128 bits
    assert(p5.lo.isZero());
    
    std::cout << "  Int256::mul128: PASSED" << std::endl;
}

void test_carry_propagation() {
    std::cout << "Testing carry propagation..." << std::endl;
    
    // Test Int128 carry
    Int128 a(int64_t(0), UINT64_MAX);
    Int128 b(int64_t(0), uint64_t(1));
    Int128 c = a + b;
    assert(c.lo == uint64_t(0));
    assert(c.hi == int64_t(1));
    
    // Test Int256 carry
    Int256 d(Int128::zero(), Int128(int64_t(0), UINT64_MAX));
    Int256 e(Int128::zero(), Int128(int64_t(0), uint64_t(1)));
    Int256 f = d + e;
    assert(f.lo.lo == uint64_t(0));
    assert(f.lo.hi == int64_t(1));
    assert(f.hi.isZero());
    
    // Test borrow
    Int128 g(int64_t(1), uint64_t(0));
    Int128 h(int64_t(0), uint64_t(1));
    Int128 i = g - h;
    assert(i.lo == UINT64_MAX);
    assert(i.hi == int64_t(0));
    
    std::cout << "  Carry propagation: PASSED" << std::endl;
}

void test_ember_budget() {
    std::cout << "Testing EMBER bit-width budget..." << std::endl;
    
    // 26-bit input coordinates
    Coord coord1 = (1 << 25) - 1;  // 33,554,431
    Coord coord2 = (1 << 25) - 1;
    
    // 2×2 minor: 53-bit (fits in int64_t)
    Minor2 minor = static_cast<Minor2>(coord1) * coord2;
    assert(minor > 0);
    
    // 3×3 determinant: 80-bit (fits in Int128)
    // Simulate a 3×3 determinant calculation
    Int128 det = Int128::mul64(static_cast<int64_t>(coord1), static_cast<int64_t>(coord2));
    det = det * static_cast<int64_t>(coord1);  // coord1^3
    assert(!det.isZero());
    
    // classify() dot product: 160-bit (fits in Int256)
    Int128 a(Int128::mul64(static_cast<int64_t>(coord1), static_cast<int64_t>(coord2)));
    Int128 b(Int128::mul64(static_cast<int64_t>(coord1), static_cast<int64_t>(coord2)));
    Int256 dot = Int256::mul128(a, b);
    assert(!dot.isZero());
    
    std::cout << "  EMBER bit-width budget: PASSED" << std::endl;
}

void test_type_sizes() {
    std::cout << "Testing type sizes..." << std::endl;
    
    assert(sizeof(Int128) == 16);
    assert(sizeof(Int256) == 32);
    assert(alignof(Int128) >= 8);
    assert(alignof(Int256) >= 8);
    
    std::cout << "  Type sizes: PASSED" << std::endl;
    std::cout << "    sizeof(Int128) = " << sizeof(Int128) << " bytes" << std::endl;
    std::cout << "    sizeof(Int256) = " << sizeof(Int256) << " bytes" << std::endl;
}

int main() {
    std::cout << "======================================" << std::endl;
    std::cout << "EMBER IntegerTypes Test Suite" << std::endl;
    std::cout << "======================================" << std::endl;
    
    test_type_sizes();
    test_int128_basic();
    test_int128_mul64();
    test_int128_bitwise();
    test_int256_basic();
    test_int256_mul128();
    test_carry_propagation();
    test_ember_budget();
    
    std::cout << "======================================" << std::endl;
    std::cout << "All tests PASSED!" << std::endl;
    std::cout << "======================================" << std::endl;
    
    return 0;
}
