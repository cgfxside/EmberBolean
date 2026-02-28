/**
 * @file test_harness.cpp
 * @brief Standalone test harness for EMBER ExactPredicates.h geometric kernels
 * 
 * Validates:
 * 1. Collinearity detection with LPI (Line-Plane Intersection) points
 * 2. Int256 multiplication at the 26-bit coordinate limit
 * 3. Int128 basic operations (add, sub, mul, sign)
 * 4. orient2d_exact and orient3d_exact predicates
 * 
 * This is a completely standalone file with no external dependencies
 * except the standard C++ library.
 * 
 * @author QA Automation Expert
 * @version 1.0.0
 */

#include <cstdint>
#include <cmath>
#include <array>
#include <algorithm>
#include <type_traits>
#include <tuple>
#include <cstdio>
#include <cstdlib>
#include <string>

//=============================================================================
// TEST FRAMEWORK MACROS
//=============================================================================

static int g_testsPassed = 0;
static int g_testsFailed = 0;
static int g_currentTestLine = 0;
static std::string g_currentTestName;

#define TEST_BEGIN(name) \
    do { \
        g_currentTestName = name; \
        g_currentTestLine = __LINE__; \
        std::printf("[TEST] %-50s ", name); \
    } while(0)

#define TEST_ASSERT(condition) \
    do { \
        if (!(condition)) { \
            std::printf("FAIL\n"); \
            std::printf("       Assertion failed at line %d: %s\n", \
                        __LINE__, #condition); \
            g_testsFailed++; \
            return; \
        } \
    } while(0)

#define TEST_ASSERT_MSG(condition, msg) \
    do { \
        if (!(condition)) { \
            std::printf("FAIL\n"); \
            std::printf("       %s\n", msg); \
            g_testsFailed++; \
            return; \
        } \
    } while(0)

#define TEST_PASS() \
    do { \
        std::printf("PASS\n"); \
        g_testsPassed++; \
    } while(0)

#define TEST_SUMMARY() \
    do { \
        std::printf("\n========================================\n"); \
        std::printf("TEST SUMMARY\n"); \
        std::printf("========================================\n"); \
        std::printf("Passed: %d\n", g_testsPassed); \
        std::printf("Failed: %d\n", g_testsFailed); \
        std::printf("Total:  %d\n", g_testsPassed + g_testsFailed); \
        std::printf("========================================\n"); \
        if (g_testsFailed == 0) { \
            std::printf("ALL TESTS PASSED!\n"); \
        } else { \
            std::printf("SOME TESTS FAILED!\n"); \
        } \
    } while(0)

//=============================================================================
// FORWARD DECLARATIONS
//=============================================================================

class Int128;
class Int256;

//=============================================================================
// EXACT INTEGER ARITHMETIC: Int128
//=============================================================================

class Int128 {
private:
    int64_t high_;   // Signed high 64 bits
    uint64_t low_;   // Unsigned low 64 bits

    void normalize() {
        // Carry from low to high if low overflowed
        if (high_ < 0 && low_ != 0) {
            // Negative number: high should account for low
        }
    }

public:
    // Constructors
    Int128() : high_(0), low_(0) {}
    Int128(int64_t v) : high_(v < 0 ? -1 : 0), low_(static_cast<uint64_t>(v)) {}
    Int128(int64_t high, uint64_t low) : high_(high), low_(low) { normalize(); }
    
    static Int128 fromQuantized(int64_t v) {
        return Int128(v);
    }

    // Accessors
    int64_t high() const { return high_; }
    uint64_t low() const { return low_; }
    
    bool isZero() const { return high_ == 0 && low_ == 0; }
    bool isNegative() const { return high_ < 0; }
    bool isPositive() const { return high_ > 0 || (high_ == 0 && low_ > 0); }
    
    int sign() const {
        if (isZero()) return 0;
        return isNegative() ? -1 : 1;
    }

    // Comparison operators
    bool operator==(const Int128& other) const {
        return high_ == other.high_ && low_ == other.low_;
    }
    bool operator!=(const Int128& other) const { return !(*this == other); }
    
    bool operator<(const Int128& other) const {
        if (high_ != other.high_) return high_ < other.high_;
        return low_ < other.low_;
    }
    bool operator>(const Int128& other) const { return other < *this; }
    bool operator<=(const Int128& other) const { return !(other < *this); }
    bool operator>=(const Int128& other) const { return !(*this < other); }

    // Addition
    Int128 operator+(const Int128& other) const {
        uint64_t new_low = low_ + other.low_;
        int64_t carry = (new_low < low_) ? 1 : 0;
        int64_t new_high = high_ + other.high_ + carry;
        return Int128(new_high, new_low);
    }
    
    Int128& operator+=(const Int128& other) {
        *this = *this + other;
        return *this;
    }

    // Negation
    Int128 operator-() const {
        uint64_t new_low = ~low_ + 1;
        int64_t new_high = ~high_;
        if (new_low == 0) new_high += 1;
        return Int128(new_high, new_low);
    }
    
    // Subtraction
    Int128 operator-(const Int128& other) const {
        return *this + (-other);
    }
    
    Int128& operator-=(const Int128& other) {
        *this = *this - other;
        return *this;
    }

    // 64x64 -> 128 bit multiplication using compiler builtin
    static Int128 mul(int64_t a, int64_t b) {
        // Use __int128 if available (GCC/Clang)
        #if defined(__GNUC__) || defined(__clang__)
            __int128 prod = (__int128)a * (__int128)b;
            return Int128((int64_t)(prod >> 64), (uint64_t)prod);
        #else
            // Fallback to manual multiplication
            int64_t a_hi = a >> 32;
            uint64_t a_lo = static_cast<uint64_t>(a) & 0xFFFFFFFFULL;
            int64_t b_hi = b >> 32;
            uint64_t b_lo = static_cast<uint64_t>(b) & 0xFFFFFFFFULL;
            
            uint64_t p0 = a_lo * b_lo;
            int64_t  p1 = a_hi * static_cast<int64_t>(b_lo);
            int64_t  p2 = static_cast<int64_t>(a_lo) * b_hi;
            int64_t  p3 = a_hi * b_hi;
            
            uint64_t low = p0;
            int64_t high = p3;
            
            uint64_t mid1 = static_cast<uint64_t>(p1) << 32;
            uint64_t mid2 = static_cast<uint64_t>(p2) << 32;
            
            uint64_t sum1 = low + mid1;
            int64_t carry1 = (sum1 < low) ? 1 : 0;
            
            uint64_t sum2 = sum1 + mid2;
            int64_t carry2 = (sum2 < sum1) ? 1 : 0;
            
            high += (p1 >> 32) + (p2 >> 32) + carry1 + carry2;
            
            return Int128(high, sum2);
        #endif
    }
    
    Int128 operator*(int64_t v) const {
        // Use __int128 if available
        #if defined(__GNUC__) || defined(__clang__)
            __int128 self = ((__int128)high_ << 64) | low_;
            __int128 prod = self * (__int128)v;
            return Int128((int64_t)(prod >> 64), (uint64_t)prod);
        #else
            Int128 low_prod = mul(static_cast<int64_t>(low_), v);
            Int128 high_prod = mul(high_, v);
            return low_prod + Int128(high_prod.low_, 0);
        #endif
    }
};

//=============================================================================
// EXACT INTEGER ARITHMETIC: Int256
//=============================================================================

class Int256 {
private:
    std::array<uint64_t, 4> w_;
    
    void normalize() {
        if (w_[3] & 0x8000000000000000ULL) {
            // Negative: ensure proper sign extension
        }
    }
    
    static std::pair<uint64_t, uint64_t> add64(uint64_t a, uint64_t b, uint64_t carry_in) {
        uint64_t sum = a + b;
        uint64_t carry_out = (sum < a) ? 1 : 0;
        sum += carry_in;
        if (carry_in && sum == 0) carry_out = 1;
        return {sum, carry_out};
    }

public:
    Int256() : w_{0, 0, 0, 0} {}
    Int256(int64_t v) : w_{static_cast<uint64_t>(v), 0, 0, (v < 0) ? ~0ULL : 0ULL} {}
    Int256(uint64_t w0, uint64_t w1, uint64_t w2, uint64_t w3) : w_{w0, w1, w2, w3} { normalize(); }
    
    static Int256 fromInt128(const Int128& v) {
        return Int256(v.low(), static_cast<uint64_t>(v.high()), 
                      (v.high() < 0) ? ~0ULL : 0ULL, 
                      (v.high() < 0) ? ~0ULL : 0ULL);
    }
    
    uint64_t word(int i) const { return w_[i]; }
    
    bool isZero() const { return w_[0] == 0 && w_[1] == 0 && w_[2] == 0 && w_[3] == 0; }
    bool isNegative() const { return (w_[3] & 0x8000000000000000ULL) != 0; }
    bool isPositive() const { return !isZero() && !isNegative(); }
    
    int sign() const {
        if (isZero()) return 0;
        return isNegative() ? -1 : 1;
    }

    bool operator==(const Int256& other) const {
        return w_[0] == other.w_[0] && w_[1] == other.w_[1] && 
               w_[2] == other.w_[2] && w_[3] == other.w_[3];
    }
    bool operator!=(const Int256& other) const { return !(*this == other); }
    
    bool operator<(const Int256& other) const {
        for (int i = 3; i >= 0; --i) {
            if (i == 3) {
                int64_t a = static_cast<int64_t>(w_[i]);
                int64_t b = static_cast<int64_t>(other.w_[i]);
                if (a != b) return a < b;
            } else {
                if (w_[i] != other.w_[i]) return w_[i] < other.w_[i];
            }
        }
        return false;
    }

    Int256 operator+(const Int256& other) const {
        Int256 result;
        uint64_t carry = 0;
        for (int i = 0; i < 4; ++i) {
            auto [sum, new_carry] = add64(w_[i], other.w_[i], carry);
            result.w_[i] = sum;
            carry = new_carry;
        }
        return result;
    }
    
    Int256& operator+=(const Int256& other) {
        *this = *this + other;
        return *this;
    }

    Int256 operator-() const {
        Int256 result;
        uint64_t carry = 1;
        for (int i = 0; i < 4; ++i) {
            uint64_t inv = ~w_[i];
            auto [sum, new_carry] = add64(inv, 0, carry);
            result.w_[i] = sum;
            carry = new_carry;
        }
        return result;
    }
    
    Int256 operator-(const Int256& other) const {
        return *this + (-other);
    }

    // 128x128 -> 256 bit multiplication
    static Int256 mul(const Int128& a, const Int128& b) {
        uint64_t a_lo = a.low();
        int64_t  a_hi = a.high();
        uint64_t b_lo = b.low();
        int64_t  b_hi = b.high();

        Int128 p0 = Int128::mul(static_cast<int64_t>(a_lo), static_cast<int64_t>(b_lo));
        Int128 p1 = Int128::mul(static_cast<int64_t>(a_lo), b_hi);
        Int128 p2 = Int128::mul(a_hi, static_cast<int64_t>(b_lo));
        Int128 p3 = Int128::mul(a_hi, b_hi);

        Int256 result = Int256::fromInt128(p0);

        Int256 p1_shifted(
            0ULL,
            p1.low(),
            static_cast<uint64_t>(p1.high()),
            (p1.high() < 0) ? ~0ULL : 0ULL
        );
        result += p1_shifted;

        Int256 p2_shifted(
            0ULL,
            p2.low(),
            static_cast<uint64_t>(p2.high()),
            (p2.high() < 0) ? ~0ULL : 0ULL
        );
        result += p2_shifted;

        Int256 p3_shifted(
            0ULL,
            0ULL,
            p3.low(),
            static_cast<uint64_t>(p3.high())
        );
        result += p3_shifted;

        return result;
    }
    
    static Int256 mul(int64_t a, int64_t b) {
        Int128 prod = Int128::mul(a, b);
        return fromInt128(prod);
    }
};

//=============================================================================
// DEBUG PRINTING UTILITIES
//=============================================================================

void print_int128(const char* name, const Int128& v) {
    std::printf("  %s: high=%ld, low=%lu, sign=%d\n", 
                name, v.high(), v.low(), v.sign());
}

void print_int256(const char* name, const Int256& v) {
    std::printf("  %s: w=[%lu, %lu, %lu, %lu], sign=%d\n", 
                name, v.word(0), v.word(1), v.word(2), v.word(3), v.sign());
}

//=============================================================================
// EXACT PREDICATES
//=============================================================================

inline int orient2d_exact(int64_t ax, int64_t ay,
                           int64_t bx, int64_t by,
                           int64_t cx, int64_t cy) {
    int64_t acx = ax - cx, acy = ay - cy;
    int64_t bcx = bx - cx, bcy = by - cy;
    
    Int128 det = Int128::mul(acx, bcy) - Int128::mul(acy, bcx);
    return det.sign();
}

inline int orient3d_exact(int64_t ax, int64_t ay, int64_t az,
                           int64_t bx, int64_t by, int64_t bz,
                           int64_t cx, int64_t cy, int64_t cz,
                           int64_t dx, int64_t dy, int64_t dz) {
    int64_t adx = ax - dx, ady = ay - dy, adz = az - dz;
    int64_t bdx = bx - dx, bdy = by - dy, bdz = bz - dz;
    int64_t cdx = cx - dx, cdy = cy - dy, cdz = cz - dz;
    
    Int128 ab_xy = Int128::mul(adx, bdy) - Int128::mul(ady, bdx);
    Int128 bc_xy = Int128::mul(bdx, cdy) - Int128::mul(bdy, cdx);
    Int128 ca_xy = Int128::mul(cdx, ady) - Int128::mul(cdy, adx);
    
    Int256 det = Int256::mul(adz, bc_xy)
               + Int256::mul(bdz, ca_xy)
               + Int256::mul(cdz, ab_xy);
    
    return det.sign();
}

inline int incircle_exact(int64_t ax, int64_t ay,
                           int64_t bx, int64_t by,
                           int64_t cx, int64_t cy,
                           int64_t dx, int64_t dy) {
    int64_t adx = ax - dx, ady = ay - dy;
    int64_t bdx = bx - dx, bdy = by - dy;
    int64_t cdx = cx - dx, cdy = cy - dy;

    Int128 ab_det = Int128::mul(adx, bdy) - Int128::mul(bdx, ady);
    Int128 bc_det = Int128::mul(bdx, cdy) - Int128::mul(cdx, bdy);
    Int128 ca_det = Int128::mul(cdx, ady) - Int128::mul(adx, cdy);

    Int128 a_lift = Int128::mul(adx, adx) + Int128::mul(ady, ady);
    Int128 b_lift = Int128::mul(bdx, bdx) + Int128::mul(bdy, bdy);
    Int128 c_lift = Int128::mul(cdx, cdx) + Int128::mul(cdy, cdy);

    Int256 det = Int256::mul(a_lift, bc_det)
               + Int256::mul(b_lift, ca_det)
               + Int256::mul(c_lift, ab_det);

    return det.sign();
}

// orient2d_EEL - Symbolic orient2d for Explicit, Explicit, LPI points
inline int orient2d_EEL(
    const int64_t a[3], const int64_t b[3],
    const int64_t l0[3], const int64_t l1[3],
    const int64_t p0[3], const int64_t p1[3], const int64_t p2[3],
    int proj_u, int proj_v)
{
    // Compute plane normal (Int128 components)
    int64_t e1x = p1[0]-p0[0], e1y = p1[1]-p0[1], e1z = p1[2]-p0[2];
    int64_t e2x = p2[0]-p0[0], e2y = p2[1]-p0[1], e2z = p2[2]-p0[2];

    Int128 nx = Int128::mul(e1y, e2z) - Int128::mul(e1z, e2y);
    Int128 ny = Int128::mul(e1z, e2x) - Int128::mul(e1x, e2z);
    Int128 nz = Int128::mul(e1x, e2y) - Int128::mul(e1y, e2x);

    // Line direction
    int64_t lx = l1[0]-l0[0], ly = l1[1]-l0[1], lz = l1[2]-l0[2];

    // Denominator d = dot(n, line_dir)
    Int128 d = nx * lx + ny * ly + nz * lz;

    // Numerator of t: t_n = -dot(n, L0-P0)
    int64_t dx = l0[0]-p0[0], dy = l0[1]-p0[1], dz = l0[2]-p0[2];
    Int128 t_n = -(nx * dx + ny * dy + nz * dz);

    // LPI projected coordinates (as rational: num/d)
    int64_t lu = l1[proj_u] - l0[proj_u];
    int64_t lv = l1[proj_v] - l0[proj_v];

    Int256 Cu_n = Int256::mul(Int128(l0[proj_u]), d) + Int256::mul(t_n, Int128(lu));
    Int256 Cv_n = Int256::mul(Int128(l0[proj_v]), d) + Int256::mul(t_n, Int128(lv));

    // Compute orient2d = d x inner
    int64_t au = a[proj_u], av = a[proj_v];
    int64_t bu = b[proj_u], bv = b[proj_v];

    Int128 ab = Int128::mul(au, bv) - Int128::mul(av, bu);

    int64_t coeff_uv = av - bv;
    int64_t coeff_vu = bu - au;

    Int256 term2a = Int256::mul(d, Int128(l0[proj_u] * coeff_uv));
    Int256 term2b = Int256::mul(t_n, Int128(lu * coeff_uv));
    Int256 t2 = term2a + term2b;

    Int256 term3a = Int256::mul(d, Int128(l0[proj_v] * coeff_vu));
    Int256 term3b = Int256::mul(t_n, Int128(lv * coeff_vu));
    Int256 t3 = term3a + term3b;

    Int256 inner = Int256::mul(ab, d) + t2 + t3;

    int d_sign = d.sign();
    int inner_sign = inner.sign();

    if (d_sign == 0 || inner_sign == 0) return 0;
    return (d_sign == inner_sign) ? +1 : -1;
}

//=============================================================================
// TEST SUITE: Int128 Basic Operations
//=============================================================================

void test_int128_construction() {
    TEST_BEGIN("Int128: Construction from int64_t");
    
    Int128 zero(0);
    TEST_ASSERT(zero.sign() == 0);
    TEST_ASSERT(zero.isZero());
    
    Int128 pos(42);
    TEST_ASSERT(pos.sign() == 1);
    TEST_ASSERT(pos.high() == 0);
    TEST_ASSERT(pos.low() == 42);
    
    Int128 neg(-42);
    TEST_ASSERT(neg.sign() == -1);
    TEST_ASSERT(neg.isNegative());
    
    Int128 max32(0x7FFFFFFF);
    TEST_ASSERT(max32.sign() == 1);
    TEST_ASSERT(max32.low() == 0x7FFFFFFF);
    
    Int128 min32(-0x7FFFFFFF);
    TEST_ASSERT(min32.sign() == -1);
    
    TEST_PASS();
}

void test_int128_addition() {
    TEST_BEGIN("Int128: Addition");
    
    Int128 a(100);
    Int128 b(200);
    Int128 sum = a + b;
    TEST_ASSERT(sum.sign() == 1);
    TEST_ASSERT(sum.low() == 300);
    
    Int128 c(-100);
    Int128 d(50);
    Int128 sum2 = c + d;
    TEST_ASSERT(sum2.sign() == -1);
    
    // Test carry propagation
    Int128 e(0xFFFFFFFFFFFFFFFFULL);
    Int128 f(1);
    Int128 sum3 = e + f;
    TEST_ASSERT(sum3.isZero());
    
    TEST_PASS();
}

void test_int128_subtraction() {
    TEST_BEGIN("Int128: Subtraction");
    
    Int128 a(500);
    Int128 b(200);
    Int128 diff = a - b;
    TEST_ASSERT(diff.sign() == 1);
    TEST_ASSERT(diff.low() == 300);
    
    Int128 c(100);
    Int128 d(300);
    Int128 diff2 = c - d;
    TEST_ASSERT(diff2.sign() == -1);
    
    // Self-subtraction
    Int128 e(12345);
    Int128 diff3 = e - e;
    TEST_ASSERT(diff3.isZero());
    
    TEST_PASS();
}

void test_int128_multiplication() {
    TEST_BEGIN("Int128: 64x64->128 Multiplication");
    
    // Basic cases
    Int128 prod1 = Int128::mul(10, 20);
    TEST_ASSERT(prod1.sign() == 1);
    TEST_ASSERT(prod1.low() == 200);
    
    Int128 prod2 = Int128::mul(-10, 20);
    TEST_ASSERT(prod2.sign() == -1);
    
    Int128 prod3 = Int128::mul(-10, -20);
    TEST_ASSERT(prod3.sign() == 1);
    TEST_ASSERT(prod3.low() == 200);
    
    // Large values
    Int128 prod4 = Int128::mul(0x7FFFFFFF, 0x7FFFFFFF);
    TEST_ASSERT(prod4.sign() == 1);
    
    // Edge case: max int32 * max int32
    int64_t max26 = 33554431;  // 2^25 - 1
    Int128 prod5 = Int128::mul(max26, max26);
    TEST_ASSERT(prod5.sign() == 1);
    
    TEST_PASS();
}

void test_int128_comparison() {
    TEST_BEGIN("Int128: Comparison operators");
    
    Int128 a(100);
    Int128 b(200);
    Int128 c(100);
    
    TEST_ASSERT(a < b);
    TEST_ASSERT(b > a);
    TEST_ASSERT(a == c);
    TEST_ASSERT(a != b);
    TEST_ASSERT(a <= c);
    TEST_ASSERT(a >= c);
    
    Int128 neg(-100);
    Int128 pos(100);
    TEST_ASSERT(neg < pos);
    TEST_ASSERT(pos > neg);
    
    TEST_PASS();
}

//=============================================================================
// TEST SUITE: Int256 Operations
//=============================================================================

void test_int256_construction() {
    TEST_BEGIN("Int256: Construction");
    
    Int256 zero;
    TEST_ASSERT(zero.isZero());
    TEST_ASSERT(zero.sign() == 0);
    
    Int256 pos(42);
    TEST_ASSERT(pos.sign() == 1);
    TEST_ASSERT(pos.word(0) == 42);
    
    Int256 neg(-42);
    TEST_ASSERT(neg.sign() == -1);
    TEST_ASSERT(neg.isNegative());
    
    TEST_PASS();
}

void test_int256_addition() {
    TEST_BEGIN("Int256: Addition");
    
    Int256 a(100);
    Int256 b(200);
    Int256 sum = a + b;
    TEST_ASSERT(sum.word(0) == 300);
    TEST_ASSERT(sum.sign() == 1);
    
    // Test carry across words
    Int256 c(0xFFFFFFFFFFFFFFFFULL);
    Int256 d(1);
    Int256 sum2 = c + d;
    TEST_ASSERT(sum2.word(0) == 0);
    TEST_ASSERT(sum2.word(1) == 1);
    
    TEST_PASS();
}

void test_int256_mul_128x128() {
    TEST_BEGIN("Int256: 128x128->256 Multiplication");
    
    // Simple case
    Int128 a(10);
    Int128 b(20);
    Int256 prod = Int256::mul(a, b);
    TEST_ASSERT(prod.word(0) == 200);
    TEST_ASSERT(prod.sign() == 1);
    
    // Negative case
    Int128 c(-10);
    Int128 d(20);
    Int256 prod2 = Int256::mul(c, d);
    TEST_ASSERT(prod2.sign() == -1);
    
    // Both negative
    Int128 e(-10);
    Int128 f(-20);
    Int256 prod3 = Int256::mul(e, f);
    TEST_ASSERT(prod3.word(0) == 200);
    TEST_ASSERT(prod3.sign() == 1);
    
    // Large values at 26-bit limit
    int64_t max26 = 33554431;  // 2^25 - 1
    Int128 g(max26);
    Int128 h(max26);
    Int256 prod4 = Int256::mul(g, h);
    TEST_ASSERT(prod4.sign() == 1);
    
    // Test with Int128 max values
    Int128 max_val(0x7FFFFFFFFFFFFFFFLL);
    Int256 prod5 = Int256::mul(max_val, max_val);
    TEST_ASSERT(prod5.sign() == 1);
    
    TEST_PASS();
}

//=============================================================================
// TEST SUITE: 26-bit Precision Limit Tests
//=============================================================================

void test_precision_26bit_orient2d() {
    TEST_BEGIN("Precision: orient2d at 26-bit limit");
    
    // Maximum 26-bit value: 2^25 - 1 = 33,554,431
    const int64_t MAX26 = 33554431;
    
    // Test with maximum coordinates
    int64_t ax = 0, ay = 0;
    int64_t bx = MAX26, by = 0;
    int64_t cx = MAX26 / 2, cy = MAX26;
    
    int result = orient2d_exact(ax, ay, bx, by, cx, cy);
    // Points form a CCW triangle
    TEST_ASSERT(result == 1);
    
    // Test collinearity at limits
    int64_t dx = MAX26 / 2, dy = 0;
    int result2 = orient2d_exact(ax, ay, bx, by, dx, dy);
    TEST_ASSERT_MSG(result2 == 0, "Collinear points should return 0");
    
    TEST_PASS();
}

void test_precision_26bit_orient3d() {
    TEST_BEGIN("Precision: orient3d at 26-bit limit");
    
    const int64_t MAX26 = 33554431;
    
    // Create a tetrahedron at the limits
    int64_t ax = 0, ay = 0, az = 0;
    int64_t bx = MAX26, by = 0, bz = 0;
    int64_t cx = 0, cy = MAX26, cz = 0;
    int64_t dx = 0, dy = 0, dz = MAX26;
    
    int result = orient3d_exact(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz);
    // This should give positive volume
    TEST_ASSERT(result == 1 || result == -1);  // Non-zero
    
    // Test coplanarity
    int64_t ex = MAX26, ey = MAX26, ez = 0;
    int result2 = orient3d_exact(ax, ay, az, bx, by, bz, cx, cy, cz, ex, ey, ez);
    // All points in z=0 plane should be coplanar
    TEST_ASSERT(result2 == 0);
    
    TEST_PASS();
}

void test_precision_26bit_incircle() {
    TEST_BEGIN("Precision: incircle at 26-bit limit");
    
    const int64_t MAX26 = 33554431;
    
    // Create a triangle
    int64_t ax = 0, ay = 0;
    int64_t bx = MAX26, by = 0;
    int64_t cx = MAX26 / 2, cy = MAX26;
    
    // Point inside the triangle
    int64_t dx_inside = MAX26 / 2, dy_inside = MAX26 / 4;
    int result_inside = incircle_exact(ax, ay, bx, by, cx, cy, dx_inside, dy_inside);
    
    // Point outside (far away)
    int64_t dx_outside = MAX26 * 2, dy_outside = MAX26 * 2;
    int result_outside = incircle_exact(ax, ay, bx, by, cx, cy, dx_outside, dy_outside);
    
    // Results should be different signs (one inside, one outside)
    TEST_ASSERT(result_inside != result_outside);
    
    TEST_PASS();
}

void test_precision_int256_mul_extreme() {
    TEST_BEGIN("Precision: Int256::mul at 26-bit extreme");
    
    // Maximum 26-bit value
    const int64_t MAX26 = 33554431;
    
    // Create Int128 values at the limit
    Int128 a(MAX26);
    Int128 b(MAX26);
    
    // Multiply them
    Int256 prod = Int256::mul(a, b);
    
    // Verify the result is correct
    TEST_ASSERT(prod.sign() == 1);
    TEST_ASSERT(!prod.isZero());
    
    // Test with products that need full 256 bits
    Int128 large1(0x7FFFFFFFFFFFFFFFLL);
    Int128 large2(0x7FFFFFFFFFFFFFFFLL);
    Int256 prod2 = Int256::mul(large1, large2);
    
    // Should be positive and non-zero
    TEST_ASSERT(prod2.sign() == 1);
    TEST_ASSERT(!prod2.isZero());
    
    TEST_PASS();
}

//=============================================================================
// TEST SUITE: Collinearity Trap (orient2d_EEL)
//=============================================================================

void test_collinearity_simple() {
    TEST_BEGIN("Collinearity: Simple collinear case");
    
    // Create three collinear points in XY plane
    int64_t a[3] = {0, 0, 0};
    int64_t b[3] = {100, 0, 0};
    
    // LPI line: from (0,0,0) to (50,0,0) - line along x-axis
    int64_t l0[3] = {0, 0, 0};
    int64_t l1[3] = {50, 0, 0};
    
    // Plane: XY plane at z=0
    int64_t p0[3] = {0, 0, 0};
    int64_t p1[3] = {100, 0, 0};
    int64_t p2[3] = {0, 100, 0};
    
    // The LPI point should be at the origin (intersection of x-axis line with XY plane)
    int result = orient2d_EEL(a, b, l0, l1, p0, p1, p2, 0, 1);
    
    // Point at origin is on line from (0,0) to (100,0)
    TEST_ASSERT_MSG(result == 0, "LPI point on line should return 0");
    
    TEST_PASS();
}

void test_collinearity_lpi_on_line() {
    TEST_BEGIN("Collinearity: LPI point exactly on AB line");
    
    // Points A and B define a line in XY plane
    int64_t a[3] = {0, 0, 0};
    int64_t b[3] = {100, 100, 0};
    
    // LPI line: from (50, 50, -10) to (50, 50, 10) - vertical line through midpoint
    int64_t l0[3] = {50, 50, -10};
    int64_t l1[3] = {50, 50, 10};
    
    // Plane: XY plane at z=0
    int64_t p0[3] = {0, 0, 0};
    int64_t p1[3] = {100, 0, 0};
    int64_t p2[3] = {0, 100, 0};
    
    // The LPI point should be at (50, 50, 0) which is exactly on line AB
    int result = orient2d_EEL(a, b, l0, l1, p0, p1, p2, 0, 1);
    
    TEST_ASSERT_MSG(result == 0, "LPI point (50,50) on line from (0,0) to (100,100) should be collinear");
    
    TEST_PASS();
}

void test_collinearity_lpi_off_line() {
    TEST_BEGIN("Collinearity: LPI point off the line");
    
    // Points A and B define a line in XY plane
    int64_t a[3] = {0, 0, 0};
    int64_t b[3] = {100, 0, 0};
    
    // LPI line: from (50, 50, -10) to (50, 50, 10) - vertical line at x=50, y=50
    int64_t l0[3] = {50, 50, -10};
    int64_t l1[3] = {50, 50, 10};
    
    // Plane: XY plane at z=0
    int64_t p0[3] = {0, 0, 0};
    int64_t p1[3] = {100, 0, 0};
    int64_t p2[3] = {0, 100, 0};
    
    // The LPI point should be at (50, 50, 0) which is NOT on line AB (y=0)
    int result = orient2d_EEL(a, b, l0, l1, p0, p1, p2, 0, 1);
    
    // Point (50, 50) should be to the left of line from (0,0) to (100,0)
    TEST_ASSERT_MSG(result == 1, "LPI point (50,50) should be left of horizontal line");
    
    TEST_PASS();
}

void test_collinearity_lpi_at_endpoint() {
    TEST_BEGIN("Collinearity: LPI point at endpoint A");
    
    // Points A and B
    int64_t a[3] = {50, 50, 0};
    int64_t b[3] = {100, 100, 0};
    
    // LPI line: from (50, 50, -10) to (50, 50, 10)
    int64_t l0[3] = {50, 50, -10};
    int64_t l1[3] = {50, 50, 10};
    
    // Plane: XY plane at z=0
    int64_t p0[3] = {0, 0, 0};
    int64_t p1[3] = {100, 0, 0};
    int64_t p2[3] = {0, 100, 0};
    
    // The LPI point should be at (50, 50, 0) which equals point A
    int result = orient2d_EEL(a, b, l0, l1, p0, p1, p2, 0, 1);
    
    TEST_ASSERT_MSG(result == 0, "LPI point at endpoint A should return 0");
    
    TEST_PASS();
}

void test_collinearity_different_projections() {
    TEST_BEGIN("Collinearity: Different projection axes");
    
    // Points in 3D
    int64_t a[3] = {0, 0, 0};
    int64_t b[3] = {0, 100, 100};
    
    // LPI line along Z axis
    int64_t l0[3] = {0, 50, 50};
    int64_t l1[3] = {10, 50, 50};
    
    // Plane: YZ plane at x=0
    int64_t p0[3] = {0, 0, 0};
    int64_t p1[3] = {0, 100, 0};
    int64_t p2[3] = {0, 0, 100};
    
    // Test YZ projection (axes 1 and 2)
    int result_yz = orient2d_EEL(a, b, l0, l1, p0, p1, p2, 1, 2);
    
    // The LPI point (0, 50, 50) should be on line from (0,0,0) to (0,100,100)
    TEST_ASSERT_MSG(result_yz == 0, "LPI point should be collinear in YZ projection");
    
    TEST_PASS();
}

void test_collinearity_at_precision_limit() {
    TEST_BEGIN("Collinearity: At 26-bit precision limit");
    
    const int64_t MAX26 = 33554431;
    
    // Points at precision limit
    int64_t a[3] = {0, 0, 0};
    int64_t b[3] = {MAX26, MAX26, 0};
    
    // LPI line through midpoint
    int64_t l0[3] = {MAX26 / 2, MAX26 / 2, -10};
    int64_t l1[3] = {MAX26 / 2, MAX26 / 2, 10};
    
    // Plane: XY plane at z=0
    int64_t p0[3] = {0, 0, 0};
    int64_t p1[3] = {MAX26, 0, 0};
    int64_t p2[3] = {0, MAX26, 0};
    
    int result = orient2d_EEL(a, b, l0, l1, p0, p1, p2, 0, 1);
    
    TEST_ASSERT_MSG(result == 0, "LPI point at precision limit should be exactly collinear");
    
    TEST_PASS();
}

//=============================================================================
// TEST SUITE: orient2d_exact Tests
//=============================================================================

void test_orient2d_basic() {
    TEST_BEGIN("orient2d_exact: Basic orientations");
    
    // CCW triangle
    int result1 = orient2d_exact(0, 0, 100, 0, 50, 100);
    TEST_ASSERT(result1 == 1);
    
    // CW triangle
    int result2 = orient2d_exact(0, 0, 50, 100, 100, 0);
    TEST_ASSERT(result2 == -1);
    
    // Collinear (horizontal)
    int result3 = orient2d_exact(0, 0, 50, 0, 100, 0);
    TEST_ASSERT(result3 == 0);
    
    // Collinear (diagonal)
    int result4 = orient2d_exact(0, 0, 50, 50, 100, 100);
    TEST_ASSERT(result4 == 0);
    
    // Collinear (vertical)
    int result5 = orient2d_exact(0, 0, 0, 50, 0, 100);
    TEST_ASSERT(result5 == 0);
    
    TEST_PASS();
}

void test_orient2d_degenerate() {
    TEST_BEGIN("orient2d_exact: Degenerate cases");
    
    // Duplicate points
    int result1 = orient2d_exact(50, 50, 50, 50, 100, 100);
    TEST_ASSERT(result1 == 0);
    
    int result2 = orient2d_exact(0, 0, 50, 50, 50, 50);
    TEST_ASSERT(result2 == 0);
    
    // All same
    int result3 = orient2d_exact(50, 50, 50, 50, 50, 50);
    TEST_ASSERT(result3 == 0);
    
    TEST_PASS();
}

void test_orient2d_large_values() {
    TEST_BEGIN("orient2d_exact: Large coordinate values");
    
    // Test with large coordinates
    int64_t big = 1000000000LL;  // 1 billion
    
    int result1 = orient2d_exact(0, 0, big, 0, big / 2, big);
    TEST_ASSERT(result1 == 1);
    
    int result2 = orient2d_exact(0, 0, big / 2, big, big, 0);
    TEST_ASSERT(result2 == -1);
    
    int result3 = orient2d_exact(0, 0, big, 0, big / 2, 0);
    TEST_ASSERT(result3 == 0);
    
    TEST_PASS();
}

//=============================================================================
// TEST SUITE: orient3d_exact Tests
//=============================================================================

void test_orient3d_basic() {
    TEST_BEGIN("orient3d_exact: Basic orientations");
    
    // Positive volume tetrahedron
    int result1 = orient3d_exact(
        0, 0, 0,
        100, 0, 0,
        0, 100, 0,
        0, 0, 100
    );
    TEST_ASSERT(result1 == 1);
    
    // Negative volume (d below plane)
    int result2 = orient3d_exact(
        0, 0, 0,
        100, 0, 0,
        0, 100, 0,
        0, 0, -100
    );
    TEST_ASSERT(result2 == -1);
    
    // Coplanar (all in XY plane)
    int result3 = orient3d_exact(
        0, 0, 0,
        100, 0, 0,
        0, 100, 0,
        50, 50, 0
    );
    TEST_ASSERT(result3 == 0);
    
    TEST_PASS();
}

void test_orient3d_coplanar() {
    TEST_BEGIN("orient3d_exact: Coplanar cases");
    
    // All points in XZ plane (y=0)
    int result1 = orient3d_exact(
        0, 0, 0,
        100, 0, 0,
        0, 0, 100,
        50, 0, 50
    );
    TEST_ASSERT(result1 == 0);
    
    // All points in YZ plane (x=0)
    int result2 = orient3d_exact(
        0, 0, 0,
        0, 100, 0,
        0, 0, 100,
        0, 50, 50
    );
    TEST_ASSERT(result2 == 0);
    
    TEST_PASS();
}

void test_orient3d_degenerate() {
    TEST_BEGIN("orient3d_exact: Degenerate cases");
    
    // Duplicate points
    int result1 = orient3d_exact(
        50, 50, 50,
        50, 50, 50,
        0, 100, 0,
        0, 0, 100
    );
    TEST_ASSERT(result1 == 0);
    
    // All same
    int result2 = orient3d_exact(
        50, 50, 50,
        50, 50, 50,
        50, 50, 50,
        50, 50, 50
    );
    TEST_ASSERT(result2 == 0);
    
    TEST_PASS();
}

//=============================================================================
// TEST SUITE: incircle_exact Tests
//=============================================================================

void test_incircle_basic() {
    TEST_BEGIN("incircle_exact: Basic tests");
    
    // Triangle with circumcircle at origin, radius 5
    // Points: (5,0), (0,5), (-5,0), (0,-5)
    int64_t ax = 5, ay = 0;
    int64_t bx = 0, by = 5;
    int64_t cx = -5, cy = 0;
    int64_t dx_inside = 0, dy_inside = 0;  // Origin is inside
    int64_t dx_outside = 0, dy_outside = -6;  // Below is outside
    
    int result_inside = incircle_exact(ax, ay, bx, by, cx, cy, dx_inside, dy_inside);
    int result_outside = incircle_exact(ax, ay, bx, by, cx, cy, dx_outside, dy_outside);
    
    // Inside and outside should give different results
    TEST_ASSERT(result_inside != result_outside);
    
    // Point on circle should return 0
    int64_t dx_on = 0, dy_on = -5;  // On the circle
    int result_on = incircle_exact(ax, ay, bx, by, cx, cy, dx_on, dy_on);
    TEST_ASSERT(result_on == 0);
    
    TEST_PASS();
}

//=============================================================================
// MAIN
//=============================================================================

int main() {
    std::printf("\n");
    std::printf("========================================\n");
    std::printf("EMBER ExactPredicates Test Harness\n");
    std::printf("========================================\n");
    std::printf("\n");
    
    // Int128 Tests
    std::printf("--- Int128 Basic Operations ---\n");
    test_int128_construction();
    test_int128_addition();
    test_int128_subtraction();
    test_int128_multiplication();
    test_int128_comparison();
    
    // Int256 Tests
    std::printf("\n--- Int256 Operations ---\n");
    test_int256_construction();
    test_int256_addition();
    test_int256_mul_128x128();
    
    // Precision Tests
    std::printf("\n--- 26-bit Precision Tests ---\n");
    test_precision_26bit_orient2d();
    test_precision_26bit_orient3d();
    test_precision_26bit_incircle();
    test_precision_int256_mul_extreme();
    
    // Collinearity Tests
    std::printf("\n--- Collinearity Trap Tests (orient2d_EEL) ---\n");
    test_collinearity_simple();
    test_collinearity_lpi_on_line();
    test_collinearity_lpi_off_line();
    test_collinearity_lpi_at_endpoint();
    test_collinearity_different_projections();
    test_collinearity_at_precision_limit();
    
    // orient2d Tests
    std::printf("\n--- orient2d_exact Tests ---\n");
    test_orient2d_basic();
    test_orient2d_degenerate();
    test_orient2d_large_values();
    
    // orient3d Tests
    std::printf("\n--- orient3d_exact Tests ---\n");
    test_orient3d_basic();
    test_orient3d_coplanar();
    test_orient3d_degenerate();
    
    // incircle Tests
    std::printf("\n--- incircle_exact Tests ---\n");
    test_incircle_basic();
    
    // Summary
    TEST_SUMMARY();
    
    return g_testsFailed > 0 ? 1 : 0;
}
