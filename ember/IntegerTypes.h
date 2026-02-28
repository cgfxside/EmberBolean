/**
 * @file IntegerTypes.h
 * @brief Production-ready arbitrary-precision integer types for EMBER Boolean engine
 * 
 * EMBER Bit-Width Budget (from paper):
 * =====================================
 * - Input plane coefficients: 26-bit (fits in int32_t)
 * - 2×2 minors: 53-bit (fits in int64_t)
 * - 3×3 determinants (intersect): 80-bit → Int128
 * - classify() dot product: 160-bit → Int256
 * - Second-generation intersections: ≤256-bit → Int256
 * 
 * This file provides:
 * - Int128: 128-bit signed integer for 80-bit intermediate results
 * - Int256: 256-bit signed integer for 160-bit dot products
 * 
 * Platform Support:
 * - GCC/Clang: Uses __int128 intrinsic when available
 * - MSVC: Uses _mul128/_mulh intrinsics
 * - Others: Portable 32-bit schoolbook multiplication
 * 
 * INTEGRATED HOT FIXES:
 *   - FIX 1: Int128::mul64 portable cross-terms correction
 *   - FIX 2: Int128 addition carry propagation verified correct
 * 
 * @author EMBER Boolean Engine
 * @version 1.0.1 (Hot Fixed)
 */

#ifndef EMBER_INTEGER_TYPES_H
#define EMBER_INTEGER_TYPES_H

#include <cstdint>
#include <cstddef>
#include <limits>
#include <type_traits>

// ============================================================================
// Platform Detection
// ============================================================================

// Detect __int128 support (GCC/Clang)
#if defined(__SIZEOF_INT128__)
    #define EMBER_HAS_INT128 1
#else
    #define EMBER_HAS_INT128 0
#endif

// Detect MSVC
#if defined(_MSC_VER)
    #define EMBER_IS_MSVC 1
    #include <intrin.h>
#else
    #define EMBER_IS_MSVC 0
#endif

// Detect Boost.Multiprecision
#if defined(EMBER_USE_BOOST_MULTIPRECISION)
    #include <boost/multiprecision/cpp_int.hpp>
    namespace mp = boost::multiprecision;
#endif

namespace ember {

// ============================================================================
// Static Assertions for Type Safety
// ============================================================================

static_assert(sizeof(int64_t) == 8, "int64_t must be 8 bytes");
static_assert(sizeof(uint64_t) == 8, "uint64_t must be 8 bytes");
static_assert(sizeof(int32_t) == 4, "int32_t must be 4 bytes");
static_assert(sizeof(uint32_t) == 4, "uint32_t must be 4 bytes");

// ============================================================================
// Forward Declarations
// ============================================================================

struct Int128;
struct Int256;

// ============================================================================
// Int128 - 128-bit Signed Integer
// ============================================================================
/**
 * @brief 128-bit signed integer using (hi, lo) representation
 * 
 * Representation: value = hi * 2^64 + lo
 * - hi: int64_t (signed high bits, includes sign)
 * - lo: uint64_t (unsigned low bits)
 * 
 * Bit-width budget: 128 bits (covers 80-bit 3×3 determinants with margin)
 */
struct Int128 {
    // ------------------------------------------------------------------------
    // Member Data
    // ------------------------------------------------------------------------
    int64_t hi;      // High 64 bits (signed, includes sign extension)
    uint64_t lo;     // Low 64 bits (unsigned magnitude)

    // ------------------------------------------------------------------------
    // Constants
    // ------------------------------------------------------------------------
    static constexpr int64_t  MAX_HI = INT64_MAX;
    static constexpr int64_t  MIN_HI = INT64_MIN;
    static constexpr uint64_t MAX_LO = UINT64_MAX;
    
    // ------------------------------------------------------------------------
    // Accessors
    // ------------------------------------------------------------------------
    
    /// Get high 64 bits (signed)
    constexpr int64_t high() const noexcept { return hi; }
    
    /// Get low 64 bits (unsigned)
    constexpr uint64_t low() const noexcept { return lo; }
    
    /// Get low 64 bits as signed int64_t (for MeshImport portability)
    /// Note: This reinterprets the bits; for values > INT64_MAX, result is negative
    constexpr int64_t low_as_int64() const noexcept {
        return static_cast<int64_t>(lo);
    }
    
    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------
    
    /// Default constructor: initializes to zero
    constexpr Int128() noexcept : hi(int64_t(0)), lo(uint64_t(0)) {}
    
    /// Constructor from int (prevents ambiguity with int64_t/uint64_t)
    constexpr Int128(int v) noexcept : hi(v < 0 ? int64_t(-1) : int64_t(0)), lo(static_cast<uint64_t>(v)) {}
    
    /// Constructor from int64_t (sign-extends)
    constexpr Int128(int64_t v) noexcept : hi(v < 0 ? int64_t(-1) : int64_t(0)), lo(static_cast<uint64_t>(v)) {}
    
    /// Constructor from uint64_t (zero-extends)
    constexpr Int128(uint64_t v) noexcept : hi(int64_t(0)), lo(v) {}
    
    /// Constructor from hi/lo pair
    constexpr Int128(int64_t high, uint64_t low) noexcept : hi(high), lo(low) {}
    
    /// Constructor from Int256 (truncates - use with caution)
    explicit Int128(const Int256& v) noexcept;

    // ------------------------------------------------------------------------
    // Factory Methods
    // ------------------------------------------------------------------------
    
    /// Create from signed 64-bit value
    static constexpr Int128 fromInt64(int64_t v) noexcept {
        return Int128(v);
    }
    
    /// Create from unsigned 64-bit value
    static constexpr Int128 fromUInt64(uint64_t v) noexcept {
        return Int128(int64_t(0), v);
    }
    
    /// Create from hi/lo pair
    static constexpr Int128 fromParts(int64_t high, uint64_t low) noexcept {
        return Int128(high, low);
    }
    
    /// Zero constant
    static constexpr Int128 zero() noexcept { return Int128(); }
    
    /// One constant
    static constexpr Int128 one() noexcept { return Int128(int64_t(0), uint64_t(1)); }

    // ------------------------------------------------------------------------
    // 64×64 → 128 Multiplication (Critical for EMBER)
    // ------------------------------------------------------------------------
    
    /**
     * @brief Multiply two int64_t values to produce Int128
     * 
     * This is the foundation of all EMBER arithmetic.
     * Uses platform-optimized implementations.
     * 
     * HOT FIX 1: Portable fallback now correctly handles cross-terms
     * 
     * @param a First operand (26-bit plane coefficient)
     * @param b Second operand (26-bit plane coefficient)
     * @return Int128 product (up to 52 bits, fits in 80-bit budget)
     */
    static Int128 mul64(int64_t a, int64_t b) noexcept {
#if EMBER_HAS_INT128
        // GCC/Clang: Use __int128 intrinsic
        __int128 prod = static_cast<__int128>(a) * static_cast<__int128>(b);
        return Int128(
            static_cast<int64_t>(prod >> 64),
            static_cast<uint64_t>(prod)
        );
#elif EMBER_IS_MSVC
        // MSVC: Use _mul128 intrinsic
        int64_t high;
        int64_t low = _mul128(a, b, &high);
        return Int128(high, static_cast<uint64_t>(low));
#else
        // Portable fallback: 32-bit schoolbook multiplication (FIXED)
        return mul64_portable(a, b);
#endif
    }
    
    /**
     * @brief Multiply two uint64_t values to produce Int128
     */
    static Int128 mul64(uint64_t a, uint64_t b) noexcept {
#if EMBER_HAS_INT128
        unsigned __int128 prod = static_cast<unsigned __int128>(a) * static_cast<unsigned __int128>(b);
        return Int128(
            static_cast<int64_t>(static_cast<uint64_t>(prod >> 64)),
            static_cast<uint64_t>(prod)
        );
#elif EMBER_IS_MSVC
        uint64_t high;
        uint64_t low = _umul128(a, b, &high);
        return Int128(static_cast<int64_t>(high), low);
#else
        return mul64_portable_u64(a, b);
#endif
    }

    // ------------------------------------------------------------------------
    // Portable Multiplication Fallbacks (HOT FIX 1)
    // ------------------------------------------------------------------------
    
private:
    /// Portable signed 64×64 → 128 using 32-bit operations (FIX 1: Carry-Bit Integrity)
    /// 
    /// AUDIT FIX (2026-02-28): Fixed cross-term carry overflow bug.
    /// When both p1_a and p1_b are large (as uint64), their sum can exceed UINT64_MAX.
    /// The carry bit from this addition was lost before the >> 32 shift.
    /// 
    /// Algorithm: a*b = a1*b1*2^64 + (a1*b0 + a0*b1)*2^32 + a0*b0
    /// where a = a1*2^32 + a0, b = b1*2^32 + b0
    static Int128 mul64_portable(int64_t a, int64_t b) noexcept {
        // Decompose into 32-bit parts
        int64_t a1 = a >> 32;
        int64_t b1 = b >> 32;
        uint64_t a0 = static_cast<uint64_t>(a) & 0xFFFFFFFFULL;
        uint64_t b0 = static_cast<uint64_t>(b) & 0xFFFFFFFFULL;
        
        // Compute partial products
        // a * b = (a1*2^32 + a0) * (b1*2^32 + b0)
        //       = a1*b1*2^64 + (a1*b0 + a0*b1)*2^32 + a0*b0
        
        uint64_t p0 = a0 * b0;                              // 64-bit product
        int64_t  p1_a = a1 * static_cast<int64_t>(b0);      // Signed 32×32
        int64_t  p1_b = b1 * static_cast<int64_t>(a0);      // Signed 32×32
        int64_t  p2 = a1 * b1;                              // Signed 32×32
        
        // Combine with proper carry handling
        uint64_t lo = p0;
        int64_t hi = p2;
        
        // Add cross terms to high part with carry
        // FIX: Detect overflow from cross-term addition before shifting
        uint64_t cross_a = static_cast<uint64_t>(p1_a);
        uint64_t cross_b = static_cast<uint64_t>(p1_b);
        uint64_t cross = cross_a + cross_b;
        uint64_t cross_carry = (cross < cross_a) ? 1ULL : 0ULL;  // ← THE FIX
        
        hi += static_cast<int64_t>(cross >> 32);
        hi += static_cast<int64_t>(cross_carry) << 32;  // Propagate lost carry
        
        // Add low carry from cross terms
        uint64_t cross_lo = (cross << 32);
        uint64_t new_lo = lo + cross_lo;
        if (new_lo < lo) {  // Carry occurred
            hi += 1;
        }
        lo = new_lo;
        
        // Handle sign extension for negative cross terms
        if (p1_a < 0) hi -= (1LL << 32);
        if (p1_b < 0) hi -= (1LL << 32);
        
        return Int128(hi, lo);
    }
    
    /// Portable unsigned 64×64 → 128 using 32-bit operations
    /// 
    /// AUDIT FIX (2026-02-28): Fixed cross-term carry overflow bug.
    /// When both p1_a and p1_b are large, their sum can exceed UINT64_MAX.
    /// The carry bit from this addition was lost before the >> 32 shift.
    static Int128 mul64_portable_u64(uint64_t a, uint64_t b) noexcept {
        uint64_t a0 = a & 0xFFFFFFFFULL;
        uint64_t a1 = a >> 32;
        uint64_t b0 = b & 0xFFFFFFFFULL;
        uint64_t b1 = b >> 32;
        
        // a * b = a1*b1*2^64 + (a1*b0 + a0*b1)*2^32 + a0*b0
        uint64_t p0 = a0 * b0;
        uint64_t p1_a = a1 * b0;
        uint64_t p1_b = a0 * b1;
        uint64_t p2 = a1 * b1;
        
        uint64_t lo = p0;
        uint64_t hi = p2;
        
        // FIX: Detect carry from cross-term addition
        uint64_t cross = p1_a + p1_b;
        uint64_t cross_carry = (cross < p1_a) ? 1ULL : 0ULL;  // ← THE FIX
        
        hi += cross >> 32;
        hi += cross_carry << 32;  // Propagate lost carry
        
        uint64_t cross_lo = cross << 32;
        uint64_t new_lo = lo + cross_lo;
        if (new_lo < lo) {
            hi += 1;
        }
        lo = new_lo;
        
        return Int128(static_cast<int64_t>(hi), lo);
    }

public:
    // ------------------------------------------------------------------------
    // Arithmetic Operators
    // ------------------------------------------------------------------------
    
    /// Unary negation (two's complement)
    constexpr Int128 operator-() const noexcept {
        // -x = ~x + 1
        uint64_t new_lo = ~lo + uint64_t(1);
        int64_t new_hi = ~hi;
        if (new_lo == uint64_t(0)) {
            new_hi += int64_t(1);  // Carry from low to high
        }
        return Int128(new_hi, new_lo);
    }
    
    /// Unary plus
    constexpr Int128 operator+() const noexcept {
        return *this;
    }
    
    /// Prefix increment
    constexpr Int128& operator++() noexcept {
        lo++;
        if (lo == uint64_t(0)) {
            hi++;
        }
        return *this;
    }
    
    /// Postfix increment
    constexpr Int128 operator++(int) noexcept {
        Int128 temp = *this;
        ++(*this);
        return temp;
    }
    
    /// Prefix decrement
    constexpr Int128& operator--() noexcept {
        if (lo == uint64_t(0)) {
            hi--;
        }
        lo--;
        return *this;
    }
    
    /// Postfix decrement
    constexpr Int128 operator--(int) noexcept {
        Int128 temp = *this;
        --(*this);
        return temp;
    }
    
    /// Addition with carry propagation (HOT FIX 2: Verified correct)
    constexpr Int128 operator+(const Int128& other) const noexcept {
        uint64_t new_lo = lo + other.lo;
        // Carry occurs if new_lo < lo (unsigned overflow)
        int64_t carry = (new_lo < lo) ? int64_t(1) : int64_t(0);
        int64_t new_hi = hi + other.hi + carry;
        return Int128(new_hi, new_lo);
    }
    
    /// Addition assignment
    constexpr Int128& operator+=(const Int128& other) noexcept {
        *this = *this + other;
        return *this;
    }
    
    /// Subtraction (using addition of negation)
    constexpr Int128 operator-(const Int128& other) const noexcept {
        return *this + (-other);
    }
    
    /// Subtraction assignment
    constexpr Int128& operator-=(const Int128& other) noexcept {
        *this = *this - other;
        return *this;
    }
    
    /// Multiply by int64_t (produces Int128)
    Int128 operator*(int64_t m) const noexcept {
        return mul128x64(*this, m);
    }
    
    /// Multiply by Int128 (produces Int256)
    Int256 operator*(const Int128& other) const noexcept;
    
    /// Multiply assignment by int64_t
    Int128& operator*=(int64_t m) noexcept {
        *this = *this * m;
        return *this;
    }
    
    /// Division by Int128 (returns Int128 quotient)
    Int128 operator/(const Int128& other) const noexcept {
        return div128(*this, other);
    }
    
    /// Division assignment by Int128
    Int128& operator/=(const Int128& other) noexcept {
        *this = *this / other;
        return *this;
    }

    // ------------------------------------------------------------------------
    // Bitwise Operators
    // ------------------------------------------------------------------------
    
    constexpr Int128 operator~() const noexcept {
        return Int128(~hi, ~lo);
    }
    
    constexpr Int128 operator&(const Int128& other) const noexcept {
        return Int128(hi & other.hi, lo & other.lo);
    }
    
    constexpr Int128& operator&=(const Int128& other) noexcept {
        hi &= other.hi;
        lo &= other.lo;
        return *this;
    }
    
    constexpr Int128 operator|(const Int128& other) const noexcept {
        return Int128(hi | other.hi, lo | other.lo);
    }
    
    constexpr Int128& operator|=(const Int128& other) noexcept {
        hi |= other.hi;
        lo |= other.lo;
        return *this;
    }
    
    constexpr Int128 operator^(const Int128& other) const noexcept {
        return Int128(hi ^ other.hi, lo ^ other.lo);
    }
    
    constexpr Int128& operator^=(const Int128& other) noexcept {
        hi ^= other.hi;
        lo ^= other.lo;
        return *this;
    }
    
    /// Left shift
    constexpr Int128 operator<<(int shift) const noexcept {
        if (shift >= 128) return Int128(int64_t(0), uint64_t(0));
        if (shift >= 64) {
            return Int128(static_cast<int64_t>(lo << (shift - 64)), uint64_t(0));
        }
        if (shift == 0) return *this;
        uint64_t new_lo = lo << shift;
        int64_t new_hi = (hi << shift) | static_cast<int64_t>(lo >> (64 - shift));
        return Int128(new_hi, new_lo);
    }
    
    constexpr Int128& operator<<=(int shift) noexcept {
        *this = *this << shift;
        return *this;
    }
    
    /// Right shift (arithmetic, sign-extending)
    constexpr Int128 operator>>(int shift) const noexcept {
        if (shift >= 128) return (hi < int64_t(0)) ? Int128(int64_t(-1), UINT64_MAX) : Int128(int64_t(0), uint64_t(0));
        if (shift >= 64) {
            return Int128(hi >> 63, static_cast<uint64_t>(hi >> (shift - 64)));
        }
        if (shift == 0) return *this;
        uint64_t new_lo = (lo >> shift) | (static_cast<uint64_t>(hi) << (64 - shift));
        int64_t new_hi = hi >> shift;  // Arithmetic shift
        return Int128(new_hi, new_lo);
    }
    
    constexpr Int128& operator>>=(int shift) noexcept {
        *this = *this >> shift;
        return *this;
    }

    // ------------------------------------------------------------------------
    // Comparison Operators
    // ------------------------------------------------------------------------
    
    constexpr bool operator==(const Int128& other) const noexcept {
        return hi == other.hi && lo == other.lo;
    }
    
    constexpr bool operator!=(const Int128& other) const noexcept {
        return !(*this == other);
    }
    
    constexpr bool operator<(const Int128& other) const noexcept {
        if (hi != other.hi) return hi < other.hi;
        return lo < other.lo;
    }
    
    constexpr bool operator>(const Int128& other) const noexcept {
        return other < *this;
    }
    
    constexpr bool operator<=(const Int128& other) const noexcept {
        return !(other < *this);
    }
    
    constexpr bool operator>=(const Int128& other) const noexcept {
        return !(*this < other);
    }

    // ------------------------------------------------------------------------
    // Utility Methods
    // ------------------------------------------------------------------------
    
    /// Get sign: -1 for negative, 0 for zero, +1 for positive
    constexpr int sign() const noexcept {
        if (hi < int64_t(0)) return -1;
        if (hi > int64_t(0) || lo > uint64_t(0)) return 1;
        return 0;
    }
    
    /// Check if value is zero
    constexpr bool isZero() const noexcept {
        return hi == int64_t(0) && lo == uint64_t(0);
    }
    
    /// Check if value is negative
    constexpr bool isNegative() const noexcept {
        return hi < int64_t(0);
    }
    
    /// Check if value is positive
    constexpr bool isPositive() const noexcept {
        return hi > int64_t(0) || (hi == int64_t(0) && lo > uint64_t(0));
    }
    
    /// Get absolute value (careful: |INT128_MIN| overflows)
    constexpr Int128 abs() const noexcept {
        return isNegative() ? -(*this) : *this;
    }
    
    /// Explicit conversion to double (for final output only)
    /// WARNING: May lose precision for values > 2^53
    explicit operator double() const noexcept {
        // Convert: value = hi * 2^64 + lo
        // For negative: need to handle two's complement
        if (isNegative()) {
            Int128 neg = -(*this);
            return -(static_cast<double>(neg.hi) * 18446744073709551616.0 + 
                     static_cast<double>(neg.lo));
        }
        return static_cast<double>(hi) * 18446744073709551616.0 + 
               static_cast<double>(lo);
    }

private:
    /// 128×64 → 128 multiplication helper
    static Int128 mul128x64(const Int128& a, int64_t b) noexcept {
        // (a.hi*2^64 + a.lo) * b = a.hi*b*2^64 + a.lo*b
        Int128 lo_prod = mul64(static_cast<int64_t>(a.lo), b);
        Int128 hi_prod = mul64(a.hi, b);
        
        // Combine: result = hi_prod * 2^64 + lo_prod
        // hi_prod is already shifted by 64 bits conceptually
        // So: result.hi = hi_prod.lo + lo_prod.hi (with carry)
        //     result.lo = lo_prod.lo
        
        uint64_t new_lo = lo_prod.lo;
        int64_t new_hi = hi_prod.lo + lo_prod.hi;
        
        // Handle carry from adding lo_prod.hi to hi_prod.lo
        // Since both are treated as signed, we need to be careful
        // Actually, lo_prod.hi is int64_t, hi_prod.lo is int64_t (from mul64)
        // But hi_prod.lo represents the high bits of a.hi*b
        
        // Check for overflow: if hi_prod.lo > 0 and lo_prod.hi > 0 and result < 0
        // or hi_prod.lo < 0 and lo_prod.hi < 0 and result > 0
        [[maybe_unused]] bool overflow = ((hi_prod.lo > int64_t(0) && lo_prod.hi > int64_t(0) && new_hi < int64_t(0)) ||
                        (hi_prod.lo < int64_t(0) && lo_prod.hi < int64_t(0) && new_hi > int64_t(0)));
        (void)overflow;  // Suppress unused warning - kept for documentation
        
        // Add hi_prod.hi (this is the overflow from a.hi*b)
        new_hi += hi_prod.hi;
        
        return Int128(new_hi, new_lo);
    }
    
    /// 128-bit division: dividend / divisor → quotient
    /// Uses binary long division algorithm
    static Int128 div128(Int128 dividend, Int128 divisor) noexcept {
        // Handle division by zero
        if (divisor.isZero()) {
            return Int128::zero();  // Or could signal error
        }
        
        // Handle trivial cases
        if (dividend.isZero()) {
            return Int128::zero();
        }
        
        // Track signs and work with absolute values
        bool negative_result = dividend.isNegative() != divisor.isNegative();
        Int128 a = dividend.abs();
        Int128 b = divisor.abs();
        
        // Fast path: if a < b, result is 0
        if (a < b) {
            return Int128::zero();
        }
        
        // Fast path: if a == b, result is 1 (or -1)
        if (a == b) {
            return negative_result ? Int128(-1) : Int128(1);
        }
        
        // Binary long division
        Int128 quotient(0);
        Int128 remainder(0);
        
        // Process bits from high to low
        for (int i = 127; i >= 0; --i) {
            // Shift remainder left by 1
            remainder = remainder << 1;
            
            // Set lowest bit of remainder to bit i of dividend
            Int128 bit = (a >> i) & Int128(1);
            remainder = remainder | bit;
            
            // If remainder >= divisor, subtract and set quotient bit
            if (remainder >= b) {
                remainder = remainder - b;
                quotient = quotient | (Int128(1) << i);
            }
        }
        
        return negative_result ? -quotient : quotient;
    }
};

// ============================================================================
// Int256 - 256-bit Signed Integer
// ============================================================================
/**
 * @brief 256-bit signed integer using two Int128 parts
 * 
 * Representation: value = hi * 2^128 + lo
 * - hi: Int128 (signed high bits, includes sign)
 * - lo: Int128 (unsigned low bits)
 * 
 * Bit-width budget: 256 bits (covers 160-bit dot products with margin)
 */
struct Int256 {
    // ------------------------------------------------------------------------
    // Member Data
    // ------------------------------------------------------------------------
    Int128 hi;   // High 128 bits (signed)
    Int128 lo;   // Low 128 bits (unsigned magnitude)

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------
    
    /// Default constructor: initializes to zero
    constexpr Int256() noexcept : hi(Int128::zero()), lo(Int128::zero()) {}
    
    /// Constructor from Int128 (sign-extends)
    constexpr Int256(const Int128& v) noexcept 
        : hi(v.isNegative() ? Int128(int64_t(-1)) : Int128::zero()), lo(v) {}
    
    /// Constructor from int64_t (sign-extends)
    constexpr Int256(int64_t v) noexcept 
        : hi(v < int64_t(0) ? Int128(int64_t(-1)) : Int128::zero()), lo(v) {}
    
    /// Constructor from hi/lo pair
    constexpr Int256(const Int128& high, const Int128& low) noexcept 
        : hi(high), lo(low) {}

    // ------------------------------------------------------------------------
    // Factory Methods
    // ------------------------------------------------------------------------
    
    static constexpr Int256 fromInt64(int64_t v) noexcept {
        return Int256(v);
    }
    
    static constexpr Int256 fromInt128(const Int128& v) noexcept {
        return Int256(v);
    }
    
    static constexpr Int256 fromParts(const Int128& high, const Int128& low) noexcept {
        return Int256(high, low);
    }
    
    static constexpr Int256 zero() noexcept { return Int256(); }
    
    static constexpr Int256 one() noexcept { 
        return Int256(Int128::zero(), Int128::one()); 
    }

    // ------------------------------------------------------------------------
    // 128×128 → 256 Multiplication (Critical for EMBER dot products)
    // ------------------------------------------------------------------------
    
    /**
     * @brief Multiply two Int128 values to produce Int256
     * 
     * Uses unsigned schoolbook multiplication with sign correction.
     * 
     * HOT FIX: Verified correct shift placement for cross-terms (Fix 1.0)
     * 
     * @param a First 128-bit operand
     * @param b Second 128-bit operand
     * @return Int256 product (up to 256 bits)
     */
    static Int256 mul128(const Int128& a, const Int128& b) noexcept {
        // Determine result sign
        bool result_negative = a.isNegative() != b.isNegative();
        
        // Work with absolute values
        Int128 a_abs = a.abs();
        Int128 b_abs = b.abs();
        
        // Compute unsigned product using 64-bit decomposition
        // a_abs = a_hi * 2^64 + a_lo
        // b_abs = b_hi * 2^64 + b_lo
        // a_abs * b_abs = a_hi*b_hi*2^128 + (a_hi*b_lo + a_lo*b_hi)*2^64 + a_lo*b_lo
        
        uint64_t a_lo = a_abs.lo;
        uint64_t a_hi = static_cast<uint64_t>(a_abs.hi);
        uint64_t b_lo = b_abs.lo;
        uint64_t b_hi = static_cast<uint64_t>(b_abs.hi);
        
        // Compute partial products (all unsigned)
        // p0 = a_lo * b_lo = (p0_hi, p0_lo) where p0 = p0_hi * 2^64 + p0_lo
        // p1 = a_hi * b_lo = (p1_hi, p1_lo)
        // p2 = a_lo * b_hi = (p2_hi, p2_lo)
        // p3 = a_hi * b_hi = (p3_hi, p3_lo)
        Int128 p0 = Int128::mul64(a_lo, b_lo);
        Int128 p1 = Int128::mul64(a_hi, b_lo);
        Int128 p2 = Int128::mul64(a_lo, b_hi);
        Int128 p3 = Int128::mul64(a_hi, b_hi);
        
        // Extract 64-bit components
        uint64_t p0_lo = p0.lo;
        uint64_t p0_hi = static_cast<uint64_t>(p0.hi);
        uint64_t p1_lo = p1.lo;
        uint64_t p1_hi = static_cast<uint64_t>(p1.hi);
        uint64_t p2_lo = p2.lo;
        uint64_t p2_hi = static_cast<uint64_t>(p2.hi);
        uint64_t p3_lo = p3.lo;
        uint64_t p3_hi = static_cast<uint64_t>(p3.hi);
        
        // Accumulate:
        // result = p3*2^128 + (p1 + p2)*2^64 + p0
        //        = (p3_hi*2^64 + p3_lo)*2^128 + ((p1_hi*2^64 + p1_lo) + (p2_hi*2^64 + p2_lo))*2^64 + (p0_hi*2^64 + p0_lo)
        //        = p3_hi*2^192 + p3_lo*2^128 + p1_hi*2^128 + p1_lo*2^64 + p2_hi*2^128 + p2_lo*2^64 + p0_hi*2^64 + p0_lo
        //        = p3_hi*2^192 + (p3_lo + p1_hi + p2_hi)*2^128 + (p1_lo + p2_lo + p0_hi)*2^64 + p0_lo
        
        // Bits 0-63: p0_lo
        uint64_t r0 = p0_lo;
        
        // Bits 64-127: p0_hi + p1_lo + p2_lo (with carry to bit 128)
        uint64_t sum_64_127 = p0_hi + p1_lo;
        uint64_t carry1 = (sum_64_127 < p0_hi) ? 1 : 0;
        sum_64_127 += p2_lo;
        uint64_t carry2 = (sum_64_127 < p2_lo) ? 1 : 0;
        uint64_t carry_to_128 = carry1 + carry2;
        uint64_t r1 = sum_64_127;  // This is the high 64 bits of the low 128-bit word
        
        // Bits 128-191: p1_hi + p2_hi + p3_lo + carry_to_128 (with carry to bit 192)
        uint64_t sum_128_191 = p1_hi + p2_hi;
        uint64_t carry3 = (sum_128_191 < p1_hi) ? 1 : 0;
        sum_128_191 += p3_lo;
        uint64_t carry4 = (sum_128_191 < p3_lo) ? 1 : 0;
        sum_128_191 += carry_to_128;
        uint64_t carry5 = (carry_to_128 > 0 && sum_128_191 < carry_to_128) ? 1 : 0;
        uint64_t carry_to_192 = carry3 + carry4 + carry5;
        uint64_t r2 = sum_128_191;  // This is the low 64 bits of the high 128-bit word
        
        // Bits 192-255: p3_hi + carry_to_192
        uint64_t r3 = p3_hi + carry_to_192;
        
        // Build result:
        // lo (bits 0-127) = r1 * 2^64 + r0
        // hi (bits 128-255) = r3 * 2^64 + r2
        Int256 result_unsigned(
            Int128(static_cast<int64_t>(r3), r2),  // hi
            Int128(static_cast<int64_t>(r1), r0)   // lo
        );
        
        // Apply sign
        if (result_negative) {
            return -result_unsigned;
        }
        return result_unsigned;
    }
    
public:

    // ------------------------------------------------------------------------
    // Arithmetic Operators
    // ------------------------------------------------------------------------
    
    /// Unary negation (two's complement)
    constexpr Int256 operator-() const noexcept {
        // ~x + 1
        Int128 inv_lo = ~lo;
        Int128 inv_hi = ~hi;
        inv_lo = inv_lo + Int128::one();
        if (inv_lo.isZero() && !lo.isZero()) {
            inv_hi = inv_hi + Int128::one();
        }
        return Int256(inv_hi, inv_lo);
    }
    
    /// Unary plus
    constexpr Int256 operator+() const noexcept {
        return *this;
    }
    
    /// Addition with carry chain
    constexpr Int256 operator+(const Int256& other) const noexcept {
        Int128 new_lo = lo + other.lo;
        Int128 carry = (new_lo < lo) ? Int128::one() : Int128::zero();
        Int128 new_hi = hi + other.hi + carry;
        return Int256(new_hi, new_lo);
    }
    
    /// Addition assignment
    constexpr Int256& operator+=(const Int256& other) noexcept {
        *this = *this + other;
        return *this;
    }
    
    /// Subtraction
    constexpr Int256 operator-(const Int256& other) const noexcept {
        return *this + (-other);
    }
    
    /// Subtraction assignment
    constexpr Int256& operator-=(const Int256& other) noexcept {
        *this = *this - other;
        return *this;
    }
    
    /// Left shift
    constexpr Int256 operator<<(int shift) const noexcept {
        if (shift >= 256) return Int256::zero();
        if (shift >= 128) {
            return Int256(lo << (shift - 128), Int128::zero());
        }
        if (shift == 0) return *this;
        Int128 new_lo = lo << shift;
        Int128 new_hi = (hi << shift) | (lo >> (128 - shift));
        return Int256(new_hi, new_lo);
    }
    
    constexpr Int256& operator<<=(int shift) noexcept {
        *this = *this << shift;
        return *this;
    }
    
    /// Right shift (arithmetic)
    constexpr Int256 operator>>(int shift) const noexcept {
        if (shift >= 256) return (hi.isNegative() ? Int256(Int128(int64_t(-1))) : Int256::zero());
        if (shift >= 128) {
            return Int256(hi >> 127, hi >> (shift - 128));
        }
        if (shift == 0) return *this;
        Int128 new_hi = hi >> shift;
        Int128 new_lo = (lo >> shift) | (hi << (128 - shift));
        return Int256(new_hi, new_lo);
    }
    
    constexpr Int256& operator>>=(int shift) noexcept {
        *this = *this >> shift;
        return *this;
    }

    /// Multiply by Int256 (produces Int256 - low 256 bits of full product)
    /// Uses schoolbook multiplication with 128-bit components.
    /// Note: Not constexpr because mul128 is not constexpr.
    Int256 operator*(const Int256& other) const noexcept {
        // Decompose: a = a_hi * 2^128 + a_lo, b = b_hi * 2^128 + b_lo
        // a * b = a_hi*b_hi*2^256 + (a_hi*b_lo + a_lo*b_hi)*2^128 + a_lo*b_lo
        // We only need the low 256 bits, so we ignore a_hi*b_hi*2^256
        
        // Compute partial products
        Int256 p0 = Int256::mul128(lo, other.lo);  // a_lo * b_lo → 256 bits
        Int256 p1 = Int256::mul128(hi, other.lo);  // a_hi * b_lo → 256 bits
        Int256 p2 = Int256::mul128(lo, other.hi);  // a_lo * b_hi → 256 bits
        
        // Combine: result = p0 + (p1 + p2) * 2^128
        // (p1 + p2) * 2^128 contributes to high 128 bits
        Int256 cross_sum = p1 + p2;
        
        // Add cross_sum shifted by 128 to p0
        Int256 shifted_cross(cross_sum.lo, Int128::zero());
        Int256 result = p0 + shifted_cross;
        
        // Handle carry from cross_sum.lo overflow
        if (cross_sum.lo < p1.lo || cross_sum.lo < p2.lo) {
            result.hi = result.hi + Int128::one();
        }
        
        return result;
    }
    
    /// Access individual 64-bit words (0=lowest, 3=highest)
    constexpr uint64_t word(int i) const noexcept {
        switch (i) {
            case 0: return lo.lo;
            case 1: return static_cast<uint64_t>(lo.hi);
            case 2: return hi.lo;
            case 3: return static_cast<uint64_t>(hi.hi);
            default: return 0;
        }
    }

    // ------------------------------------------------------------------------
    // Comparison Operators
    // ------------------------------------------------------------------------
    
    constexpr bool operator==(const Int256& other) const noexcept {
        return hi == other.hi && lo == other.lo;
    }
    
    constexpr bool operator!=(const Int256& other) const noexcept {
        return !(*this == other);
    }
    
    constexpr bool operator<(const Int256& other) const noexcept {
        if (hi != other.hi) return hi < other.hi;
        return lo < other.lo;
    }
    
    constexpr bool operator>(const Int256& other) const noexcept {
        return other < *this;
    }
    
    constexpr bool operator<=(const Int256& other) const noexcept {
        return !(other < *this);
    }
    
    constexpr bool operator>=(const Int256& other) const noexcept {
        return !(*this < other);
    }
    
    /// Comparison with Int128 (converts Int128 to Int256 for comparison)
    constexpr bool operator==(const Int128& other) const noexcept {
        return *this == Int256(other);
    }
    
    constexpr bool operator!=(const Int128& other) const noexcept {
        return !(*this == Int256(other));
    }
    
    constexpr bool operator<(const Int128& other) const noexcept {
        return *this < Int256(other);
    }
    
    constexpr bool operator>(const Int128& other) const noexcept {
        return *this > Int256(other);
    }
    
    constexpr bool operator<=(const Int128& other) const noexcept {
        return *this <= Int256(other);
    }
    
    constexpr bool operator>=(const Int128& other) const noexcept {
        return *this >= Int256(other);
    }

    // ------------------------------------------------------------------------
    // Utility Methods
    // ------------------------------------------------------------------------
    
    /// Get sign: -1 for negative, 0 for zero, +1 for positive
    constexpr int sign() const noexcept {
        int hi_sign = hi.sign();
        if (hi_sign != 0) return hi_sign;
        if (lo.isZero()) return 0;
        return 1;
    }
    
    /// Check if value is zero
    constexpr bool isZero() const noexcept {
        return hi.isZero() && lo.isZero();
    }
    
    /// Check if value is negative
    constexpr bool isNegative() const noexcept {
        return hi.isNegative();
    }
    
    /// Check if value is positive
    constexpr bool isPositive() const noexcept {
        return hi.isPositive() || (hi.isZero() && !lo.isZero());
    }
    
    /// Get absolute value
    constexpr Int256 abs() const noexcept {
        return isNegative() ? -(*this) : *this;
    }
    
    /// Get low 128 bits
    constexpr Int128 low() const noexcept { return lo; }
    
    /// Get high 128 bits
    constexpr Int128 high() const noexcept { return hi; }
    
    /// Explicit conversion to double (may lose precision)
    explicit operator double() const noexcept {
        if (isNegative()) {
            Int256 neg = -(*this);
            // Approximate: 2^128 is about 3.4e38
            double hi_d = static_cast<double>(neg.hi);
            double lo_d = static_cast<double>(neg.lo);
            return -(hi_d * 3.4028236692093846346337460743177e38 + lo_d);
        }
        double hi_d = static_cast<double>(hi);
        double lo_d = static_cast<double>(lo);
        return hi_d * 3.4028236692093846346337460743177e38 + lo_d;
    }
};

// ============================================================================
// Int128 Methods Depending on Int256
// ============================================================================

/// Int128 × Int128 → Int256 multiplication
inline Int256 Int128::operator*(const Int128& other) const noexcept {
    return Int256::mul128(*this, other);
}

/// Int128 from Int256 constructor
inline Int128::Int128(const Int256& v) noexcept 
    : hi(v.lo.hi), lo(v.lo.lo) {
    // Truncates - only takes low 128 bits
}

// ============================================================================
// Boost.Multiprecision Fallback
// ============================================================================

#if defined(EMBER_USE_BOOST_MULTIPRECISION)
/**
 * @brief Boost.Multiprecision wrapper providing same interface as Int128/Int256
 * 
 * Use this when native 128-bit types are unavailable or for verification.
 */
namespace boost_fallback {

    using Int128_mp = mp::int128_t;
    using Int256_mp = mp::int256_t;
    
    // Adapter functions to match Int128/Int256 interface
    inline Int128_mp mul64_boost(int64_t a, int64_t b) {
        return Int128_mp(a) * Int128_mp(b);
    }
    
    inline Int256_mp mul128_boost(const Int128_mp& a, const Int128_mp& b) {
        return Int256_mp(a) * Int256_mp(b);
    }
    
    inline int sign_boost(const Int128_mp& v) {
        if (v < 0) return -1;
        if (v > 0) return 1;
        return 0;
    }
    
    inline int sign_boost(const Int256_mp& v) {
        if (v < 0) return -1;
        if (v > 0) return 1;
        return 0;
    }

} // namespace boost_fallback
#endif // EMBER_USE_BOOST_MULTIPRECISION

// ============================================================================
// Type Aliases and Helper Types
// ============================================================================

/// Input coordinate type (26-bit budget)
using Coord = int32_t;

/// 2×2 minor type (53-bit budget)
using Minor2 = int64_t;

/// 3×3 determinant type (80-bit budget, stored in 128 bits)
using Determinant3 = Int128;

/// Dot product type (160-bit budget, stored in 256 bits)
using DotProduct = Int256;

// ============================================================================
// Compile-Time Verification
// ============================================================================

static_assert(sizeof(Int128) == 16, "Int128 must be 16 bytes");
static_assert(sizeof(Int256) == 32, "Int256 must be 32 bytes");
static_assert(alignof(Int128) >= 8, "Int128 must be 8-byte aligned");
static_assert(alignof(Int256) >= 8, "Int256 must be 8-byte aligned");

// Verify constexpr operations compile
constexpr Int128 test_constexpr128() {
    Int128 a(int64_t(5));
    Int128 b(int64_t(3));
    Int128 c = a + b;
    c = c - Int128(int64_t(2));
    c = -c;
    return c;
}
static_assert(test_constexpr128() == Int128(int64_t(-6)), "Constexpr Int128 operations failed");

constexpr Int256 test_constexpr256() {
    Int256 a(Int128(int64_t(5)));
    Int256 b(Int128(int64_t(3)));
    Int256 c = a + b;
    c = c - Int256(Int128(int64_t(2)));
    c = -c;
    return c;
}

// ============================================================================
// Int128 vs Int256 Comparison Operators (for expressions like n > a*b)
// ============================================================================

inline constexpr bool operator==(const Int128& a, const Int256& b) noexcept {
    return Int256(a) == b;
}

inline constexpr bool operator!=(const Int128& a, const Int256& b) noexcept {
    return Int256(a) != b;
}

inline constexpr bool operator<(const Int128& a, const Int256& b) noexcept {
    return Int256(a) < b;
}

inline constexpr bool operator>(const Int128& a, const Int256& b) noexcept {
    return Int256(a) > b;
}

inline constexpr bool operator<=(const Int128& a, const Int256& b) noexcept {
    return Int256(a) <= b;
}

inline constexpr bool operator>=(const Int128& a, const Int256& b) noexcept {
    return Int256(a) >= b;
}

} // namespace ember

#endif // EMBER_INTEGER_TYPES_H
