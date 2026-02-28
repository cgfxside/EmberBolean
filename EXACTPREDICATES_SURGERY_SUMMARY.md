# ExactPredicates.h Unification Surgery — Integration Summary

**Date:** 2026-02-28  
**Task:** Implement ExactPredicates_Surgery_Changelog.md fixes

---

## Changes Applied

### 1. Type Unification (Phase A)

| Metric | Before | After |
|--------|--------|-------|
| Total lines | 1,233 | 944 |
| Int128 class definitions | 1 (duplicate) | 0 (uses IntegerTypes.h) |
| Int256 class definitions | 1 (duplicate) | 0 (uses IntegerTypes.h) |
| `Int128::mul()` calls | 18 | 0 (replaced with `mul64`) |
| `Int128::mul64()` calls | 0 | 20 |
| `Int256::mul()` calls | 8 | 0 (replaced with `mul128`) |
| `Int256::mul128()` calls | 0 | 7 |
| FP divisions in orient2d_indirect | 3 | 0 |
| Epsilon in pointCompare_indirect | 3 | 0 |
| `#include "IntegerTypes.h"` | 0 | 1 |
| `#include "ExactPredicates_Indirect.h"` | 0 | 1 |

**Deleted:** Lines 74-417 of original (331 lines)
- `class Int128` (152 lines) — replaced by `ember::Int128` from IntegerTypes.h
- `class Int256` (163 lines) — replaced by `ember::Int256` from IntegerTypes.h
- `Int128::fromDouble()` — removed (quantize before constructing)
- `Int256::fromDouble()` — removed (same reason)

**Added to ExactPredicates.h:**
```cpp
#include "IntegerTypes.h"
#include "ExactPredicates_Indirect.h"

using ember::Int128;
using ember::Int256;
```

**API Remapping:**
| Old Call | New Call | Count |
|----------|----------|-------|
| `Int128::mul(a, b)` | `Int128::mul64(a, b)` | 18→20 |
| `Int256::mul(Int128(x), y)` | `Int256::mul128(Int128(x), y)` | 8→7 |

---

### 2. New File: ExactPredicates_Indirect.h

**Created:** `/mnt/okcomputer/output/ember/ExactPredicates_Indirect.h` (355 lines)

**Provides:**
- `LPI_Exact` struct — exact LPI representation with Int128 homogeneous coordinates
- `HomogPoint2D` struct — 2D point in homogeneous coordinates
- `compute_LPI_exact()` — computes LPI without FP division
- `HomogPoint2D::fromExplicit()` — creates from explicit coordinates
- `HomogPoint2D::fromLPI()` — creates from LPI with projection
- `orient2d_homogeneous_exact()` — ZERO FP division orient2d
- `pointCompare_exact()` — ZERO epsilon LPI comparison
- `pointCompare_explicit_vs_lpi()` — mixed explicit/LPI comparison

---

### 3. orient2d_indirect Fix (Phase B — CRITICAL)

**Deleted:** Lines 1015-1097 of original (82 lines)
- Three identical blocks of FP-division LPI computation
- `orient2d_filtered(ax, ay, bx, by, cx, cy)` call on FP results

**Added:** 78 lines of exact homogeneous orient2d
- `makeHomogPoint` lambda: builds `HomogPoint2D` from explicit or LPI
- For explicit points: quantizes to integers, wraps with `denom=1`
- For LPI points: calls `compute_LPI_exact()` → `HomogPoint2D::fromLPI()`
- Final call: `orient2d_homogeneous_exact(A, B, C)` — ZERO FP divisions

---

### 4. pointCompare_indirect Fix (Phase B)

**Deleted:** Lines 1100-1176 of original (76 lines)
- FP-division LPI coordinate computation
- Epsilon-based lexicographic comparison:
  ```cpp
  const double eps = EPSILON * (std::abs(ax) + std::abs(bx) + 1.0);
  if (ax < bx - eps) return -1;
  ```

**Added:** 73 lines of exact comparison
- Both-explicit path: direct double comparison (bit-exact for identical coords)
- Both-implicit path: `compute_LPI_exact()` + `pointCompare_exact()` (Int256 products)
- Mixed path: `pointCompare_explicit_vs_lpi()` (Int128 cross-multiplication)

---

## Verification Results

```bash
# Verify no duplicate type definitions
grep -c "^class Int128\|^class Int256" ExactPredicates.h
# Result: 0 ✓

# Verify no old API calls
grep "Int128::mul(" ExactPredicates.h | grep -v "mul64\|//"
# Result: (empty) ✓

# Verify no FP division in code
grep -v "//" ExactPredicates.h | grep "num / den"
# Result: (empty) ✓

# Verify no epsilon in point comparison
grep "eps.*abs" ExactPredicates.h | grep -v "//\|ERRBOUND\|\*"
# Result: (empty) ✓
```

---

## Dependencies

The unified file requires:
1. `IntegerTypes.h` — canonical Int128/Int256 (carry fix already applied)
2. `ExactPredicates_Indirect.h` — LPI_Exact, HomogPoint2D, exact predicates

---

## Files Modified/Created

| File | Action | Lines |
|------|--------|-------|
| `ember/ExactPredicates.h` | Replaced | 944 (was 1,233) |
| `ember/ExactPredicates_Indirect.h` | Created | 355 |

---

## Status: ✅ COMPLETE

All changes from ExactPredicates_Surgery_Changelog.md have been implemented:
- ✓ Type unification (deleted duplicate Int128/Int256)
- ✓ API remapping (mul → mul64, mul → mul128)
- ✓ Zero FP divisions in orient2d_indirect
- ✓ Zero epsilon in pointCompare_indirect
- ✓ New ExactPredicates_Indirect.h created
- ✓ All verification checks pass
