# EMBER Boolean Engine — Windows Build Guide for Houdini

This guide explains how to compile the EMBER Boolean Engine on Windows for Houdini.

---

## Prerequisites

### Required Software

1. **Houdini 20.0+** (tested with 20.0.625)
   - Download from [SideFX.com](https://www.sidefx.com/download/)
   - Install with "Houdini Command Line Tools" option

2. **Visual Studio 2022** (Community edition is free)
   - Download from [Visual Studio](https://visualstudio.microsoft.com/)
   - Required components:
     - "Desktop development with C++" workload
     - Windows 10/11 SDK
     - CMake tools for Windows

3. **CMake 3.28+**
   - Included with Visual Studio 2022, or download from [cmake.org](https://cmake.org/download/)

---

## Quick Start (One Command)

### Option 1: Using Houdini Command Line Tools

1. Open **"Houdini Command Line Tools"** from the Start Menu
2. Navigate to the EMBER directory:
   ```cmd
   cd C:\path\to\EMBER
   ```
3. Run the build script:
   ```cmd
   build_windows.bat install
   ```

### Option 2: Manual Build

1. Open **"x64 Native Tools Command Prompt for VS 2022"** from Start Menu
2. Set Houdini environment:
   ```cmd
   set HFS=C:\Program Files\Side Effects Software\Houdini 20.0.XXX
   set PATH=%HFS%\bin;%PATH%
   ```
3. Build:
   ```cmd
   cd C:\path\to\EMBER
   mkdir build && cd build
   cmake .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_BUILD_TYPE=Release
   cmake --build . --config Release --parallel
   ```

---

## Detailed Build Instructions

### Step 1: Verify Prerequisites

Open a command prompt and verify:

```cmd
cmake --version
REM Should show 3.28 or higher

cl
REM Should show Microsoft C/C++ Optimizing Compiler
```

### Step 2: Set Houdini Environment

#### Method A: Using Houdini Command Line Tools (Recommended)

Simply open **"Houdini Command Line Tools"** from the Start Menu. This sets up all required environment variables automatically.

#### Method B: Manual Setup

```cmd
set HFS=C:\Program Files\Side Effects Software\Houdini 20.0.XXX
set PATH=%HFS%\bin;%PATH%
set HOUDINI_DSO_PATH=%USERPROFILE%\Documents\houdini20.0\dso
```

### Step 3: Configure Build

```cmd
cd C:\path\to\EMBER
mkdir build
cd build

cmake .. -G "Visual Studio 17 2022" -A x64 ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DEMBER_BUILD_TESTS=ON
```

**CMake Options:**
- `-G "Visual Studio 17 2022"` — Use VS 2022 generator
- `-A x64` — Build for 64-bit Windows
- `-DCMAKE_BUILD_TYPE=Release` — Release build (optimized)
- `-DEMBER_BUILD_TESTS=ON` — Build unit tests
- `-DEMBER_ENABLE_CUDA=ON` — Enable CUDA (requires CUDA toolkit)

### Step 4: Build

```cmd
cmake --build . --config Release --parallel
```

This creates:
- `build\Release\SOP_EmberBoolean.dll` — The Houdini plugin
- `build\Release\test_predicates.exe` — Unit tests
- `build\Release\test_cdt.exe` — CDT tests

### Step 5: Install

#### Automatic Install

```cmd
cmake --install . --config Release
```

This installs to `%USERPROFILE%\Documents\houdini20.0\dso\sop\`

#### Manual Install

```cmd
copy build\Release\SOP_EmberBoolean.dll %USERPROFILE%\Documents\houdini20.0\dso\sop\
```

---

## Verification

### Test 1: Check Plugin Loads

1. Start Houdini
2. Create a geometry object (e.g., a sphere)
3. Press **TAB** and type **"ember"**
4. Select **"Ember Boolean"** from the list

### Test 2: Verify with Houdini Console

In Houdini's Python Shell:

```python
import hou
node = hou.node('/obj/geo1').createNode('sop_ember_boolean')
print(node)
# Should print: <SopNode of type sop_ember_boolean at /obj/geo1/sop_ember_boolean1>
```

### Test 3: Run Unit Tests

```cmd
cd C:\path\to\EMBER\build\Release
test_predicates.exe
test_cdt.exe
```

---

## Troubleshooting

### Error: "Could not find Houdini"

**Solution:** Set `HFS` environment variable:

```cmd
set HFS=C:\Program Files\Side Effects Software\Houdini 20.0.XXX
```

### Error: "CMake Error: Could not find CMAKE_ROOT"

**Solution:** Add CMake to PATH:

```cmd
set PATH=C:\Program Files\CMake\bin;%PATH%
```

### Error: "The C++ compiler does not support C++17"

**Solution:** Update Visual Studio 2022 or use a newer version.

### Error: "Cannot open include file: 'HoudiniAPI.h'"

**Solution:** Houdini environment not set up. Run from "Houdini Command Line Tools" or set:

```cmd
set HFS=C:\Program Files\Side Effects Software\Houdini 20.0.XXX
set INCLUDE=%HFS%\toolkit\include;%INCLUDE%
```

### Plugin Not Appearing in Houdini

1. Check the DSO path:
   ```python
   import hou
   print(hou.expandString('$HOUDINI_DSO_PATH'))
   ```

2. Verify the plugin file exists:
   ```cmd
   dir %USERPROFILE%\Documents\houdini20.0\dso\sop\SOP_EmberBoolean.dll
   ```

3. Check Houdini's console for errors (Windows → Console)

### Build Warnings About Integer Types

These are expected and harmless:
- C4244: conversion from 'int64_t' to 'int32_t'
- C4267: conversion from 'size_t' to 'int'

The code handles these conversions safely.

---

## Platform-Specific Notes

### MSVC Compiler Flags

The CMakeLists.txt sets these flags for MSVC:

- `/W4` — Warning level 4
- `/WX-` — Don't treat warnings as errors
- `/permissive-` — Strict conformance mode
- `/O2` — Optimize for speed (Release)
- `/fp:strict` — Strict floating-point model

### __int128 vs MSVC

The code **does NOT use `__int128`** on Windows. Instead, it uses:

- MSVC intrinsics: `_mul128()`, `_umul128()`
- Portable fallback: 32-bit schoolbook multiplication

This is handled automatically in `IntegerTypes.h`:

```cpp
#if EMBER_IS_MSVC
    // MSVC: Use _mul128 intrinsic
    int64_t high;
    int64_t low = _mul128(a, b, &high);
    return Int128(high, static_cast<uint64_t>(low));
#endif
```

---

## Quick Reference

| Task | Command |
|------|---------|
| Configure | `cmake .. -G "Visual Studio 17 2022" -A x64` |
| Build | `cmake --build . --config Release --parallel` |
| Install | `cmake --install . --config Release` |
| Test | `Release\test_predicates.exe` |
| Clean | `rmdir /s /q build` |

---

## Support

For issues specific to Windows builds:
1. Check Houdini version compatibility (20.0+ required)
2. Verify Visual Studio 2022 is installed with C++ workload
3. Ensure Houdini Command Line Tools are used
4. Check Windows Event Viewer for DLL load errors
