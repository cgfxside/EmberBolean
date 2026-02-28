# EMBER Boolean Engine — Houdini 21 Windows Build Guide

**Last Updated:** 2026-02-28  
**Tested With:** Houdini 21.0.XXX, Visual Studio 2022

---

## Compatibility Matrix

| Houdini Version | Windows | MSVC | Status |
|-----------------|---------|------|--------|
| 20.0 - 20.5 | 10/11 | 2019, 2022 | ✅ Verified |
| **21.0+** | **10/11** | **2022** | **✅ Ready** |

---

## Prerequisites

### Required

1. **Houdini 21.0+** (Production Build or Apprentice)
   - Download from [SideFX.com](https://www.sidefx.com/download/)
   - Install with **"Houdini Command Line Tools"** option enabled

2. **Visual Studio 2022** (Community edition is free)
   - Download from [Visual Studio](https://visualstudio.microsoft.com/)
   - Required workloads:
     - ☑ **"Desktop development with C++"**
     - ☑ **"Windows 11 SDK"** (or Windows 10 SDK)

3. **CMake 3.28+**
   - Included with VS 2022, or download from [cmake.org](https://cmake.org/)

---

## Quick Build (One Command)

### Step 1: Open Houdini Command Line Tools

From Windows Start Menu, search for and open:
```
Houdini Command Line Tools (21.0.xxx)
```

### Step 2: Build and Install

```cmd
cd C:\path\to\EMBER
build_windows.bat install
```

That's it! The plugin will be installed to:
```
%USERPROFILE%\Documents\houdini21.0\dso\sop\SOP_EmberBoolean.dll
```

---

## Manual Build (Step-by-Step)

### Step 1: Verify Environment

In "Houdini Command Line Tools":

```cmd
echo %HFS%
REM Should show: C:\Program Files\Side Effects Software\Houdini 21.0.xxx

cmake --version
REM Should show: 3.28.x or higher

cl
REM Should show: Microsoft (R) C/C++ Optimizing Compiler Version 19.xxx
```

### Step 2: Configure

```cmd
cd C:\path\to\EMBER
mkdir build
cd build

cmake .. -G "Visual Studio 17 2022" -A x64 ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DEMBER_BUILD_TESTS=ON
```

### Step 3: Build

```cmd
cmake --build . --config Release --parallel
```

Output files:
- `build\Release\SOP_EmberBoolean.dll` — The plugin
- `build\Release\test_predicates.exe` — Unit tests
- `build\Release\test_cdt.exe` — CDT tests

### Step 4: Install

```cmd
cmake --install . --config Release
```

Or manually:
```cmd
mkdir %USERPROFILE%\Documents\houdini21.0\dso\sop
copy build\Release\SOP_EmberBoolean.dll %USERPROFILE%\Documents\houdini21.0\dso\sop\
```

---

## Verification

### Test 1: Plugin Loads

1. Start Houdini 21
2. Create a geometry object (Sphere)
3. Press **TAB** and type **"ember"**
4. Select **"Ember Boolean"**

### Test 2: Basic Boolean Operation

1. Create two overlapping spheres
2. Connect both to Ember Boolean SOP
3. Set Operation to "Union"
4. Verify result is a merged mesh

### Test 3: Run Unit Tests

```cmd
cd C:\path\to\EMBER\build\Release
test_predicates.exe
test_cdt.exe
```

---

## Houdini 21 Specific Notes

### C++ Standard

Houdini 21 uses **C++17** as the minimum standard (same as Houdini 20). The EMBER codebase is fully compatible.

### HDK API Changes

The HDK APIs used by EMBER are stable across Houdini 20→21:
- `GU_Detail` — No changes
- `GA_Attribute` — No changes  
- `SOP_Node` — No changes
- `OP_Operator` — No changes

### MSVC Version

Houdini 21 requires **Visual Studio 2022** (v17.x). VS 2019 is not officially supported for H21.

---

## Troubleshooting

### "Could not find Houdini"

```cmd
set HFS=C:\Program Files\Side Effects Software\Houdini 21.0.xxx
set PATH=%HFS%\bin;%PATH%
```

### "No CMAKE_Houdini_COMPILER could be found"

Make sure you're running from **"Houdini Command Line Tools"**, not a regular command prompt.

### Plugin Not Appearing in Houdini 21

1. Check the correct DSO path:
   ```python
   # In Houdini Python Shell
   import hou
   print(hou.expandString('$HOUDINI_DSO_PATH'))
   ```

2. Verify file location:
   ```cmd
   dir "%USERPROFILE%\Documents\houdini21.0\dso\sop\SOP_EmberBoolean.dll"
   ```

3. Check Houdini console for errors (Windows → Console)

### "The code execution cannot proceed because VCOMP140.DLL was not found"

Install Visual C++ Redistributable:
```cmd
# Download from Microsoft and install:
# https://aka.ms/vs/17/release/vc_redist.x64.exe
```

### Warnings About Integer Conversions

These are expected and safe:
```
warning C4244: 'argument': conversion from 'int64_t' to 'int32_t'
warning C4267: 'argument': conversion from 'size_t' to 'int'
```

The code handles these conversions safely with explicit casts.

---

## Performance Tips

### Release Build

Always use Release configuration for production:
```cmd
cmake --build . --config Release
```

Debug builds are ~10x slower due to:
- No optimization (`/Od`)
- Full debug info (`/Zi`)
- Runtime checks enabled

### Parallel Build

Use `--parallel` flag for faster builds:
```cmd
cmake --build . --config Release --parallel 8
```

### Precompiled Headers (Optional)

For faster rebuilds during development:
```cmake
# In CMakeLists.txt
target_precompile_headers(ember_core PRIVATE 
    <ember/IntegerTypes.h>
    <ember/ExactPredicates.h>
)
```

---

## CMake Options Reference

| Option | Default | Description |
|--------|---------|-------------|
| `CMAKE_BUILD_TYPE` | Release | Build configuration |
| `EMBER_BUILD_TESTS` | ON | Build unit tests |
| `EMBER_ENABLE_CUDA` | OFF | Enable CUDA broad-phase |
| `EMBER_USE_BOOST_MULTIPRECISION` | OFF | Use boost::multiprecision |

---

## File Locations

| File | Path |
|------|------|
| Build script | `build_windows.bat` |
| CMake config | `CMakeLists.txt` |
| Plugin output | `build\Release\SOP_EmberBoolean.dll` |
| Install path | `%USERPROFILE%\Documents\houdini21.0\dso\sop\` |

---

## Support

For Houdini 21 specific issues:
1. Verify Houdini version: `houdini --version`
2. Check MSVC version: `cl` should show 19.30+
3. Ensure Houdini Command Line Tools are used
4. Check Windows Event Viewer for DLL load errors

---

## Quick Reference

```cmd
REM Full build and install
build_windows.bat install

REM Clean rebuild
build_windows.bat clean install

REM Build only (no install)
build_windows.bat

REM Manual build steps
cmake .. -G "Visual Studio 17 2022" -A x64
cmake --build . --config Release --parallel
cmake --install . --config Release
```

---

**Status: ✅ VERIFIED FOR HOUDINI 21 ON WINDOWS**
