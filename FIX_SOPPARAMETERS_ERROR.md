# Fix: SOPParameters.h Include Error

**Error:**
```
E:\Tools\EMBRA\EMBER_04\ember\CutterDispatch.cpp(1,10): fatal error : 'SOPParameters.h' file not found
```

**Root Cause:**
The `ember/` core library was trying to include `SOPParameters.h` from `src/`, but `src/` is NOT in the include path for `ember_core`. This created an architectural dependency issue - the core library should NOT depend on SOP-specific code.

**Solution:**
Created a Houdini-agnostic parameter struct in `ember/Parameters.h` that the core library can use.

---

## Files Changed

### 1. NEW: `ember/Parameters.h`
- Core parameter struct (`ember::Parameters`) 
- Houdini-agnostic (no HDK dependencies)
- Used by `ember/` library and `backend/` library

### 2. MODIFIED: `ember/CutterDispatch.h`
- Now includes `ember/Parameters.h`
- Function signature changed from `const struct SOPParameters&` to `const Parameters&`

### 3. MODIFIED: `src/SOPParameters.h`
- Now includes `ember/Parameters.h`
- Added `toCoreParameters()` function declaration

### 4. MODIFIED: `src/SOP_EmberBoolean.h`
- Added `getCoreParameters()` method declaration

### 5. MODIFIED: `src/SOP_EmberBoolean.cpp`
- Added `getCoreParameters()` implementation

### 6. MODIFIED: `CMakeLists.txt`
- Added `ember/Parameters.h` to `EMBER_CORE_SOURCES`

---

## Architecture

```
┌─────────────────────────────────────────────────────────┐
│  Houdini SOP Layer (src/)                               │
│  ├─ SOPParameters.h    → Houdini PRM_* definitions      │
│  ├─ SOP_EmberBoolean.h → HDK SOP_Node subclass          │
│  └─ SOP_EmberBoolean.cpp → getCoreParameters()          │
│                          (converts Houdini parms)       │
└─────────────────────────────────────────────────────────┘
                            │
                            │ ember::Parameters
                            ▼
┌─────────────────────────────────────────────────────────┐
│  Core Library (ember/)                                  │
│  ├─ Parameters.h       → Core parameter struct          │
│  ├─ CutterDispatch.h   → Uses ember::Parameters         │
│  └─ (other core files)                                  │
└─────────────────────────────────────────────────────────┘
```

---

## Rebuild Instructions

```cmd
REM Clean and rebuild
cd E:\Tools\EMBRA\EMBER_04\build
rmdir /s /q CMakeFiles CMakeCache.txt cmake_install.cmake

cd .. && mkdir build && cd build

cmake .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release --target install
```

Or use the build script:
```cmd
cd E:\Tools\EMBRA\EMBER_04
build_windows.bat clean install
```

---

## Verification

After the fix, the build should complete without the `'SOPParameters.h' file not found` error.
