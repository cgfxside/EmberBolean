# EMBER Boolean Engine — Build Instructions

Quick reference for building EMBER on all platforms.

---

## Platform Quick Reference

| Platform | Compiler | Houdini | Command |
|----------|----------|---------|---------|
| **Windows** | VS 2022 | 21.0 | `build_windows.bat install` |
| **Windows** | VS 2019/2022 | 20.0 | `build_windows.bat install` |
| **Linux** | GCC 9+ | 20/21 | `./build_for_houdini.sh install` |
| **macOS** | Clang 12+ | 20/21 | `./build_for_houdini.sh install` |

---

## Windows (Houdini 20/21)

### Prerequisites
- Houdini 20.0+ or 21.0+
- Visual Studio 2022 (required for H21)
- CMake 3.28+

### Build
```cmd
REM Open "Houdini Command Line Tools" from Start Menu
cd C:\path\to\EMBER
build_windows.bat install
```

### Output
- Plugin: `%USERPROFILE%\Documents\houdini21.0\dso\sop\SOP_EmberBoolean.dll`

---

## Linux (Houdini 20/21)

### Prerequisites
- Houdini 20.0+ or 21.0+
- GCC 9+ or Clang 12+
- CMake 3.28+

### Build
```bash
# Source Houdini environment
source /opt/hfs21.0.xxx/houdini_setup

# Build
cd /path/to/EMBER
./build_for_houdini.sh install
```

### Output
- Plugin: `$HOME/houdini21.0/dso/sop/SOP_EmberBoolean.so`

---

## macOS (Houdini 20/21)

### Prerequisites
- Houdini 20.0+ or 21.0+
- Xcode 12+ or Command Line Tools
- CMake 3.28+

### Build
```bash
# Source Houdini environment
source /Applications/Houdini/Houdini21.0.xxx/Frameworks/Houdini.framework/Versions/Current/Resources/houdini_setup

# Build
cd /path/to/EMBER
./build_for_houdini.sh install
```

### Output
- Plugin: `$HOME/Library/Preferences/houdini/21.0/dso/sop/SOP_EmberBoolean.dylib`

---

## Manual Build (All Platforms)

```bash
# 1. Set Houdini environment
export HFS=/opt/hfs21.0.xxx  # Linux/macOS
# set HFS=C:\Program Files\Side Effects Software\Houdini 21.0.xxx  # Windows

# 2. Create build directory
mkdir build && cd build

# 3. Configure
cmake .. -DCMAKE_BUILD_TYPE=Release

# 4. Build
cmake --build . --config Release --parallel

# 5. Install
cmake --install . --config Release
```

---

## CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `CMAKE_BUILD_TYPE` | Release | Build configuration |
| `EMBER_BUILD_TESTS` | ON | Build unit tests |
| `EMBER_ENABLE_CUDA` | OFF | Enable CUDA broad-phase |

---

## Verification

### Test 1: Plugin Loads
In Houdini:
1. Create geometry (Sphere)
2. Press TAB, type "ember"
3. Select "Ember Boolean"

### Test 2: Basic Boolean
1. Create two overlapping spheres
2. Connect to Ember Boolean
3. Set Operation to "Union"
4. Verify merged output

### Test 3: Unit Tests
```bash
# Linux/macOS
./build/test_predicates
./build/test_cdt

# Windows
.\build\Release\test_predicates.exe
.\build\Release\test_cdt.exe
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "HFS not set" | Run from Houdini Command Line Tools |
| "CMake not found" | Install CMake 3.28+ |
| "Compiler not found" | Install VS 2022 (Windows) or GCC 9+ (Linux) |
| "Plugin not appearing" | Check `%USERPROFILE%\Documents\houdini21.0\dso\sop\` |

---

## Documentation

- `BUILD_WINDOWS.md` — Detailed Windows guide
- `BUILD_HOUDINI21_WINDOWS.md` — Houdini 21 specific
- `HOUDINI21_COMPATIBILITY.md` — Compatibility matrix
