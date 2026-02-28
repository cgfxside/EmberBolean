# EMBER Boolean Engine — Build Instructions for Houdini

This guide explains how to compile the EMBER Boolean Engine as a Houdini HDK plugin.

---

## Prerequisites

### Required
- **Houdini 20.0+** (tested with 20.0.625)
- **CMake 3.28+**
- **C++17 compatible compiler**:
  - Linux: GCC 9+ or Clang 12+
  - Windows: MSVC 2019+
  - macOS: Xcode 12+

### Optional (but recommended)
- **TBB** (Threading Building Blocks) — for parallel processing
- **Eigen3** — for matrix operations (if using CDT features)

---

## Step 1: Set Up Houdini Environment

### Linux/macOS

```bash
# Source Houdini environment (adjust path for your installation)
source /opt/hfs20.0.625/houdini_setup

# Verify Houdini is in PATH
which houdini
# Output: /opt/hfs20.0.625/bin/houdini
```

### Windows

```powershell
# Run Houdini Command Line Tools from Start Menu
# Or manually set environment:
set HFS=C:\Program Files\Side Effects Software\Houdini 20.0.625
set PATH=%HFS%\bin;%PATH%
set HOUDINI_DSO_PATH=%HFS%\dso
```

---

## Step 2: Configure Build

### Option A: Quick Build (Recommended)

```bash
cd /mnt/okcomputer/output

# Create build directory
mkdir -p build && cd build

# Configure with Houdini HDK
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
make -j$(nproc)  # Linux/macOS
# or
cmake --build . --config Release --parallel  # Windows
```

### Option B: Custom Build with Options

```bash
cd /mnt/okcomputer/output
mkdir -p build && cd build

cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DEMBER_BUILD_TESTS=ON \
    -DCMAKE_INSTALL_PREFIX=$HOME/houdini20.0/dso

make -j$(nproc)
```

---

## Step 3: Install the Plugin

### Automatic Install (CMake)

```bash
cd /mnt/okcomputer/output/build
make install
```

This installs `SOP_EmberBoolean.so` (Linux) or `SOP_EmberBoolean.dll` (Windows) to:
- Linux: `$HFS/dso/sop/`
- Windows: `%HFS%/dso/sop/`
- macOS: `$HFS/dso/sop/`

### Manual Install

```bash
# Linux/macOS
cp build/SOP_EmberBoolean.so $HOME/houdini20.0/dso/sop/

# Windows
copy build\Release\SOP_EmberBoolean.dll %USERPROFILE%\Documents\houdini20.0\dso\sop\
```

---

## Step 4: Verify Installation

### Check Plugin Loads

```bash
# Start Houdini and check for errors
houdini

# Or test from command line:
houdini -c "opcf /obj; opadd -n sop_ember_boolean test_node; opinfo test_node"
```

### Expected Output

In Houdini, you should see:
- **SOP node**: "Ember Boolean" in the TAB menu under "Polygon"
- **Parameters**: Operation (Union/Intersect/Difference), Backend selection

---

## Troubleshooting

### Error: "Could not find Houdini"

```bash
# Set HFS environment variable explicitly
export HFS=/opt/hfs20.0.625
cmake .. -DHoudini_ROOT=$HFS
```

### Error: "__int128 is not supported"

This is expected on MSVC. The code now uses `ember::Int128` which is portable.

```bash
# Ensure you're using the unified ExactPredicates.h
grep -c "using ember::Int128" ember/ExactPredicates.h
# Should output: 1
```

### Error: "undefined reference to `isqrt128`"

```bash
# Check MeshImport.cpp uses Int128
grep "isqrt128" ember/MeshImport.cpp
# Should show Int128 version, not __int128
```

### Plugin Not Appearing in Houdini

1. Check DSO path:
   ```bash
   houdini -c "echo $HOUDINI_DSO_PATH"
   ```

2. Verify plugin loads without errors:
   ```bash
   houdini -c "opadd -n sop_ember_boolean test; opdestroy test"
   ```

3. Check Houdini console for errors (Windows: Alt+Shift+C, Linux/macOS: Terminal)

---

## Build Verification Script

Run this to verify all components are correctly built:

```bash
cd /mnt/okcomputer/output

# Test 1: Check headers compile
./tests/verify_audit.sh

# Test 2: Build and run unit tests
mkdir -p build && cd build
cmake .. -DEMBER_BUILD_TESTS=ON
make -j$(nproc)
./test_predicates
./test_cdt

# Test 3: Verify plugin symbols
nm -D SOP_EmberBoolean.so | grep newSopOperator
```

---

## Platform-Specific Notes

### Linux (CentOS 7/8, Rocky Linux, Ubuntu)

```bash
# Install dependencies
sudo yum install gcc-c++ cmake  # CentOS/RHEL
sudo apt install g++ cmake      # Ubuntu/Debian

# Build
source /opt/hfs20.0.625/houdini_setup
cd /mnt/okcomputer/output
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### Windows (Visual Studio 2019+)

```powershell
# Open "x64 Native Tools Command Prompt for VS 2019"
# Run Houdini Command Line Tools

# Build
mkdir build
cd build
cmake .. -G "Visual Studio 16 2019" -A x64
cmake --build . --config Release --parallel

# Install
copy Release\SOP_EmberBoolean.dll %USERPROFILE%\Documents\houdini20.0\dso\sop\
```

### macOS (Intel/Apple Silicon)

```bash
# Install dependencies
brew install cmake

# Build
source /Applications/Houdini/Houdini20.0.625/Frameworks/Houdini.framework/Versions/Current/Resources/houdini_setup
cd /mnt/okcomputer/output
mkdir build && cd build
cmake ..
make -j$(sysctl -n hw.ncpu)
```

---

## Quick Reference

| Command | Purpose |
|---------|---------|
| `cmake ..` | Configure build |
| `make -j4` | Compile with 4 threads |
| `make install` | Install to Houdini DSO path |
| `nm -D SOP_EmberBoolean.so \| grep newSopOperator` | Verify plugin exports |
| `houdini -c "opadd -n sop_ember_boolean test"` | Test plugin loads |

---

## Next Steps

1. **Test basic boolean operations**: Create two spheres, apply Ember Boolean SOP
2. **Check performance**: Compare with Houdini's native Boolean SOP
3. **Report issues**: File bugs with test geometry and Houdini version
