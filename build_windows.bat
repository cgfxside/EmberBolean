@echo off
REM EMBER Boolean Engine — Windows Build Script for Houdini
REM Supports: Houdini 20.0+, Houdini 21.0+
REM Usage: build_windows.bat [clean] [install]

setlocal EnableDelayedExpansion

echo ============================================
echo EMBER Boolean Engine — Windows Build
echo ============================================
echo.

REM Check for Houdini environment
if not defined HFS (
    echo ERROR: HFS environment variable not set.
    echo.
    echo Please run this script from "Houdini Command Line Tools"
    echo or set HFS manually:
    echo.
    echo   set HFS=C:\Program Files\Side Effects Software\Houdini 21.0.xxx
    echo.
    exit /b 1
)

REM Detect Houdini version
for /f "tokens=3 delims=\." %%a in ("%HFS%") do set HOUDINI_MAJOR=%%a
for /f "tokens=4 delims=\." %%a in ("%HFS%") do set HOUDINI_MINOR=%%a

echo Houdini: %HFS%
echo Version: %HOUDINI_MAJOR%.%HOUDINI_MINOR%
echo.

REM Check minimum version
if %HOUDINI_MAJOR% LSS 20 (
    echo WARNING: Houdini %HOUDINI_MAJOR%.%HOUDINI_MINOR% may not be fully supported.
    echo Minimum recommended: Houdini 20.0 or 21.0
    echo.
)

REM Check for Visual Studio
where cl >nul 2>&1
if errorlevel 1 (
    echo ERROR: Visual Studio compiler not found in PATH.
    echo Please run from "Houdini Command Line Tools" or "x64 Native Tools Command Prompt"
    exit /b 1
)

for /f "tokens=1-4" %%a in ('cl 2^>^&1 ^| findstr /C:"Version"') do (
    set MSVC_VERSION=%%d
)
echo MSVC: %MSVC_VERSION%
echo.

REM Parse arguments
set DO_CLEAN=0
set DO_INSTALL=0

for %%a in (%*) do (
    if "%%a"=="clean" set DO_CLEAN=1
    if "%%a"=="install" set DO_INSTALL=1
)

REM Clean if requested
if %DO_CLEAN%==1 (
    echo Cleaning build directory...
    if exist build rmdir /s /q build
    echo Clean complete.
    echo.
)

REM Create build directory
if not exist build mkdir build
cd build

REM Configure with CMake
echo Configuring with CMake...
echo.

cmake .. -G "Visual Studio 17 2022" -A x64 ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DEMBER_BUILD_TESTS=ON ^
    -DHoudini_ROOT="%HFS%" ^
    2>&1

if errorlevel 1 (
    echo.
    echo ERROR: CMake configuration failed!
    echo.
    echo Common causes:
    echo   - Houdini not installed or HFS not set
    echo   - Visual Studio 2022 not installed
    echo   - CMake not in PATH
    echo.
    exit /b 1
)

echo.
echo Configuration complete.
echo.

REM Build
echo Building EMBER Boolean Engine...
echo.

cmake --build . --config Release --parallel

if errorlevel 1 (
    echo.
    echo ERROR: Build failed!
    echo.
    exit /b 1
)

echo.
echo Build complete.
echo.

REM Run tests
echo Running tests...
echo.

if exist Release\test_predicates.exe (
    echo [test_predicates]
    Release\test_predicates.exe
    echo.
)
if exist Release\test_cdt.exe (
    echo [test_cdt]
    Release\test_cdt.exe
    echo.
)

REM Install if requested
if %DO_INSTALL%==1 (
    echo.
    echo Installing plugin...
    echo.
    
    cmake --install . --config Release
    
    if errorlevel 1 (
        echo.
        echo ERROR: Installation failed!
        echo.
        exit /b 1
    )
    
    echo.
    echo Plugin installed successfully.
)

echo.
echo ============================================
echo Build Complete!
echo ============================================
echo.
echo Houdini: %HOUDINI_MAJOR%.%HOUDINI_MINOR%
echo MSVC: %MSVC_VERSION%
echo.

if %DO_INSTALL%==1 (
    echo Install path: %USERPROFILE%\Documents\houdini%HOUDINI_MAJOR%.%HOUDINI_MINOR%\dso\sop\
    echo.
)

echo To test in Houdini:
echo   1. Start Houdini %HOUDINI_MAJOR%.%HOUDINI_MINOR%
echo   2. Create a geometry object
echo   3. Press TAB and type 'ember'
echo   4. Select 'Ember Boolean' SOP
echo.

cd ..
endlocal
