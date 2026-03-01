# ═══════════════════════════════════════════════════════════════════════════════
# ADDITIONS TO EXISTING CMakeLists.txt FOR MANIFOLD BACKEND
# ═══════════════════════════════════════════════════════════════════════════════
#
# Add these blocks to your existing E:\Tools\EMBER\CMakeLists.txt
#
# ═══════════════════════════════════════════════════════════════════════════════

# ┌─────────────────────────────────────────────────────────────────────────────
# │ ADD NEAR THE TOP (after project() and option declarations):
# └─────────────────────────────────────────────────────────────────────────────

option(EMBER_HAS_MANIFOLD "Build with Manifold library backend" ON)

# ┌─────────────────────────────────────────────────────────────────────────────
# │ ADD AFTER OPTIONS, BEFORE source file lists:
# │ Fetches elalish/manifold v3.0.1 (MIT license, no TBB needed)
# └─────────────────────────────────────────────────────────────────────────────

if(EMBER_HAS_MANIFOLD)
    message(STATUS "EMBER: Manifold backend ENABLED — fetching library...")

    include(FetchContent)
    FetchContent_Declare(
        manifold
        GIT_REPOSITORY https://github.com/elalish/manifold.git
        GIT_TAG        v3.0.1
        GIT_SHALLOW    TRUE
    )

    # Disable Manifold extras we don't need
    set(MANIFOLD_TEST   OFF CACHE BOOL "" FORCE)
    set(MANIFOLD_PYBIND OFF CACHE BOOL "" FORCE)
    set(MANIFOLD_CBIND  OFF CACHE BOOL "" FORCE)
    set(MANIFOLD_JSBIND OFF CACHE BOOL "" FORCE)
    set(MANIFOLD_EXPORT OFF CACHE BOOL "" FORCE)
    set(MANIFOLD_DEBUG  OFF CACHE BOOL "" FORCE)

    # Serial backend = no TBB dependency from Manifold
    if(NOT DEFINED MANIFOLD_PAR)
        set(MANIFOLD_PAR "NONE" CACHE STRING "" FORCE)
    endif()

    FetchContent_MakeAvailable(manifold)
endif()

# ┌─────────────────────────────────────────────────────────────────────────────
# │ ADD ManifoldBackend.cpp TO EMBER_BACKEND_SOURCES conditionally:
# │ (find your existing EMBER_BACKEND_SOURCES list and add this)
# └─────────────────────────────────────────────────────────────────────────────

# In your existing set(EMBER_BACKEND_SOURCES ...) block:
#   backend/BackendFactory.cpp
#   backend/CherchiBackend_stub.cpp
# ADD after it:
if(EMBER_HAS_MANIFOLD)
    list(APPEND EMBER_BACKEND_SOURCES backend/ManifoldBackend.cpp)
endif()

# ┌─────────────────────────────────────────────────────────────────────────────
# │ ADD compile definition to ember_core:
# │ (find your target_compile_definitions for ember_core, or add new one)
# └─────────────────────────────────────────────────────────────────────────────

target_compile_definitions(ember_core PUBLIC
    $<$<BOOL:${EMBER_HAS_MANIFOLD}>:EMBER_HAS_MANIFOLD>
)

# ┌─────────────────────────────────────────────────────────────────────────────
# │ ADD link to ember_core:
# │ (find your existing target_link_libraries(ember_core ...) and add)
# └─────────────────────────────────────────────────────────────────────────────

if(EMBER_HAS_MANIFOLD)
    target_link_libraries(ember_core PUBLIC manifold)
endif()

# ═══════════════════════════════════════════════════════════════════════════════
# BUILD COMMAND:
#
#   cd E:\Tools\EMBER\build
#   cmake .. -G "Visual Studio 17 2022" -A x64 ^
#       -DCMAKE_BUILD_TYPE=Release ^
#       -DEMBER_HAS_MANIFOLD=ON
#   cmake --build . --config Release
#
# To disable Manifold (if fetch fails or not needed):
#   -DEMBER_HAS_MANIFOLD=OFF
# ═══════════════════════════════════════════════════════════════════════════════
