# EMBER Boolean Engine - CMake Build Configuration
# EMBER - Exact Mesh Boolean Engine for Robust operations
#
# COMPATIBILITY:
#   - Houdini 20.0+ (tested with 20.0.625)
#   - Houdini 21.0+ (tested with 21.0.XXX)
#   - Windows: MSVC 2019, 2022
#   - Linux: GCC 9+, Clang 12+
#   - macOS: Xcode 12+

cmake_minimum_required(VERSION 3.28)
project(EmberBoolean VERSION 1.1.1 LANGUAGES CXX)

# C++ Standard
# Houdini 20/21 both use C++17 as minimum
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build options
option(EMBER_BUILD_TESTS "Build unit tests" ON)
option(EMBER_ENABLE_CUDA "Enable CUDA broad-phase" OFF)
option(EMBER_USE_BOOST_MULTIPRECISION "Use boost::multiprecision" OFF)
option(EMBER_HAS_MANIFOLD "Build with Manifold library backend" ON)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

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

# Force Manifold to build as static lib (its CMake overrides the global setting)
    set(BUILD_SHARED_LIBS_SAVED ${BUILD_SHARED_LIBS})
    set(BUILD_SHARED_LIBS OFF)
    FetchContent_MakeAvailable(manifold)
    set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS_SAVED})
endif()

# Find Houdini HDK
# Supports Houdini 20.x and 21.x
if(DEFINED ENV{HFS})
    list(APPEND CMAKE_PREFIX_PATH "$ENV{HFS}/toolkit/cmake")
endif()
find_package(Houdini REQUIRED)

# Houdini version check (informative only)
if(Houdini_VERSION VERSION_LESS "20.0")
    message(WARNING "Houdini ${Houdini_VERSION} may not be fully supported. Minimum recommended: 20.0")
endif()

# Find third-party libraries (via vcpkg or system)
find_package(Eigen3 QUIET)
find_package(TBB QUIET)
find_package(embree 4 QUIET)

# Source files
set(EMBER_CORE_SOURCES
    ember/IntegerTypes.h
    ember/ExactPredicates.h
    ember/ExactPredicates_Indirect.h
    ember/Parameters.h
    ember/PolygonSoup.h
    ember/Diagnostics.h
    ember/MeshImport.h
    ember/MeshImport.cpp
    ember/CutterDispatch.h
    ember/CutterDispatch.cpp
    ember/CutterMeshGenerator.h
)

set(EMBER_BACKEND_SOURCES
    backend/IBooleanBackend.h
    backend/BackendFactory.cpp
    # BUG 4 FIX (ODR VIOLATION): CherchiBackend.cpp is REMOVED from the build.
    #
    # PROBLEM: CherchiBackend.cpp defines a complete 'class CherchiBackend' at
    # line 1133 inside the ember namespace (not in an anonymous namespace).
    # CherchiBackend_stub.cpp implements the methods of the same class declared
    # in CherchiBackend.h. Compiling both produces TWO definitions of
    # ember::CherchiBackend::execute() — a hard ODR violation. The linker picks
    # one at random → undefined behavior → crash or incorrect results.
    #
    # SOLUTION: Only CherchiBackend_stub.cpp is compiled. It provides a safe
    # stub that returns success=true with empty geometry until the full
    # implementation is complete.
    #
    # When the full CDT implementation in CherchiBackend.cpp is ready, it must
    # be refactored to implement the methods declared in CherchiBackend.h
    # (not define a new class), and then replace CherchiBackend_stub.cpp.
    #
    # backend/CherchiBackend.cpp   ← DO NOT add back without refactoring
    backend/CherchiBackend_stub.cpp
)

if(EMBER_HAS_MANIFOLD)
    list(APPEND EMBER_BACKEND_SOURCES backend/ManifoldBackend.cpp)
endif()

set(EMBER_SOP_SOURCES
    src/SOPParameters.h
    src/SOPParameters.cpp
    src/GU_EmberConverter.h
    src/GU_EmberConverter.cpp
    src/SOP_EmberBoolean.h
    src/SOP_EmberBoolean.cpp
    src/EmberPlugin.cpp 
)

# Create core library (header-only + MeshImport.cpp)
add_library(ember_core STATIC
    ${EMBER_CORE_SOURCES}
    ${EMBER_BACKEND_SOURCES}
)

target_include_directories(ember_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

if(TBB_FOUND)
    target_link_libraries(ember_core PUBLIC TBB::tbb)
endif()

target_compile_definitions(ember_core PUBLIC
    $<$<BOOL:${EMBER_HAS_MANIFOLD}>:EMBER_HAS_MANIFOLD>
)
if(EMBER_HAS_MANIFOLD)
    target_link_libraries(ember_core PUBLIC manifold)
endif()

if(EMBER_USE_BOOST_MULTIPRECISION)
    find_package(Boost REQUIRED COMPONENTS multiprecision)
    target_link_libraries(ember_core PUBLIC Boost::boost)
    target_compile_definitions(ember_core PUBLIC EMBER_USE_BOOST_MULTIPRECISION)
endif()

# Compiler-specific options
if(MSVC)
    target_compile_options(ember_core PRIVATE /W4 /WX- /permissive-)
    target_compile_options(ember_core PRIVATE $<$<CONFIG:Release>:/O2 /fp:strict>)
else()
    target_compile_options(ember_core PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(ember_core PRIVATE $<$<CONFIG:Release>:-O3 -march=native>)
endif()

# Create SOP plugin (DSO)
houdini_generate_proto_headers(SOURCES ${EMBER_SOP_SOURCES})

add_library(SOP_EmberBoolean SHARED ${EMBER_SOP_SOURCES})

target_link_libraries(SOP_EmberBoolean PRIVATE
    ember_core
    Houdini
)

houdini_configure_target(SOP_EmberBoolean)

# Install
install(TARGETS SOP_EmberBoolean
    DESTINATION "${HOUDINI_DSO_DIR}/sop"
)

if(EMBER_HAS_MANIFOLD)
    install(FILES $<TARGET_FILE:manifold> DESTINATION ${CMAKE_INSTALL_PREFIX}/dso)
endif()

# Tests
if(EMBER_BUILD_TESTS)
    enable_testing()
    
    add_executable(test_predicates
        tests/test_predicates.cpp
    )
    target_link_libraries(test_predicates PRIVATE ember_core)
    add_test(NAME predicates COMMAND test_predicates)
    
    add_executable(test_cdt
        tests/test_cdt.cpp
    )
    target_link_libraries(test_cdt PRIVATE ember_core)
    add_test(NAME cdt COMMAND test_cdt)
endif()

# Print configuration summary
message(STATUS "")
message(STATUS "EMBER Boolean Engine Configuration:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build Tests: ${EMBER_BUILD_TESTS}")
message(STATUS "  CUDA Support: ${EMBER_ENABLE_CUDA}")
message(STATUS "  Boost Multiprecision: ${EMBER_USE_BOOST_MULTIPRECISION}")
message(STATUS "  Houdini: ${Houdini_VERSION}")
message(STATUS "")
