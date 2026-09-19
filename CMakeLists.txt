cmake_minimum_required(VERSION 3.16)

project(CFDX
    VERSION 0.7
    LANGUAGES CXX C
    DESCRIPTION "CFDX - High-Performance General-Purpose CFD Framework"
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------------------------------------------------------------------
# Build options
# ---------------------------------------------------------------------------
option(CFDX_BUILD_TESTS "Build CFDX tests" ON)
option(CFDX_BUILD_PYTHON "Build CFDX Python bindings" OFF)
option(CFDX_ENABLE_GPU "Enable CUDA backend (requires full CUDA toolkit)" OFF)
option(CFDX_ENABLE_MPI "Enable MPI support" ON)
option(CFDX_USE_HDF5_SERIAL "Use locally-built HDF5 serial" ON)

# ---------------------------------------------------------------------------
# Compiler flags
# ---------------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-O2 -DNDEBUG)
    endif()
endif()

# ---------------------------------------------------------------------------
# HDF5 detection (before any target uses it)
# ---------------------------------------------------------------------------
if(CFDX_USE_HDF5_SERIAL)
    set(CFDX_HDF5_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/hdf5")
endif()

find_package(HDF5 QUIET)
if(HDF5_FOUND)
    message(STATUS "Using system HDF5 ${HDF5_VERSION}")
elseif(EXISTS "${CFDX_HDF5_DIR}/include/H5public.h")
    message(STATUS "Using locally-built HDF5 (static) at ${CFDX_HDF5_DIR}")
    set(HDF5_INCLUDE_DIRS "${CFDX_HDF5_DIR}/include")
    set(HDF5_C_LIBRARIES "${CFDX_HDF5_DIR}/lib/libhdf5.a")
    set(HDF5_HL_LIBRARIES "${CFDX_HDF5_DIR}/lib/libhdf5_hl.a")
else()
    message(WARNING "HDF5 not found. HDF5-dependent features disabled.")
endif()

# ---------------------------------------------------------------------------
# CFDX core library
# ---------------------------------------------------------------------------
set(CFDX_CORE_SOURCES
    src/cfdx/core/mesh/point.cpp
)

set(CFDX_CORE_HEADERS
    src/cfdx/core/mesh/point.h
)

add_library(cfdx_core STATIC ${CFDX_CORE_SOURCES} ${CFDX_CORE_HEADERS})
target_include_directories(cfdx_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

if(HDF5_FOUND AND HDF5_C_LIBRARIES)
    target_include_directories(cfdx_core PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(cfdx_core PUBLIC ${HDF5_C_LIBRARIES} ${HDF5_HL_LIBRARIES})
endif()

if(CFDX_ENABLE_MPI)
    find_package(MPI QUIET)
    if(MPI_FOUND)
        target_link_libraries(cfdx_core PUBLIC MPI::MPI_CXX)
    endif()
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(CFDX_BUILD_TESTS)
    enable_testing()

    add_executable(test_point
        tests/unit/test_point.cpp
    )
    target_link_libraries(test_point PRIVATE cfdx_core)
    target_include_directories(test_point PRIVATE tests)
    add_test(NAME point COMMAND test_point)
endif()