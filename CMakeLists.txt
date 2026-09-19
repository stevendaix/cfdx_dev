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
    src/cfdx/core/mesh/face.cpp
    src/cfdx/core/mesh/ownership.cpp
    src/cfdx/core/mesh/cell.cpp
    src/cfdx/core/mesh/mesh.cpp
    src/cfdx/core/geometry/face_geometry.cpp
    src/cfdx/core/geometry/cell_geometry.cpp
    src/cfdx/core/geometry/mesh_quality.cpp
    src/cfdx/core/geometry/mesh_validator.cpp
    src/cfdx/core/field/field.cpp
    src/cfdx/core/field/storage.cpp
)

set(CFDX_CORE_HEADERS
    src/cfdx/core/mesh/point.h
    src/cfdx/core/mesh/face.h
    src/cfdx/core/mesh/ownership.h
    src/cfdx/core/mesh/cell.h
    src/cfdx/core/mesh/boundary.h
    src/cfdx/core/mesh/mesh.h
    src/cfdx/core/geometry/face_geometry.h
    src/cfdx/core/geometry/cell_geometry.h
    src/cfdx/core/geometry/mesh_quality.h
    src/cfdx/core/geometry/mesh_validator.h
src/cfdx/core/field/field.h
    src/cfdx/core/field/storage.h
    src/cfdx/core/boundary/boundary_field.h
    src/cfdx/core/numerics/interpolation.h
    src/cfdx/core/numerics/gradient.h
    src/cfdx/core/numerics/divergence.h
    src/cfdx/core/numerics/laplacian.h
    src/cfdx/core/numerics/flux.h
    src/cfdx/core/numerics/integrate.h
    src/cfdx/core/linalg/sparse_matrix.h
    src/cfdx/core/linalg/vector.h
    src/cfdx/core/linalg/linear_system.h
    src/cfdx/core/linalg/cg_solver.h
    src/cfdx/core/linalg/bicgstab_solver.h
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

    add_executable(test_face
        tests/unit/test_face.cpp
    )
    target_link_libraries(test_face PRIVATE cfdx_core)
    target_include_directories(test_face PRIVATE tests)
    add_test(NAME face COMMAND test_face)

    add_executable(test_ownership
        tests/unit/test_ownership.cpp
    )
    target_link_libraries(test_ownership PRIVATE cfdx_core)
    target_include_directories(test_ownership PRIVATE tests)
    add_test(NAME ownership COMMAND test_ownership)

    add_executable(test_cell
        tests/unit/test_cell.cpp
    )
    target_link_libraries(test_cell PRIVATE cfdx_core)
    target_include_directories(test_cell PRIVATE tests)
    add_test(NAME cell COMMAND test_cell)

    add_executable(test_boundary
        tests/unit/test_boundary.cpp
    )
    target_link_libraries(test_boundary PRIVATE cfdx_core)
    target_include_directories(test_boundary PRIVATE tests)
    add_test(NAME boundary COMMAND test_boundary)

    add_executable(test_mesh
        tests/unit/test_mesh.cpp
    )
    target_link_libraries(test_mesh PRIVATE cfdx_core)
    target_include_directories(test_mesh PRIVATE tests)
    add_test(NAME mesh COMMAND test_mesh)

    add_executable(test_face_geometry
        tests/numerical/test_face_geometry.cpp
    )
    target_link_libraries(test_face_geometry PRIVATE cfdx_core)
    target_include_directories(test_face_geometry PRIVATE tests)
    add_test(NAME face_geometry COMMAND test_face_geometry)

    add_executable(test_cell_geometry
        tests/numerical/test_cell_geometry.cpp
    )
    target_link_libraries(test_cell_geometry PRIVATE cfdx_core)
    target_include_directories(test_cell_geometry PRIVATE tests)
    add_test(NAME cell_geometry COMMAND test_cell_geometry)

    add_executable(test_mesh_quality
        tests/numerical/test_mesh_quality.cpp
    )
    target_link_libraries(test_mesh_quality PRIVATE cfdx_core)
    target_include_directories(test_mesh_quality PRIVATE tests)
    add_test(NAME mesh_quality COMMAND test_mesh_quality)

    add_executable(test_field
        tests/unit/test_field.cpp
    )
    target_link_libraries(test_field PRIVATE cfdx_core)
    target_include_directories(test_field PRIVATE tests)
    add_test(NAME field COMMAND test_field)

    add_executable(test_field_storage
        tests/unit/test_field_storage.cpp
    )
    target_link_libraries(test_field_storage PRIVATE cfdx_core)
    target_include_directories(test_field_storage PRIVATE tests)
    add_test(NAME field_storage COMMAND test_field_storage)

    add_executable(test_boundary_field
        tests/unit/test_boundary_field.cpp
    )
    target_link_libraries(test_boundary_field PRIVATE cfdx_core)
    target_include_directories(test_boundary_field PRIVATE tests)
    add_test(NAME boundary_field COMMAND test_boundary_field)

    add_executable(test_interpolation
        tests/unit/test_interpolation.cpp
    )
    target_link_libraries(test_interpolation PRIVATE cfdx_core)
    target_include_directories(test_interpolation PRIVATE tests)
    add_test(NAME interpolation COMMAND test_interpolation)

    add_executable(test_gradient
        tests/unit/test_gradient.cpp
    )
    target_link_libraries(test_gradient PRIVATE cfdx_core)
    target_include_directories(test_gradient PRIVATE tests)
    add_test(NAME gradient COMMAND test_gradient)

    add_executable(test_divergence
        tests/unit/test_divergence.cpp
    )
    target_link_libraries(test_divergence PRIVATE cfdx_core)
    target_include_directories(test_divergence PRIVATE tests)
    add_test(NAME divergence COMMAND test_divergence)

    add_executable(test_laplacian
        tests/unit/test_laplacian.cpp
    )
    target_link_libraries(test_laplacian PRIVATE cfdx_core)
    target_include_directories(test_laplacian PRIVATE tests)
    add_test(NAME laplacian COMMAND test_laplacian)

    add_executable(test_flux
        tests/unit/test_flux.cpp
    )
    target_link_libraries(test_flux PRIVATE cfdx_core)
    target_include_directories(test_flux PRIVATE tests)
    add_test(NAME flux COMMAND test_flux)

    add_executable(test_integrate
        tests/unit/test_integrate.cpp
    )
    target_link_libraries(test_integrate PRIVATE cfdx_core)
    target_include_directories(test_integrate PRIVATE tests)
    add_test(NAME integrate COMMAND test_integrate)

    add_executable(test_sparse_matrix
        tests/unit/test_sparse_matrix.cpp
    )
    target_link_libraries(test_sparse_matrix PRIVATE cfdx_core)
    target_include_directories(test_sparse_matrix PRIVATE tests)
    add_test(NAME sparse_matrix COMMAND test_sparse_matrix)

    add_executable(test_vector
        tests/unit/test_vector.cpp
    )
    target_link_libraries(test_vector PRIVATE cfdx_core)
    target_include_directories(test_vector PRIVATE tests)
    add_test(NAME vector COMMAND test_vector)

    add_executable(test_cg_solver
        tests/unit/test_cg_solver.cpp
    )
    target_link_libraries(test_cg_solver PRIVATE cfdx_core)
    target_include_directories(test_cg_solver PRIVATE tests)
    add_test(NAME cg_solver COMMAND test_cg_solver)

    add_executable(test_bicgstab_solver
        tests/unit/test_bicgstab_solver.cpp
    )
    target_link_libraries(test_bicgstab_solver PRIVATE cfdx_core)
    target_include_directories(test_bicgstab_solver PRIVATE tests)
    add_test(NAME bicgstab_solver COMMAND test_bicgstab_solver)
endif()