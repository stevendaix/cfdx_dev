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
# Sanitizer build type
# ---------------------------------------------------------------------------
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG")

# Sanitizer build type (DebugSanitizers)
if(NOT CMAKE_BUILD_TYPE STREQUAL "DebugSanitizers")
    set(CMAKE_CXX_FLAGS_DEBUGSANITIZERS "-g -O1 -fsanitize=address,undefined -fno-omit-frame-pointer")
    set(CMAKE_C_FLAGS_DEBUGSANITIZERS "-g -O1 -fsanitize=address,undefined -fno-omit-frame-pointer")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUGSANITIZERS "-fsanitize=address,undefined")
    set(CMAKE_SHARED_LINKER_FLAGS_DEBUGSANITIZERS "-fsanitize=address,undefined")
    set(CMAKE_MODULE_LINKER_FLAGS_DEBUGSANITIZERS "-fsanitize=address,undefined")
endif()

# ---------------------------------------------------------------------------
# Compiler flags
# ---------------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion)
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug" AND NOT CMAKE_BUILD_TYPE STREQUAL "DebugSanitizers")
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
    # System HDF5 provides these via find_package
    if(NOT HDF5_INCLUDE_DIRS)
        set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
    endif()
    if(NOT HDF5_C_LIBRARIES)
        set(HDF5_C_LIBRARIES ${HDF5_C_LIBRARIES})
    endif()
    if(NOT HDF5_HL_LIBRARIES)
        set(HDF5_HL_LIBRARIES ${HDF5_HL_LIBRARIES})
    endif()
elseif(EXISTS "${CFDX_HDF5_DIR}/include/H5public.h")
    message(STATUS "Using locally-built HDF5 (static) at ${CFDX_HDF5_DIR}")
    set(HDF5_INCLUDE_DIRS "${CFDX_HDF5_DIR}/include")
    set(HDF5_C_LIBRARIES "${CFDX_HDF5_DIR}/lib/libhdf5.a")
    set(HDF5_HL_LIBRARIES "${CFDX_HDF5_DIR}/lib/libhdf5_hl.a")
else()
    message(WARNING "HDF5 not found. HDF5-dependent features disabled.")
endif()

# ---------------------------------------------------------------------------
# CFDX core library - Sources organized by module
# ---------------------------------------------------------------------------

# Mesh module
set(CFDX_MESH_SOURCES
    src/cfdx/core/mesh/point.cpp
    src/cfdx/core/mesh/face.cpp
    src/cfdx/core/mesh/ownership.cpp
    src/cfdx/core/mesh/cell.cpp
    src/cfdx/core/mesh/mesh.cpp
    src/cfdx/core/mesh/boundary.cpp
)

# Importers module
set(CFDX_IMPORTER_SOURCES
    src/cfdx/core/importers/opensfoam_importer.cpp
    src/cfdx/core/importers/gmsh_importer.cpp
)

add_executable(opensfoam_importer src/cfdx/core/importers/opensfoam_importer.cpp)
add_executable(gmsh_importer src/cfdx/core/importers/gmsh_importer.cpp)

# Geometry module
set(CFDX_GEOMETRY_SOURCES
    src/cfdx/core/geometry/face_geometry.cpp
    src/cfdx/core/geometry/cell_geometry.cpp
    src/cfdx/core/geometry/mesh_quality.cpp
    src/cfdx/core/geometry/mesh_validator.cpp
)

# Field module
set(CFDX_FIELD_SOURCES
    src/cfdx/core/field/field.cpp
    src/cfdx/core/field/storage.cpp
)

# IO module
set(CFDX_IO_SOURCES
    src/cfdx/io/hdf5/hdf5_writer.cpp
    src/cfdx/io/hdf5/hdf5_reader.cpp
    src/cfdx/io/vtu/vtu_writer.cpp
)

set(CFDX_CORE_SOURCES
    ${CFDX_MESH_SOURCES}
    ${CFDX_GEOMETRY_SOURCES}
    ${CFDX_FIELD_SOURCES}
    ${CFDX_IO_SOURCES}
)

# Headers (for IDE visibility, not compiled)
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
    src/cfdx/core/numerics/source_term.h
    src/cfdx/core/numerics/temporal.h
    src/cfdx/physics/equation_of_state.h
    src/cfdx/physics/transport_models.h
    src/cfdx/physics/multi_component.h
    src/cfdx/core/linalg/sparse_matrix.h
    src/cfdx/core/linalg/vector.h
    src/cfdx/core/linalg/linear_system.h
    src/cfdx/core/linalg/cg_solver.h
    src/cfdx/core/linalg/bicgstab_solver.h
    src/cfdx/io/hdf5/hdf5_writer.h
    src/cfdx/io/hdf5/hdf5_reader.h
    src/cfdx/io/vtu/vtu_writer.h
)

add_library(cfdx_core STATIC ${CFDX_CORE_SOURCES})
target_include_directories(cfdx_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

if(HDF5_FOUND OR EXISTS "${CFDX_HDF5_DIR}/include/H5public.h")
    target_include_directories(cfdx_core PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(cfdx_core PUBLIC ${HDF5_C_LIBRARIES} ${HDF5_HL_LIBRARIES} z)
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

    # Helper function for tests
    macro(add_cfdx_test test_name source_file)
        add_executable(${test_name} ${source_file})
        target_link_libraries(${test_name} PRIVATE cfdx_core)
        target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        add_test(NAME ${test_name} COMMAND ${test_name})
    endmacro()

    # Unit tests - Core
    add_cfdx_test(test_point tests/unit/test_point.cpp)
    add_cfdx_test(test_face tests/unit/test_face.cpp)
    add_cfdx_test(test_ownership tests/unit/test_ownership.cpp)
    add_cfdx_test(test_cell tests/unit/test_cell.cpp)
    add_cfdx_test(test_boundary tests/unit/test_boundary.cpp)
    add_cfdx_test(test_mesh tests/unit/test_mesh.cpp)
    add_cfdx_test(test_field tests/unit/test_field.cpp)
    add_cfdx_test(test_field_storage tests/unit/test_field_storage.cpp)
    add_cfdx_test(test_boundary_field tests/unit/test_boundary_field.cpp)
    add_cfdx_test(test_interpolation tests/unit/test_interpolation.cpp)
    add_cfdx_test(test_gradient tests/unit/test_gradient.cpp)
    add_cfdx_test(test_divergence tests/unit/test_divergence.cpp)
    add_cfdx_test(test_flux tests/unit/test_flux.cpp)
    add_cfdx_test(test_integrate tests/unit/test_integrate.cpp)

    # Unit tests - Linear Algebra
    add_cfdx_test(test_sparse_matrix tests/unit/test_sparse_matrix.cpp)
    add_cfdx_test(test_vector tests/unit/test_vector.cpp)
    add_cfdx_test(test_linear_system tests/unit/test_linear_system.cpp)
    add_cfdx_test(test_cg_solver tests/unit/test_cg_solver.cpp)
    add_cfdx_test(test_bicgstab_solver tests/unit/test_bicgstab_solver.cpp)

    # Unit tests - Numerical
    add_cfdx_test(test_laplacian tests/unit/test_laplacian.cpp)

    # Numerical tests - Geometry
    add_cfdx_test(test_face_geometry tests/numerical/test_face_geometry.cpp)
    add_cfdx_test(test_cell_geometry tests/numerical/test_cell_geometry.cpp)
    add_cfdx_test(test_mesh_quality tests/numerical/test_mesh_quality.cpp)
    add_cfdx_test(test_source_term tests/numerical/test_source_term.cpp)
    add_cfdx_test(test_temporal tests/numerical/test_temporal.cpp)
    # add_cfdx_test(test_geometry_cache tests/unit/test_geometry_cache.cpp) 
    # add_cfdx_test(test_mpi_partition tests/unit/test_mpi_partition.cpp)  

    # HDF5 tests
    add_cfdx_test(test_hdf5_writer tests/unit/test_hdf5_writer.cpp)
    add_cfdx_test(test_hdf5_roundtrip tests/unit/test_hdf5_roundtrip.cpp)
    add_cfdx_test(test_field_roundtrip tests/unit/test_field_roundtrip.cpp)

endif()
