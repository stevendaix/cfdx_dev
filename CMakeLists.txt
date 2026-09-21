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
option(CFDX_BUILD_RUNTIME_TESTS "Build execution/GPU-OOC foundation tests" ON)

if(CFDX_ENABLE_GPU)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

find_package(Python3 COMPONENTS Interpreter QUIET)

# ---------------------------------------------------------------------------
# Sanitizer build type
# ---------------------------------------------------------------------------
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG")

# Sanitizer build type (DebugSanitizers)
if(CMAKE_BUILD_TYPE STREQUAL "DebugSanitizers")
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
set(CFDX_LINALG_SOURCES
    src/cfdx/core/linalg/hypre_amg.cpp
)

# Memory sources
set(CFDX_MEMORY_SOURCES
    src/cfdx/core/memory/memory_planner.cpp
    src/cfdx/runtime/gpu/gpu_kernels.cpp
)

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
    src/cfdx/io/openfoam/openfoam_importer.cpp
    src/cfdx/io/gmsh/gmsh_importer.cpp
    src/cfdx/io/mesh/mesh_importer.cpp
)

# Physics sources
set(CFDX_PHYSICS_SOURCES
    src/cfdx/physics/equation_of_state.h
    src/cfdx/physics/multi_component.h
    src/cfdx/physics/source_term_linearization.h
    src/cfdx/physics/transport_models.h
    src/cfdx/physics/temporal.h
    src/cfdx/physics/low_storage_time_integration.h
    src/cfdx/physics/adaptive_cfl.h
    src/cfdx/physics/advanced_convergence.h
    src/cfdx/physics/local_time_stepping.h
    src/cfdx/physics/low_mach.h
    src/cfdx/physics/compressible_flux.h
    src/cfdx/physics/radiation_models.h
    src/cfdx/physics/boussinesq.h
    src/cfdx/physics/spalart_allmaras.h
    src/cfdx/physics/radiation_streaming.h
    src/cfdx/thermodynamics/thermo_state.h
    src/cfdx/thermodynamics/thermo_cache.h
)

# Thermodynamics sources (P2 reorganization)
set(CFDX_THERMODYNAMICS_SOURCES
    src/cfdx/thermodynamics/equation_of_state.h
)

# Chemistry sources (P2 reorganization)
set(CFDX_CHEMISTRY_SOURCES
    src/cfdx/chemistry/source_terms.h
)

# Transport sources (P2 reorganizization)
set(CFDX_TRANSPORT_SOURCES
    src/cfdx/transport/viscosity/viscosity.h
    src/cfdx/transport/conductivity/conductivity.h
    src/cfdx/transport/diffusivity/diffusivity.h
)

set(CFDX_CORE_SOURCES
    ${CFDX_MESH_SOURCES}
    ${CFDX_GEOMETRY_SOURCES}
    ${CFDX_FIELD_SOURCES}
    ${CFDX_IO_SOURCES}
    ${CFDX_PHYSICS_SOURCES}
    ${CFDX_THERMODYNAMICS_SOURCES}
    ${CFDX_CHEMISTRY_SOURCES}
    ${CFDX_TRANSPORT_SOURCES}
    ${CFDX_LINALG_SOURCES}
    ${CFDX_MEMORY_SOURCES}
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
    src/cfdx/core/numerics/convection.h
    src/cfdx/physics/incompressible.h
    src/cfdx/physics/pressure_velocity.h
    src/cfdx/physics/pressure_velocity_algorithms.h
    src/cfdx/physics/turbulence.h
    src/cfdx/physics/thermal.h
    src/cfdx/physics/radiation.h
    src/cfdx/physics/m1_m4_models.h
    src/cfdx/physics/solver_control.h
    src/cfdx/physics/finite_volume_transport.h
    src/cfdx/physics/steady_incompressible_solver.h
    src/cfdx/physics/turbulence_transport.h
    src/cfdx/physics/turbulence_solver.h
    src/cfdx/physics/sst_solver.h
    src/cfdx/physics/energy_solver.h
    src/cfdx/physics/cht_solver.h
    src/cfdx/physics/radiation_solver.h
    src/cfdx/physics/equation_of_state.h
    src/cfdx/physics/transport_models.h
    src/cfdx/physics/multi_component.h
    src/cfdx/core/linalg/sparse_matrix.h
    src/cfdx/core/linalg/vector.h
    src/cfdx/core/linalg/linear_system.h
    src/cfdx/core/linalg/cg_solver.h
    src/cfdx/core/linalg/bicgstab_solver.h
    src/cfdx/core/linalg/gmres_solver.h
    src/cfdx/core/linalg/linear_operator.h
    src/cfdx/core/linalg/solver_workspace.h
    src/cfdx/core/linalg/matrix_free_fv_operator.h
    src/cfdx/core/linalg/amg_preconditioner.h
    src/cfdx/core/linalg/mixed_precision.h
    src/cfdx/core/linalg/chebyshev_smoother.h
    src/cfdx/core/linalg/fv_operator.h
    src/cfdx/core/linalg/krylov_controls.h
    src/cfdx/core/linalg/communication_avoiding.h
    src/cfdx/core/memory/reuse_pool.h
    src/cfdx/core/memory/field_lifetime.h
    src/cfdx/core/linalg/mpi_overlap.h
    src/cfdx/core/linalg/block_schur.h
    src/cfdx/core/mesh/solver_mesh.h
    src/cfdx/core/mesh/sfc_ordering.h
    src/cfdx/thermodynamics/equation_of_state.h
    src/cfdx/chemistry/source_terms.h
    src/cfdx/transport/viscosity/viscosity.h
    src/cfdx/transport/conductivity/conductivity.h
    src/cfdx/transport/diffusivity/diffusivity.h
    src/cfdx/io/hdf5/hdf5_writer.h
    src/cfdx/io/hdf5/hdf5_reader.h
    src/cfdx/io/vtu/vtu_writer.h
)

add_library(cfdx_core STATIC ${CFDX_CORE_SOURCES})
target_include_directories(cfdx_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_compile_definitions(cfdx_core PRIVATE CFDX_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
if(Python3_Interpreter_FOUND)
    target_compile_definitions(cfdx_core PRIVATE CFDX_PYTHON_EXECUTABLE="${Python3_EXECUTABLE}")
endif()

if(HDF5_FOUND OR EXISTS "${CFDX_HDF5_DIR}/include/H5public.h")
    target_include_directories(cfdx_core PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(cfdx_core PUBLIC ${HDF5_C_LIBRARIES} ${HDF5_HL_LIBRARIES} z)
endif()

if(CFDX_ENABLE_GPU)
    target_sources(cfdx_core PRIVATE src/cfdx/runtime/gpu/cuda_kernels.cu)
    target_compile_definitions(cfdx_core PUBLIC CFDX_ENABLE_GPU=1)
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
    add_cfdx_test(test_oriented_numerics tests/unit/test_oriented_numerics.cpp)

    # Unit tests - Linear Algebra
    add_cfdx_test(test_sparse_matrix tests/unit/test_sparse_matrix.cpp)
    add_cfdx_test(test_vector tests/unit/test_vector.cpp)
    add_cfdx_test(test_linear_system tests/unit/test_linear_system.cpp)
    add_cfdx_test(test_cg_solver tests/unit/test_cg_solver.cpp)
    add_cfdx_test(test_bicgstab_solver tests/unit/test_bicgstab_solver.cpp)
    add_cfdx_test(test_gmres_solver tests/unit/test_gmres_solver.cpp)
    add_cfdx_test(test_performance_foundations tests/unit/test_performance_foundations.cpp)
    add_cfdx_test(test_performance_runtime tests/unit/test_performance_runtime.cpp)
    add_cfdx_test(test_amg_chebyshev tests/unit/test_amg_chebyshev.cpp)

    # Unit tests - Numerical
    add_cfdx_test(test_laplacian tests/unit/test_laplacian.cpp)
    add_cfdx_test(test_phase123_regression tests/unit/test_phase123_regression.cpp)
    add_cfdx_test(test_scalar_diffusion tests/unit/test_scalar_diffusion.cpp)
    add_cfdx_test(test_advanced_preconditioners tests/unit/test_advanced_preconditioners.cpp)

    # Numerical tests - Geometry
    add_cfdx_test(test_face_geometry tests/numerical/test_face_geometry.cpp)
    add_cfdx_test(test_cell_geometry tests/numerical/test_cell_geometry.cpp)
    add_cfdx_test(test_mesh_quality tests/numerical/test_mesh_quality.cpp)
    add_cfdx_test(test_geometry_invariants tests/unit/test_geometry_invariants.cpp)
    add_cfdx_test(test_mesh_family_volumes tests/numerical/test_mesh_family_volumes.cpp)
    add_cfdx_test(test_source_term tests/numerical/test_source_term.cpp)
    add_cfdx_test(test_temporal_physics tests/unit/test_temporal.cpp)
    add_cfdx_test(test_core_temporal tests/unit/test_core_temporal.cpp)
    add_cfdx_test(test_incompressible tests/unit/test_incompressible.cpp)
    add_cfdx_test(test_m1_m4_physics tests/unit/test_m1_m4_physics.cpp)
    add_cfdx_test(test_m1_m4_validation tests/validation/test_m1_m4_validation.cpp)
    add_cfdx_test(test_solver_control tests/unit/test_solver_control.cpp)
    add_cfdx_test(test_finite_volume_transport tests/unit/test_finite_volume_transport.cpp)
    add_cfdx_test(test_steady_incompressible_solver tests/validation/test_steady_incompressible_solver.cpp)
    add_cfdx_test(test_ghia_cavity tests/validation/test_ghia_cavity.cpp)
    add_cfdx_test(test_m2_m4_solver_validation tests/validation/test_m2_m4_solver_validation.cpp)
    add_cfdx_test(test_analytical_benchmarks tests/validation/test_analytical_benchmarks.cpp)
    add_cfdx_test(test_fluent_vmfl_reference tests/validation/test_fluent_vmfl_reference.cpp)
    add_cfdx_test(test_benchmark_matrix tests/validation/test_benchmark_matrix.cpp)
    add_cfdx_test(test_level_b_reference_benchmarks tests/validation/test_level_b_reference_benchmarks.cpp)
    add_cfdx_test(test_level_c_coupled_verification tests/validation/test_level_c_coupled_verification.cpp)
    add_cfdx_test(test_cht_validation tests/validation/test_cht_validation.cpp)
    add_cfdx_test(test_geometry_cache tests/unit/test_geometry_cache.cpp)
    if(CFDX_ENABLE_MPI AND MPI_FOUND)
        add_executable(test_mpi_partition tests/unit/test_mpi_partition.cpp)
        target_link_libraries(test_mpi_partition PRIVATE cfdx_core MPI::MPI_CXX)
        target_include_directories(test_mpi_partition PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        add_test(NAME test_mpi_partition
                 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_partition> ${MPIEXEC_POSTFLAGS})
    endif()
    add_cfdx_test(test_transport_models tests/unit/test_transport_models.cpp)
    add_cfdx_test(test_m5_m7_foundations tests/unit/test_m5_m7_foundations.cpp)
    if(CFDX_BUILD_RUNTIME_TESTS)
        add_cfdx_test(test_runtime_execution_ooc tests/unit/test_runtime_execution_ooc.cpp)
    endif()
    if(CFDX_ENABLE_GPU)
        add_cfdx_test(test_cuda_runtime tests/unit/test_cuda_runtime.cpp)
    endif()

    # HDF5 tests
    add_cfdx_test(test_hdf5_writer tests/unit/test_hdf5_writer.cpp)
    add_cfdx_test(test_hdf5_roundtrip tests/unit/test_hdf5_roundtrip.cpp)
    add_cfdx_test(test_field_roundtrip tests/unit/test_field_roundtrip.cpp)

    # Memory Planner tests
    add_cfdx_test(test_memory_planner_validation tests/unit/test_memory_planner_validation.cpp)

    # Plan completion tests
    add_cfdx_test(test_plan_complet tests/unit/test_plan_complet.cpp)

    # Module validation tests (M0.8-T06 to T10, M0.9-T05, M0.10-T01)
    add_cfdx_test(test_validation_simple tests/unit/test_validation_simple.cpp)
    add_cfdx_test(test_openfoam_importer tests/unit/test_openfoam_importer.cpp)
    if(Python3_Interpreter_FOUND)
        add_test(NAME test_meshio_bridge
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_SOURCE_DIR}/tests/mesh/test_meshio_bridge.py
                         ${CMAKE_SOURCE_DIR}/scripts/meshio_import.py)
        add_test(NAME test_mesh_geometry_matrix
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_SOURCE_DIR}/tests/mesh/test_mesh_geometry_matrix.py
                         ${CMAKE_SOURCE_DIR}/scripts/meshio_import.py)
    endif()

endif()
