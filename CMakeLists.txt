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
option(CFDX_ENABLE_PARALLEL_HDF5 "Enable MPI-enabled HDF5 file access" OFF)
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
if(CFDX_ENABLE_PARALLEL_HDF5 AND NOT CFDX_ENABLE_MPI)
    message(FATAL_ERROR "CFDX_ENABLE_PARALLEL_HDF5 requires CFDX_ENABLE_MPI=ON")
endif()
if(HDF5_FOUND)
    if(CFDX_ENABLE_PARALLEL_HDF5 AND NOT HDF5_IS_PARALLEL)
        message(FATAL_ERROR "CFDX_ENABLE_PARALLEL_HDF5=ON but detected HDF5 is not parallel-enabled")
    endif()
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
set(CFDX_MESH_SOURCES
    src/cfdx/core/mesh/point.cpp
    src/cfdx/core/mesh/face.cpp
    src/cfdx/core/mesh/ownership.cpp
    src/cfdx/core/mesh/cell.cpp
    src/cfdx/core/mesh/mesh.cpp
    src/cfdx/core/mesh/boundary.cpp)
set(CFDX_LINALG_SOURCES src/cfdx/core/linalg/hypre_amg.cpp)
set(CFDX_MEMORY_SOURCES
    src/cfdx/core/memory/memory_planner.cpp
    src/cfdx/runtime/gpu/gpu_kernels.cpp)
set(CFDX_GEOMETRY_SOURCES
    src/cfdx/core/geometry/face_geometry.cpp
    src/cfdx/core/geometry/cell_geometry.cpp
    src/cfdx/core/geometry/mesh_quality.cpp
    src/cfdx/core/geometry/mesh_validator.cpp)
set(CFDX_FIELD_SOURCES
    src/cfdx/core/field/field.cpp
    src/cfdx/core/field/storage.cpp)
set(CFDX_IO_SOURCES
    src/cfdx/io/hdf5/hdf5_writer.cpp
    src/cfdx/io/hdf5/hdf5_reader.cpp
    src/cfdx/io/vtu/vtu_writer.cpp
    src/cfdx/io/openfoam/openfoam_importer.cpp
    src/cfdx/io/gmsh/gmsh_importer.cpp
    src/cfdx/io/mesh/mesh_importer.cpp)
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
    src/cfdx/thermodynamics/thermo_cache.h)
set(CFDX_THERMODYNAMICS_SOURCES src/cfdx/thermodynamics/equation_of_state.h)
set(CFDX_CHEMISTRY_SOURCES src/cfdx/chemistry/source_terms.h)
set(CFDX_TRANSPORT_SOURCES
    src/cfdx/transport/viscosity/viscosity.h
    src/cfdx/transport/conductivity/conductivity.h
    src/cfdx/transport/diffusivity/diffusivity.h)

set(CFDX_CORE_SOURCES
    src/cfdx/application/simulation_controller.cpp
    ${CFDX_MESH_SOURCES} ${CFDX_GEOMETRY_SOURCES} ${CFDX_FIELD_SOURCES}
    ${CFDX_IO_SOURCES} ${CFDX_PHYSICS_SOURCES}
    ${CFDX_THERMODYNAMICS_SOURCES} ${CFDX_CHEMISTRY_SOURCES}
    ${CFDX_TRANSPORT_SOURCES} ${CFDX_LINALG_SOURCES} ${CFDX_MEMORY_SOURCES})

add_library(cfdx_core STATIC ${CFDX_CORE_SOURCES})
target_include_directories(cfdx_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
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

if(CFDX_ENABLE_PARALLEL_HDF5)
    target_compile_definitions(cfdx_core PUBLIC CFDX_ENABLE_PARALLEL_HDF5=1)
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(CFDX_BUILD_TESTS)
    enable_testing()

    macro(add_cfdx_test test_name source_file)
        add_executable(${test_name} ${source_file})
        target_link_libraries(${test_name} PRIVATE cfdx_core)
        target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        add_test(NAME ${test_name} COMMAND ${test_name})
    endmacro()

    add_cfdx_test(test_point tests/unit/test_point.cpp)
    add_cfdx_test(test_face tests/unit/test_face.cpp)
    add_cfdx_test(test_ownership tests/unit/test_ownership.cpp)
    add_cfdx_test(test_cell tests/unit/test_cell.cpp)
    add_cfdx_test(test_boundary tests/unit/test_boundary.cpp)
    add_cfdx_test(test_mesh tests/unit/test_mesh.cpp)
    add_cfdx_test(test_field tests/unit/test_field.cpp)
    add_cfdx_test(test_field_storage tests/unit/test_field_storage.cpp)
    add_cfdx_test(test_field_wall_hardening tests/unit/test_field_wall_hardening.cpp)
    add_cfdx_test(test_boundary_field tests/unit/test_boundary_field.cpp)
    add_cfdx_test(test_interpolation tests/unit/test_interpolation.cpp)
    add_cfdx_test(test_gradient tests/unit/test_gradient.cpp)
    add_cfdx_test(test_oriented_numerics tests/unit/test_oriented_numerics.cpp)
    add_cfdx_test(test_divergence tests/unit/test_divergence.cpp)
    add_cfdx_test(test_flux tests/unit/test_flux.cpp)
    add_cfdx_test(test_integrate tests/unit/test_integrate.cpp)
    add_cfdx_test(test_sparse_matrix tests/unit/test_sparse_matrix.cpp)
    add_cfdx_test(test_vector tests/unit/test_vector.cpp)
    add_cfdx_test(test_linear_system tests/unit/test_linear_system.cpp)
    add_cfdx_test(test_cg_solver tests/unit/test_cg_solver.cpp)
    add_cfdx_test(test_bicgstab_solver tests/unit/test_bicgstab_solver.cpp)
    add_cfdx_test(test_gmres_solver tests/unit/test_gmres_solver.cpp)
    add_cfdx_test(test_krylov_preconditioning tests/unit/test_krylov_preconditioning.cpp)
    add_cfdx_test(test_stationary_solvers tests/unit/test_stationary_solvers.cpp)
    add_cfdx_test(test_performance_foundations tests/unit/test_performance_foundations.cpp)
    add_cfdx_test(test_performance_runtime tests/unit/test_performance_runtime.cpp)
    add_cfdx_test(test_amg_chebyshev tests/unit/test_amg_chebyshev.cpp)
    add_cfdx_test(test_matrix_free_fv_operator tests/unit/test_matrix_free_fv_operator.cpp)
    add_cfdx_test(test_laplacian tests/unit/test_laplacian.cpp)
    add_cfdx_test(test_phase123_regression tests/unit/test_phase123_regression.cpp)
    add_cfdx_test(test_scalar_diffusion tests/unit/test_scalar_diffusion.cpp)
    add_cfdx_test(test_advanced_preconditioners tests/unit/test_advanced_preconditioners.cpp)
    add_cfdx_test(test_schur_preconditioner tests/unit/test_schur_preconditioner.cpp)
    add_cfdx_test(test_face_geometry tests/numerical/test_face_geometry.cpp)
    add_cfdx_test(test_cell_geometry tests/numerical/test_cell_geometry.cpp)
    add_cfdx_test(test_mesh_quality tests/numerical/test_mesh_quality.cpp)
    add_cfdx_test(test_geometry_invariants tests/unit/test_geometry_invariants.cpp)
    add_cfdx_test(test_mesh_family_volumes tests/numerical/test_mesh_family_volumes.cpp)
    add_cfdx_test(test_source_term tests/numerical/test_source_term.cpp)
    add_cfdx_test(test_temporal_physics tests/unit/test_temporal.cpp)
    add_cfdx_test(test_core_temporal tests/unit/test_core_temporal.cpp)
    add_cfdx_test(test_simulation_controller tests/unit/test_simulation_controller.cpp)
    add_cfdx_test(test_adaptive_cfl tests/unit/test_adaptive_cfl.cpp)
    add_cfdx_test(test_local_time_stepping tests/unit/test_local_time_stepping.cpp)
    add_cfdx_test(test_incompressible tests/unit/test_incompressible.cpp)
    add_cfdx_test(test_m1_m4_physics tests/unit/test_m1_m4_physics.cpp)
    add_cfdx_test(test_m1_m4_validation tests/validation/test_m1_m4_validation.cpp)
    add_cfdx_test(test_solver_control tests/unit/test_solver_control.cpp)
    add_cfdx_test(test_steady_incompressible_solver tests/validation/test_steady_incompressible_solver.cpp)
    add_executable(test_ghia_cavity tests/validation/test_ghia_cavity.cpp)
    target_link_libraries(test_ghia_cavity PRIVATE cfdx_core)
    add_cfdx_test(test_m2_m4_solver_validation tests/validation/test_m2_m4_solver_validation.cpp)
    add_cfdx_test(test_analytical_benchmarks tests/validation/test_analytical_benchmarks.cpp)
    add_cfdx_test(test_mms_scalar_diffusion tests/validation/test_mms_scalar_diffusion.cpp)
    add_cfdx_test(test_fluent_vmfl_reference tests/validation/test_fluent_vmfl_reference.cpp)
    add_cfdx_test(test_benchmark_matrix tests/validation/test_benchmark_matrix.cpp)
    add_cfdx_test(test_level_b_reference_benchmarks tests/validation/test_level_b_reference_benchmarks.cpp)
    add_cfdx_test(test_level_c_coupled_verification tests/validation/test_level_c_coupled_verification.cpp)
    add_cfdx_test(test_cht_validation tests/validation/test_cht_validation.cpp)
    add_cfdx_test(test_geometry_cache tests/unit/test_geometry_cache.cpp)

    if(CFDX_ENABLE_MPI AND MPI_FOUND)
        add_executable(test_mpi_partition tests/unit/test_mpi_partition.cpp)
        add_executable(test_mpi_deterministic tests/unit/test_mpi_deterministic.cpp)
        add_executable(test_mpi_fused_reduction tests/unit/test_mpi_fused_reduction.cpp)
        add_executable(test_mpi_deterministic_fused tests/unit/test_mpi_deterministic_fused.cpp)
        target_link_libraries(test_mpi_deterministic PRIVATE cfdx_core MPI::MPI_CXX)
        target_link_libraries(test_mpi_partition PRIVATE cfdx_core MPI::MPI_CXX)
        target_link_libraries(test_mpi_fused_reduction PRIVATE cfdx_core MPI::MPI_CXX)
        target_link_libraries(test_mpi_deterministic_fused PRIVATE cfdx_core MPI::MPI_CXX)
        target_include_directories(test_mpi_deterministic PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        target_include_directories(test_mpi_partition PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        target_include_directories(test_mpi_fused_reduction PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        target_include_directories(test_mpi_deterministic_fused PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        set(CFDX_MPIEXEC_NUMPROCS_FLAG "${MPIEXEC_NUMPROCS_FLAG}")
        if(CFDX_MPIEXEC_NUMPROCS_FLAG STREQUAL "")
            set(CFDX_MPIEXEC_NUMPROCS_FLAG "-n")
        endif()
        add_test(NAME test_mpi_deterministic
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_deterministic> ${MPIEXEC_POSTFLAGS})
        add_test(NAME test_mpi_partition
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_partition> ${MPIEXEC_POSTFLAGS})
        add_test(NAME test_mpi_fused_reduction
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_fused_reduction> ${MPIEXEC_POSTFLAGS})
        add_test(NAME test_mpi_deterministic_fused
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_deterministic_fused> ${MPIEXEC_POSTFLAGS})
        if(CMAKE_BUILD_TYPE STREQUAL "DebugSanitizers")
            set_tests_properties(test_mpi_deterministic test_mpi_partition test_mpi_fused_reduction test_mpi_deterministic_fused PROPERTIES
                ENVIRONMENT "LSAN_OPTIONS=detect_leaks=0")
        endif()
    endif()

    add_cfdx_test(test_restart_mapping tests/unit/test_restart_mapping.cpp)
    add_cfdx_test(test_transport_models tests/unit/test_transport_models.cpp)
    add_cfdx_test(test_equation_of_state_consistency tests/unit/test_equation_of_state_consistency.cpp)
    add_cfdx_test(test_sst_blending tests/unit/test_sst_blending.cpp)
    add_cfdx_test(test_thermo_transport_contracts tests/unit/test_thermo_transport_contracts.cpp)
    add_cfdx_test(test_m5_m7_foundations tests/unit/test_m5_m7_foundations.cpp)
    if(CFDX_BUILD_RUNTIME_TESTS)
        add_cfdx_test(test_runtime_execution_ooc tests/unit/test_runtime_execution_ooc.cpp)
        add_cfdx_test(test_memory_planner tests/unit/test_memory_planner.cpp)
        add_cfdx_test(test_ooc_pipeline tests/unit/test_ooc_pipeline.cpp)
    endif()

    # Parallel HDF5 regression uses real MPI. Legacy HDF5 tests also open
    # files through the same parallel HDF5 library, so execute them under
    # a single MPI rank to satisfy the MPI initialization contract without
    # changing their serial test logic.
    if(CFDX_ENABLE_PARALLEL_HDF5 AND CFDX_ENABLE_MPI AND MPI_FOUND)
        add_executable(test_parallel_hdf5 tests/unit/test_parallel_hdf5.cpp)
        target_link_libraries(test_parallel_hdf5 PRIVATE cfdx_core MPI::MPI_CXX)
        add_test(NAME test_parallel_hdf5
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_parallel_hdf5> ${MPIEXEC_POSTFLAGS})
    endif()

    add_executable(test_hdf5_writer tests/unit/test_hdf5_writer.cpp)
    add_executable(test_hdf5_roundtrip tests/unit/test_hdf5_roundtrip.cpp)
    add_executable(test_field_roundtrip tests/unit/test_field_roundtrip.cpp)
    target_link_libraries(test_hdf5_writer PRIVATE cfdx_core)
    target_link_libraries(test_hdf5_roundtrip PRIVATE cfdx_core)
    target_link_libraries(test_field_roundtrip PRIVATE cfdx_core)
    target_include_directories(test_hdf5_writer PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    target_include_directories(test_hdf5_roundtrip PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    target_include_directories(test_field_roundtrip PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    if(CFDX_ENABLE_PARALLEL_HDF5 AND CFDX_ENABLE_MPI AND MPI_FOUND)
        add_test(NAME test_hdf5_writer
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 1
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_hdf5_writer> ${MPIEXEC_POSTFLAGS})
        add_test(NAME test_hdf5_roundtrip
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 1
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_hdf5_roundtrip> ${MPIEXEC_POSTFLAGS})
        add_test(NAME test_field_roundtrip
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 1
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_field_roundtrip> ${MPIEXEC_POSTFLAGS})
    else()
        add_test(NAME test_hdf5_writer COMMAND test_hdf5_writer)
        add_test(NAME test_hdf5_roundtrip COMMAND test_hdf5_roundtrip)
        add_test(NAME test_field_roundtrip COMMAND test_field_roundtrip)
    endif()

    add_cfdx_test(test_memory_planner_validation tests/unit/test_memory_planner_validation.cpp)
    add_cfdx_test(test_plan_complet tests/unit/test_plan_complet.cpp)
    add_cfdx_test(test_validation_simple tests/unit/test_validation_simple.cpp)
    add_cfdx_test(test_openfoam_importer tests/unit/test_openfoam_importer.cpp)
    add_executable(test_cross_import_equivalence tests/unit/test_cross_import_equivalence.cpp)
    target_link_libraries(test_cross_import_equivalence PRIVATE cfdx_core)
    target_include_directories(test_cross_import_equivalence PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    if(CFDX_ENABLE_PARALLEL_HDF5 AND CFDX_ENABLE_MPI AND MPI_FOUND)
        add_test(NAME test_cross_import_equivalence
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 1
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_cross_import_equivalence> ${MPIEXEC_POSTFLAGS})
    else()
        add_test(NAME test_cross_import_equivalence COMMAND test_cross_import_equivalence)
    endif()

    if(Python3_Interpreter_FOUND)
        add_test(NAME test_meshio_bridge
                 COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/mesh/test_meshio_bridge.py
                         ${CMAKE_SOURCE_DIR}/scripts/meshio_import.py)
        add_test(NAME test_mesh_geometry_matrix
                 COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/mesh/test_mesh_geometry_matrix.py
                         ${CMAKE_SOURCE_DIR}/scripts/meshio_import.py)
    endif()
endif()
