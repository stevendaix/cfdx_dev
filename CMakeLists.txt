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
option(CFDX_ENABLE_PARALLEL_HDF5 "Enable MPI-enabled HDF5 distributed checkpoints" OFF)
option(CFDX_BUILD_RUNTIME_TESTS "Build execution/GPU-OOC foundation tests" ON)
option(CFDX_ENABLE_LONG_VALIDATION "Register full long-running validation campaigns" OFF)

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
    set(CFDX_HAS_HDF5 ON)
    message(STATUS "Using system HDF5 ${HDF5_VERSION}")
    if(CFDX_ENABLE_PARALLEL_HDF5 AND NOT HDF5_IS_PARALLEL)
        message(FATAL_ERROR "CFDX_ENABLE_PARALLEL_HDF5=ON but detected HDF5 is not parallel-enabled")
    endif()
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
    set(CFDX_HAS_HDF5 ON)
    message(STATUS "Using locally-built HDF5 (static) at ${CFDX_HDF5_DIR}")
    set(HDF5_INCLUDE_DIRS "${CFDX_HDF5_DIR}/include")
    set(HDF5_C_LIBRARIES "${CFDX_HDF5_DIR}/lib/libhdf5.a")
    set(HDF5_HL_LIBRARIES "${CFDX_HDF5_DIR}/lib/libhdf5_hl.a")
else()
    set(CFDX_HAS_HDF5 OFF)
    message(WARNING "HDF5 not found. HDF5-dependent features disabled.")
endif()

if(CFDX_ENABLE_PARALLEL_HDF5 AND NOT HDF5_FOUND)
    message(FATAL_ERROR "CFDX_ENABLE_PARALLEL_HDF5=ON requires a parallel HDF5 installation")
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
    src/cfdx/io/vtu/vtu_writer.cpp
    src/cfdx/io/openfoam/openfoam_importer.cpp
    src/cfdx/io/gmsh/gmsh_importer.cpp
    src/cfdx/io/mesh/mesh_importer.cpp
    src/cfdx/io/cfdx_io/gap_analysis.cpp
    src/cfdx/io/cfdx_io/mapping_rules.cpp
)
if(CFDX_HAS_HDF5)
    list(APPEND CFDX_IO_SOURCES
        src/cfdx/io/hdf5/hdf5_writer.cpp
        src/cfdx/io/hdf5/hdf5_reader.cpp
    )
    list(APPEND CFDX_IO_SOURCES
        src/cfdx/io/hdf5/case_hdf5_io.cpp
    )
else()
    list(APPEND CFDX_IO_SOURCES src/cfdx/io/hdf5/hdf5_unavailable.cpp)
endif()

# Physics sources
set(CFDX_PHYSICS_SOURCES
    src/cfdx/physics/equation_of_state.h
    src/cfdx/physics/multi_component.h
    src/cfdx/physics/source_term_linearization.h
    src/cfdx/physics/transport_models.h
    src/cfdx/physics/temporal.h
    src/cfdx/physics/low_storage_time_integration.h
    src/cfdx/physics/adaptive_cfl.h
    src/cfdx/physics/timestep_control.h
    src/cfdx/physics/dual_time_physical_lifecycle.h
    src/cfdx/physics/dual_time_navier_stokes.h
    src/cfdx/physics/advanced_convergence.h
    src/cfdx/physics/continuation.h
    src/cfdx/physics/forces.h
    src/cfdx/physics/local_time_stepping.h
    src/cfdx/physics/low_mach.h
    src/cfdx/physics/compressible_flux.h
    src/cfdx/physics/radiation_models.h
    src/cfdx/physics/radiation_s2s.h
    src/cfdx/physics/boussinesq.h
    src/cfdx/physics/spalart_allmaras.h
    src/cfdx/physics/wall_distance.h
    src/cfdx/physics/wall_distance_turbulence.h
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
    src/cfdx/application/simulation_controller.cpp
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
    src/cfdx/core/boundary/boundary.h
    src/cfdx/core/boundary/boundary_role.h
    src/cfdx/core/boundary/boundary_condition.h
    src/cfdx/core/boundary/boundary_constraint.h
    src/cfdx/core/boundary/boundary_validator.h
    src/cfdx/core/boundary/mathematical_condition.h
    src/cfdx/core/boundary/value_provider.h
    src/cfdx/core/boundary/legacy_patch_field_adapter.h
    src/cfdx/core/boundary/flow_boundary_conditions.h
    src/cfdx/core/numerics/interpolation.h
    src/cfdx/core/numerics/gradient.h
    src/cfdx/core/numerics/gradient_stencil.h
    src/cfdx/core/numerics/divergence.h
    src/cfdx/core/numerics/laplacian.h
    src/cfdx/core/numerics/flux.h
    src/cfdx/core/numerics/integrate.h
    src/cfdx/core/numerics/source_term.h
    src/cfdx/core/numerics/temporal.h
    src/cfdx/core/numerics/convection.h
    src/cfdx/core/numerics/backend_equivalence.h
    src/cfdx/physics/incompressible.h
    src/cfdx/physics/pressure_velocity.h
    src/cfdx/physics/pressure_velocity_algorithms.h
    src/cfdx/physics/turbulence.h
    src/cfdx/physics/thermal.h
    src/cfdx/physics/thermophysical_models.h
    src/cfdx/physics/turbulence_models.h
    src/cfdx/physics/radiation.h
    src/cfdx/physics/m1_m4_models.h
    src/cfdx/physics/solver_control.h
    src/cfdx/physics/finite_volume_transport.h
    src/cfdx/physics/boundary_constraint_fvm.h
    src/cfdx/physics/steady_incompressible_solver.h
    src/cfdx/physics/dual_time_driver.h
    src/cfdx/physics/dual_time_adaptive.h
    src/cfdx/physics/dual_time_field_driver.h
    src/cfdx/physics/dual_time_field_integrator.h
    src/cfdx/physics/dual_time_navier_stokes.h
    src/cfdx/physics/turbulence_transport.h
    src/cfdx/physics/turbulence_solver.h
    src/cfdx/physics/sst_solver.h
    src/cfdx/physics/energy_solver.h
    src/cfdx/physics/cht_solver.h
    src/cfdx/physics/radiation_solver.h
    src/cfdx/physics/radiation_s2s.h
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
    src/cfdx/core/linalg/linear_solver_context.h
    src/cfdx/core/linalg/linear_solver_models.h
    src/cfdx/core/linalg/linear_solver_dispatch.h
    src/cfdx/core/linalg/null_space.h
    src/cfdx/core/linalg/matrix_free_fv_operator.h
    src/cfdx/core/linalg/amg_preconditioner.h
    src/cfdx/core/linalg/hypre_amg.h
    src/cfdx/core/linalg/coupled_amg_schur.h
    src/cfdx/core/linalg/native_fieldsplit.h
    src/cfdx/core/linalg/mixed_precision.h
    src/cfdx/core/linalg/chebyshev_smoother.h
    src/cfdx/core/linalg/fv_operator.h
    src/cfdx/core/linalg/krylov_controls.h
    src/cfdx/core/linalg/communication_avoiding.h
    src/cfdx/core/parallel/distributed_execution.h
    src/cfdx/core/parallel/distributed_poisson.h
    src/cfdx/runtime/gpu/cuda_poisson.h
    src/cfdx/core/memory/reuse_pool.h
    src/cfdx/core/memory/field_lifetime.h
    src/cfdx/core/linalg/mpi_overlap.h
    src/cfdx/core/linalg/block_schur.h
    src/cfdx/core/mesh/solver_mesh.h
    src/cfdx/core/mesh/sfc_ordering.h
    src/cfdx/core/parallel/restart_mapping.h
    src/cfdx/thermodynamics/equation_of_state.h
    src/cfdx/chemistry/source_terms.h
    src/cfdx/transport/viscosity/viscosity.h
    src/cfdx/transport/conductivity/conductivity.h
    src/cfdx/transport/diffusivity/diffusivity.h
    src/cfdx/io/hdf5/hdf5_writer.h
    src/cfdx/io/hdf5/hdf5_reader.h
    src/cfdx/io/hdf5/case_hdf5_io.h
    src/cfdx/io/hdf5/mini_json.h
    src/cfdx/io/vtu/vtu_writer.h
)

add_executable(cfdx_production_solver apps/production_solver/main.cpp)
target_link_libraries(cfdx_production_solver PRIVATE cfdx_core)

add_executable(cfdx_convert
    apps/cfdx_convert/main.cpp
    apps/cfdx_convert/PythonBridge.cpp
)
target_link_libraries(cfdx_convert PRIVATE cfdx_core)
if(CFDX_HAS_HDF5)
    target_include_directories(cfdx_convert PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(cfdx_convert PRIVATE ${HDF5_C_LIBRARIES} ${HDF5_HL_LIBRARIES} z)
    target_compile_definitions(cfdx_convert PRIVATE CFDX_HAS_HDF5=1)
else()
    target_compile_definitions(cfdx_convert PRIVATE CFDX_HAS_HDF5=0)
endif()

add_executable(cfdx_wall_distance_complex_benchmark apps/verification/wall_distance_complex_benchmark.cpp)
target_link_libraries(cfdx_wall_distance_complex_benchmark PRIVATE cfdx_core)

add_executable(cfdx_wall_distance_controlled_benchmark apps/verification/wall_distance_controlled_benchmark.cpp)
target_link_libraries(cfdx_wall_distance_controlled_benchmark PRIVATE cfdx_core)

add_executable(cfdx_wall_distance_poisson_manufactured apps/verification/wall_distance_poisson_manufactured.cpp)
target_link_libraries(cfdx_wall_distance_poisson_manufactured PRIVATE cfdx_core)

add_executable(cfdx_wall_distance_poisson_decomposition apps/verification/wall_distance_poisson_decomposition.cpp)
target_link_libraries(cfdx_wall_distance_poisson_decomposition PRIVATE cfdx_core)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(cfdx_wall_distance_poisson_manufactured PRIVATE
        -Werror=multichar -Werror=unused-parameter)
endif()

add_library(cfdx_core STATIC ${CFDX_CORE_SOURCES})
target_include_directories(cfdx_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_compile_definitions(cfdx_core PRIVATE CFDX_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
if(Python3_Interpreter_FOUND)
    target_compile_definitions(cfdx_core PRIVATE CFDX_PYTHON_EXECUTABLE="${Python3_EXECUTABLE}")
endif()

if(CFDX_HAS_HDF5)
    target_include_directories(cfdx_core PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(cfdx_core PUBLIC ${HDF5_C_LIBRARIES} ${HDF5_HL_LIBRARIES} z)
    target_compile_definitions(cfdx_core PUBLIC CFDX_HAS_HDF5=1)
else()
    target_compile_definitions(cfdx_core PUBLIC CFDX_HAS_HDF5=0)
endif()

if(CFDX_ENABLE_GPU)
    target_sources(cfdx_core PRIVATE
        src/cfdx/runtime/gpu/cuda_kernels.cu
        src/cfdx/runtime/gpu/cuda_poisson.cu
    )
    target_compile_definitions(cfdx_core PUBLIC CFDX_ENABLE_GPU=1)
endif()

if(CFDX_ENABLE_MPI)
    find_package(MPI QUIET)
    if(MPI_FOUND)
        target_link_libraries(cfdx_core PUBLIC MPI::MPI_CXX)
        target_compile_definitions(cfdx_core PUBLIC CFDX_HAS_MPI=1)
    else()
        target_compile_definitions(cfdx_core PUBLIC CFDX_HAS_MPI=0)
    endif()
else()
    target_compile_definitions(cfdx_core PUBLIC CFDX_HAS_MPI=0)
endif()

if(CFDX_ENABLE_PARALLEL_HDF5)
    target_compile_definitions(cfdx_core PUBLIC CFDX_ENABLE_PARALLEL_HDF5=1)
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

    # N8 physical/production qualification reuses the existing validation
    # executables and records a machine-readable campaign report. It is
    # intentionally a separate CTest entry so the individual physical or
    # algebraic tests remain independently visible in CI.
    if(Python3_Interpreter_FOUND)
        add_test(NAME test_n8_physical_qualification
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_SOURCE_DIR}/scripts/n8_physical_qualification.py
                         --build-dir ${CMAKE_BINARY_DIR})
        set_tests_properties(test_n8_physical_qualification PROPERTIES
            LABELS "n8;qualification;validation;sanitizer-slow"
            RUN_SERIAL TRUE)
    endif()

    # Unit tests - Core
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
    add_cfdx_test(test_boundary_contract tests/unit/test_boundary_contract.cpp)
    add_cfdx_test(test_boundary_fvm_constraints tests/unit/test_boundary_fvm_constraints.cpp)
    add_cfdx_test(test_boundary_validator tests/unit/test_boundary_validator.cpp)
    add_cfdx_test(test_interpolation tests/unit/test_interpolation.cpp)
    add_cfdx_test(test_gradient tests/unit/test_gradient.cpp)
    add_cfdx_test(test_pressure_velocity_system tests/unit/test_pressure_velocity_system.cpp)
    add_cfdx_test(test_oriented_numerics tests/unit/test_oriented_numerics.cpp)
    add_cfdx_test(test_divergence tests/unit/test_divergence.cpp)
    add_cfdx_test(test_flux tests/unit/test_flux.cpp)
    add_cfdx_test(test_spatial_scheme_contract tests/unit/test_spatial_scheme_contract.cpp)
    add_cfdx_test(test_conservation_boundedness tests/unit/test_conservation_boundedness.cpp)
    add_cfdx_test(test_transport_conservation tests/unit/test_transport_conservation.cpp)
    add_cfdx_test(test_conservation_assembly tests/unit/test_conservation_assembly.cpp)
    add_cfdx_test(test_integrate tests/unit/test_integrate.cpp)
    add_cfdx_test(test_forces tests/unit/test_forces.cpp)
    add_cfdx_test(test_probe_csv tests/unit/test_probe_csv.cpp)

    # Unit tests - Linear Algebra
    add_cfdx_test(test_sparse_matrix tests/unit/test_sparse_matrix.cpp)
    add_cfdx_test(test_vector tests/unit/test_vector.cpp)
    add_cfdx_test(test_linear_system tests/unit/test_linear_system.cpp)
    add_cfdx_test(test_cg_solver tests/unit/test_cg_solver.cpp)
    add_cfdx_test(test_null_space tests/unit/test_null_space.cpp)
    add_cfdx_test(test_bicgstab_solver tests/unit/test_bicgstab_solver.cpp)
    add_cfdx_test(test_gmres_solver tests/unit/test_gmres_solver.cpp)
    add_cfdx_test(test_fgmres_solver tests/unit/test_fgmres_solver.cpp)
    add_cfdx_test(test_krylov_preconditioning tests/unit/test_krylov_preconditioning.cpp)
    add_cfdx_test(test_stationary_solvers tests/unit/test_stationary_solvers.cpp)
    add_cfdx_test(test_performance_foundations tests/unit/test_performance_foundations.cpp)
    add_cfdx_test(test_performance_runtime tests/unit/test_performance_runtime.cpp)
    add_cfdx_test(test_phase6_memory_benchmark tests/benchmark/benchmark_phase6_memory.cpp)
    add_cfdx_test(test_amg_chebyshev tests/unit/test_amg_chebyshev.cpp)
    add_cfdx_test(test_native_external_methods tests/unit/test_native_external_methods.cpp)
    add_cfdx_test(test_amg_preconditioner_qualification tests/validation/test_amg_preconditioner_qualification.cpp)
    add_cfdx_test(test_matrix_free_fv_operator tests/unit/test_matrix_free_fv_operator.cpp)
    add_cfdx_test(test_n11_scaling_equivalence tests/validation/test_n11_scaling_equivalence.cpp)
    add_cfdx_test(test_n11_production_scaling tests/validation/test_n11_production_scaling.cpp)
    add_cfdx_test(test_backend_equivalence tests/unit/test_backend_equivalence.cpp)
    add_executable(test_n13_cuda_equivalence tests/validation/test_n13_cuda_equivalence.cpp)
    target_link_libraries(test_n13_cuda_equivalence PRIVATE cfdx_core)
    add_cfdx_test(test_n13_restart_equivalence tests/validation/test_n13_restart_equivalence.cpp)

    # Unit tests - Numerical
    add_cfdx_test(test_laplacian tests/unit/test_laplacian.cpp)
    add_cfdx_test(test_diffusion_diagnostics tests/unit/test_diffusion_diagnostics.cpp)
    add_cfdx_test(test_phase123_regression tests/unit/test_phase123_regression.cpp)
add_cfdx_test(test_phase12_radiation_vv tests/validation/test_phase12_radiation_vv.cpp)
set_tests_properties(test_phase12_radiation_vv PROPERTIES LABELS "phase12;validation")
    add_cfdx_test(test_scalar_diffusion tests/unit/test_scalar_diffusion.cpp)
    add_cfdx_test(test_advanced_preconditioners tests/unit/test_advanced_preconditioners.cpp)
    add_cfdx_test(test_schur_preconditioner tests/unit/test_schur_preconditioner.cpp)

    # Numerical tests - Geometry
    add_cfdx_test(test_face_geometry tests/numerical/test_face_geometry.cpp)
    add_cfdx_test(test_cell_geometry tests/numerical/test_cell_geometry.cpp)
    add_cfdx_test(test_mesh_quality tests/numerical/test_mesh_quality.cpp)
    add_cfdx_test(test_geometry_invariants tests/unit/test_geometry_invariants.cpp)
    add_cfdx_test(test_dat_restart_identity tests/unit/test_dat_restart_identity.cpp)
    add_cfdx_test(test_mesh_family_volumes tests/numerical/test_mesh_family_volumes.cpp)
    add_cfdx_test(test_source_term tests/numerical/test_source_term.cpp)
    add_cfdx_test(test_temporal_physics tests/unit/test_temporal.cpp)
    add_cfdx_test(test_dual_time_stepping tests/unit/test_dual_time_stepping.cpp)
    add_cfdx_test(test_dual_time_adaptive tests/unit/test_dual_time_adaptive.cpp)
    add_cfdx_test(test_dual_time_field_driver tests/unit/test_dual_time_field_driver.cpp)
    add_cfdx_test(test_dual_time_field_integrator tests/unit/test_dual_time_field_integrator.cpp)
    add_cfdx_test(test_dual_time_physical_lifecycle tests/unit/test_dual_time_physical_lifecycle.cpp)
    add_cfdx_test(test_dual_time_navier_stokes tests/unit/test_dual_time_navier_stokes.cpp)
    add_cfdx_test(test_case_numerics_report tests/unit/test_case_numerics_report.cpp)
    add_cfdx_test(test_temporal_history tests/unit/test_temporal_history.cpp)
    add_cfdx_test(test_core_temporal tests/unit/test_core_temporal.cpp)
    add_cfdx_test(test_newton_field tests/unit/test_newton_field.cpp)
    add_cfdx_test(test_simulation_controller tests/unit/test_simulation_controller.cpp)
    add_cfdx_test(test_adaptive_cfl tests/unit/test_adaptive_cfl.cpp)
    add_cfdx_test(test_timestep_control tests/unit/test_timestep_control.cpp)
    add_cfdx_test(test_continuation tests/unit/test_continuation.cpp)
    add_cfdx_test(test_local_time_stepping tests/unit/test_local_time_stepping.cpp)
    add_cfdx_test(test_incompressible tests/unit/test_incompressible.cpp)
    add_cfdx_test(test_m1_m4_physics tests/unit/test_m1_m4_physics.cpp)
    add_cfdx_test(test_m1_m4_validation tests/validation/test_m1_m4_validation.cpp)
    add_cfdx_test(test_solver_control tests/unit/test_solver_control.cpp)
    add_cfdx_test(test_convergence_monitor tests/unit/test_convergence_monitor.cpp)
    add_cfdx_test(test_initialization tests/unit/test_initialization.cpp)
    add_executable(test_n7_convergence_qualification tests/validation/test_n7_convergence_qualification.cpp)
    target_link_libraries(test_n7_convergence_qualification PRIVATE cfdx_core)
    target_include_directories(test_n7_convergence_qualification PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    add_test(NAME test_n7_convergence_qualification COMMAND test_n7_convergence_qualification)
    set_tests_properties(test_n7_convergence_qualification PROPERTIES
        LABELS "n7;qualification;validation;sanitizer-slow")
    add_executable(test_n7_adaptive_relaxation_re400
        tests/validation/test_n7_adaptive_relaxation_re400.cpp)
    target_link_libraries(test_n7_adaptive_relaxation_re400 PRIVATE cfdx_core)
    target_include_directories(test_n7_adaptive_relaxation_re400 PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    add_test(NAME test_n7_adaptive_relaxation_re400
             COMMAND test_n7_adaptive_relaxation_re400)
    set_tests_properties(test_n7_adaptive_relaxation_re400 PROPERTIES
        LABELS "n7;qualification;validation;sanitizer-slow"
        TIMEOUT 1800)
    add_cfdx_test(test_n7_causal_overshoot_regression tests/unit/test_n7_causal_overshoot_regression.cpp)
    set_tests_properties(test_n7_causal_overshoot_regression PROPERTIES
        LABELS "n7;regression;unit")
    add_cfdx_test(test_finite_volume_transport tests/unit/test_finite_volume_transport.cpp)
    add_cfdx_test(test_steady_incompressible_solver tests/validation/test_steady_incompressible_solver.cpp)
    add_cfdx_test(test_vtu_writer_validation tests/validation/test_vtu_writer_validation.cpp)
    add_cfdx_test(test_solver_vtu_provenance tests/validation/test_solver_vtu_provenance.cpp)
    add_executable(test_ghia_cavity tests/validation/test_ghia_cavity.cpp)
    target_link_libraries(test_ghia_cavity PRIVATE cfdx_core)
    target_include_directories(test_ghia_cavity PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    add_test(NAME test_ghia_cavity_quick COMMAND test_ghia_cavity --quick)
    set_tests_properties(test_ghia_cavity_quick PROPERTIES
        LABELS "phase13;validation")
    # Production PCD is an explicit N8 qualification test. It must never be
    # pulled into the quick smoke path because its coupled Krylov/preconditioner
    # cost is materially different from the SIMPLE Ghia regression.
    if(CFDX_ENABLE_LONG_VALIDATION)
        # Run a bounded low-resolution preflight first. The dependency prevents
        # the physical PCD qualification from starting when the nested solver
        # already exhibits runaway behaviour on the smaller case.
        add_test(NAME test_ghia_cavity_n8_pcd_preflight COMMAND test_ghia_cavity --pcd-preflight)
        set_tests_properties(test_ghia_cavity_n8_pcd_preflight PROPERTIES
            LABELS "n8;qualification;validation;long;preflight"
            TIMEOUT 120)
        add_test(NAME test_ghia_cavity_n8_pcd COMMAND test_ghia_cavity --pcd)
        set_tests_properties(test_ghia_cavity_n8_pcd PROPERTIES
            LABELS "n8;qualification;validation;long"
            TIMEOUT 720
            DEPENDS test_ghia_cavity_n8_pcd_preflight)
    endif()
    # Merge-gated observed-order tripwire. Three resolutions at 16/32/64 give a
    # Richardson-style order in about 90s, so observed-order evidence now gates
    # every merge instead of running only in the skipped long-validation
    # workflow. The 32/64/128 campaign stays behind CFDX_ENABLE_LONG_VALIDATION.
    add_test(NAME test_ghia_cavity_order COMMAND test_ghia_cavity --gate)
    set_tests_properties(test_ghia_cavity_order PROPERTIES
        LABELS "phase13;validation;campaign"
        TIMEOUT 900)
    if(CFDX_ENABLE_LONG_VALIDATION)
        add_test(NAME test_ghia_cavity COMMAND test_ghia_cavity)
        set_tests_properties(test_ghia_cavity PROPERTIES
            LABELS "phase13;validation;long")
    endif()
    add_executable(test_phase9_acceptance tests/validation/test_phase9_acceptance.cpp)
    target_link_libraries(test_phase9_acceptance PRIVATE cfdx_core)
    target_include_directories(test_phase9_acceptance PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    # Full N8 pressure-velocity matrix: deliberately not the --quick smoke test.
    # This exercises SIMPLE, SIMPLEC, PISO, PIMPLE, fractional-step and COUPLED.
    add_test(NAME test_n8_pressure_velocity_matrix COMMAND test_phase9_acceptance)
    set_tests_properties(test_n8_pressure_velocity_matrix PROPERTIES
        LABELS "n8;qualification;validation;sanitizer-slow"
        RUN_SERIAL TRUE)
    add_test(NAME test_couette_quick COMMAND test_phase9_acceptance --quick)
    # Phase 9 diagnostics must be deterministic under CTest. The validation
    # workflow may set these variables at the shell level, but CTest test
    # properties are the authoritative per-test environment and avoid losing
    # the microscope when the test is launched through CTest.
    set_tests_properties(test_couette_quick PROPERTIES
        LABELS "phase9;validation;coupled-ns"
        ENVIRONMENT "CFDX_DEBUG_CELL=33;CFDX_FREEZE_STATE=1;CFDX_DEBUG_COUPLED=1")
    if(CFDX_ENABLE_LONG_VALIDATION)
        add_test(NAME test_phase9_acceptance COMMAND test_phase9_acceptance)
        set_tests_properties(test_phase9_acceptance PROPERTIES
            LABELS "phase9;validation;coupled-ns;long"
            ENVIRONMENT "CFDX_DEBUG_CELL=33;CFDX_FREEZE_STATE=1;CFDX_DEBUG_COUPLED=1")
    endif()
    add_cfdx_test(test_m2_m4_solver_validation tests/validation/test_m2_m4_solver_validation.cpp)
    add_cfdx_test(test_phase10_turbulence_hardening tests/validation/test_phase10_turbulence_hardening.cpp)
    add_cfdx_test(test_turbulence_qualification_equations tests/unit/test_turbulence_qualification_equations.cpp)
    add_cfdx_test(test_phase10_model_registry tests/validation/test_phase10_model_registry.cpp)
    add_cfdx_test(test_wall_distance_methods tests/unit/test_wall_distance_methods.cpp)
    add_test(NAME benchmark_wall_distance_complex
             COMMAND cfdx_wall_distance_complex_benchmark ${CMAKE_BINARY_DIR}/wall_distance_complex_benchmark.csv)
    set_tests_properties(benchmark_wall_distance_complex PROPERTIES LABELS "phase13;validation;benchmark")
    add_cfdx_test(test_phase10_sst_limiter tests/validation/test_phase10_sst_limiter.cpp)
    add_cfdx_test(test_analytical_benchmarks tests/validation/test_analytical_benchmarks.cpp)
    add_cfdx_test(test_vmfl036_reference_case tests/validation/test_vmfl036_reference_case.cpp)
    add_cfdx_test(test_axisymmetric_vmfl036_kernel tests/validation/test_axisymmetric_vmfl036_kernel.cpp)
    add_executable(test_vmfl036_axisymmetric tests/validation/test_vmfl036_axisymmetric.cpp)
    target_link_libraries(test_vmfl036_axisymmetric PRIVATE cfdx_core)
    target_include_directories(test_vmfl036_axisymmetric PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    target_compile_definitions(test_vmfl036_axisymmetric PRIVATE CFDX_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
    if(Python3_Interpreter_FOUND)
        target_compile_definitions(test_vmfl036_axisymmetric PRIVATE CFDX_PYTHON_EXECUTABLE="${Python3_EXECUTABLE}")
    endif()
    add_test(NAME test_vmfl036_axisymmetric_quick COMMAND test_vmfl036_axisymmetric --quick)
    set_tests_properties(test_vmfl036_axisymmetric_quick PROPERTIES LABELS "phase14;validation;axisymmetric")
    if(CFDX_ENABLE_LONG_VALIDATION)
        add_test(NAME test_vmfl036_axisymmetric COMMAND test_vmfl036_axisymmetric)
        set_tests_properties(test_vmfl036_axisymmetric PROPERTIES LABELS "phase14;validation;axisymmetric;long" SKIP_RETURN_CODE 1)
    endif()
    set_tests_properties(test_axisymmetric_vmfl036_kernel PROPERTIES LABELS "phase14;validation;axisymmetric")
    target_compile_definitions(test_vmfl036_reference_case PRIVATE
        CFDX_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    )
    if(Python3_Interpreter_FOUND)
        target_compile_definitions(test_vmfl036_reference_case PRIVATE
            CFDX_PYTHON_EXECUTABLE="${Python3_EXECUTABLE}"
        )
    endif()
    set_tests_properties(test_vmfl036_reference_case PROPERTIES LABELS "phase14;validation;reference")
    add_executable(test_bfs_qualification tests/validation/test_bfs_driver_smoke.cpp)
    target_link_libraries(test_bfs_qualification PRIVATE cfdx_core)
    add_test(NAME test_bfs_quick COMMAND test_bfs_qualification --quick ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/bfs_re200_n16.h5)
    set_tests_properties(test_bfs_quick PROPERTIES LABELS "qualification;validation")
    add_executable(test_naca0012_qualification tests/validation/test_naca0012_geometry.cpp)
    target_link_libraries(test_naca0012_qualification PRIVATE cfdx_core)
    add_test(NAME test_naca0012_quick COMMAND test_naca0012_qualification --quick ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/naca0012_laminar_n64.h5)
    set_tests_properties(test_naca0012_quick PROPERTIES LABELS "qualification;validation")
    if(Python3_Interpreter_FOUND)
        set_tests_properties(test_bfs_quick test_naca0012_quick PROPERTIES DEPENDS test_generate_2d_validation_meshes)
    endif()
    if(CFDX_ENABLE_LONG_VALIDATION)
        add_test(NAME test_bfs_qualification COMMAND test_bfs_qualification ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/bfs_re200_n16.h5 ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/bfs_re200_n32.h5 ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/bfs_re200_n64.h5)
        set_tests_properties(test_bfs_qualification PROPERTIES LABELS "qualification;validation;long" SKIP_RETURN_CODE 1)
        add_test(NAME test_naca0012_qualification COMMAND test_naca0012_qualification ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/naca0012_laminar_n64.h5 ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/naca0012_laminar_n128.h5 ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/naca0012_laminar_n256.h5)
        set_tests_properties(test_naca0012_qualification PROPERTIES LABELS "qualification;validation;long" SKIP_RETURN_CODE 1)
        if(Python3_Interpreter_FOUND)
            set_tests_properties(test_bfs_qualification test_naca0012_qualification PROPERTIES DEPENDS test_generate_2d_validation_meshes)
        endif()
    endif()
    if(Python3_Interpreter_FOUND)
        add_test(NAME test_generate_2d_validation_meshes COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/generate_2d_validation_meshes.py all --output-dir ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes)
        set_tests_properties(test_generate_2d_validation_meshes PROPERTIES LABELS "qualification;mesh;validation")
    endif()
    add_executable(test_channel_qualification tests/validation/test_channel_qualification.cpp)
    target_link_libraries(test_channel_qualification PRIVATE cfdx_core)
    target_include_directories(test_channel_qualification PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    add_test(NAME test_channel_quick COMMAND test_channel_qualification --quick ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/channel_n16.h5)
    set_tests_properties(test_channel_quick PROPERTIES LABELS "qualification;validation;sanitizer-slow")
    if(Python3_Interpreter_FOUND)
        set_tests_properties(test_channel_quick PROPERTIES DEPENDS test_generate_2d_validation_meshes)
    endif()
    if(CFDX_ENABLE_LONG_VALIDATION)
        add_test(NAME test_channel_qualification COMMAND test_channel_qualification ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/channel_n16.h5 ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/channel_n32.h5 ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/channel_n64.h5)
        set_tests_properties(test_channel_qualification PROPERTIES LABELS "qualification;validation;long" SKIP_RETURN_CODE 1)
        if(Python3_Interpreter_FOUND)
            set_tests_properties(test_channel_qualification PROPERTIES DEPENDS test_generate_2d_validation_meshes)
        endif()
    endif()
    add_cfdx_test(test_numerical_model_verification tests/validation/test_numerical_model_verification.cpp)
    add_executable(test_poiseuille_diagnostics tests/validation/test_poiseuille_diagnostics.cpp)
    target_link_libraries(test_poiseuille_diagnostics PRIVATE cfdx_core)
    target_include_directories(test_poiseuille_diagnostics PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    add_test(NAME test_poiseuille_quick COMMAND test_poiseuille_diagnostics --quick)
    set_tests_properties(test_poiseuille_quick PROPERTIES
        LABELS "phase13;validation")
    if(CFDX_ENABLE_LONG_VALIDATION)
        add_test(NAME test_poiseuille_diagnostics COMMAND test_poiseuille_diagnostics)
        set_tests_properties(test_poiseuille_diagnostics PROPERTIES
            LABELS "phase13;validation;long")
    endif()
    add_cfdx_test(test_mms_scalar_diffusion tests/validation/test_mms_scalar_diffusion.cpp)
    add_executable(test_mesh_refinement_order tests/validation/test_mesh_refinement_order.cpp)
    target_link_libraries(test_mesh_refinement_order PRIVATE cfdx_core)
    target_include_directories(test_mesh_refinement_order PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    if(CFDX_ENABLE_LONG_VALIDATION)
        add_test(NAME test_mesh_refinement_order COMMAND test_mesh_refinement_order)
        set_tests_properties(test_mesh_refinement_order PROPERTIES
            LABELS "phase13;validation;long")
    endif()
    add_cfdx_test(test_nonorthogonal_skew_campaign tests/validation/test_nonorthogonal_skew_campaign.cpp)
    add_cfdx_test(test_nonorthogonal_laplacian_campaign tests/validation/test_nonorthogonal_laplacian_campaign.cpp)
    add_cfdx_test(test_transient_mms tests/validation/test_transient_mms.cpp)
    add_cfdx_test(test_transient_mms_closure tests/validation/test_transient_mms_closure.cpp)
    add_cfdx_test(test_transient_mms_3d_stiff tests/validation/test_transient_mms_3d_stiff.cpp)
    set_tests_properties(
        test_mms_scalar_diffusion
        test_transient_mms
        test_transient_mms_closure
        test_transient_mms_3d_stiff
        PROPERTIES LABELS "n12;mms;validation")
    add_cfdx_test(test_polyhedral_laplacian_campaign tests/validation/test_polyhedral_laplacian_campaign.cpp)
    add_cfdx_test(test_convection_scheme_verification tests/validation/test_convection_scheme_verification.cpp)
    add_cfdx_test(test_convection_3d_verification tests/validation/test_convection_3d_verification.cpp)
    add_cfdx_test(test_convection_blending tests/validation/test_convection_blending.cpp)
    add_cfdx_test(test_n16_weno5z tests/validation/test_n16_weno5z.cpp)
    set_tests_properties(test_n16_weno5z PROPERTIES LABELS "n16;validation;experimental")
    add_cfdx_test(test_temporal_order_matrix tests/validation/test_temporal_order_matrix.cpp)
    add_cfdx_test(test_gradient_verification tests/validation/test_gradient_verification.cpp)
    add_cfdx_test(test_polyhedral_gradient_campaign tests/validation/test_polyhedral_gradient_campaign.cpp)
    add_cfdx_test(test_n10_difficult_mesh_robustness tests/validation/test_n10_difficult_mesh_robustness.cpp)
    add_executable(test_n10_external_fixture_import tests/validation/test_n10_external_fixture_import.cpp)
    target_link_libraries(test_n10_external_fixture_import PRIVATE cfdx_core)
    target_include_directories(test_n10_external_fixture_import PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    if(Python3_Interpreter_FOUND)
        add_test(NAME test_n10_qualification
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/n10_qualification.py
                    --build-dir ${CMAKE_BINARY_DIR})
        set_tests_properties(test_n10_qualification PROPERTIES
            LABELS "n10;qualification;validation"
            RUN_SERIAL TRUE)
        add_test(NAME test_n10_mesh_fixture_manifest
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/validate_mesh_fixture_manifest.py)
        set_tests_properties(test_n10_mesh_fixture_manifest PROPERTIES
            LABELS "n10;mesh;provenance;validation"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endif()
    set_tests_properties(
        test_n10_difficult_mesh_robustness
        test_nonorthogonal_skew_campaign
        test_nonorthogonal_laplacian_campaign
        test_polyhedral_laplacian_campaign
        test_gradient_verification
        test_polyhedral_gradient_campaign
        PROPERTIES LABELS "validation;qualification;n10")
    add_cfdx_test(test_convection_polyhedral_campaign tests/validation/test_convection_polyhedral_campaign.cpp)
    add_cfdx_test(test_convection_solver_qualification tests/validation/test_convection_solver_qualification.cpp)
    set_tests_properties(test_convection_solver_qualification PROPERTIES LABELS "phase13;validation;qualification")
    set_tests_properties(test_convection_polyhedral_campaign PROPERTIES LABELS "phase13;validation;qualification")
    add_cfdx_test(test_fluent_vmfl_reference tests/validation/test_fluent_vmfl_reference.cpp)
    add_cfdx_test(test_benchmark_matrix tests/validation/test_benchmark_matrix.cpp)
    add_cfdx_test(test_level_b_reference_benchmarks tests/validation/test_level_b_reference_benchmarks.cpp)
    add_cfdx_test(test_level_c_coupled_verification tests/validation/test_level_c_coupled_verification.cpp)
    set_tests_properties(test_level_c_coupled_verification PROPERTIES LABELS "validation;coupled-ns")

    # Canonical engineering qualification cases. These labels make the
    # qualification subset independently selectable from ordinary unit tests.
    add_cfdx_test(test_cht_validation tests/validation/test_cht_validation.cpp)
    add_cfdx_test(test_thermal_vv tests/validation/test_thermal_vv.cpp)
    add_cfdx_test(test_radiation_vv tests/validation/test_radiation_vv.cpp)
    add_cfdx_test(test_s2s_radiation_vv tests/validation/test_s2s_radiation_vv.cpp)
    add_cfdx_test(test_thermal_radiation_model_matrix tests/validation/test_thermal_radiation_model_matrix.cpp)
    add_cfdx_test(test_thermal_radiation_regression tests/validation/test_thermal_radiation_regression.cpp)
    add_cfdx_test(test_thermal_challenge_vv tests/validation/test_thermal_challenge_vv.cpp)
    add_cfdx_test(test_thermophysical_models_vv tests/validation/test_thermophysical_models_vv.cpp)
    set_tests_properties(
        test_ghia_cavity_quick
        test_couette_quick
        test_poiseuille_quick
        test_mms_scalar_diffusion
        test_cht_validation
        test_thermal_vv
        test_radiation_vv
        test_s2s_radiation_vv
        PROPERTIES LABELS "qualification;validation")
    set_tests_properties(test_couette_quick test_ghia_cavity_quick
        PROPERTIES LABELS "qualification;validation;sanitizer-slow")

    add_test(NAME test_turbulence_qualification_matrix
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/validate_turbulence_qualification_matrix.py
                ${CMAKE_SOURCE_DIR}/docs/validation/TURBULENCE_QUALIFICATION_MATRIX.json)
    set_tests_properties(test_turbulence_qualification_matrix PROPERTIES
        LABELS "qualification;validation;turbulence")

    # Engineering qualification registry contract. This is intentionally a
    # metadata gate: numerical qualification remains the responsibility of
    # the executable case and its retained calculation evidence.
    if(Python3_Interpreter_FOUND)
        add_test(NAME test_cfdx_qualification_registry
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/validate_qualification_registry.py
                    ${CMAKE_SOURCE_DIR}/docs/validation/CFDX_QUALIFICATION_REGISTRY.json)
        set_tests_properties(test_cfdx_qualification_registry PROPERTIES
            LABELS "qualification;validation")
    endif()

    # Phase 13 verification labels. The campaign driver selects these tests
    # without conflating unit/regression coverage with numerical verification.
    add_test(NAME test_n13_cuda_equivalence COMMAND test_n13_cuda_equivalence)
    set_tests_properties(test_n13_cuda_equivalence PROPERTIES LABELS "phase13;validation;n13;hardware" SKIP_RETURN_CODE 2)
    set_tests_properties(test_backend_equivalence test_n13_restart_equivalence
        PROPERTIES LABELS "phase13;validation;n13")
    set_tests_properties(test_mms_scalar_diffusion
        test_m1_m4_validation test_m2_m4_solver_validation test_cht_validation
        test_thermal_vv test_radiation_vv test_s2s_radiation_vv test_thermal_radiation_model_matrix
        test_thermal_challenge_vv
        test_phase10_turbulence_hardening test_thermophysical_models_vv
        test_analytical_benchmarks test_benchmark_matrix
        test_level_b_reference_benchmarks test_level_c_coupled_verification
        test_steady_incompressible_solver test_matrix_free_fv_operator
        test_temporal_physics test_core_temporal
        PROPERTIES LABELS "phase13;validation")
    add_cfdx_test(test_schur_infrastructure tests/unit/test_schur_infrastructure.cpp)
    add_cfdx_test(test_exact_schur tests/unit/test_exact_schur.cpp)
    add_cfdx_test(test_n9_schur_oracle tests/unit/test_n9_schur_oracle.cpp)
    add_cfdx_test(test_simplerc_schur tests/unit/test_simplerc_schur.cpp)
    add_cfdx_test(test_lsc_bfbt_schur tests/unit/test_lsc_bfbt_schur.cpp)
    add_cfdx_test(test_pcd_schur tests/unit/test_pcd_schur.cpp)
    add_cfdx_test(test_coupled_schur_selection tests/unit/test_coupled_schur_selection.cpp)
    add_cfdx_test(test_lsc_bfbt_schur_null_space tests/unit/test_lsc_bfbt_schur_null_space.cpp)
    add_cfdx_test(test_schur_approximation_comparison tests/unit/test_schur_approximation_comparison.cpp)
    add_cfdx_test(test_schur_quantitative_qualification tests/unit/test_schur_quantitative_qualification.cpp)
    add_cfdx_test(test_mgr_preconditioner tests/unit/test_mgr_preconditioner.cpp)
    add_cfdx_test(test_wall_distance_turbulence tests/unit/test_wall_distance_turbulence.cpp)
    add_cfdx_test(test_geometry_cache tests/unit/test_geometry_cache.cpp)
    add_cfdx_test(test_coupled_block_schur_amg tests/validation/test_coupled_block_schur_amg.cpp)
    add_cfdx_test(test_n8_schur_production_benchmark tests/validation/test_n8_schur_production_benchmark.cpp)
    add_cfdx_test(test_n8_schur_representative_cfd tests/validation/test_n8_schur_representative_cfd.cpp)
    add_cfdx_test(test_n8_schur_complete_krylov tests/validation/test_n8_schur_complete_krylov.cpp)
    add_cfdx_test(test_n8_amg_anisotropy_scaling tests/validation/test_n8_amg_anisotropy_scaling.cpp)
    add_cfdx_test(test_n9_common_schur_qualification tests/validation/test_n9_common_schur_qualification.cpp)
    add_cfdx_test(test_n9_s9_rhie_chow_consistency tests/validation/test_n9_s9_rhie_chow_consistency.cpp)
    set_tests_properties(test_n9_s9_rhie_chow_consistency PROPERTIES
        LABELS "n9;qualification;validation;long;validation-total"
        TIMEOUT 300
        RUN_SERIAL TRUE)
    if(CFDX_ENABLE_MPI AND MPI_FOUND)
        add_executable(test_mpi_partition tests/unit/test_mpi_partition.cpp)
        add_executable(test_n13_mpi_equivalence tests/validation/test_n13_mpi_equivalence.cpp)
        add_executable(test_mpi_deterministic tests/unit/test_mpi_deterministic.cpp)
        add_executable(test_mpi_fused_reduction tests/unit/test_mpi_fused_reduction.cpp)
        add_executable(test_mpi_krylov_reduction tests/unit/test_mpi_krylov_reduction.cpp)
        add_executable(test_phase8_mpi_poisson tests/validation/test_phase8_mpi_poisson.cpp)
        if(CFDX_ENABLE_PARALLEL_HDF5 AND HDF5_FOUND)
            add_executable(cfdx_mpi_poisson_restart apps/mpi_poisson_restart/main.cpp)
            target_link_libraries(cfdx_mpi_poisson_restart PRIVATE cfdx_core MPI::MPI_CXX)
            add_executable(test_phase5_distributed tests/unit/test_phase5_distributed.cpp)
            add_executable(test_phase5_restart_writer tests/unit/test_phase5_restart_writer.cpp)
            add_executable(test_phase5_restart_reader tests/unit/test_phase5_restart_reader.cpp)
            target_link_libraries(test_phase5_distributed PRIVATE cfdx_core MPI::MPI_CXX)
            target_link_libraries(test_phase5_restart_writer PRIVATE cfdx_core MPI::MPI_CXX)
            target_link_libraries(test_phase5_restart_reader PRIVATE cfdx_core MPI::MPI_CXX)
            target_include_directories(test_phase5_distributed PRIVATE ${CMAKE_SOURCE_DIR}/tests)
            target_include_directories(test_phase5_restart_writer PRIVATE ${CMAKE_SOURCE_DIR}/tests)
            target_include_directories(test_phase5_restart_reader PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        endif()
        target_link_libraries(test_mpi_deterministic PRIVATE cfdx_core MPI::MPI_CXX)
        target_link_libraries(test_mpi_partition PRIVATE cfdx_core MPI::MPI_CXX)
        target_link_libraries(test_n13_mpi_equivalence PRIVATE cfdx_core MPI::MPI_CXX)
        target_link_libraries(test_mpi_fused_reduction PRIVATE cfdx_core MPI::MPI_CXX)
        target_link_libraries(test_mpi_krylov_reduction PRIVATE cfdx_core MPI::MPI_CXX)
        target_link_libraries(test_phase8_mpi_poisson PRIVATE cfdx_core MPI::MPI_CXX)
        target_include_directories(test_mpi_partition PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        target_include_directories(test_n13_mpi_equivalence PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        target_include_directories(test_mpi_fused_reduction PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        target_include_directories(test_mpi_krylov_reduction PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        set(CFDX_MPIEXEC_NUMPROCS_FLAG "${MPIEXEC_NUMPROCS_FLAG}")
        if(CFDX_MPIEXEC_NUMPROCS_FLAG STREQUAL "")
            set(CFDX_MPIEXEC_NUMPROCS_FLAG "-n")
        endif()
        add_test(NAME test_mpi_deterministic
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_deterministic> ${MPIEXEC_POSTFLAGS})
        add_test(NAME test_n13_mpi_equivalence
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_n13_mpi_equivalence> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(test_n13_mpi_equivalence PROPERTIES
            LABELS "phase13;validation;n13"
            ENVIRONMENT "LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/tests/sanitizers/lsan_mpi.supp")
        add_test(NAME test_mpi_partition
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_partition> ${MPIEXEC_POSTFLAGS})
        add_test(NAME test_mpi_fused_reduction
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_fused_reduction> ${MPIEXEC_POSTFLAGS})
        add_test(NAME test_mpi_krylov_reduction
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_krylov_reduction> ${MPIEXEC_POSTFLAGS})
        add_test(NAME test_phase8_mpi_poisson
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_phase8_mpi_poisson> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(test_phase8_mpi_poisson PROPERTIES LABELS "phase13;validation")
        if(CFDX_ENABLE_PARALLEL_HDF5 AND HDF5_FOUND)
            add_test(NAME test_phase5_distributed
                     COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_phase5_distributed> ${MPIEXEC_POSTFLAGS})
            add_test(NAME test_production_mpi_poisson_writer
                     COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:cfdx_mpi_poisson_restart>
                             --mode writer
                             --checkpoint ${CMAKE_BINARY_DIR}/production_poisson_checkpoint.h5
                             --reference ${CMAKE_BINARY_DIR}/production_poisson_reference.h5
                             ${MPIEXEC_POSTFLAGS})
            add_test(NAME test_production_mpi_poisson_reader
                     COMMAND ${MPIEXEC_EXECUTABLE} --oversubscribe ${CFDX_MPIEXEC_NUMPROCS_FLAG} 3
                             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:cfdx_mpi_poisson_restart>
                             --mode reader
                             --checkpoint ${CMAKE_BINARY_DIR}/production_poisson_checkpoint.h5
                             --reference ${CMAKE_BINARY_DIR}/production_poisson_reference.h5
                             ${MPIEXEC_POSTFLAGS})
            set_tests_properties(test_production_mpi_poisson_reader PROPERTIES
                DEPENDS test_production_mpi_poisson_writer)
            add_test(NAME test_phase5_restart_writer
                     COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_phase5_restart_writer> ${MPIEXEC_POSTFLAGS})
            add_test(NAME test_phase5_restart_reader
                     COMMAND ${MPIEXEC_EXECUTABLE} --oversubscribe ${CFDX_MPIEXEC_NUMPROCS_FLAG} 3
                             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_phase5_restart_reader> ${MPIEXEC_POSTFLAGS})
            set_tests_properties(test_phase5_restart_reader PROPERTIES DEPENDS test_phase5_restart_writer)
        endif()
        if(CMAKE_BUILD_TYPE STREQUAL "DebugSanitizers")
            # OpenMPI/hwloc allocates process-global state during MPI_Init which
            # is reported as a leak after MPI_Finalize by LeakSanitizer. Keep
            # ASan/UBSan active for the MPI regression, but exclude those
            # external-library teardown allocations from the leak check.
            # N13 has the same external OpenMPI/hwloc teardown behavior as the
            # other MPI sanitizer tests. Its numerical equivalence and reduction
            # assertions remain active; only LeakSanitizer's process-global
            # external-library teardown check is disabled for these MPI tests.
            set_tests_properties(test_n13_mpi_equivalence test_mpi_deterministic
                test_mpi_partition test_mpi_fused_reduction test_mpi_krylov_reduction
                test_phase8_mpi_poisson PROPERTIES
                ENVIRONMENT "LSAN_OPTIONS=detect_leaks=0")
        endif()
    endif()
    add_cfdx_test(test_convergence_history tests/unit/test_convergence_history.cpp)
    add_cfdx_test(test_restart_mapping tests/unit/test_restart_mapping.cpp)
    add_cfdx_test(test_transport_models tests/unit/test_transport_models.cpp)
    add_cfdx_test(test_equation_of_state_consistency tests/unit/test_equation_of_state_consistency.cpp)
    add_cfdx_test(test_sst_blending tests/unit/test_sst_blending.cpp)
    add_cfdx_test(test_thermo_transport_contracts tests/unit/test_thermo_transport_contracts.cpp)
    add_cfdx_test(test_thermophysical_models tests/unit/test_thermophysical_models.cpp)
    add_cfdx_test(test_m5_m7_foundations tests/unit/test_m5_m7_foundations.cpp)
    if(CFDX_BUILD_RUNTIME_TESTS)
        add_cfdx_test(test_runtime_execution_ooc tests/unit/test_runtime_execution_ooc.cpp)
        add_cfdx_test(test_runtime_pipeline tests/unit/test_runtime_pipeline.cpp)
        add_cfdx_test(test_phase6_runtime_acceptance tests/validation/test_phase6_runtime_acceptance.cpp)
        add_cfdx_test(test_ooc_pipeline tests/unit/test_ooc_pipeline.cpp)
        add_cfdx_test(test_phase7_ooc_validation tests/validation/test_phase7_ooc_validation.cpp)
    endif()
    if(CFDX_ENABLE_GPU)
        add_cfdx_test(test_cuda_runtime tests/unit/test_cuda_runtime.cpp)
        add_cfdx_test(test_phase8_cuda_poisson tests/validation/test_phase8_cuda_poisson.cpp)
    endif()

    # HDF5 tests require the native reader/writer implementation.
    if(CFDX_HAS_HDF5)
        add_cfdx_test(test_hdf5_writer tests/unit/test_hdf5_writer.cpp)
        add_cfdx_test(test_hdf5_roundtrip tests/unit/test_hdf5_roundtrip.cpp)
        add_cfdx_test(test_field_roundtrip tests/unit/test_field_roundtrip.cpp)
        add_cfdx_test(test_case_hdf5_io tests/unit/test_case_hdf5_io.cpp)
    else()
        add_cfdx_test(test_hdf5_unavailable tests/unit/test_hdf5_unavailable.cpp)
    endif()

    # Memory Planner tests
    add_cfdx_test(test_memory_planner_validation tests/unit/test_memory_planner_validation.cpp)

    # Plan completion tests
    add_cfdx_test(test_plan_complet tests/unit/test_plan_complet.cpp)

    # Module validation tests (M0.8-T06 to T10, M0.9-T05, M0.10-T01)
    add_cfdx_test(test_validation_simple tests/unit/test_validation_simple.cpp)
    add_cfdx_test(test_openfoam_importer tests/unit/test_openfoam_importer.cpp)
    if(CFDX_HAS_HDF5)
        add_cfdx_test(test_cross_import_equivalence tests/unit/test_cross_import_equivalence.cpp)
    endif()
    if(Python3_Interpreter_FOUND)
        add_test(NAME test_n9_s8_naca0012_contract
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_SOURCE_DIR}/tests/python/test_n9_s8_naca0012_contract.py)
        set_tests_properties(test_n9_s8_naca0012_contract PROPERTIES
            LABELS "n9;qualification;validation;total-validation"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endif()

    if(Python3_Interpreter_FOUND)
        add_test(NAME test_missing_qualification_case_contracts COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/test_missing_qualification_case_contracts.py)
        set_tests_properties(test_missing_qualification_case_contracts PROPERTIES LABELS "qualification;validation")
    endif()

    if(Python3_Interpreter_FOUND)
        if(CFDX_HAS_HDF5)
            add_test(NAME test_production_solver_e2e
                     COMMAND ${Python3_EXECUTABLE}
                             ${CMAKE_SOURCE_DIR}/tests/python/test_production_solver_e2e.py)
            set_tests_properties(test_production_solver_e2e PROPERTIES
                ENVIRONMENT "CFDX_PRODUCTION_SOLVER=$<TARGET_FILE:cfdx_production_solver>;CFDX_PRODUCTION_MESH=${CMAKE_SOURCE_DIR}/tests/data/cavity_with_patches.msh;PYTHONPATH=${CMAKE_SOURCE_DIR}/src/cfdx/python")
            add_test(NAME test_production_solver_numerics_e2e
                     COMMAND ${Python3_EXECUTABLE}
                             ${CMAKE_SOURCE_DIR}/tests/python/test_production_solver_numerics_e2e.py)
            set_tests_properties(test_production_solver_numerics_e2e PROPERTIES
                ENVIRONMENT "CFDX_PRODUCTION_SOLVER=$<TARGET_FILE:cfdx_production_solver>;PYTHONPATH=${CMAKE_SOURCE_DIR}/src/cfdx/python")
        endif()
        add_test(NAME test_meshio_bridge
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_SOURCE_DIR}/tests/mesh/test_meshio_bridge.py
                         ${CMAKE_SOURCE_DIR}/scripts/meshio_import.py)
        add_test(NAME test_mesh_geometry_matrix
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_SOURCE_DIR}/tests/mesh/test_mesh_geometry_matrix.py
                         ${CMAKE_SOURCE_DIR}/scripts/meshio_import.py)
    endif()

    # I/O adapter tests (M0.10-T03)
    # C++ solver adapters are now thin stubs (format parsing moved to Python).
    add_cfdx_test(test_gap_analysis tests/unit/test_gap_analysis.cpp)

endif()
