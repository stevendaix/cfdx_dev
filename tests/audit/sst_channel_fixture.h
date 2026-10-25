#pragma once

#include "chain_fixture.h"
#include "cfdx/physics/steady_incompressible_solver.h"

using namespace cfdx::physics;

// Three-cell production lifecycle fixture with an accepted SST transport step.
// It is not a physically validated turbulent channel benchmark.
struct SSTChannelCase {
    Mesh mesh = make_1d_chain(3);
    Field<double,Location::CELL> U{3,"U","m/s",3};
    Field<double,Location::CELL> p{3,"p","Pa",1};
    SSTProductionState state{3};
    SSTProductionRuntime runtime;
    IncompressibleSolverControls controls;
    VelocityBoundaryConditions ubc;
    ScalarBoundaryConditions pbc;

    SSTChannelCase() {
        mesh.boundary().patch(0).type = PatchType::INLET;
        mesh.boundary().patch(1).type = PatchType::OUTLET;
        U.fill(-7.0); p.fill(4.0);
        state.k.fill(0.01); state.omega.fill(2.0); state.nu_t.fill(9.0);
        state.wall_distance.fill(17.0); state.strain_rate.fill(11.0);
        controls.density = 1.2;
        controls.kinematic_viscosity = 0.01;
        controls.convergence.max_iterations = 1;
        controls.linear_tolerance = 1e-12;
        controls.convergence.continuity_tolerance = 1e-12;
        controls.initialization.mode = InitializationMode::Uniform;
        controls.initialization.uniform_velocity = {1,0,0};
        controls.initialization.uniform_pressure = 0.0;
        ubc["inlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{1,0,0}};
        ubc["outlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
        // Moving walls make the carrier a uniform exact solution; this is a
        // production lifecycle test, not a channel-flow validation benchmark.
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{1,0,0}};
        pbc["inlet"] = {ScalarBoundaryType::ZERO_GRADIENT,0,0};
        pbc["outlet"] = {ScalarBoundaryType::FIXED_VALUE,0,0};
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0,0};
        runtime.state = &state;
        runtime.transport.model = TurbulenceModel::SST;
        runtime.transport.density = controls.density;
        runtime.transport.molecular_viscosity = controls.kinematic_viscosity;
        runtime.wall_patch_names = {"wall"};
        runtime.wall_omega_sst_patches = {"wall"};
        runtime.k_boundary_conditions = {
            {"inlet",{ScalarBoundaryType::FIXED_VALUE,0.01,0}},
            {"outlet",{ScalarBoundaryType::ZERO_GRADIENT,0,0}},
            {"wall",{ScalarBoundaryType::FIXED_VALUE,0,0}}};
        runtime.omega_boundary_conditions = {
            {"inlet",{ScalarBoundaryType::FIXED_VALUE,2,0}},
            {"outlet",{ScalarBoundaryType::ZERO_GRADIENT,0,0}}};
        runtime.omega_face_values = make_sst_wall_omega_face_values(
            mesh, build_fv_geometry(mesh), runtime.wall_omega_sst_patches, runtime.transport);
        controls.sst_runtime = &runtime;
    }
    SSTChannelCase(const SSTChannelCase&) = delete;
};
