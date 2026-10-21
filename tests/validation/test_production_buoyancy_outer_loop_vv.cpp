// Production buoyancy coupling - external flow-energy iteration through
// the Boussinesq momentum source (capability #885 step 4, issue #885).
//
// Verified properties:
//   * Channel identity: a uniform cell_body_force channel assembles
//     BIT-IDENTICALLY to the uniform body_force value (same body_x/y/z
//     arithmetic, same iterations, same authoritative mass flux).
//   * Source builder: T == T_ref gives exactly zero, a hot cell pushes
//     against gravity (upward for g = (0,0,-9.81)), magnitude is
//     beta dT |g| exactly, and the strict gravity parser rejects empty,
//     wrong-arity, non-numeric, nan/inf and zero declarations.
//   * Coupled fixed point on a uniform response: a buoyant case whose
//     walls force a uniform temperature converges at the second pass
//     (k = 1 can never converge), and its final U, p reproduce a single
//     uniform-body-force run of the same discrete problem.
//   * Real coupling: differential heating (305 K / 295 K) on a closed
//     cavity drives a nonzero circulation, every inner gate is green,
//     and the loop reports stationary T/U/p responses within 1e-10.
//     The verification regime is deliberately weak-coupling (small beta,
//     Stokes-dominated carrier) so the relaxed coupled map contracts
//     geometrically well inside the 200-iteration budget: this step
//     verifies the exact coupled machinery, not physical-strength
//     natural convection.
//   * Replay: from the converged state, one hand-rebuilt pass (source
//     rebuilt from the final temperature, carrier solve, energy pass)
//     reproduces the loop's responses.
//   * Nested radiation: a buoyant case declaring S2S runs the step-3
//     energy-radiation iteration inside every buoyancy pass, and both
//     loop verdicts are green (the step-3 gates are unchanged, only the
//     routing was widened).
//   * Fail-closed catalogue: missing/malformed gravity vector, missing /
//     non-positive / non-numeric thermal expansion coefficient, zero or
//     non-finite reference temperature, uniform body_force conflict,
//     missing material, incomplete energy material, running the loop
//     under a non-buoyant model, and the solver channel misuse (wrong
//     size, both channels, non-finite) are all rejected.
#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/physics/production_energy_stage.h"
#include "cfdx/physics/production_buoyancy_stage.h"
#include "cfdx/physics/production_radiation_contract.h"
#include "common/test_harness.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::io;
using namespace cfdx::physics;
using namespace cfdx::testing;

// ---------------------------------------------------------------------------
// Structured box mesh: nx x ny x 1 hexahedra of size dx x dy x 1 along x,
// y, z. Two point columns per (i, j) node (l = 0, 1 for z = 0, z = 1).
// Patches in order: left, right, bottom, top, front (z = 0), back (z = 1).
// Same construction as the stage-2 / stage-3 VVs so all production-path
// tests exercise the identical production geometry path.
// ---------------------------------------------------------------------------
static Mesh box_mesh(std::size_t nx, std::size_t ny, double dx, double dy)
{
    const auto pid = [&](std::size_t i, std::size_t j, std::size_t l) {
        return (i * (ny + 1) + j) * 2 + l;
    };
    const std::size_t n_points = 2 * (nx + 1) * (ny + 1);
    const std::size_t base_fy = (nx + 1) * ny;
    const std::size_t base_fz = base_fy + nx * (ny + 1);
    const std::size_t n_faces = base_fz + 2 * nx * ny;
    const auto fx = [&](std::size_t i, std::size_t j) { return i * ny + j; };
    const auto fy = [&](std::size_t i, std::size_t j) {
        return base_fy + i * (ny + 1) + j;
    };
    const auto fz = [&](std::size_t i, std::size_t j, std::size_t l) {
        return base_fz + 2 * (i * ny + j) + l;
    };
    const auto cell = [&](std::size_t i, std::size_t j) {
        return i * ny + j;
    };

    Mesh m;
    m.points().resize(n_points);
    for (std::size_t i = 0; i <= nx; ++i)
        for (std::size_t j = 0; j <= ny; ++j)
            for (std::size_t l = 0; l < 2; ++l)
                m.points().set(pid(i, j, l),
                               static_cast<double>(i) * dx,
                               static_cast<double>(j) * dy,
                               static_cast<double>(l));

    for (std::size_t i = 0; i <= nx; ++i)
        for (std::size_t j = 0; j < ny; ++j) {
            // Face normals point owner-to-neighbour for internal faces and
            // outward from the owner for boundary faces.
            if (i == 0)
                m.faces().push_face({pid(i, j, 0), pid(i, j, 1),
                                     pid(i, j + 1, 1), pid(i, j + 1, 0)});
            else
                m.faces().push_face({pid(i, j, 0), pid(i, j + 1, 0),
                                     pid(i, j + 1, 1), pid(i, j, 1)});
        }
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j <= ny; ++j) {
            if (j == 0)
                m.faces().push_face({pid(i, j, 0), pid(i + 1, j, 0),
                                     pid(i + 1, j, 1), pid(i, j, 1)});
            else
                m.faces().push_face({pid(i, j, 0), pid(i, j, 1),
                                     pid(i + 1, j, 1), pid(i + 1, j, 0)});
        }
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t l = 0; l < 2; ++l) {
                if (l == 0)
                    m.faces().push_face({pid(i, j, l), pid(i, j + 1, l),
                                         pid(i + 1, j + 1, l), pid(i + 1, j, l)});
                else
                    m.faces().push_face({pid(i, j, l), pid(i + 1, j, l),
                                         pid(i + 1, j + 1, l), pid(i, j + 1, l)});
            }

    m.ownership().resize(n_faces);
    for (std::size_t j = 0; j < ny; ++j) {
        m.ownership().set_owner(fx(0, j), cell(0, j));
        m.ownership().set_neighbour(fx(0, j), FaceOwnership::BOUNDARY);
        m.ownership().set_owner(fx(nx, j), cell(nx - 1, j));
        m.ownership().set_neighbour(fx(nx, j), FaceOwnership::BOUNDARY);
    }
    for (std::size_t i = 1; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j) {
            m.ownership().set_owner(fx(i, j), cell(i - 1, j));
            m.ownership().set_neighbour(fx(i, j),
                                        static_cast<int>(cell(i, j)));
        }
    for (std::size_t i = 0; i < nx; ++i) {
        m.ownership().set_owner(fy(i, 0), cell(i, 0));
        m.ownership().set_neighbour(fy(i, 0), FaceOwnership::BOUNDARY);
        m.ownership().set_owner(fy(i, ny), cell(i, ny - 1));
        m.ownership().set_neighbour(fy(i, ny), FaceOwnership::BOUNDARY);
    }
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 1; j < ny; ++j) {
            m.ownership().set_owner(fy(i, j), cell(i, j - 1));
            m.ownership().set_neighbour(fy(i, j),
                                        static_cast<int>(cell(i, j)));
        }
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t l = 0; l < 2; ++l) {
            m.ownership().set_owner(fz(i, j, l), cell(i, j));
            m.ownership().set_neighbour(fz(i, j, l),
                                        FaceOwnership::BOUNDARY);
        }

    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
            m.cells().push_cell({fx(i, j), fx(i + 1, j),
                                 fy(i, j), fy(i, j + 1),
                                 fz(i, j, 0), fz(i, j, 1)});

    Patch left;
    left.name = "left";
    left.type = PatchType::WALL;
    for (std::size_t j = 0; j < ny; ++j) left.face_ids.push_back(fx(0, j));
    Patch right;
    right.name = "right";
    right.type = PatchType::WALL;
    for (std::size_t j = 0; j < ny; ++j) right.face_ids.push_back(fx(nx, j));
    Patch bottom;
    bottom.name = "bottom";
    bottom.type = PatchType::WALL;
    for (std::size_t i = 0; i < nx; ++i) bottom.face_ids.push_back(fy(i, 0));
    Patch top;
    top.name = "top";
    top.type = PatchType::WALL;
    for (std::size_t i = 0; i < nx; ++i) top.face_ids.push_back(fy(i, ny));
    Patch front;
    front.name = "front";
    front.type = PatchType::WALL;
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
            front.face_ids.push_back(fz(i, j, 0));
    Patch back;
    back.name = "back";
    back.type = PatchType::WALL;
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
            back.face_ids.push_back(fz(i, j, 1));
    m.boundary().add_patch(left);
    m.boundary().add_patch(right);
    m.boundary().add_patch(bottom);
    m.boundary().add_patch(top);
    m.boundary().add_patch(front);
    m.boundary().add_patch(back);
    return m;
}

// No-slip walls everywhere, zero-gradient pressure: a closed cavity.
static void add_closed_cavity_conditions(
    const Mesh& mesh, VelocityBoundaryConditions& ubc,
    ScalarBoundaryConditions& pbc)
{
    for (std::size_t i = 0; i < mesh.boundary().n_patches(); ++i) {
        const std::string name = mesh.boundary().patch(i).name;
        ubc[name] = {VelocityBoundaryCondition::Type::FIXED_VALUE,
                     {0.0, 0.0, 0.0}};
        pbc[name] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    }
}

// Stokes-dominated verification carrier (mu = 10 Pa.s, rho = 1000): the
// same declared fluid the energy stage audits. Weak buoyancy (small beta)
// keeps the relaxed coupled map geometrically contractive well inside the
// 200-iteration budget.
static IncompressibleSolverControls make_flow_controls()
{
    IncompressibleSolverControls controls;
    controls.algorithm = PressureVelocityAlgorithm::SIMPLE;
    controls.convection_scheme = ConvectionScheme::UPWIND;
    controls.convergence.max_iterations = 400;
    controls.convergence.continuity_tolerance = 1e-12;
    controls.linear_max_iterations = 1000;
    controls.linear_tolerance = 1e-11;
    controls.pressure_reference_cell = 0;
    controls.pressure_reference_value = 0.0;
    controls.density = 1000.0;
    controls.kinematic_viscosity = 1.0e-2;
    return controls;
}

// Buoyant case, uniform temperature forced by all-Dirichlet walls
// (310 K). beta = 2e-4 is physical-strength: the response stays uniform
// whatever the source, so the loop must converge at the second pass.
static CaseSetup make_uniform_buoyant_setup()
{
    CaseSetup setup;
    setup.energy_model = "buoyant";
    setup.gravity_vector = "0 0 -9.81";
    setup.ref_temperature = 300.0;

    MaterialSpec fluid;
    fluid.name = "fluid";
    fluid.density = 1000.0;
    fluid.dynamic_viscosity = 10.0;
    fluid.thermal_conductivity = 1.0;
    fluid.specific_heat = 1000.0;
    fluid.extra_properties["thermal_expansion_coefficient"] = "2e-4";
    setup.materials.push_back(fluid);

    setup.has_initial_condition = true;
    setup.initial_condition.temperature = 310.0;

    const char* walls[6] = {"left", "right", "bottom", "top", "front",
                            "back"};
    for (const char* wall : walls) {
        BoundarySpec spec;
        spec.patch_name = wall;
        spec.type = BCType::WALL;
        spec.value_type = BCValueType::WALL_THERMAL;
        spec.wall_temperature = 310.0;
        setup.boundary_conditions.push_back(spec);
    }
    return setup;
}

// Buoyant case, differential heating: left 305 K, right 295 K, other
// patches adiabatic (no thermal spec). Weak coupling (beta = 1e-11,
// Ra ~ 0.1) so the coupled map contracts geometrically inside the budget.
static CaseSetup make_buoyant_setup()
{
    CaseSetup setup;
    setup.energy_model = "buoyant";
    setup.gravity_vector = "0 0 -9.81";
    setup.ref_temperature = 300.0;

    MaterialSpec fluid;
    fluid.name = "fluid";
    fluid.density = 1000.0;
    fluid.dynamic_viscosity = 10.0;
    fluid.thermal_conductivity = 1.0;
    fluid.specific_heat = 1000.0;
    fluid.extra_properties["thermal_expansion_coefficient"] = "1e-11";
    setup.materials.push_back(fluid);

    // Seed away from T_ref so the first carrier solve carries a non-zero
    // source (lecon #663: a zero first source risks an honestly
    // stagnated first solve).
    setup.has_initial_condition = true;
    setup.initial_condition.temperature = 302.0;

    BoundarySpec left;
    left.patch_name = "left";
    left.type = BCType::WALL;
    left.value_type = BCValueType::WALL_THERMAL;
    left.wall_temperature = 305.0;
    setup.boundary_conditions.push_back(left);

    BoundarySpec right;
    right.patch_name = "right";
    right.type = BCType::WALL;
    right.value_type = BCValueType::WALL_THERMAL;
    right.wall_temperature = 295.0;
    setup.boundary_conditions.push_back(right);

    // bottom / top / front / back: deliberately no spec (adiabatic).
    return setup;
}

int main()
{
    // -----------------------------------------------------------------
    // Channel identity: the cellwise body-force channel is bit-identical
    // to the uniform value when fed a uniform field. Same discrete
    // problem, same iterations, same authoritative flux.
    // -----------------------------------------------------------------
    run_case("buoyancy_cellwise_body_force_channel_bit_identical", [] {
        const std::size_t nx = 4;
        const std::size_t ny = 4;
        auto mesh = box_mesh(nx, ny, 0.25, 0.25);
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        add_closed_cavity_conditions(mesh, ubc, pbc);
        const cfdx::core::Vec3 force{0.0, 0.0, 0.01962};

        Field<double, Location::CELL> u_uniform(
            mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p_uniform(
            mesh.n_cells(), "p", "Pa", 1);
        u_uniform.fill(0.0);
        p_uniform.fill(0.0);
        auto uniform_controls = make_flow_controls();
        uniform_controls.body_force = force;
        const auto uniform_result = solve_steady_incompressible(
            mesh, u_uniform, p_uniform, ubc, pbc, uniform_controls);
        EXPECT_TRUE(uniform_result.converged);

        Field<double, Location::CELL> u_cellwise(
            mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p_cellwise(
            mesh.n_cells(), "p", "Pa", 1);
        u_cellwise.fill(0.0);
        p_cellwise.fill(0.0);
        auto cellwise_controls = make_flow_controls();
        cellwise_controls.cell_body_force.assign(mesh.n_cells(), force);
        const auto cellwise_result = solve_steady_incompressible(
            mesh, u_cellwise, p_cellwise, ubc, pbc, cellwise_controls);
        EXPECT_TRUE(cellwise_result.converged);
        EXPECT_TRUE(uniform_result.iterations ==
                    cellwise_result.iterations);

        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            for (std::size_t d = 0; d < 3; ++d)
                EXPECT_NEAR(u_cellwise.component_data(d)[c],
                            u_uniform.component_data(d)[c], 0.0);
            EXPECT_NEAR(p_cellwise(c), p_uniform(c), 0.0);
        }
        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
            EXPECT_NEAR(cellwise_result.authoritative_mass_flux(f),
                        uniform_result.authoritative_mass_flux(f), 0.0);
    });

    // -----------------------------------------------------------------
    // Source builder: identity at T_ref, signs against gravity, exact
    // magnitude, and the strict gravity parser catalogue.
    // -----------------------------------------------------------------
    run_case("buoyancy_source_builder_identity_signs_magnitude", [] {
        Field<double, Location::CELL> temperature(3, "T", "K", 1);
        temperature(0) = 300.0;
        temperature(1) = 310.0;
        temperature(2) = 290.0;
        const cfdx::core::Vec3 gravity{0.0, 0.0, -9.81};
        const auto acceleration = build_cell_buoyancy_acceleration(
            temperature, 2e-4, 300.0, gravity);

        // Identity: T == T_ref gives exactly zero.
        EXPECT_TRUE(acceleration[0].x == 0.0);
        EXPECT_TRUE(acceleration[0].y == 0.0);
        EXPECT_TRUE(acceleration[0].z == 0.0);
        // Hot cell pushes against gravity: upward for g = (0,0,-9.81).
        EXPECT_TRUE(acceleration[1].z > 0.0);
        EXPECT_NEAR(acceleration[1].z, 2e-4 * 10.0 * 9.81, 1e-15);
        EXPECT_TRUE(acceleration[1].x == 0.0);
        EXPECT_TRUE(acceleration[1].y == 0.0);
        // Cold cell sinks.
        EXPECT_TRUE(acceleration[2].z < 0.0);
        EXPECT_NEAR(acceleration[2].z, -2e-4 * 10.0 * 9.81, 1e-15);
        // Arbitrary gravity direction: exact componentwise product.
        const cfdx::core::Vec3 tilted{1.0, -2.0, 0.5};
        const auto tilted_acceleration = build_cell_buoyancy_acceleration(
            temperature, 3e-3, 293.15, tilted);
        for (std::size_t c = 0; c < 3; ++c) {
            const double scaled =
                3e-3 * (293.15 - temperature(c));
            EXPECT_NEAR(tilted_acceleration[c].x, scaled * 1.0, 1e-18);
            EXPECT_NEAR(tilted_acceleration[c].y, scaled * -2.0, 1e-18);
            EXPECT_NEAR(tilted_acceleration[c].z, scaled * 0.5, 1e-18);
        }

        // Non-positive / non-finite temperatures are rejected.
        Field<double, Location::CELL> bad(1, "T", "K", 1);
        bad(0) = 0.0;
        EXPECT_THROW(
            build_cell_buoyancy_acceleration(bad, 2e-4, 300.0, gravity),
            std::invalid_argument);
        bad(0) = -5.0;
        EXPECT_THROW(
            build_cell_buoyancy_acceleration(bad, 2e-4, 300.0, gravity),
            std::invalid_argument);
        bad(0) = std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(
            build_cell_buoyancy_acceleration(bad, 2e-4, 300.0, gravity),
            std::invalid_argument);

        // Strict gravity parsing: whole-string contract.
        EXPECT_THROW(parse_gravity_vector(""), std::invalid_argument);
        EXPECT_THROW(parse_gravity_vector("   "), std::invalid_argument);
        EXPECT_THROW(parse_gravity_vector("0 0"), std::invalid_argument);
        EXPECT_THROW(parse_gravity_vector("0 0 -9.81 1"),
                     std::invalid_argument);
        EXPECT_THROW(parse_gravity_vector("0 0 -9.81x"),
                     std::invalid_argument);
        EXPECT_THROW(parse_gravity_vector("a b c"), std::invalid_argument);
        EXPECT_THROW(parse_gravity_vector("nan 0 0"),
                     std::invalid_argument);
        EXPECT_THROW(parse_gravity_vector("0 0 inf"),
                     std::invalid_argument);
        EXPECT_THROW(parse_gravity_vector("0 0 0"), std::invalid_argument);
        const auto parsed = parse_gravity_vector("  0 0 -9.81  ");
        EXPECT_NEAR(parsed.x, 0.0, 0.0);
        EXPECT_NEAR(parsed.y, 0.0, 0.0);
        EXPECT_NEAR(parsed.z, -9.81, 0.0);
    });

    // -----------------------------------------------------------------
    // Coupled fixed point on a uniform response: all-Dirichlet walls at
    // 310 K force a uniform temperature, so the loop must converge at the
    // second pass (k = 1 has no previous response) and its final U, p
    // reproduce a single uniform-body-force run of the same discrete
    // problem.
    // -----------------------------------------------------------------
    run_case("buoyant_uniform_temperature_converges_at_second_pass", [] {
        const std::size_t nx = 4;
        const std::size_t ny = 4;
        auto mesh = box_mesh(nx, ny, 0.25, 0.25);
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        add_closed_cavity_conditions(mesh, ubc, pbc);
        auto setup = make_uniform_buoyant_setup();
        const auto controls = make_flow_controls();

        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        U.fill(0.0);
        p.fill(0.0);
        const auto result = run_production_buoyant_outer_loop(
            mesh, ubc, pbc, controls, setup, U, p);

        EXPECT_TRUE(result.converged);
        EXPECT_TRUE(result.outer_iterations == 2);
        EXPECT_TRUE(result.energy_ran);
        EXPECT_TRUE(!result.s2s_ran);
        EXPECT_TRUE(result.max_temperature_change_normalized <= 1e-10);
        EXPECT_TRUE(result.max_velocity_change_normalized <= 1e-10);
        EXPECT_TRUE(result.max_pressure_change_normalized <= 1e-10);
        EXPECT_TRUE(result.energy.converged);
        EXPECT_TRUE(result.energy.global_balance_normalized <= 1e-12);
        EXPECT_TRUE(result.energy.local_balance_normalized <= 1e-12);
        EXPECT_TRUE(result.energy.max_cell_residual_normalized <= 1e-12);
        for (const auto& entry : result.history) {
            EXPECT_TRUE(entry.flow_converged);
            EXPECT_TRUE(entry.energy_converged);
        }
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            EXPECT_NEAR(result.temperature(c), 310.0, 1e-9);

        // Reference: a single uniform-body-force run of the same discrete
        // problem (a uniform cellwise field assembles identically, and
        // T - T_ref = 10 K everywhere gives exactly this source).
        Field<double, Location::CELL> u_reference(
            mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p_reference(
            mesh.n_cells(), "p", "Pa", 1);
        u_reference.fill(0.0);
        p_reference.fill(0.0);
        auto reference_controls = make_flow_controls();
        reference_controls.body_force = cfdx::core::Vec3{
            0.0, 0.0, 2e-4 * 10.0 * 9.81};
        const auto reference = solve_steady_incompressible(
            mesh, u_reference, p_reference, ubc, pbc, reference_controls);
        EXPECT_TRUE(reference.converged);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            for (std::size_t d = 0; d < 3; ++d)
                EXPECT_NEAR(U.component_data(d)[c],
                            u_reference.component_data(d)[c], 1e-9);
            EXPECT_NEAR(p(c), p_reference(c), 1e-9);
        }
    });

    // -----------------------------------------------------------------
    // Real coupling: differential heating on a closed cavity drives a
    // circulation; every inner gate is green and the loop reports
    // stationary T/U/p responses. The converged state is replayed
    // independently: one hand-rebuilt pass reproduces the responses.
    // -----------------------------------------------------------------
    run_case("buoyant_differential_heating_coupled_fixed_point", [] {
        const std::size_t nx = 4;
        const std::size_t ny = 4;
        auto mesh = box_mesh(nx, ny, 0.25, 0.25);
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        add_closed_cavity_conditions(mesh, ubc, pbc);
        auto setup = make_buoyant_setup();
        const auto controls = make_flow_controls();

        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        U.fill(0.0);
        p.fill(0.0);
        const auto result = run_production_buoyant_outer_loop(
            mesh, ubc, pbc, controls, setup, U, p);

        EXPECT_TRUE(result.converged);
        EXPECT_TRUE(result.outer_iterations >= 2);
        EXPECT_TRUE(result.outer_iterations <= 200);
        EXPECT_TRUE(result.energy_ran);
        EXPECT_TRUE(result.flow.converged);
        EXPECT_TRUE(result.max_temperature_change_normalized <= 1e-10);
        EXPECT_TRUE(result.max_velocity_change_normalized <= 1e-10);
        EXPECT_TRUE(result.max_pressure_change_normalized <= 1e-10);
        EXPECT_TRUE(result.energy.converged);
        EXPECT_TRUE(result.energy.global_balance_normalized <= 1e-12);
        EXPECT_TRUE(result.energy.local_balance_normalized <= 1e-12);
        EXPECT_TRUE(result.energy.max_cell_residual_normalized <= 1e-12);
        for (const auto& entry : result.history) {
            EXPECT_TRUE(entry.flow_converged);
            EXPECT_TRUE(entry.energy_converged);
            EXPECT_TRUE(std::isfinite(
                entry.max_temperature_change_normalized));
        }

        // The differential heating exists in the final field (discrete
        // maximum principle keeps it within the Dirichlet bounds).
        double t_min = result.temperature(0);
        double t_max = result.temperature(0);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            t_min = std::min(t_min, result.temperature(c));
            t_max = std::max(t_max, result.temperature(c));
            EXPECT_TRUE(result.temperature(c) >= 295.0 - 1e-6);
            EXPECT_TRUE(result.temperature(c) <= 305.0 + 1e-6);
        }
        EXPECT_TRUE(t_max - t_min > 9.0);

        // The coupling exists: the buoyancy source drives a nonzero
        // circulation and a nonzero authoritative mass flux.
        double max_speed = 0.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            for (std::size_t d = 0; d < 3; ++d)
                max_speed = std::max(
                    max_speed,
                    std::abs(U.component_data(d)[c]));
        EXPECT_TRUE(max_speed > 1e-12);
        double max_flux = 0.0;
        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
            max_flux = std::max(
                max_flux,
                std::abs(result.flow.authoritative_mass_flux(f)));
        EXPECT_TRUE(max_flux > 1e-12);

        // -------------------------------------------------------------
        // Independent replay from the converged state: rebuild the
        // source from the loop's final temperature, run ONE carrier
        // solve warm-started from the final fields and ONE energy pass,
        // and require the responses to reproduce the loop's responses
        // (the coupled fixed point). Tolerances bound the stationarity
        // tolerance (1e-10 normalized) amplified through the response
        // map, with orders of margin.
        // -------------------------------------------------------------
        const auto geometry = build_fv_geometry(mesh);
        const auto replay_gravity =
            parse_gravity_vector(setup.gravity_vector);
        const double replay_beta = production_buoyancy_beta(setup);
        auto replay_controls = make_flow_controls();
        replay_controls.cell_body_force = build_cell_buoyancy_acceleration(
            result.temperature, replay_beta, setup.ref_temperature,
            replay_gravity);
        replay_controls.initialization.mode = InitializationMode::Provided;
        Field<double, Location::CELL> u_replay = U;
        Field<double, Location::CELL> p_replay = p;
        const auto replay_flow = solve_steady_incompressible(
            mesh, u_replay, p_replay, ubc, pbc, replay_controls);
        EXPECT_TRUE(replay_flow.converged);
        Field<double, Location::CELL> t_replay = result.temperature;
        const auto replay_energy = run_production_energy_stage(
            mesh, geometry, setup, replay_flow.authoritative_mass_flux,
            t_replay);
        EXPECT_TRUE(replay_energy.converged);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            EXPECT_NEAR(t_replay(c), result.temperature(c), 1e-6);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            for (std::size_t d = 0; d < 3; ++d)
                EXPECT_NEAR(u_replay.component_data(d)[c],
                            U.component_data(d)[c], 1e-10);
            EXPECT_NEAR(p_replay(c), p(c), 1e-8);
        }
    });

    // -----------------------------------------------------------------
    // Nested radiation: a buoyant case declaring S2S runs the step-3
    // energy-radiation iteration inside every buoyancy pass. Both loop
    // verdicts are green; the step-3 gates are unchanged, only the
    // routing was widened.
    // -----------------------------------------------------------------
    run_case("buoyant_s2s_nested_coupled_fixed_point", [] {
        const std::size_t nx = 4;
        const std::size_t ny = 4;
        auto mesh = box_mesh(nx, ny, 0.25, 0.25);
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        add_closed_cavity_conditions(mesh, ubc, pbc);

        auto setup = make_buoyant_setup();
        setup.radiation_model = "S2S";
        // Two-surface parallel-plate enclosure, same as the step-3 VV;
        // weak beta keeps the buoyancy map contractive.
        setup.materials.front()
            .extra_properties["thermal_expansion_coefficient"] = "1e-12";
        setup.initial_condition.temperature = 350.0;

        BoundarySpec bottom;
        bottom.patch_name = "bottom";
        bottom.type = BCType::WALL;
        bottom.raw_params["emissivity"] = "0.8";
        bottom.raw_params["heat_flux"] = "50";
        bottom.raw_params["s2s_view_factor_to_top"] = "1";
        setup.boundary_conditions.push_back(bottom);

        BoundarySpec top;
        top.patch_name = "top";
        top.type = BCType::WALL;
        top.raw_params["emissivity"] = "0.9";
        top.raw_params["s2s_view_factor_to_bottom"] = "1";
        setup.boundary_conditions.push_back(top);

        // The widened routing accepts S2S on a buoyant case.
        validate_production_radiation_selection(setup);

        const auto controls = make_flow_controls();
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        U.fill(0.0);
        p.fill(0.0);
        const auto result = run_production_buoyant_outer_loop(
            mesh, ubc, pbc, controls, setup, U, p);

        EXPECT_TRUE(result.converged);
        EXPECT_TRUE(result.outer_iterations >= 2);
        EXPECT_TRUE(result.outer_iterations <= 200);
        EXPECT_TRUE(result.energy_ran);
        EXPECT_TRUE(result.s2s_ran);
        EXPECT_TRUE(result.flow.converged);
        EXPECT_TRUE(result.max_temperature_change_normalized <= 1e-10);
        EXPECT_TRUE(result.max_velocity_change_normalized <= 1e-10);
        EXPECT_TRUE(result.max_pressure_change_normalized <= 1e-10);
        // The nested step-3 loop reached its own coupled fixed point on
        // the final pass, with its own separate criteria green.
        EXPECT_TRUE(result.s2s_outer_result.converged);
        EXPECT_TRUE(result.s2s_outer_result.fixed_point_flux_passed);
        EXPECT_TRUE(
            result.s2s_outer_result.fixed_point_wall_temperature_passed);
        EXPECT_TRUE(result.energy.converged);
        EXPECT_TRUE(result.energy.global_balance_normalized <= 1e-12);
        EXPECT_TRUE(result.energy.local_balance_normalized <= 1e-12);
        EXPECT_TRUE(
            result.energy.max_cell_residual_normalized <= 1e-12);
        for (const auto& entry : result.history) {
            EXPECT_TRUE(entry.flow_converged);
            EXPECT_TRUE(entry.energy_converged);
            EXPECT_TRUE(entry.s2s_converged);
        }
        // The final temperature is finite, positive and within the
        // declared thermal bounds (discrete maximum principle).
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            EXPECT_TRUE(std::isfinite(result.temperature(c)));
            EXPECT_TRUE(result.temperature(c) > 0.0);
            EXPECT_TRUE(result.temperature(c) >= 295.0 - 1e-6);
            EXPECT_TRUE(result.temperature(c) <= 400.0 + 1e-6);
        }
    });

    // -----------------------------------------------------------------
    // Fail-closed catalogue.
    // -----------------------------------------------------------------
    run_case("buoyancy_fail_closed_catalogue", [] {
        const std::size_t nx = 2;
        const std::size_t ny = 2;
        auto mesh = box_mesh(nx, ny, 0.25, 0.25);
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        add_closed_cavity_conditions(mesh, ubc, pbc);
        const auto controls = make_flow_controls();

        auto fresh_fields = [&](Field<double, Location::CELL>& U,
                                Field<double, Location::CELL>& p) {
            U = Field<double, Location::CELL>(
                mesh.n_cells(), "U", "m/s", 3);
            p = Field<double, Location::CELL>(
                mesh.n_cells(), "p", "Pa", 1);
            U.fill(0.0);
            p.fill(0.0);
        };
        Field<double, Location::CELL> U;
        Field<double, Location::CELL> p;

        // Sanity: the reference declaration validates.
        { auto setup = make_buoyant_setup();
          validate_production_buoyancy_model(setup); }

        // Running the loop under a non-buoyant declared model.
        {
            auto setup = make_buoyant_setup();
            setup.energy_model = "conjugate";
            fresh_fields(U, p);
            EXPECT_THROW(
                run_production_buoyant_outer_loop(
                    mesh, ubc, pbc, controls, setup, U, p),
                std::invalid_argument);
        }

        // Missing / blank gravity declaration.
        {
            auto setup = make_buoyant_setup();
            setup.gravity_vector = "";
            EXPECT_THROW(
                validate_production_buoyancy_model(setup),
                std::invalid_argument);
            fresh_fields(U, p);
            EXPECT_THROW(
                run_production_buoyant_outer_loop(
                    mesh, ubc, pbc, controls, setup, U, p),
                std::invalid_argument);
        }
        {
            auto setup = make_buoyant_setup();
            setup.gravity_vector = "   ";
            EXPECT_THROW(
                validate_production_buoyancy_model(setup),
                std::invalid_argument);
        }

        // Malformed gravity: wrong arity, trailing junk, non-numeric,
        // nan, zero vector.
        for (const char* bad_gravity :
             {"0 0", "0 0 -9.81 1", "0 0 -9.81x", "a b c", "nan 0 0",
              "0 0 0"}) {
            auto setup = make_buoyant_setup();
            setup.gravity_vector = bad_gravity;
            EXPECT_THROW(
                validate_production_buoyancy_model(setup),
                std::invalid_argument);
        }

        // Missing / invalid thermal expansion coefficient.
        {
            auto setup = make_buoyant_setup();
            setup.materials.front().extra_properties.erase(
                "thermal_expansion_coefficient");
            EXPECT_THROW(
                validate_production_buoyancy_model(setup),
                std::invalid_argument);
        }
        for (const char* bad_beta :
             {"", "0", "-1e-4", "abc", "500abc", "inf"}) {
            auto setup = make_buoyant_setup();
            setup.materials.front()
                .extra_properties["thermal_expansion_coefficient"] =
                bad_beta;
            EXPECT_THROW(
                validate_production_buoyancy_model(setup),
                std::invalid_argument);
        }

        // No material at all.
        {
            auto setup = make_buoyant_setup();
            setup.materials.clear();
            EXPECT_THROW(
                validate_production_buoyancy_model(setup),
                std::invalid_argument);
        }

        // Incomplete energy material is pre-flighted before any solve.
        {
            auto setup = make_buoyant_setup();
            setup.materials.front().specific_heat = 0.0;
            EXPECT_THROW(
                validate_production_buoyancy_model(setup),
                std::invalid_argument);
        }

        // Reference temperature must be finite and positive.
        {
            auto setup = make_buoyant_setup();
            setup.ref_temperature = 0.0;
            EXPECT_THROW(
                validate_production_buoyancy_model(setup),
                std::invalid_argument);
        }
        {
            auto setup = make_buoyant_setup();
            setup.ref_temperature =
                std::numeric_limits<double>::quiet_NaN();
            EXPECT_THROW(
                validate_production_buoyancy_model(setup),
                std::invalid_argument);
        }

        // A declared uniform body force conflicts with the cellwise
        // buoyancy source.
        {
            auto setup = make_buoyant_setup();
            setup.body_force = {0.0, 0.0, -9.81};
            EXPECT_THROW(
                validate_production_buoyancy_model(setup),
                std::invalid_argument);
        }

        // Radiation routing: S2S stays rejected on a non-coupled energy
        // model and is accepted on a buoyant case (widened routing).
        {
            auto setup = make_buoyant_setup();
            setup.radiation_model = "S2S";
            validate_production_radiation_selection(setup);
            auto rejected = make_buoyant_setup();
            rejected.radiation_model = "S2S";
            rejected.energy_model = "isothermal";
            EXPECT_THROW(
                validate_production_radiation_selection(rejected),
                std::invalid_argument);
        }

        // Solver channel: wrong size, both channels, non-finite values.
        {
            auto bad_controls = make_flow_controls();
            bad_controls.cell_body_force.assign(
                mesh.n_cells() - 1, cfdx::core::Vec3{0.0, 0.0, 1e-3});
            fresh_fields(U, p);
            EXPECT_THROW(
                solve_steady_incompressible(
                    mesh, U, p, ubc, pbc, bad_controls),
                std::invalid_argument);
        }
        {
            auto bad_controls = make_flow_controls();
            bad_controls.cell_body_force.assign(
                mesh.n_cells(), cfdx::core::Vec3{0.0, 0.0, 1e-3});
            bad_controls.body_force = cfdx::core::Vec3{0.0, 0.0, -9.81};
            fresh_fields(U, p);
            EXPECT_THROW(
                solve_steady_incompressible(
                    mesh, U, p, ubc, pbc, bad_controls),
                std::invalid_argument);
        }
        {
            auto bad_controls = make_flow_controls();
            bad_controls.cell_body_force.assign(
                mesh.n_cells(),
                cfdx::core::Vec3{0.0, 0.0,
                                 std::numeric_limits<double>::quiet_NaN()});
            fresh_fields(U, p);
            EXPECT_THROW(
                solve_steady_incompressible(
                    mesh, U, p, ubc, pbc, bad_controls),
                std::invalid_argument);
        }
    });

    return run_all();
}
