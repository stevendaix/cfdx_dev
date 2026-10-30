// Production buoyancy stage - coupled inverse energy->flow iteration
// (capability #885 step 4, issue #885): the flow-energy Boussinesq
// outer loop.
//
// Verified properties (fiche methode docs/audits/10-thermique-cht.md
// APPROVED Boussinesq contract; fiche 04 for the pressure-velocity
// carrier; fiche 14 for the reused energy stage; fiche 11 for the
// nested S2S fixed point):
//   * UNIFORM-SOURCE EQUIVALENCE: a buoyant run whose temperature is
//     uniform in space (Dirichlet walls at the initial temperature)
//     builds a constant cellwise Boussinesq source. The loop's final
//     carrier state is BIT-IDENTICAL to the equivalent uniform
//     body_force run repeated with the same two-solve sequence
//     (cold start + warm restart), proving the cellwise channel
//     reaches the momentum assembly and the conservation audit with
//     the exact uniform semantics.
//   * REFERENCE-TEMPERATURE IDENTITY: at T = T_ref everywhere the
//     Boussinesq source is identically zero, and the buoyant loop is
//     bit-identical to the same carrier solved WITHOUT any body force
//     (a wrong zero point of the source would show up immediately).
//   * NATURAL CONVECTION FIXED POINT: a differentially heated cavity
//     (hot left / cold right walls, gravity in -y) converges to the
//     coupled fixed point with the three SEPARATE stationarity gates
//     (unrelaxed T/U/p, gross-scale normalized) all <= 1e-10, a
//     nonzero circulation with the physical sign (hot wall rises),
//     a temperature field different from the frozen-flow conduction
//     solution, and a replayed fixed point: rebuilding the source
//     from the final unrelaxed response, re-solving the carrier warm
//     and replaying one energy pass reproduces the final field within
//     the documented relaxation-lag bound.
//   * NESTED S2S: with declared S2S radiation the loop's inner
//     response operator is the stage-3 energy-radiation fixed point;
//     its gates stay green on every outer iteration.
//   * Fail-closed: the buoyant coupling rejects a missing/short/long/
//     non-numeric/zero gravity_vector, a missing/zero/negative/
//     trailing-garbage/non-finite thermal expansion coefficient, a
//     declared uniform case body_force, wrong material counts and a
//     non-positive ref_temperature; the loop itself refuses a
//     non-buoyant case; the solver's cellwise body-force channel
//     rejects a size mismatch, a uniform+cellwise combination and
//     non-finite components. Positives: the strict gravity parse, the
//     one-shot energy stage and the S2S stages all accept a valid
//     buoyant case.
#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/physics/production_energy_stage.h"
#include "cfdx/physics/production_radiation_contract.h"
#include "cfdx/physics/production_radiation_s2s_stage.h"
#include "cfdx/physics/production_buoyancy_stage.h"
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
// Same construction as the stage-2 / stage-3 VVs so all tests exercise
// the identical production geometry path. Cell (i, j) has index
// i * ny + j.
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

// ---------------------------------------------------------------------------
// Helpers. All flow solves use the production solver controls of the
// production driver (SIMPLE + upwind, strict continuity 1e-12, linear
// 1e-11 / 1000, pressure reference cell 0): the VV exercises the same
// numerical path, never a loosened test-only variant.
// ---------------------------------------------------------------------------

// Single complete material carrying the Boussinesq beta through the
// established extra_properties convention (strict case-double parse).
static MaterialSpec buoyant_material(const std::string& beta_text)
{
    MaterialSpec fluid;
    fluid.name = "fluid";
    fluid.density = 1000.0;
    fluid.dynamic_viscosity = 10.0;      // nu = 1e-2 m^2/s: weakly coupled
    fluid.thermal_conductivity = 100.0;  // conduction-dominated response
    fluid.specific_heat = 1000.0;
    fluid.extra_properties["thermal_expansion_coefficient"] = beta_text;
    return fluid;
}

static BoundarySpec wall_thermal(const std::string& patch, double temperature)
{
    BoundarySpec bc;
    bc.patch_name = patch;
    bc.type = BCType::WALL;
    bc.value_type = BCValueType::WALL_THERMAL;
    bc.wall_temperature = temperature;
    return bc;
}

static VelocityBoundaryConditions no_slip_velocity_bcs(const Mesh& mesh)
{
    VelocityBoundaryConditions ubc;
    for (std::size_t i = 0; i < mesh.boundary().n_patches(); ++i)
        ubc[mesh.boundary().patch(i).name] = {
            VelocityBoundaryCondition::Type::FIXED_VALUE,
            {0.0, 0.0, 0.0}};
    return ubc;
}

static ScalarBoundaryConditions zero_gradient_pressure_bcs(const Mesh& mesh)
{
    ScalarBoundaryConditions pbc;
    for (std::size_t i = 0; i < mesh.boundary().n_patches(); ++i)
        pbc[mesh.boundary().patch(i).name] = {
            ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    return pbc;
}

static IncompressibleSolverControls vv_controls(double density,
                                                double kinematic_viscosity)
{
    IncompressibleSolverControls controls;
    controls.algorithm = PressureVelocityAlgorithm::SIMPLE;
    controls.convection_scheme = ConvectionScheme::UPWIND;
    controls.convergence.max_iterations = 2000;
    controls.convergence.continuity_tolerance = 1e-12;
    controls.linear_max_iterations = 1000;
    controls.linear_tolerance = 1e-11;
    controls.pressure_reference_cell = 0;
    controls.pressure_reference_value = 0.0;
    controls.density = density;
    controls.kinematic_viscosity = kinematic_viscosity;
    return controls;
}

// Bit-exact field comparison (a plain == : any bit difference, including
// a NaN on either side, fails).
static void expect_identical_fields(
    const Field<double, Location::CELL>& a,
    const Field<double, Location::CELL>& b,
    const std::string& what)
{
    if (a.size() != b.size())
        throw std::runtime_error("size mismatch comparing " + what);
    for (std::size_t c = 0; c < a.size(); ++c) {
        if (!(a(c) == b(c)))
            throw std::runtime_error(
                "fields differ (" + what + ") at cell " +
                std::to_string(c) + ": " + std::to_string(a(c)) +
                " vs " + std::to_string(b(c)));
    }
}

static void expect_identical_velocity(
    const Field<double, Location::CELL>& a,
    const Field<double, Location::CELL>& b,
    const std::string& what)
{
    for (std::size_t d = 0; d < 3; ++d) {
        const auto& ca = a.component_data(d);
        const auto& cb = b.component_data(d);
        for (std::size_t c = 0; c < a.size(); ++c) {
            if (!(ca[c] == cb[c]))
                throw std::runtime_error(
                    "velocity fields differ (" + what + ") component " +
                    std::to_string(d) + " cell " + std::to_string(c));
        }
    }
}

static std::size_t surface_index(const ProductionS2SRadiationResult& result,
                                 const std::string& name)
{
    for (std::size_t i = 0; i < result.surface_names.size(); ++i)
        if (result.surface_names[i] == name)
            return i;
    throw std::runtime_error("surface not found: " + name);
}

int main()
{
    // -----------------------------------------------------------------
    // A. Uniform-source equivalence: walls and initial field at 320 K,
    //    T_ref 300 K, beta 0.001, gravity "0 0 -9.81". The Boussinesq
    //    source is the constant -beta (T - T_ref) g_vec; the loop must
    //    reproduce the equivalent uniform body_force run BIT-EXACTLY
    //    through the same cold-start + warm-restart sequence.
    // -----------------------------------------------------------------
    run_case("production_buoyant_outer_loop_uniform_source_equivalence", [] {
        const std::size_t nx = 4;
        const std::size_t ny = 4;
        auto mesh = box_mesh(nx, ny, 0.25, 0.25);
        const auto ubc = no_slip_velocity_bcs(mesh);
        const auto pbc = zero_gradient_pressure_bcs(mesh);

        CaseSetup setup;
        setup.energy_model = "buoyant";
        setup.ref_temperature = 300.0;
        setup.gravity_vector = "0 0 -9.81";
        setup.materials.push_back(buoyant_material("0.001"));
        for (std::size_t i = 0; i < mesh.boundary().n_patches(); ++i)
            setup.boundary_conditions.push_back(
                wall_thermal(mesh.boundary().patch(i).name, 320.0));
        setup.has_initial_condition = true;
        setup.initial_condition.temperature = 320.0;

        EXPECT_TRUE(production_buoyancy_requested(setup));
        EXPECT_TRUE(production_energy_requested(setup));

        Field<double, Location::CELL> U(
            mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(
            mesh.n_cells(), "p", "Pa", 1);
        auto controls = vv_controls(1000.0, 1e-6);
        const auto loop = run_production_buoyant_outer_loop(
            mesh, ubc, pbc, controls, setup, U, p);

        // Coupled fixed point reached on the third outer iteration:
        // the first can never prove it, the second is a full warm
        // carrier solve whose response still moves, and the third is
        // accepted at entry (diag run 39633631822: the warm carrier
        // re-solve ends STAGNATED with every residual at the machine
        // floor, the entry state is restored and the previously
        // accepted carrier reused), with all three separate
        // stationarity gates green and every inner stage green.
        EXPECT_TRUE(loop.converged);
        EXPECT_TRUE(loop.fixed_point_temperature_passed);
        EXPECT_TRUE(loop.fixed_point_velocity_passed);
        EXPECT_TRUE(loop.fixed_point_pressure_passed);
        EXPECT_TRUE(loop.max_temperature_change_normalized <= 1e-10);
        EXPECT_TRUE(loop.max_velocity_change_normalized <= 1e-10);
        EXPECT_TRUE(loop.max_pressure_change_normalized <= 1e-10);
        EXPECT_TRUE(loop.outer_iterations == 3);
        EXPECT_TRUE(loop.energy_ran);
        EXPECT_FALSE(loop.s2s_ran);
        EXPECT_TRUE(loop.energy.converged);
        EXPECT_TRUE(loop.flow.converged);
        EXPECT_TRUE(loop.history.size() == 3);
        for (const auto& entry : loop.history) {
            EXPECT_TRUE(entry.flow_converged);
            EXPECT_TRUE(entry.energy_converged);
            EXPECT_FALSE(entry.s2s_ran);
        }
        // Convergence-at-entry bookkeeping (diag run 39633631822):
        // only the final (third) iteration is accepted at entry - the
        // warm carrier re-solve ended STAGNATED at machine-precision
        // residuals, the entry state was restored and the previously
        // accepted carrier reused; the first two iterations are full
        // carrier solves.
        EXPECT_TRUE(loop.flow_converged_at_entry ==
                    loop.history[2].flow_converged_at_entry);
        EXPECT_FALSE(loop.history[1].flow_converged_at_entry);
        EXPECT_FALSE(loop.history[0].flow_converged_at_entry);
        EXPECT_TRUE(loop.flow_converged_at_entry);
        // The first iteration has no previous response: its
        // stationarity diagnostics must remain infinite (a single
        // Picard pass is never a fixed point).
        EXPECT_FALSE(std::isfinite(
            loop.history[0].max_temperature_change_normalized));

        // Energy response: uniform 320 K field.
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            EXPECT_TRUE(std::isfinite(loop.temperature(c)));
            EXPECT_NEAR(loop.temperature(c), 320.0, 1e-8);
        }

        // Equivalent uniform body force computed with the SAME
        // expression as build_cell_buoyancy_acceleration (including
        // the negative-zero components), then the SAME two-solve
        // sequence as the loop (cold start, warm restart).
        const Vec3 gravity = parse_gravity_vector(setup.gravity_vector);
        const double perturbation = 0.001 * (320.0 - 300.0);
        const Vec3 equivalent{-perturbation * gravity.x,
                               -perturbation * gravity.y,
                               -perturbation * gravity.z};

        Field<double, Location::CELL> U_ref(
            mesh.n_cells(), "U_ref", "m/s", 3);
        Field<double, Location::CELL> p_ref(
            mesh.n_cells(), "p_ref", "Pa", 1);
        auto ref_controls = vv_controls(1000.0, 1e-6);
        ref_controls.body_force = equivalent;
        const auto ref1 = solve_steady_incompressible(
            mesh, U_ref, p_ref, ubc, pbc, ref_controls);
        EXPECT_TRUE(ref1.converged);
        ref_controls.initialization.mode = InitializationMode::Provided;
        const auto ref2 = solve_steady_incompressible(
            mesh, U_ref, p_ref, ubc, pbc, ref_controls);
        EXPECT_TRUE(ref2.converged);

        // BIT-IDENTICAL carrier state: the cellwise channel must
        // reproduce the uniform semantics exactly.
        expect_identical_velocity(U, U_ref, "uniform equivalence U");
        expect_identical_fields(p, p_ref, "uniform equivalence p");
        EXPECT_TRUE(loop.flow.converged == ref2.converged);
    });

    // -----------------------------------------------------------------
    // B. Reference-temperature identity: at T = T_ref everywhere the
    //    Boussinesq source is identically zero and the buoyant loop is
    //    bit-identical to the same carrier solved with NO body force.
    // -----------------------------------------------------------------
    run_case("production_buoyant_outer_loop_reference_temperature_identity", [] {
        const std::size_t nx = 4;
        const std::size_t ny = 4;
        auto mesh = box_mesh(nx, ny, 0.25, 0.25);
        const auto ubc = no_slip_velocity_bcs(mesh);
        const auto pbc = zero_gradient_pressure_bcs(mesh);

        CaseSetup setup;
        setup.energy_model = "buoyant";
        setup.ref_temperature = 300.0;
        setup.gravity_vector = "0 0 -9.81";
        setup.materials.push_back(buoyant_material("0.001"));
        for (std::size_t i = 0; i < mesh.boundary().n_patches(); ++i)
            setup.boundary_conditions.push_back(
                wall_thermal(mesh.boundary().patch(i).name, 300.0));
        setup.has_initial_condition = true;
        setup.initial_condition.temperature = 300.0;

        Field<double, Location::CELL> U(
            mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(
            mesh.n_cells(), "p", "Pa", 1);
        auto controls = vv_controls(1000.0, 1e-6);
        const auto loop = run_production_buoyant_outer_loop(
            mesh, ubc, pbc, controls, setup, U, p);

        EXPECT_TRUE(loop.converged);
        EXPECT_TRUE(loop.outer_iterations == 2);
        EXPECT_TRUE(loop.energy.converged);
        EXPECT_TRUE(loop.energy_ran);

        Field<double, Location::CELL> U_ref(
            mesh.n_cells(), "U_ref", "m/s", 3);
        Field<double, Location::CELL> p_ref(
            mesh.n_cells(), "p_ref", "Pa", 1);
        auto ref_controls = vv_controls(1000.0, 1e-6);
        const auto ref1 = solve_steady_incompressible(
            mesh, U_ref, p_ref, ubc, pbc, ref_controls);
        EXPECT_TRUE(ref1.converged);
        ref_controls.initialization.mode = InitializationMode::Provided;
        const auto ref2 = solve_steady_incompressible(
            mesh, U_ref, p_ref, ubc, pbc, ref_controls);
        EXPECT_TRUE(ref2.converged);

        expect_identical_velocity(U, U_ref, "identity U");
        expect_identical_fields(p, p_ref, "identity p");
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            EXPECT_TRUE(std::isfinite(loop.temperature(c)));
            EXPECT_NEAR(loop.temperature(c), 300.0, 1e-8);
        }
    });

    // -----------------------------------------------------------------
    // C. Natural convection fixed point: differentially heated cavity
    //    (left 320 K, right 280 K, other walls adiabatic), gravity
    //    "0 -9.81 0" in the mesh plane. The loop must converge the
    //    COUPLED fixed point, not a single Picard pass.
    // -----------------------------------------------------------------
    run_case("production_buoyant_outer_loop_natural_convection_fixed_point", [] {
        const std::size_t nx = 4;
        const std::size_t ny = 4;
        const double dx = 0.0625;
        const double dy = 0.0625;
        auto mesh = box_mesh(nx, ny, dx, dy);
        const auto ubc = no_slip_velocity_bcs(mesh);
        const auto pbc = zero_gradient_pressure_bcs(mesh);
        const auto geometry = build_fv_geometry(mesh);

        CaseSetup setup;
        setup.energy_model = "buoyant";
        setup.ref_temperature = 300.0;
        setup.gravity_vector = "0 -9.81 0";
        setup.materials.push_back(buoyant_material("0.0001"));
        setup.boundary_conditions.push_back(wall_thermal("left", 320.0));
        setup.boundary_conditions.push_back(wall_thermal("right", 280.0));
        // bottom / top / front / back: deliberately no spec (adiabatic).
        setup.has_initial_condition = true;
        setup.initial_condition.temperature = 300.0;

        Field<double, Location::CELL> U(
            mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(
            mesh.n_cells(), "p", "Pa", 1);
        auto controls = vv_controls(1000.0, 0.01);
        const auto loop = run_production_buoyant_outer_loop(
            mesh, ubc, pbc, controls, setup, U, p);

        EXPECT_TRUE(loop.converged);
        EXPECT_TRUE(loop.fixed_point_temperature_passed);
        EXPECT_TRUE(loop.fixed_point_velocity_passed);
        EXPECT_TRUE(loop.fixed_point_pressure_passed);
        EXPECT_TRUE(loop.max_temperature_change_normalized <= 1e-10);
        EXPECT_TRUE(loop.max_velocity_change_normalized <= 1e-10);
        EXPECT_TRUE(loop.max_pressure_change_normalized <= 1e-10);
        EXPECT_TRUE(loop.outer_iterations >= 3);
        EXPECT_TRUE(loop.outer_iterations <= 200);
        EXPECT_TRUE(loop.energy.converged);
        EXPECT_TRUE(loop.energy.global_balance_normalized <= 1e-12);
        // Calibrated local audit floor (diag run 39633631822): the
        // local and per-cell normalized balances bottom out at the
        // floating-point summation floor (~1.4e-12 measured); the
        // gate is 1e-11, a margin of more than 7x above the
        // measured floor, matching the calibrated stage control
        // (fiche 14).
        EXPECT_TRUE(loop.energy.local_balance_normalized <= 1e-11);
        EXPECT_TRUE(loop.energy.max_cell_residual_normalized <= 1e-11);
        for (const auto& entry : loop.history) {
            EXPECT_TRUE(entry.flow_converged);
            EXPECT_TRUE(entry.energy_converged);
        }

        // Nonzero circulation: the frozen carrier of this case is
        // motionless (no body force), so any velocity is buoyancy.
        double max_speed = 0.0;
        for (std::size_t d = 0; d < 3; ++d)
            for (std::size_t c = 0; c < mesh.n_cells(); ++c)
                max_speed = std::max(
                    max_speed,
                    std::abs(U.component_data(d)[c]));
        EXPECT_TRUE(max_speed > 1e-12);

        // Physical sign (fiche 10): the hot left wall rises (+y),
        // the cold right wall sinks (-y).
        double hot_mean_vy = 0.0;
        double cold_mean_vy = 0.0;
        for (std::size_t j = 0; j < ny; ++j) {
            hot_mean_vy += U.component_data(1)[0 * ny + j];
            cold_mean_vy += U.component_data(1)[(nx - 1) * ny + j];
        }
        hot_mean_vy /= static_cast<double>(ny);
        cold_mean_vy /= static_cast<double>(ny);
        EXPECT_TRUE(hot_mean_vy > 0.0);
        EXPECT_TRUE(cold_mean_vy < 0.0);

        // The coupled temperature is NOT the frozen-flow conduction
        // solution of the same case.
        Field<double, Location::FACE> zero_flux(
            mesh.n_faces(), "mass_flux", "kg/s", 1);
        zero_flux.fill(0.0);
        Field<double, Location::CELL> conduction(
            mesh.n_cells(), "conduction", "K", 1);
        conduction.fill(300.0);
        const auto frozen = run_production_energy_stage(
            mesh, geometry, setup, zero_flux, conduction);
        EXPECT_TRUE(frozen.converged);
        double max_temperature_difference = 0.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            max_temperature_difference = std::max(
                max_temperature_difference,
                std::abs(loop.temperature(c) - conduction(c)));
        EXPECT_TRUE(max_temperature_difference > 1e-10);

        // -------------------------------------------------------------
        // Fixed point replayed INDEPENDENTLY: rebuild the source from
        // the final UNRELAXED response, re-solve the carrier warm and
        // replay one energy pass. The replayed field must reproduce
        // the loop's final field. Tolerance note: the applied source
        // was built from the relaxed field, whose lag from the final
        // response is bounded by ((1 - omega)/omega) * tolerance *
        // scale ~ 1e-7 K at exit; 1e-6 K bounds the replayed
        // difference with orders of margin while still excluding any
        // intermediate Picard pass (differences of order K).
        // -------------------------------------------------------------
        auto replay_controls = vv_controls(1000.0, 0.01);
        build_cell_buoyancy_acceleration(
            replay_controls.cell_body_force, loop.temperature,
            loop.coupling);
        Field<double, Location::CELL> U_replay(
            mesh.n_cells(), "U_replay", "m/s", 3);
        Field<double, Location::CELL> p_replay(
            mesh.n_cells(), "p_replay", "Pa", 1);
        U_replay = U;
        p_replay = p;
        const auto replay_flow = solve_steady_incompressible(
            mesh, U_replay, p_replay, ubc, pbc, replay_controls);
        EXPECT_TRUE(replay_flow.converged);
        Field<double, Location::CELL> replay_temperature(
            mesh.n_cells(), "replay", "K", 1);
        replay_temperature = loop.temperature;
        const auto replay_energy = run_production_energy_stage(
            mesh, geometry, setup, replay_flow.authoritative_mass_flux,
            replay_temperature);
        EXPECT_TRUE(replay_energy.converged);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            EXPECT_NEAR(replay_temperature(c), loop.temperature(c), 1e-6);
    });

    // -----------------------------------------------------------------
    // D. Nested S2S: with declared S2S radiation the inner response
    //    operator is the stage-3 energy-radiation fixed point, and its
    //    gates stay green on every outer iteration.
    // -----------------------------------------------------------------
    run_case("production_buoyant_outer_loop_nested_s2s_radiation", [] {
        const std::size_t nx = 4;
        const std::size_t ny = 4;
        auto mesh = box_mesh(nx, ny, 0.0625, 0.0625);
        const auto ubc = no_slip_velocity_bcs(mesh);
        const auto pbc = zero_gradient_pressure_bcs(mesh);

        CaseSetup setup;
        setup.energy_model = "buoyant";
        setup.radiation_model = "S2S";
        setup.ref_temperature = 300.0;
        setup.gravity_vector = "0 -9.81 0";
        setup.materials.push_back(buoyant_material("0.0001"));
        setup.boundary_conditions.push_back(wall_thermal("left", 320.0));
        setup.boundary_conditions.push_back(wall_thermal("right", 280.0));

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
        // front / back: deliberately no spec (adiabatic walls).

        setup.has_initial_condition = true;
        setup.initial_condition.temperature = 300.0;

        Field<double, Location::CELL> U(
            mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(
            mesh.n_cells(), "p", "Pa", 1);
        auto controls = vv_controls(1000.0, 0.01);
        const auto loop = run_production_buoyant_outer_loop(
            mesh, ubc, pbc, controls, setup, U, p);

        EXPECT_TRUE(loop.converged);
        EXPECT_TRUE(loop.s2s_ran);
        EXPECT_TRUE(loop.outer_iterations >= 2);
        EXPECT_TRUE(loop.max_temperature_change_normalized <= 1e-10);
        EXPECT_TRUE(loop.max_velocity_change_normalized <= 1e-10);
        EXPECT_TRUE(loop.max_pressure_change_normalized <= 1e-10);
        // Inner gates of the nested stage-3 fixed point (fiche 11)
        // and the energy stage (fiche 14) on the final carrier.
        EXPECT_TRUE(loop.s2s.converged);
        EXPECT_TRUE(loop.s2s.radiation.balance_normalized <= 1e-10);
        EXPECT_TRUE(loop.s2s.radiation.max_surface_identity_normalized <=
                     1e-10);
        EXPECT_TRUE(loop.energy.converged);
        EXPECT_TRUE(loop.energy.global_balance_normalized <= 1e-12);
        // Calibrated local audit floor (diag run 39633631822): the
        // local and per-cell normalized balances bottom out at the
        // floating-point summation floor (~1.4e-12 measured); the
        // gate is 1e-11, a margin of more than 7x above the
        // measured floor, matching the calibrated stage control
        // (fiche 14).
        EXPECT_TRUE(loop.energy.local_balance_normalized <= 1e-11);
        EXPECT_TRUE(loop.energy.max_cell_residual_normalized <= 1e-11);
        for (const auto& entry : loop.history) {
            EXPECT_TRUE(entry.flow_converged);
            EXPECT_TRUE(entry.energy_converged);
            EXPECT_TRUE(entry.s2s_ran);
            EXPECT_TRUE(entry.s2s_converged);
        }

        // Radiative pair closure (fiche 11 B3): the exchange between the
        // two participants is internal, the hot bottom loses, the cold
        // top gains, the pair closes to machine zero.
        const auto bi = surface_index(loop.s2s.radiation, "bottom");
        const auto ti = surface_index(loop.s2s.radiation, "top");
        EXPECT_TRUE(loop.s2s.radiation.surface_count == 2);
        EXPECT_NEAR(loop.s2s.radiation.net_flux[bi] +
                        loop.s2s.radiation.net_flux[ti],
                    0.0, 1e-11);

        // The carrier still circulates (the radiation adds to the
        // differential heating drive).
        double max_speed = 0.0;
        for (std::size_t d = 0; d < 3; ++d)
            for (std::size_t c = 0; c < mesh.n_cells(); ++c)
                max_speed = std::max(
                    max_speed,
                    std::abs(U.component_data(d)[c]));
        EXPECT_TRUE(max_speed > 1e-12);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            EXPECT_TRUE(std::isfinite(loop.temperature(c)));
    });

    // -----------------------------------------------------------------
    // Fail-closed catalogue on a 2x2x1 box, plus the accepted positives.
    // -----------------------------------------------------------------
    run_case("production_buoyant_coupling_fail_closed", [] {
        const std::size_t nx = 2;
        const std::size_t ny = 2;
        auto mesh = box_mesh(nx, ny, 0.25, 0.25);
        const auto ubc = no_slip_velocity_bcs(mesh);
        const auto pbc = zero_gradient_pressure_bcs(mesh);
        const auto geometry = build_fv_geometry(mesh);
        const std::size_t n_cells = mesh.n_cells();

        auto valid_setup = [] {
            CaseSetup setup;
            setup.energy_model = "buoyant";
            setup.ref_temperature = 300.0;
            setup.gravity_vector = "0 0 -9.81";
            setup.materials.push_back(buoyant_material("0.001"));
            setup.boundary_conditions.push_back(
                wall_thermal("left", 320.0));
            setup.boundary_conditions.push_back(
                wall_thermal("right", 280.0));
            setup.has_initial_condition = true;
            setup.initial_condition.temperature = 300.0;
            return setup;
        };

        // gravity_vector: absent / short / long / non-numeric / zero.
        {
            auto s = valid_setup();
            s.gravity_vector = "";
            EXPECT_THROW_WITH(make_production_buoyancy_coupling(s),
                              std::invalid_argument,
                              "declared gravity_vector");
        }
        {
            auto s = valid_setup();
            s.gravity_vector = "0 -9.81";
            EXPECT_THROW_WITH(make_production_buoyancy_coupling(s),
                              std::invalid_argument,
                              "exactly 3");
        }
        {
            auto s = valid_setup();
            s.gravity_vector = "0 0 -9.81 1";
            EXPECT_THROW_WITH(make_production_buoyancy_coupling(s),
                              std::invalid_argument,
                              "exactly 3");
        }
        {
            auto s = valid_setup();
            s.gravity_vector = "0 zero -9.81";
            EXPECT_THROW_WITH(make_production_buoyancy_coupling(s),
                              std::invalid_argument,
                              "exactly 3");
        }
        {
            auto s = valid_setup();
            s.gravity_vector = "0 0 0";
            EXPECT_THROW_WITH(make_production_buoyancy_coupling(s),
                              std::invalid_argument,
                              "zero gravity_vector");
        }

        // thermal_expansion_coefficient: missing / zero / negative /
        // trailing garbage / non-finite.
        {
            auto s = valid_setup();
            s.materials.front().extra_properties.clear();
            EXPECT_THROW_WITH(make_production_buoyancy_coupling(s),
                              std::invalid_argument,
                              "no default is invented");
        }
        {
            auto s = valid_setup();
            s.materials.front().extra_properties[
                "thermal_expansion_coefficient"] = "0";
            EXPECT_THROW_WITH(make_production_buoyancy_coupling(s),
                              std::invalid_argument,
                              "strictly positive");
        }
        {
            auto s = valid_setup();
            s.materials.front().extra_properties[
                "thermal_expansion_coefficient"] = "-0.001";
            EXPECT_THROW_WITH(make_production_buoyancy_coupling(s),
                              std::invalid_argument,
                              "strictly positive");
        }
        {
            auto s = valid_setup();
            s.materials.front().extra_properties[
                "thermal_expansion_coefficient"] = "0.001abc";
            EXPECT_THROW_WITH(make_production_buoyancy_coupling(s),
                              std::invalid_argument,
                              "trailing characters");
        }
        {
            auto s = valid_setup();
            s.materials.front().extra_properties[
                "thermal_expansion_coefficient"] = "inf";
            EXPECT_THROW_WITH(make_production_buoyancy_coupling(s),
                              std::invalid_argument,
                              "must be finite");
        }

        // A buoyant case must not declare a uniform case body_force.
        {
            auto s = valid_setup();
            s.body_force = {0.0, 0.0, -9.81};
            EXPECT_THROW_WITH(make_production_buoyancy_coupling(s),
                              std::invalid_argument,
                              "must not declare a uniform body_force");
        }

        // Material count: zero and two materials are both rejected.
        {
            auto s = valid_setup();
            s.materials.clear();
            EXPECT_THROW_WITH(make_production_buoyancy_coupling(s),
                              std::invalid_argument,
                              "exactly one material");
        }
        {
            auto s = valid_setup();
            s.materials.push_back(buoyant_material("0.001"));
            EXPECT_THROW_WITH(make_production_buoyancy_coupling(s),
                              std::invalid_argument,
                              "exactly one material");
        }

        // ref_temperature must be finite and positive.
        {
            auto s = valid_setup();
            s.ref_temperature = 0.0;
            EXPECT_THROW_WITH(make_production_buoyancy_coupling(s),
                              std::invalid_argument,
                              "ref_temperature must be finite and "
                              "positive");
        }

        // The case-load validator is a no-op for any other energy
        // model and accepts the valid buoyant case without throwing.
        {
            auto s = valid_setup();
            s.energy_model = "conjugate";
            validate_production_buoyancy_model(s);
        }
        validate_production_buoyancy_model(valid_setup());

        // The loop itself refuses a non-buoyant case.
        {
            auto s = valid_setup();
            s.energy_model = "conjugate";
            Field<double, Location::CELL> U(
                n_cells, "U", "m/s", 3);
            Field<double, Location::CELL> p(
                n_cells, "p", "Pa", 1);
            auto controls = vv_controls(1000.0, 1e-6);
            EXPECT_THROW_WITH(
                run_production_buoyant_outer_loop(
                    mesh, ubc, pbc, controls, s, U, p),
                std::invalid_argument, "is not declared");
        }

        // Strict gravity parse positive: surrounding whitespace and
        // exact components.
        {
            const auto g = parse_gravity_vector("  0   0  -9.81  ");
            EXPECT_TRUE(g.x == 0.0);
            EXPECT_TRUE(g.y == 0.0);
            EXPECT_TRUE(g.z == -9.81);
        }

        // Solver channel fail-closed (direct carrier solves): size
        // mismatch, uniform+cellwise combination, non-finite
        // components.
        {
            Field<double, Location::CELL> U(n_cells, "U", "m/s", 3);
            Field<double, Location::CELL> p(n_cells, "p", "Pa", 1);
            auto controls = vv_controls(1000.0, 1e-6);
            controls.cell_body_force.assign(
                n_cells + 1, Vec3{0.0, 0.0, 0.0});
            EXPECT_THROW_WITH(
                solve_steady_incompressible(mesh, U, p, ubc, pbc, controls),
                std::invalid_argument,
                "cellwise body force size does not match mesh");
        }
        {
            Field<double, Location::CELL> U(n_cells, "U", "m/s", 3);
            Field<double, Location::CELL> p(n_cells, "p", "Pa", 1);
            auto controls = vv_controls(1000.0, 1e-6);
            controls.body_force = Vec3{1.0, 0.0, 0.0};
            controls.cell_body_force.assign(
                n_cells, Vec3{0.0, 0.0, 0.0});
            EXPECT_THROW_WITH(
                solve_steady_incompressible(mesh, U, p, ubc, pbc, controls),
                std::invalid_argument,
                "not both");
        }
        {
            Field<double, Location::CELL> U(n_cells, "U", "m/s", 3);
            Field<double, Location::CELL> p(n_cells, "p", "Pa", 1);
            auto controls = vv_controls(1000.0, 1e-6);
            controls.cell_body_force.assign(
                n_cells, Vec3{0.0, 0.0, 0.0});
            controls.cell_body_force[0].x =
                std::numeric_limits<double>::quiet_NaN();
            EXPECT_THROW_WITH(
                solve_steady_incompressible(mesh, U, p, ubc, pbc, controls),
                std::invalid_argument,
                "cellwise body force components must be finite");
        }

        // Positives: the one-shot energy stage and the S2S stages
        // accept a valid buoyant case (widened gates, unchanged
        // criteria).
        {
            Field<double, Location::FACE> zero_flux(
                mesh.n_faces(), "mass_flux", "kg/s", 1);
            zero_flux.fill(0.0);
            Field<double, Location::CELL> temperature(
                n_cells, "temperature", "K", 1);
            temperature.fill(300.0);
            const auto energy = run_production_energy_stage(
                mesh, geometry, valid_setup(), zero_flux, temperature);
            EXPECT_TRUE(energy.converged);
            EXPECT_TRUE(energy.global_balance_normalized <= 1e-12);
        }
        {
            auto s = valid_setup();
            s.radiation_model = "S2S";
            BoundarySpec bottom;
            bottom.patch_name = "bottom";
            bottom.type = BCType::WALL;
            bottom.raw_params["emissivity"] = "0.8";
            bottom.raw_params["heat_flux"] = "50";
            bottom.raw_params["s2s_view_factor_to_top"] = "1";
            s.boundary_conditions.push_back(bottom);
            BoundarySpec top;
            top.patch_name = "top";
            top.type = BCType::WALL;
            top.raw_params["emissivity"] = "0.9";
            top.raw_params["s2s_view_factor_to_bottom"] = "1";
            s.boundary_conditions.push_back(top);

            Field<double, Location::FACE> zero_flux(
                mesh.n_faces(), "mass_flux", "kg/s", 1);
            zero_flux.fill(0.0);
            Field<double, Location::CELL> temperature(
                n_cells, "temperature", "K", 1);
            temperature.fill(300.0);
            const auto s2s = run_production_s2s_energy_stage(
                mesh, geometry, s, zero_flux, temperature);
            EXPECT_TRUE(s2s.converged);
            EXPECT_TRUE(s2s.energy.converged);
        }
    });

    return run_all();
}
