// Production Boussinesq buoyancy coupling (capability #885 step 4, issue
// #885): the energy field couples back into the incompressible carrier
// through the cellwise buoyancy source a_i = -beta (T_i - T_ref) g_vec
// (fiche methode 10, APPROVED; carrier fiche 04; energy stage fiche 14;
// S2S stage fiche 11 when radiation is declared).
//
// Verified properties:
//   * Unit gates: the strict gravity_vector parser accepts exactly three
//     whitespace-separated finite tokens and rejects missing components,
//     trailing text and non-numeric/non-finite tokens; the source builder
//     carries the exact Boussinesq sign convention (a hot cell accelerates
//     AGAINST gravity, a cell at T_ref carries exactly zero source) with
//     bit-identical expressions; the fail-closed declaration catalogue
//     (missing/non-positive/non-numeric beta, missing/zero/malformed
//     gravity, uniform body_force together with buoyant, two materials,
//     non-positive reference temperature) is rejected instead of silently
//     running without buoyancy.
//   * Solver channel: the cellwise body force replaces the uniform vector
//     per cell; a wrong size, a non-finite entry and the two-channel mix
//     are rejected by the carrier.
//   * Identity at T = T_ref: the source is exactly zero, so the buoyant
//     loop reproduces the buoyancy-free carrier (the uniform temperature
//     field stays at T_ref and U/p agree with a plain solve).
//   * Uniform equivalence: for a uniform temperature field away from
//     T_ref the source is constant, and a constant cellwise field
//     reproduces the uniform controls.body_force assembly BIT FOR BIT
//     (direct solves agree exactly); the loop converges to the coupled
//     fixed point whose carrier state agrees with the equivalent uniform
//     body-force solve.
//   * Coupled fixed point: a differentially heated closed cavity converges
//     with ALL FOUR blocking criteria (temperature, velocity, pressure
//     stationarity AND the source-input lag) at 1e-10 normalized, with
//     green inner stages recorded at every outer iteration and a first
//     pass whose stationarity is infinite (a single pass can never claim
//     the fixed point). The buoyancy source drives a non-zero flow where
//     the source-free carrier is exactly at rest. Replaying ONE outer
//     iteration by hand from the final field (source rebuilt from the
//     final response, warm-started carrier, energy stage on the new mass
//     flux) reproduces the loop's final responses.
//   * Nested S2S: with radiation declared the loop reuses the stage-3
//     energy-radiation outer loop UNCHANGED as its inner stage (fiche 11
//     gates widened to accept the buoyant model): the coupled fixed point
//     converges, the radiosity machine reports its own green gates, the
//     final radiative state satisfies the independent hand Cramer oracle
//     and the parallel-plate identity, the surface pair closes to machine
//     zero with the fiche 11 B3 sign convention, and the walls carry the
//     declared 50 W/m^2.
//   * Fail-closed: running the buoyancy loop under another declared
//     energy model, with a mismatched temperature field, or without any
//     thermal Dirichlet patch is rejected, never silently degraded.
#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/physics/production_energy_stage.h"
#include "cfdx/physics/production_buoyancy_stage.h"
#include "cfdx/physics/production_radiation_contract.h"
#include "cfdx/physics/production_radiation_s2s_stage.h"
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
// Same construction as the stage-2 / stage-3 VVs so all production tests
// exercise the identical geometry path.
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
            m.ownership().set_neighbour(fz(i, j, l), FaceOwnership::BOUNDARY);
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
// Shared case/controls builders.
// ---------------------------------------------------------------------------

// One viscous-oil material carrying the Boussinesq beta declaration
// (extra_properties convention shared with volumetric_power). nu = mu/rho
// = 1e-3 m^2/s, matching the carrier controls below.
static MaterialSpec make_buoyant_material(const std::string& beta_text)
{
    MaterialSpec material;
    material.name = "viscous_oil";
    material.density = 1000.0;
    material.dynamic_viscosity = 1.0;
    material.thermal_conductivity = 1.0;
    material.specific_heat = 1000.0;
    material.extra_properties["thermal_expansion_coefficient"] = beta_text;
    return material;
}

static BoundarySpec thermal_dirichlet(const std::string& patch,
                                      double wall_temperature)
{
    BoundarySpec spec;
    spec.patch_name = patch;
    spec.type = BCType::WALL;
    spec.value_type = BCValueType::WALL_THERMAL;
    spec.wall_temperature = wall_temperature;
    return spec;
}

// Carrier controls: strict per-criterion tolerances (1e-12) so a converged
// warm-started solve changes U/p by far less than the loop's 1e-10
// stationarity criteria, with a generous iteration budget.
static IncompressibleSolverControls make_flow_controls()
{
    IncompressibleSolverControls controls;
    controls.algorithm = PressureVelocityAlgorithm::SIMPLE;
    controls.initialization.mode = InitializationMode::Uniform;
    controls.convergence.max_iterations = 500;
    controls.convergence.relative_tolerance = 1e-12;
    controls.convergence.continuity_tolerance = 1e-12;
    controls.linear_tolerance = 1e-12;
    controls.pressure_reference_cell = 0;
    controls.pressure_reference_value = 0.0;
    controls.density = 1000.0;
    controls.kinematic_viscosity = 1e-3;
    return controls;
}

static void add_closed_cavity_bcs(VelocityBoundaryConditions& ubc,
                                   ScalarBoundaryConditions& pbc)
{
    for (const char* name :
         {"left", "right", "bottom", "top", "front", "back"}) {
        ubc[name] = {VelocityBoundaryCondition::Type::FIXED_VALUE,
                     {0.0, 0.0, 0.0}};
        pbc[name] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    }
}

// Base buoyant case: one material with beta, a gravity vector, a
// reference temperature and one thermal Dirichlet patch (the minimum the
// stages require); every case specializes the walls.
static CaseSetup make_buoyant_base()
{
    CaseSetup setup;
    setup.energy_model = "buoyant";
    setup.radiation_model = "none";
    setup.gravity_vector = "0 0 -9.81";
    setup.ref_temperature = 300.0;
    setup.materials.push_back(make_buoyant_material("2e-4"));
    setup.boundary_conditions.push_back(thermal_dirichlet("left", 300.0));
    return setup;
}

static double max_abs_scalar(const Field<double, Location::CELL>& field)
{
    double value = 0.0;
    for (std::size_t c = 0; c < field.size(); ++c)
        value = std::max(value, std::abs(field(c)));
    return value;
}

static double max_abs_velocity(const Field<double, Location::CELL>& field)
{
    double value = 0.0;
    for (std::size_t d = 0; d < 3; ++d)
        for (std::size_t c = 0; c < field.size(); ++c)
            value = std::max(value, std::abs(field.component_data(d)[c]));
    return value;
}

// Independent |Sf|-weighted wall temperature of a patch (the only link
// between the resolved thermal field and the radiosity system), recomputed
// here without any production helper. Same oracle as the stage-3 VV.
static double patch_wall_temperature(
    const Mesh& mesh, const FvGeometry& geometry, const std::string& name,
    const Field<double, Location::CELL>& temperature)
{
    double weighted = 0.0;
    double total = 0.0;
    for (std::size_t p = 0; p < mesh.boundary().n_patches(); ++p) {
        const auto& patch = mesh.boundary().patch(p);
        if (patch.name != name)
            continue;
        for (const auto face : patch.face_ids) {
            const auto& area = geometry.face_area_vectors[face];
            const double weight = std::sqrt(
                area.x * area.x + area.y * area.y + area.z * area.z);
            const auto owner = static_cast<std::size_t>(
                mesh.ownership().owner(face));
            weighted += weight * temperature(owner);
            total += weight;
        }
    }
    return weighted / total;
}

static std::size_t surface_index(const ProductionS2SRadiationResult& result,
                                  const std::string& name)
{
    for (std::size_t i = 0; i < result.surface_names.size(); ++i)
        if (result.surface_names[i] == name)
            return i;
    throw std::runtime_error("surface not found: " + name);
}

// Common verdict of a converged buoyancy loop: all four blocking criteria
// green, honest bookkeeping, green inner stages at every iteration and a
// first pass that can never claim the fixed point.
static void expect_converged_loop(const ProductionBuoyancyOuterLoopResult& result)
{
    EXPECT_TRUE(result.converged);
    EXPECT_TRUE(result.failure_reason.empty());
    EXPECT_TRUE(result.outer_iterations >= 2);
    EXPECT_TRUE(result.max_temperature_change_normalized <= 1e-10);
    EXPECT_TRUE(result.max_velocity_change_normalized <= 1e-10);
    EXPECT_TRUE(result.max_pressure_change_normalized <= 1e-10);
    EXPECT_TRUE(result.max_source_input_lag_normalized <= 1e-10);
    EXPECT_TRUE(result.energy_ran);
    EXPECT_TRUE(result.outer_history.size() == result.outer_iterations);
    for (const auto& entry : result.outer_history) {
        EXPECT_TRUE(entry.carrier_converged);
        EXPECT_TRUE(entry.energy_converged);
        EXPECT_TRUE(std::isfinite(entry.max_source_input_lag_normalized));
        EXPECT_TRUE(entry.max_source_input_lag_normalized >= 0.0);
    }
    EXPECT_TRUE(std::isinf(
        result.outer_history.front().max_temperature_change_normalized));
    EXPECT_TRUE(std::isinf(
        result.outer_history.front().max_velocity_change_normalized));
    EXPECT_TRUE(std::isinf(
        result.outer_history.front().max_pressure_change_normalized));
}

int main()
{
    // -----------------------------------------------------------------
    // Unit gates: strict gravity parsing, exact Boussinesq signs, the
    // fail-closed declaration catalogue and the carrier channel
    // validation.
    // -----------------------------------------------------------------
    run_case("production_buoyancy_unit_gates", [] {
        // Strict gravity_vector parser: exactly three finite tokens.
        const Vec3 canonical = parse_gravity_vector("0 0 -9.81");
        EXPECT_TRUE(canonical.x == 0.0);
        EXPECT_TRUE(canonical.y == 0.0);
        EXPECT_TRUE(canonical.z == -9.81);
        const Vec3 padded = parse_gravity_vector("  0.5  -1.0   9.81  ");
        EXPECT_TRUE(padded.x == 0.5);
        EXPECT_TRUE(padded.y == -1.0);
        EXPECT_TRUE(padded.z == 9.81);
        EXPECT_THROW(parse_gravity_vector(""), std::invalid_argument);
        EXPECT_THROW(parse_gravity_vector("0 0"), std::invalid_argument);
        EXPECT_THROW(parse_gravity_vector("0"), std::invalid_argument);
        EXPECT_THROW(parse_gravity_vector("a b c"), std::invalid_argument);
        EXPECT_THROW(parse_gravity_vector("0 0 -9.81 trailing"),
                     std::invalid_argument);
        EXPECT_THROW(parse_gravity_vector("0 0 inf"), std::invalid_argument);
        EXPECT_THROW(parse_gravity_vector("0 0 nan"), std::invalid_argument);
        EXPECT_THROW(parse_gravity_vector("1e309 0 0"),
                     std::invalid_argument);

        // Exact source builder signs (fiche 10 Boussinesq contract):
        // a_i = -beta (T_i - T_ref) g_vec - a hot cell accelerates against
        // gravity, a cell at T_ref carries exactly zero source.
        ProductionBuoyancyCoupling coupling;
        coupling.beta = 2e-4;
        coupling.reference_temperature = 300.0;
        coupling.gravity = parse_gravity_vector("0.5 -1.0 -9.81");
        Field<double, Location::CELL> temperature(3, "T", "K", 1);
        temperature(0) = 310.0;  // hot
        temperature(1) = 290.0;  // cold
        temperature(2) = 300.0;  // at the reference
        std::vector<Vec3> acceleration;
        build_cell_buoyancy_acceleration(acceleration, temperature, coupling);
        EXPECT_TRUE(acceleration.size() == 3);
        const double delta_hot = 310.0 - 300.0;
        const double delta_cold = 290.0 - 300.0;
        EXPECT_NEAR(acceleration[0].x,
                    -coupling.beta * delta_hot * coupling.gravity.x, 1e-18);
        EXPECT_NEAR(acceleration[0].y,
                    -coupling.beta * delta_hot * coupling.gravity.y, 1e-18);
        EXPECT_NEAR(acceleration[0].z,
                    -coupling.beta * delta_hot * coupling.gravity.z, 1e-18);
        EXPECT_NEAR(acceleration[1].z,
                    -coupling.beta * delta_cold * coupling.gravity.z, 1e-18);
        EXPECT_TRUE(acceleration[0].z > 0.0);  // hot rises against gravity
        EXPECT_TRUE(acceleration[1].z < 0.0);  // cold sinks along gravity
        EXPECT_TRUE(acceleration[2].x == 0.0);  // T = T_ref: exactly zero
        EXPECT_TRUE(acceleration[2].y == 0.0);
        EXPECT_TRUE(acceleration[2].z == 0.0);

        // Selection + coupling parse of the declared case.
        auto good = make_buoyant_base();
        EXPECT_TRUE(production_buoyancy_requested(good));
        validate_production_buoyancy_model(good);  // no-op throw-free pass
        const ProductionBuoyancyCoupling parsed =
            production_buoyancy_coupling(good);
        EXPECT_TRUE(parsed.beta == 2e-4);
        EXPECT_TRUE(parsed.reference_temperature == 300.0);
        EXPECT_TRUE(parsed.gravity.x == 0.0);
        EXPECT_TRUE(parsed.gravity.y == 0.0);
        EXPECT_TRUE(parsed.gravity.z == -9.81);

        // A non-buoyant case is not owned by this gate: no declarations
        // required, no rejection.
        auto other = good;
        other.energy_model = "conjugate";
        EXPECT_TRUE(!production_buoyancy_requested(other));
        validate_production_buoyancy_model(other);  // silent no-op

        // Fail-closed declaration catalogue: every broken buoyant case is
        // rejected instead of silently running without buoyancy.
        const auto expect_rejected = [](const CaseSetup& broken) {
            EXPECT_THROW(validate_production_buoyancy_model(broken),
                         std::invalid_argument);
        };
        {
            auto broken = good;
            broken.body_force = {0.0, 0.0, -9.81};
            expect_rejected(broken);  // buoyant + uniform body_force
        }
        {
            auto broken = good;
            broken.materials.push_back(make_buoyant_material("2e-4"));
            expect_rejected(broken);  // two materials
        }
        {
            auto broken = good;
            broken.materials.front().extra_properties.erase(
                "thermal_expansion_coefficient");
            expect_rejected(broken);  // beta not declared
        }
        {
            auto broken = good;
            broken.materials.front().extra_properties[
                "thermal_expansion_coefficient"] = "0";
            expect_rejected(broken);  // beta not strictly positive
        }
        {
            auto broken = good;
            broken.materials.front().extra_properties[
                "thermal_expansion_coefficient"] = "-2e-4";
            expect_rejected(broken);
        }
        {
            auto broken = good;
            broken.materials.front().extra_properties[
                "thermal_expansion_coefficient"] = "two";
            expect_rejected(broken);  // strict case-double parse
        }
        {
            auto broken = good;
            broken.ref_temperature = 0.0;
            expect_rejected(broken);
        }
        {
            auto broken = good;
            broken.ref_temperature = -5.0;
            expect_rejected(broken);
        }
        {
            auto broken = good;
            broken.ref_temperature =
                std::numeric_limits<double>::quiet_NaN();
            expect_rejected(broken);
        }
        {
            auto broken = good;
            broken.gravity_vector = "";
            expect_rejected(broken);  // gravity not declared
        }
        {
            auto broken = good;
            broken.gravity_vector = "0 0 0";
            expect_rejected(broken);  // zero gravity magnitude
        }
        {
            auto broken = good;
            broken.gravity_vector = "0 0";
            expect_rejected(broken);  // malformed
        }
        {
            auto broken = good;
            broken.gravity_vector = "0 0 -9.81 extra";
            expect_rejected(broken);
        }

        // Carrier channel validation: wrong size, non-finite entries and
        // the two-channel mix are rejected by the solve itself.
        const std::size_t nx = 2;
        const std::size_t ny = 2;
        auto mesh = box_mesh(nx, ny, 0.25, 0.25);
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        add_closed_cavity_bcs(ubc, pbc);
        const std::size_t n_cells = mesh.n_cells();
        {
            Field<double, Location::CELL> U(n_cells, "U", "m/s", 3);
            Field<double, Location::CELL> p(n_cells, "p", "Pa", 1);
            auto controls = make_flow_controls();
            controls.cell_body_force.assign(n_cells + 1,
                                            {0.0, 0.0, 0.0});
            EXPECT_THROW(
                solve_steady_incompressible(mesh, U, p, ubc, pbc, controls),
                std::invalid_argument);
        }
        {
            Field<double, Location::CELL> U(n_cells, "U", "m/s", 3);
            Field<double, Location::CELL> p(n_cells, "p", "Pa", 1);
            auto controls = make_flow_controls();
            controls.cell_body_force.assign(n_cells, {0.0, 0.0, 0.0});
            controls.body_force = {0.0, 0.0, -9.81};
            EXPECT_THROW(
                solve_steady_incompressible(mesh, U, p, ubc, pbc, controls),
                std::invalid_argument);
        }
        {
            Field<double, Location::CELL> U(n_cells, "U", "m/s", 3);
            Field<double, Location::CELL> p(n_cells, "p", "Pa", 1);
            auto controls = make_flow_controls();
            controls.cell_body_force.assign(
                n_cells, {std::numeric_limits<double>::quiet_NaN(),
                          0.0, 0.0});
            EXPECT_THROW(
                solve_steady_incompressible(mesh, U, p, ubc, pbc, controls),
                std::invalid_argument);
        }
    });

    // -----------------------------------------------------------------
    // Identity at T = T_ref: every wall at the reference temperature, so
    // the Boussinesq source is exactly zero everywhere. The buoyant loop
    // must reproduce the buoyancy-free carrier and the uniform reference
    // temperature field.
    // -----------------------------------------------------------------
    run_case("production_buoyancy_identity_at_reference_temperature", [] {
        const std::size_t nx = 4;
        const std::size_t ny = 4;
        auto mesh = box_mesh(nx, ny, 0.25, 0.25);
        const std::size_t n_cells = mesh.n_cells();
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        add_closed_cavity_bcs(ubc, pbc);
        auto controls = make_flow_controls();

        auto setup = make_buoyant_base();
        setup.boundary_conditions.clear();
        for (const char* name :
             {"left", "right", "bottom", "top", "front", "back"})
            setup.boundary_conditions.push_back(
                thermal_dirichlet(name, 300.0));

        Field<double, Location::CELL> U(n_cells, "U", "m/s", 3);
        U.fill(0.0);
        Field<double, Location::CELL> p(n_cells, "p", "Pa", 1);
        p.fill(0.0);
        Field<double, Location::CELL> temperature(
            n_cells, "temperature", "K", 1);
        temperature.fill(300.0);
        const auto result = run_production_buoyant_outer_loop(
            mesh, ubc, pbc, controls, setup, U, p, temperature);
        expect_converged_loop(result);
        EXPECT_TRUE(!result.s2s_ran);
        EXPECT_TRUE(result.flow.converged);
        EXPECT_TRUE(result.energy.converged);
        EXPECT_TRUE(result.energy.global_balance_normalized <= 1e-12);
        EXPECT_TRUE(result.energy.local_balance_normalized <= 1e-12);
        EXPECT_TRUE(result.energy.max_cell_residual_normalized <= 1e-12);
        // The response reproduces its input from the very first pass
        // (the input is already the fixed point), and the field stays
        // exactly at the reference temperature.
        EXPECT_TRUE(result.outer_history.front().max_source_input_lag_normalized
                    <= 1e-10);
        for (std::size_t c = 0; c < n_cells; ++c)
            EXPECT_NEAR(temperature(c), 300.0, 1e-9);
        // The carrier is at rest, like a plain solve without any source.
        EXPECT_TRUE(max_abs_velocity(U) <= 1e-9);
        EXPECT_TRUE(max_abs_scalar(p) <= 1e-9);

        // Plain carrier solve without any body force: the buoyant run
        // must agree with it (the source is exactly zero).
        Field<double, Location::CELL> U_ref(n_cells, "U_ref", "m/s", 3);
        U_ref.fill(0.0);
        Field<double, Location::CELL> p_ref(n_cells, "p_ref", "Pa", 1);
        p_ref.fill(0.0);
        const auto plain = solve_steady_incompressible(
            mesh, U_ref, p_ref, ubc, pbc, controls);
        EXPECT_TRUE(plain.converged);
        for (std::size_t d = 0; d < 3; ++d)
            for (std::size_t c = 0; c < n_cells; ++c)
                EXPECT_NEAR(U.component_data(d)[c],
                            U_ref.component_data(d)[c], 1e-9);
        for (std::size_t c = 0; c < n_cells; ++c)
            EXPECT_NEAR(p(c), p_ref(c), 1e-9);
    });

    // -----------------------------------------------------------------
    // Uniform equivalence: every wall at 310 K (T_ref = 300 K), so the
    // source is the constant a = -beta (310 - 300) g_vec. A constant
    // cellwise field must reproduce the uniform controls.body_force
    // assembly BIT FOR BIT, and the loop's converged carrier state must
    // agree with the equivalent uniform body-force solve.
    // -----------------------------------------------------------------
    run_case("production_buoyancy_uniform_body_force_equivalence", [] {
        const std::size_t nx = 4;
        const std::size_t ny = 4;
        auto mesh = box_mesh(nx, ny, 0.25, 0.25);
        const std::size_t n_cells = mesh.n_cells();
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        add_closed_cavity_bcs(ubc, pbc);
        auto controls = make_flow_controls();

        auto setup = make_buoyant_base();
        setup.gravity_vector = "0 -9.81 0";  // non-trivial pressure profile
        setup.boundary_conditions.clear();
        for (const char* name :
             {"left", "right", "bottom", "top", "front", "back"})
            setup.boundary_conditions.push_back(
                thermal_dirichlet(name, 310.0));

        const double beta = 2e-4;
        const double delta = 310.0 - 300.0;
        const Vec3 gravity = parse_gravity_vector("0 -9.81 0");

        // Bit-level oracle: a direct solve through the cellwise channel
        // with the constant source and a direct solve through the uniform
        // channel with the same doubles assemble the same system, so the
        // two solves agree exactly.
        Field<double, Location::CELL> uniform_temperature(
            n_cells, "T_uniform", "K", 1);
        uniform_temperature.fill(310.0);
        ProductionBuoyancyCoupling coupling;
        coupling.beta = beta;
        coupling.reference_temperature = 300.0;
        coupling.gravity = gravity;
        std::vector<Vec3> constant_acceleration;
        build_cell_buoyancy_acceleration(
            constant_acceleration, uniform_temperature, coupling);
        Field<double, Location::CELL> U_cellwise(
            n_cells, "U_cellwise", "m/s", 3);
        U_cellwise.fill(0.0);
        Field<double, Location::CELL> p_cellwise(
            n_cells, "p_cellwise", "Pa", 1);
        p_cellwise.fill(0.0);
        auto cellwise_controls = controls;
        cellwise_controls.cell_body_force = constant_acceleration;
        const auto cellwise_result = solve_steady_incompressible(
            mesh, U_cellwise, p_cellwise, ubc, pbc, cellwise_controls);
        EXPECT_TRUE(cellwise_result.converged);

        Field<double, Location::CELL> U_uniform(
            n_cells, "U_uniform", "m/s", 3);
        U_uniform.fill(0.0);
        Field<double, Location::CELL> p_uniform(
            n_cells, "p_uniform", "Pa", 1);
        p_uniform.fill(0.0);
        auto uniform_controls = controls;
        uniform_controls.body_force = Vec3{
            -beta * delta * gravity.x,
            -beta * delta * gravity.y,
            -beta * delta * gravity.z};
        const auto uniform_result = solve_steady_incompressible(
            mesh, U_uniform, p_uniform, ubc, pbc, uniform_controls);
        EXPECT_TRUE(uniform_result.converged);

        for (std::size_t d = 0; d < 3; ++d)
            for (std::size_t c = 0; c < n_cells; ++c)
                EXPECT_NEAR(U_cellwise.component_data(d)[c],
                            U_uniform.component_data(d)[c], 0.0);
        for (std::size_t c = 0; c < n_cells; ++c)
            EXPECT_NEAR(p_cellwise(c), p_uniform(c), 0.0);
        // The buoyancy source is real: the driven carrier is not at rest.
        EXPECT_TRUE(max_abs_velocity(U_cellwise) > 0.0 ||
                    max_abs_scalar(p_cellwise) > 0.0);

        // The loop converges to the coupled fixed point (the response of
        // uniformly heated walls is the uniform 310 K field) and its
        // final carrier state agrees with the equivalent uniform
        // body-force solve (both solve the same discrete saddle point,
        // from different paths; the linear solve tolerance bounds the
        // difference with orders of margin).
        Field<double, Location::CELL> U(n_cells, "U", "m/s", 3);
        U.fill(0.0);
        Field<double, Location::CELL> p(n_cells, "p", "Pa", 1);
        p.fill(0.0);
        Field<double, Location::CELL> temperature(
            n_cells, "temperature", "K", 1);
        temperature.fill(300.0);
        const auto result = run_production_buoyant_outer_loop(
            mesh, ubc, pbc, controls, setup, U, p, temperature);
        expect_converged_loop(result);
        for (std::size_t c = 0; c < n_cells; ++c)
            EXPECT_NEAR(temperature(c), 310.0, 1e-9);
        for (std::size_t d = 0; d < 3; ++d)
            for (std::size_t c = 0; c < n_cells; ++c)
                EXPECT_NEAR(
                    U.component_data(d)[c],
                    U_uniform.component_data(d)[c],
                    1e-8 * std::max(
                                1.0,
                                std::abs(U_uniform.component_data(d)[c])));
        for (std::size_t c = 0; c < n_cells; ++c)
            EXPECT_NEAR(p(c), p_uniform(c),
                        1e-8 * std::max(1.0, std::abs(p_uniform(c))));
    });

    // -----------------------------------------------------------------
    // Coupled fixed point: differentially heated closed cavity (left
    // 305 K, right 295 K, everything else adiabatic), gravity down. The
    // temperature gradient drives a buoyant circulation; the loop must
    // converge to the exact coupled fixed point with all four criteria,
    // and replaying one outer iteration by hand must reproduce the final
    // responses.
    // -----------------------------------------------------------------
    run_case("production_buoyancy_coupled_fixed_point", [] {
        const std::size_t nx = 4;
        const std::size_t ny = 4;
        auto mesh = box_mesh(nx, ny, 0.25, 0.25);
        const std::size_t n_cells = mesh.n_cells();
        const auto geometry = build_fv_geometry(mesh);
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        add_closed_cavity_bcs(ubc, pbc);
        auto controls = make_flow_controls();

        auto setup = make_buoyant_base();
        setup.gravity_vector = "0 -9.81 0";
        setup.materials.front().extra_properties[
            "thermal_expansion_coefficient"] = "2e-5";  // beta dT << 1
        setup.boundary_conditions.clear();
        setup.boundary_conditions.push_back(thermal_dirichlet("left", 305.0));
        setup.boundary_conditions.push_back(
            thermal_dirichlet("right", 295.0));

        Field<double, Location::CELL> U(n_cells, "U", "m/s", 3);
        U.fill(0.0);
        Field<double, Location::CELL> p(n_cells, "p", "Pa", 1);
        p.fill(0.0);
        Field<double, Location::CELL> temperature(
            n_cells, "temperature", "K", 1);
        temperature.fill(300.0);
        const auto result = run_production_buoyant_outer_loop(
            mesh, ubc, pbc, controls, setup, U, p, temperature);
        expect_converged_loop(result);
        EXPECT_TRUE(!result.s2s_ran);
        EXPECT_TRUE(result.energy.converged);
        EXPECT_TRUE(result.energy.global_balance_normalized <= 1e-12);
        EXPECT_TRUE(result.energy.local_balance_normalized <= 1e-12);
        EXPECT_TRUE(result.energy.max_cell_residual_normalized <= 1e-12);

        // The buoyancy source drives the flow: the carrier is not at
        // rest, while the same cavity without any source is.
        EXPECT_TRUE(max_abs_velocity(U) > 1e-6);
        {
            Field<double, Location::CELL> U_rest(n_cells, "U_rest", "m/s", 3);
            U_rest.fill(0.0);
            Field<double, Location::CELL> p_rest(n_cells, "p_rest", "Pa", 1);
            p_rest.fill(0.0);
            const auto rest = solve_steady_incompressible(
                mesh, U_rest, p_rest, ubc, pbc, controls);
            EXPECT_TRUE(rest.converged);
            EXPECT_TRUE(max_abs_velocity(U_rest) <= 1e-12);
        }

        // Bounded thermal response between the two Dirichlet walls.
        for (std::size_t c = 0; c < n_cells; ++c) {
            EXPECT_TRUE(temperature(c) >= 295.0 - 1e-9);
            EXPECT_TRUE(temperature(c) <= 305.0 + 1e-9);
        }

        // Fixed-point replay: rebuild the source from the FINAL response,
        // warm-start the carrier from the final state and run ONE energy
        // pass on the new mass flux. At the coupled fixed point the
        // responses reproduce themselves.
        const auto coupling = production_buoyancy_coupling(setup);
        std::vector<Vec3> replay_acceleration;
        build_cell_buoyancy_acceleration(
            replay_acceleration, temperature, coupling);
        Field<double, Location::CELL> U_replay(
            n_cells, "U_replay", "m/s", 3);
        U_replay = U;
        Field<double, Location::CELL> p_replay(
            n_cells, "p_replay", "Pa", 1);
        p_replay = p;
        auto replay_controls = controls;
        replay_controls.cell_body_force = replay_acceleration;
        replay_controls.initialization.mode = InitializationMode::Provided;
        const auto replay_flow = solve_steady_incompressible(
            mesh, U_replay, p_replay, ubc, pbc, replay_controls);
        EXPECT_TRUE(replay_flow.converged);
        Field<double, Location::CELL> T_replay(n_cells, "T_replay", "K", 1);
        T_replay = temperature;
        const auto replay_energy = run_production_energy_stage(
            mesh, geometry, setup, replay_flow.authoritative_mass_flux,
            T_replay);
        EXPECT_TRUE(replay_energy.converged);
        for (std::size_t c = 0; c < n_cells; ++c)
            EXPECT_NEAR(T_replay(c), temperature(c), 1e-6);
        for (std::size_t d = 0; d < 3; ++d)
            for (std::size_t c = 0; c < n_cells; ++c)
                EXPECT_NEAR(
                    U_replay.component_data(d)[c],
                    U.component_data(d)[c],
                    1e-8 * std::max(
                                1.0, std::abs(U.component_data(d)[c])));
        for (std::size_t c = 0; c < n_cells; ++c)
            EXPECT_NEAR(p_replay(c), p(c),
                        1e-8 * std::max(1.0, std::abs(p(c))));
    });

    // -----------------------------------------------------------------
    // Nested S2S: the stage-3 energy-radiation outer loop (fiche 11,
    // gates widened to accept the buoyant model) runs UNCHANGED as the
    // inner stage of the buoyancy loop. Same enclosure as the stage-3
    // VV (left 400 K / right 300 K, bottom/top radiative participants
    // with view factor 1, declared 50 W/m^2 on the bottom), now coupled
    // to a buoyant carrier.
    // -----------------------------------------------------------------
    run_case("production_buoyancy_s2s_nested_fixed_point", [] {
        const std::size_t nx = 4;
        const std::size_t ny = 4;
        auto mesh = box_mesh(nx, ny, 0.25, 0.25);
        const std::size_t n_cells = mesh.n_cells();
        const auto geometry = build_fv_geometry(mesh);
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        add_closed_cavity_bcs(ubc, pbc);
        auto controls = make_flow_controls();

        auto setup = make_buoyant_base();
        setup.energy_model = "buoyant";
        setup.radiation_model = "S2S";
        setup.gravity_vector = "0 -9.81 0";
        setup.ref_temperature = 350.0;
        setup.materials.front().extra_properties[
            "thermal_expansion_coefficient"] = "2e-5";
        setup.boundary_conditions.clear();
        setup.boundary_conditions.push_back(thermal_dirichlet("left", 400.0));
        setup.boundary_conditions.push_back(
            thermal_dirichlet("right", 300.0));

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
        Field<double, Location::CELL> U(n_cells, "U", "m/s", 3);
        U.fill(0.0);
        Field<double, Location::CELL> p(n_cells, "p", "Pa", 1);
        p.fill(0.0);
        Field<double, Location::CELL> temperature(
            n_cells, "temperature", "K", 1);
        temperature.fill(350.0);
        const auto result = run_production_buoyant_outer_loop(
            mesh, ubc, pbc, controls, setup, U, p, temperature);
        expect_converged_loop(result);
        EXPECT_TRUE(result.s2s_ran);
        EXPECT_TRUE(result.s2s.converged);
        EXPECT_TRUE(result.s2s.energy.converged);
        EXPECT_TRUE(result.energy.converged);
        EXPECT_TRUE(result.energy.global_balance_normalized <= 1e-12);
        EXPECT_TRUE(result.energy.local_balance_normalized <= 1e-12);
        EXPECT_TRUE(result.energy.max_cell_residual_normalized <= 1e-12);
        for (const auto& entry : result.outer_history)
            EXPECT_TRUE(entry.s2s_converged);

        // The radiosity machine reports its own green gates on the final
        // coupled state (fiche 11).
        const auto& radiation = result.s2s.radiation;
        EXPECT_TRUE(radiation.converged);
        EXPECT_TRUE(radiation.balance_normalized <= 1e-10);
        EXPECT_TRUE(radiation.max_surface_identity_normalized <= 1e-10);
        EXPECT_TRUE(radiation.surface_count == 2);
        const auto bi = surface_index(radiation, "bottom");
        const auto ti = surface_index(radiation, "top");

        // Independent wall temperatures of the FINAL field.
        const double tb =
            patch_wall_temperature(mesh, geometry, "bottom", temperature);
        const double tt =
            patch_wall_temperature(mesh, geometry, "top", temperature);
        EXPECT_NEAR(radiation.wall_temperatures[bi], tb, 1e-12);
        EXPECT_NEAR(radiation.wall_temperatures[ti], tt, 1e-12);

        // Hand Cramer oracle of the final radiative response, plus the
        // parallel-plate identity as an independent second formula.
        const double sigma = STEFAN_BOLTZMANN;
        const double eps_b = 0.8;
        const double eps_t = 0.9;
        const double eb = sigma * std::pow(tb, 4);
        const double et = sigma * std::pow(tt, 4);
        const double det = 1.0 - (1.0 - eps_b) * (1.0 - eps_t);
        const double jb = (eps_b * eb + (1.0 - eps_b) * eps_t * et) / det;
        const double jt = (eps_t * et + (1.0 - eps_t) * eps_b * eb) / det;
        const double qb = jb - jt;
        const double qt = jt - jb;
        const double q_parallel =
            sigma * (std::pow(tb, 4) - std::pow(tt, 4)) /
            (1.0 / eps_b + 1.0 / eps_t - 1.0);
        EXPECT_NEAR(qb, q_parallel, 1e-9);
        EXPECT_NEAR(radiation.net_flux[bi], qb, 1e-11);
        EXPECT_NEAR(radiation.net_flux[ti], qt, 1e-11);
        // Fiche 11 B3 sign convention + machine closure of the pair.
        EXPECT_TRUE(radiation.net_flux[bi] > 0.0);
        EXPECT_TRUE(radiation.net_flux[ti] < 0.0);
        EXPECT_NEAR(radiation.net_flux[bi] + radiation.net_flux[ti],
                    0.0, 1e-11);

        // Physical closure (areas: bottom = top = nx * dx = 1 m^2): the
        // radiative exchange is internal to the enclosure, so the walls
        // carry the declared 50 W/m^2 in total and the adiabatic
        // front/back patches carry ~0. Tolerance note: unlike the
        // frozen-flow stage-3 VV, the carrier now convects, so the
        // discrete conservation identity closes against a convective
        // flux scale orders of magnitude above the 50 W signal; 1e-3
        // bounds the linear-solve residual with two orders of margin.
        const double area = static_cast<double>(nx) * 0.25;
        EXPECT_NEAR(result.energy.patch_flux_out.at("bottom"),
                    -area * (50.0 - qb), 1e-3);
        EXPECT_NEAR(result.energy.patch_flux_out.at("top"),
                    area * qt, 1e-3);
        EXPECT_NEAR(result.energy.patch_flux_out.at("left") +
                        result.energy.patch_flux_out.at("right"),
                    50.0 * area, 1e-3);
        EXPECT_NEAR(result.energy.patch_flux_out.at("front"), 0.0, 1e-3);
        EXPECT_NEAR(result.energy.patch_flux_out.at("back"), 0.0, 1e-3);

        // The buoyant carrier is not at rest: the cavity circulates.
        EXPECT_TRUE(max_abs_velocity(U) > 1e-6);
    });

    // -----------------------------------------------------------------
    // Fail-closed catalogue of the loop itself.
    // -----------------------------------------------------------------
    run_case("production_buoyancy_outer_loop_fail_closed", [] {
        const std::size_t nx = 2;
        const std::size_t ny = 2;
        auto mesh = box_mesh(nx, ny, 0.25, 0.25);
        const std::size_t n_cells = mesh.n_cells();
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        add_closed_cavity_bcs(ubc, pbc);
        auto controls = make_flow_controls();

        // Running the buoyancy loop under another declared energy model
        // is a contract violation, not a silent fallback.
        {
            auto setup = make_buoyant_base();
            setup.energy_model = "conjugate";
            Field<double, Location::CELL> U(n_cells, "U", "m/s", 3);
            U.fill(0.0);
            Field<double, Location::CELL> p(n_cells, "p", "Pa", 1);
            p.fill(0.0);
            Field<double, Location::CELL> temperature(
                n_cells, "T", "K", 1);
            temperature.fill(300.0);
            EXPECT_THROW(
                run_production_buoyant_outer_loop(
                    mesh, ubc, pbc, controls, setup, U, p, temperature),
                std::invalid_argument);
        }

        // A temperature field that does not match the mesh.
        {
            auto setup = make_buoyant_base();
            Field<double, Location::CELL> U(n_cells, "U", "m/s", 3);
            U.fill(0.0);
            Field<double, Location::CELL> p(n_cells, "p", "Pa", 1);
            p.fill(0.0);
            Field<double, Location::CELL> wrong(
                n_cells + 1, "T", "K", 1);
            wrong.fill(300.0);
            EXPECT_THROW(
                run_production_buoyant_outer_loop(
                    mesh, ubc, pbc, controls, setup, U, p, wrong),
                std::invalid_argument);
        }

        // A buoyant case without any thermal Dirichlet patch: the energy
        // stage inside the loop rejects it (no silent isothermal run).
        {
            auto setup = make_buoyant_base();
            setup.boundary_conditions.clear();
            Field<double, Location::CELL> U(n_cells, "U", "m/s", 3);
            U.fill(0.0);
            Field<double, Location::CELL> p(n_cells, "p", "Pa", 1);
            p.fill(0.0);
            Field<double, Location::CELL> temperature(
                n_cells, "T", "K", 1);
            temperature.fill(300.0);
            EXPECT_THROW(
                run_production_buoyant_outer_loop(
                    mesh, ubc, pbc, controls, setup, U, p, temperature),
                std::invalid_argument);
        }
    });

    return run_all();
}
