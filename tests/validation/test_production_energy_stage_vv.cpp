
// T2 qualification: frozen-flow fluid energy transport through the
// production solver's authoritative conservative face mass flux
// (issues #880, #885 stage 2).
//
// Verified properties:
//   * Plug-flow carrier (slip walls, uniform inlet) converged by
//     solve_steady_incompressible: |u - u_in| < 1e-9 and a uniform
//     authoritative face mass flux, signed outward (inlet negative,
//     outlet positive, walls zero). The gauge anchor is the outlet cell:
//     the reference cell's continuity row is eliminated from the pressure
//     system, so anchoring on the inlet cell makes the plug unreachable
//     (see the comment in solve_plug_flow and issue #663).
//   * The frozen-flow energy stage (energy_model="conjugate", one-way
//     fluid -> energy coupling on the authoritative mass flux) reproduces
//     an independent hand-assembled discrete oracle exactly: the
//     tridiagonal upwind transport system is rebuilt mechanically from the
//     fiche-14 formulas (never through assemble_scalar_equation) and solved
//     with the Thomas algorithm; agreement is < 1e-9 per cell.
//   * Independent energy conservation (reconstruct_energy_balance, not the
//     assembled matrix) is a BLOCKING verdict in three normalized forms:
//     the global closure and the local per-cell balance, both normalized by
//     the gross audit power scale (boundary throughput + source, floor
//     1 W), plus the per-cell balance normalized by each cell's OWN gross
//     power scale (its face throughputs + its source, floor 1 W) so a
//     local defect cannot hide behind a large domain scale; all three
//     (< 1e-12) must pass for the stage to be converged; the net
//     boundary flux equals the source power; per-patch identities hold
//     (adiabatic walls exactly zero, imposed wall flux = q'' * area).
//   * Discrete maximum principle with Dirichlet walls; bulk temperature
//     rise q'''*V/(mdot*cp) as a physical sanity gate.
//   * Fail-closed: unsupported energy models (including "buoyant"), wrong
//     material counts, invalid properties (density, dynamic viscosity,
//     conductivity, specific heat - the carrier flow solves the same
//     declared fluid), non-numeric / partially numeric ("500abc") / empty /
//     non-finite volumetric power and heat flux, quotient overflow of
//     heat_flux / conductivity, heat_flux on a non-WALL patch, a patch
//     declaring both a Dirichlet temperature and heat_flux, non-positive
//     declared or initial temperatures, Dirichlet-free (singular) systems,
//     non-finite mass flux, mismatched field sizes, and a geometry that
//     does not match the mesh (sizes, non-finite centres or area vectors,
//     invalid patch indices on boundary faces) are all rejected.
#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/physics/production_energy_stage.h"
#include "common/test_harness.h"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

// ---------------------------------------------------------------------------
// 1D chain mesh: n cells of width dx along x, unit cross-section. Columns of
// 4 points at each x boundary; one x-face per boundary, four wall faces per
// cell. Patches: "left", "right", "walls" (mirrors test_thermal_zones_vv).
// ---------------------------------------------------------------------------
static Mesh chain_mesh(std::size_t n, double dx)
{
    Mesh m;
    m.points().resize(4 * (n + 1));
    for (std::size_t i = 0; i <= n; ++i) {
        const double x = static_cast<double>(i) * dx;
        m.points().set(4 * i + 0, x, 0.0, 0.0);
        m.points().set(4 * i + 1, x, 1.0, 0.0);
        m.points().set(4 * i + 2, x, 1.0, 1.0);
        m.points().set(4 * i + 3, x, 0.0, 1.0);
    }

    // x-normal faces: boundary i uses column i's four points.
    for (std::size_t i = 0; i <= n; ++i)
        m.faces().push_face({4 * i + 0, 4 * i + 1, 4 * i + 2, 4 * i + 3});
    // wall faces for cell i (columns i and i+1).
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t a = 4 * i;
        const std::size_t b = 4 * (i + 1);
        m.faces().push_face({a + 0, b + 0, b + 1, a + 1});  // y = 0
        m.faces().push_face({a + 3, a + 2, b + 2, b + 3});  // y = 1
        m.faces().push_face({a + 0, a + 3, b + 3, b + 0});  // z = 0
        m.faces().push_face({a + 1, b + 1, b + 2, a + 2});  // z = 1
    }

    const std::size_t nf = (n + 1) + 4 * n;
    m.ownership().resize(nf);
    // x-face 0: boundary of cell 0; x-face i (1..n-1): cells i-1 | i;
    // x-face n: boundary of cell n-1.
    m.ownership().set_owner(0, 0);
    m.ownership().set_neighbour(0, FaceOwnership::BOUNDARY);
    for (std::size_t i = 1; i < n; ++i) {
        m.ownership().set_owner(i, i - 1);
        m.ownership().set_neighbour(i, static_cast<int>(i));
    }
    m.ownership().set_owner(n, n - 1);
    m.ownership().set_neighbour(n, FaceOwnership::BOUNDARY);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            const std::size_t f = (n + 1) + 4 * i + j;
            m.ownership().set_owner(f, i);
            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        }
    }

    for (std::size_t i = 0; i < n; ++i)
        m.cells().push_cell({i, i + 1,
                             (n + 1) + 4 * i + 0, (n + 1) + 4 * i + 1,
                             (n + 1) + 4 * i + 2, (n + 1) + 4 * i + 3});

    Patch left;
    left.name = "left";
    left.type = PatchType::WALL;
    left.face_ids = {0};
    Patch right;
    right.name = "right";
    right.type = PatchType::WALL;
    right.face_ids = {n};
    Patch walls;
    walls.name = "walls";
    walls.type = PatchType::WALL;
    walls.face_ids.reserve(4 * n);
    for (std::size_t f = n + 1; f < nf; ++f) walls.face_ids.push_back(f);
    m.boundary().add_patch(left);
    m.boundary().add_patch(right);
    m.boundary().add_patch(walls);
    return m;
}

// ---------------------------------------------------------------------------
// Plug-flow carrier: slip walls (zero-gradient velocity), uniform inlet,
// zero-gradient outlet. Initialized from rest (U = 0) with the default
// SIMPLE under-relaxation, the carrier converges geometrically onto the
// uniform plug, which is the exact discrete solution of this configuration.
// ---------------------------------------------------------------------------
struct PlugFlow {
    Mesh mesh;
    FvGeometry geometry;
    Field<double, Location::FACE> mass_flux;
    bool converged = false;
    double velocity_error = 0.0;
};

static PlugFlow solve_plug_flow(std::size_t n, double dx, double rho,
                                double mu, double u_in)
{
    PlugFlow plug;
    plug.mesh = chain_mesh(n, dx);
    Field<double, Location::CELL> U(n, "U", "m/s", 3);
    U.fill(0.0);
    Field<double, Location::CELL> p(n, "p", "Pa", 1);
    p.fill(0.0);

    VelocityBoundaryConditions ubc;
    ubc["left"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,
                   {u_in, 0.0, 0.0}};
    ubc["right"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,
                    {0.0, 0.0, 0.0}};
    ubc["walls"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,
                    {0.0, 0.0, 0.0}};
    ScalarBoundaryConditions pbc;
    pbc["left"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    pbc["right"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    pbc["walls"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

    IncompressibleSolverControls controls;
    // Carrier initialization: from REST (U = 0), never at the exact plug
    // solution. Seeding the fixed point makes the first momentum residual
    // vanish, so the relative residual used by the convergence gates stays
    // at ~1 forever and the monitor reports STAGNATED after 10 iterations
    // (the documented failure mode of issue #663). The green pattern from
    // test_channel_qualification is rest + default under-relaxation
    // (alpha_u = 0.7, alpha_p = 0.3). No tolerance is changed.
    //
    // Gauge anchor: the reference cell's continuity row is ELIMINATED from
    // the pressure system (pure-Neumann gauge handling), so that cell's mass
    // defect can never be corrected. On this chain the inlet face is the
    // only boundary flux not driven by an interior pressure gradient (the
    // outlet is zero-gradient: its flux follows the interior solution), so
    // anchoring on cell 0 removes the only row that could balance the
    // injected inlet mass: the plug becomes an UNREACHABLE fixed point and
    // the pressure CG fails with MAX_ITER_COLLAPSE (iterations ~509,
    // residual ~1e-4 - the CI failure of PR #888). Anchoring on the OUTLET
    // cell n-1 keeps every continuity row that carries an uncompensated
    // boundary flux: the predictor balances the global mass and the run
    // converges geometrically at the rate (1 - alpha_p) per iteration,
    // from rest in ~90 iterations, well inside the 200-iteration budget.
    controls.initialization.mode = InitializationMode::Uniform;
    controls.initialization.uniform_velocity = {0.0, 0.0, 0.0};
    controls.coupling.alpha_u = 0.7;
    controls.coupling.alpha_p = 0.3;
    controls.algorithm = PressureVelocityAlgorithm::SIMPLE;
    controls.convection_scheme = ConvectionScheme::UPWIND;
    controls.use_bounded_convection = true;
    controls.density = rho;
    controls.kinematic_viscosity = mu / rho;
    controls.convergence.max_iterations = 200;
    controls.convergence.relative_tolerance = 1e-9;
    controls.convergence.continuity_tolerance = 1e-11;
    controls.linear_tolerance = 1e-11;
    controls.pressure_reference_cell = n - 1;

    const auto flow = solve_steady_incompressible(plug.mesh, U, p, ubc, pbc,
                                                   controls);
    plug.converged = flow.converged;
    plug.mass_flux = flow.authoritative_mass_flux;
    for (std::size_t c = 0; c < n; ++c) {
        plug.velocity_error =
            std::max(plug.velocity_error, std::abs(U(c, 0) - u_in));
        plug.velocity_error = std::max(plug.velocity_error, std::abs(U(c, 1)));
        plug.velocity_error = std::max(plug.velocity_error, std::abs(U(c, 2)));
    }
    plug.geometry = build_fv_geometry(plug.mesh);
    return plug;
}

// ---------------------------------------------------------------------------
// Case setup for the energy stage: conjugate model, exactly one fluid
// material, volumetric power in extra_properties, inlet Dirichlet
// temperature, adiabatic/heated/imposed-flux walls.
// ---------------------------------------------------------------------------
static cfdx::io::CaseSetup make_setup(
    double rho, double mu, double k, double cp,
    const std::string& volumetric_power,
    cfdx::io::BCValueType wall_value_type,
    double wall_temperature,
    const std::string& wall_heat_flux,
    double inlet_temperature)
{
    cfdx::io::CaseSetup setup;
    setup.energy_model = "conjugate";

    cfdx::io::MaterialSpec air;
    air.name = "air";
    air.density = rho;
    air.dynamic_viscosity = mu;
    air.thermal_conductivity = k;
    air.specific_heat = cp;
    if (!volumetric_power.empty())
        air.extra_properties["volumetric_power"] = volumetric_power;
    setup.materials.push_back(air);

    cfdx::io::BoundarySpec left;
    left.patch_name = "left";
    left.type = cfdx::io::BCType::INLET;
    left.value_type = cfdx::io::BCValueType::FIXED;
    left.velocity_magnitude = 1.5;
    left.temperature = inlet_temperature;
    setup.boundary_conditions.push_back(left);

    cfdx::io::BoundarySpec right;
    right.patch_name = "right";
    right.type = cfdx::io::BCType::OUTLET;
    right.value_type = cfdx::io::BCValueType::ZERO_GRADIENT;
    setup.boundary_conditions.push_back(right);

    cfdx::io::BoundarySpec walls;
    walls.patch_name = "walls";
    walls.type = cfdx::io::BCType::WALL;
    walls.value_type = wall_value_type;
    walls.wall_temperature = wall_temperature;
    if (!wall_heat_flux.empty())
        walls.raw_params["heat_flux"] = wall_heat_flux;
    setup.boundary_conditions.push_back(walls);
    return setup;
}

// ---------------------------------------------------------------------------
// Independent hand-assembled oracle for the frozen-flow energy system on the
// chain mesh. The discrete upwind transport operator is written out
// mechanically from the fiche-14 sections 3-4 (internal upwind conductances,
// half-cell Dirichlet, fixed-gradient and bounded zero-gradient faces, the
// bounded -div(F)*T diagonal correction and the per-cell source), never by
// calling the production assembly. Geometry (areas, distances, volumes and
// face -> patch mapping) is derived analytically from the chain
// construction. Solved with the Thomas algorithm.
// ---------------------------------------------------------------------------
static std::vector<double> chain_energy_oracle(
    const Mesh& mesh,
    const FvGeometry& g,
    const Field<double, Location::FACE>& mass_flux,
    const ScalarBoundaryConditions& bcs,
    double k,
    double cp,
    double volumetric_power,
    double dx)
{
    const std::size_t nc = mesh.n_cells();
    const std::size_t n = nc;
    std::vector<double> diag(nc, 0.0);
    std::vector<double> sub(nc, 0.0);
    std::vector<double> sup(nc, 0.0);
    std::vector<double> rhs(nc, 0.0);
    std::vector<double> div(nc, 0.0);

    const auto& own = mesh.ownership();
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const std::size_t o = own.owner(f);
        const auto nraw = own.neighbour(f);
        const double F = mass_flux(f) * cp;  // enthalpy flux, W/K
        if (nraw >= 0) {
            // internal x-face: unit area, centre-to-centre distance dx
            const std::size_t nn = static_cast<std::size_t>(nraw);
            const double D = k * 1.0 / dx;
            const double a_on = D + std::max(F, 0.0);
            const double a_no = D + std::max(-F, 0.0);
            diag[o] += a_on;
            diag[nn] += a_no;
            sup[o] -= a_no;
            sub[nn] -= a_on;
            div[o] += F;
            div[nn] -= F;
        } else {
            const std::size_t patch = g.face_patch[f];
            const std::string& name = mesh.boundary().patch(patch).name;
            ScalarBoundaryCondition bc;  // defaults to zero-gradient
            const auto it = bcs.find(name);
            if (it != bcs.end()) bc = it->second;
            // wall faces come after the n+1 x-faces
            const bool lateral = f > n;
            const double area = lateral ? dx : 1.0;
            const double d_b = lateral ? 0.5 : 0.5 * dx;
            if (bc.type == ScalarBoundaryType::FIXED_VALUE) {
                const double D_b = k * area / d_b;
                diag[o] += D_b + std::max(F, 0.0);
                rhs[o] += (D_b + std::max(-F, 0.0)) * bc.value;
                div[o] += F;
            } else if (bc.type == ScalarBoundaryType::FIXED_GRADIENT) {
                diag[o] += F;
                rhs[o] -= F * bc.gradient * d_b;
                rhs[o] += k * area * bc.gradient;
                div[o] += F;
            } else {
                // bounded zero-gradient: convection only
                diag[o] += F;
                div[o] += F;
            }
        }
    }
    for (std::size_t c = 0; c < nc; ++c) {
        diag[c] -= div[c];
        rhs[c] += volumetric_power * dx;  // cell volume = dx
        if (!(diag[c] > 0.0) || !std::isfinite(diag[c]))
            throw std::runtime_error("chain_energy_oracle: non-positive diagonal");
    }

    // Thomas algorithm on the tridiagonal system.
    std::vector<double> cprime(nc, 0.0);
    std::vector<double> dprime(nc, 0.0);
    cprime[0] = sup[0] / diag[0];
    dprime[0] = rhs[0] / diag[0];
    for (std::size_t i = 1; i < nc; ++i) {
        const double m = diag[i] - sub[i] * cprime[i - 1];
        cprime[i] = sup[i] / m;
        dprime[i] = (rhs[i] - sub[i] * dprime[i - 1]) / m;
    }
    std::vector<double> t(nc, 0.0);
    t[nc - 1] = dprime[nc - 1];
    for (std::size_t i = nc - 1; i-- > 0;)
        t[i] = dprime[i] - cprime[i] * t[i + 1];
    return t;
}

int main()
{
    const std::size_t n = 16;
    const double dx = 0.25;
    const double rho = 1.2;
    const double mu = 1.8e-5;
    const double k = 0.026;
    const double cp = 1005.0;
    const double u_in = 1.5;
    const double t_in = 293.15;
    const double mdot = rho * u_in;  // kg/s through the unit cross-section
    const double length = static_cast<double>(n) * dx;

    run_case("production_energy_plug_flow_volumetric_power_oracle", [&] {
        const double q = 5000.0;  // W/m^3
        auto plug = solve_plug_flow(n, dx, rho, mu, u_in);
        EXPECT_TRUE(plug.converged);
        EXPECT_TRUE(plug.velocity_error < 1e-9);
        EXPECT_TRUE(plug.mass_flux(0) < 0.0);
        EXPECT_TRUE(plug.mass_flux(n) > 0.0);
        for (std::size_t f = 0; f <= n; ++f)
            EXPECT_NEAR(plug.mass_flux(f),
                        (f == 0 ? -1.0 : 1.0) * rho * u_in, 1e-9);
        for (std::size_t f = n + 1; f < plug.mesh.n_faces(); ++f)
            EXPECT_NEAR(plug.mass_flux(f), 0.0, 1e-12);

        const auto setup = make_setup(rho, mu, k, cp, "5000",
                                      cfdx::io::BCValueType::WALL_NO_SLIP,
                                      0.0, "", t_in);
        EXPECT_TRUE(production_energy_requested(setup));
        const auto energy_case = build_production_energy_case(setup);
        EXPECT_NEAR(energy_case.volumetric_power, q, 1e-12);
        EXPECT_TRUE(energy_case.controls.dt == 0.0);
        EXPECT_NEAR(energy_case.controls.relaxation, 1.0, 0.0);
        EXPECT_NEAR(energy_case.controls.linear_tolerance, 1e-14, 0.0);

        ScalarBoundaryConditions bcs;
        bcs["left"] = {ScalarBoundaryType::FIXED_VALUE, t_in, 0.0};
        bcs["right"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        bcs["walls"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        const auto stage_bcs =
            build_production_temperature_bcs(plug.mesh, setup, energy_case);
        EXPECT_TRUE(stage_bcs.at("left").type == ScalarBoundaryType::FIXED_VALUE);
        EXPECT_NEAR(stage_bcs.at("left").value, t_in, 0.0);
        EXPECT_TRUE(stage_bcs.at("right").type == ScalarBoundaryType::ZERO_GRADIENT);
        EXPECT_TRUE(stage_bcs.at("walls").type == ScalarBoundaryType::ZERO_GRADIENT);

        Field<double, Location::CELL> temperature(n, "temperature", "K", 1);
        temperature.fill(t_in);
        const auto stage = run_production_energy_stage(
            plug.mesh, plug.geometry, setup, plug.mass_flux, temperature);
        EXPECT_TRUE(stage.converged);

        const auto oracle = chain_energy_oracle(
            plug.mesh, plug.geometry, plug.mass_flux, bcs, k, cp, q, dx);
        double oracle_error = 0.0;
        for (std::size_t c = 0; c < n; ++c)
            oracle_error = std::max(oracle_error,
                                    std::abs(temperature(c) - oracle[c]));
        EXPECT_TRUE(oracle_error < 1e-9);

        EXPECT_TRUE(stage.nonfinite_temperature == 0);
        EXPECT_TRUE(stage.linear_solver_converged);
        EXPECT_TRUE(stage.global_energy_balance_passed);
        EXPECT_TRUE(stage.local_energy_balance_passed);
        EXPECT_TRUE(stage.nonfinite_faces == 0);
        EXPECT_TRUE(stage.nonfinite_cells == 0);
        // Blocking gates, normalized by the gross audit power scale
        // (boundary throughput + source, floor 1 W): the global closure
        // and the strictly stronger local one.
        EXPECT_TRUE(stage.audit_power_scale >= 1.0);
        EXPECT_TRUE(stage.global_balance_normalized < 1e-12);
        EXPECT_TRUE(stage.local_balance_normalized < 1e-12);
        // Blocking local gate in its per-cell normalization: each cell's
        // residual measured against its own gross power scale (its face
        // throughputs + its source, floor 1 W).
        EXPECT_TRUE(stage.max_cell_residual_normalized < 1e-12);
        EXPECT_TRUE(std::abs(stage.independent_conservation_residual) <
                    1e-9 * stage.audit_power_scale);
        const double total_power = q * length;  // 20000 W
        EXPECT_NEAR(stage.total_source_power, total_power, 1e-9);
        EXPECT_NEAR(stage.boundary_flux_out, total_power, 1e-6);
        EXPECT_NEAR(stage.patch_flux_out.at("walls"), 0.0, 1e-9);

        // Discrete maximum principle: no temperature below the inlet value.
        EXPECT_TRUE(stage.min_temperature >= t_in - 1e-9);
        // Physical sanity: bulk rise q'''V/(mdot*cp) = 11.06 K.
        const double bulk_rise = total_power / (mdot * cp);
        EXPECT_TRUE(stage.max_temperature <= t_in + bulk_rise + 1.0);
        EXPECT_TRUE(stage.max_temperature >= t_in + 0.5 * bulk_rise);
        std::cout << "case1 oracle_error=" << oracle_error
                  << " bulk_rise=" << bulk_rise
                  << " max_T=" << stage.max_temperature << "\n";
    });

    run_case("production_energy_heated_walls_bounded", [&] {
        const double t_w = 350.0;
        auto plug = solve_plug_flow(n, dx, rho, mu, u_in);
        EXPECT_TRUE(plug.converged);
        EXPECT_TRUE(plug.velocity_error < 1e-9);

        const auto setup = make_setup(rho, mu, k, cp, "",
                                      cfdx::io::BCValueType::WALL_THERMAL,
                                      t_w, "", t_in);
        const auto energy_case = build_production_energy_case(setup);
        EXPECT_NEAR(energy_case.volumetric_power, 0.0, 0.0);

        ScalarBoundaryConditions bcs;
        bcs["left"] = {ScalarBoundaryType::FIXED_VALUE, t_in, 0.0};
        bcs["right"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        bcs["walls"] = {ScalarBoundaryType::FIXED_VALUE, t_w, 0.0};
        const auto stage_bcs =
            build_production_temperature_bcs(plug.mesh, setup, energy_case);
        EXPECT_TRUE(stage_bcs.at("walls").type == ScalarBoundaryType::FIXED_VALUE);
        EXPECT_NEAR(stage_bcs.at("walls").value, t_w, 0.0);

        Field<double, Location::CELL> temperature(n, "temperature", "K", 1);
        temperature.fill(t_in);
        const auto stage = run_production_energy_stage(
            plug.mesh, plug.geometry, setup, plug.mass_flux, temperature);
        EXPECT_TRUE(stage.converged);

        const auto oracle = chain_energy_oracle(
            plug.mesh, plug.geometry, plug.mass_flux, bcs, k, cp, 0.0, dx);
        double oracle_error = 0.0;
        for (std::size_t c = 0; c < n; ++c)
            oracle_error = std::max(oracle_error,
                                    std::abs(temperature(c) - oracle[c]));
        EXPECT_TRUE(oracle_error < 1e-9);

        EXPECT_TRUE(stage.linear_solver_converged);
        EXPECT_TRUE(stage.global_energy_balance_passed);
        EXPECT_TRUE(stage.local_energy_balance_passed);
        EXPECT_TRUE(stage.global_balance_normalized < 1e-12);
        EXPECT_TRUE(stage.local_balance_normalized < 1e-12);
        EXPECT_TRUE(stage.max_cell_residual_normalized < 1e-12);
        // No source: the net boundary flux must close at zero.
        EXPECT_NEAR(stage.boundary_flux_out, 0.0, 1e-6);
        // Discrete maximum principle: bounded by the Dirichlet values.
        EXPECT_TRUE(stage.min_temperature >= t_in - 1e-9);
        EXPECT_TRUE(stage.max_temperature <= t_w + 1e-9);
        EXPECT_TRUE(stage.min_temperature <= stage.max_temperature + 1e-9);
        // Heated walls absorb nothing: heat enters through the walls and
        // leaves through the outlet; the inlet enthalpy enters.
        EXPECT_TRUE(stage.patch_flux_out.at("walls") < 0.0);
        EXPECT_TRUE(stage.patch_flux_out.at("right") > 0.0);
        EXPECT_TRUE(stage.patch_flux_out.at("left") < 0.0);
        std::cout << "case2 oracle_error=" << oracle_error
                  << " min_T=" << stage.min_temperature
                  << " max_T=" << stage.max_temperature << "\n";
    });

    run_case("production_energy_wall_heat_flux_identity", [&] {
        const double wall_flux = 500.0;  // W/m^2, positive into the domain
        auto plug = solve_plug_flow(n, dx, rho, mu, u_in);
        EXPECT_TRUE(plug.converged);
        EXPECT_TRUE(plug.velocity_error < 1e-9);

        const auto setup = make_setup(rho, mu, k, cp, "",
                                      cfdx::io::BCValueType::WALL_NO_SLIP,
                                      0.0, "500", t_in);
        const auto energy_case = build_production_energy_case(setup);
        ScalarBoundaryConditions bcs;
        bcs["left"] = {ScalarBoundaryType::FIXED_VALUE, t_in, 0.0};
        bcs["right"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        bcs["walls"] = {ScalarBoundaryType::FIXED_GRADIENT, 0.0,
                        wall_flux / k};
        const auto stage_bcs =
            build_production_temperature_bcs(plug.mesh, setup, energy_case);
        EXPECT_TRUE(stage_bcs.at("walls").type == ScalarBoundaryType::FIXED_GRADIENT);
        EXPECT_NEAR(stage_bcs.at("walls").gradient, wall_flux / k, 1e-15);

        Field<double, Location::CELL> temperature(n, "temperature", "K", 1);
        temperature.fill(t_in);
        const auto stage = run_production_energy_stage(
            plug.mesh, plug.geometry, setup, plug.mass_flux, temperature);
        EXPECT_TRUE(stage.converged);

        const auto oracle = chain_energy_oracle(
            plug.mesh, plug.geometry, plug.mass_flux, bcs, k, cp, 0.0, dx);
        double oracle_error = 0.0;
        for (std::size_t c = 0; c < n; ++c)
            oracle_error = std::max(oracle_error,
                                    std::abs(temperature(c) - oracle[c]));
        EXPECT_TRUE(oracle_error < 1e-9);

        EXPECT_TRUE(stage.linear_solver_converged);
        EXPECT_TRUE(stage.global_energy_balance_passed);
        EXPECT_TRUE(stage.local_energy_balance_passed);
        EXPECT_TRUE(stage.global_balance_normalized < 1e-12);
        EXPECT_TRUE(stage.local_balance_normalized < 1e-12);
        EXPECT_TRUE(stage.max_cell_residual_normalized < 1e-12);
        // Wall power = q'' * (4 lateral faces * L * unit depth) = 8000 W in.
        const double wall_power = wall_flux * 4.0 * length;
        EXPECT_NEAR(stage.patch_flux_out.at("walls"), -wall_power, 1e-6);
        // Steady with no volumetric source: net boundary flux closes at 0.
        EXPECT_NEAR(stage.boundary_flux_out, 0.0, 1e-6);
        // Physical sanity: bulk rise q''*4L/(mdot*cp) = 4.42 K.
        const double bulk_rise = wall_power / (mdot * cp);
        EXPECT_TRUE(stage.min_temperature >= t_in - 1e-9);
        EXPECT_TRUE(stage.max_temperature <= t_in + bulk_rise + 1.0);
        EXPECT_TRUE(stage.max_temperature >= t_in + 0.5 * bulk_rise);
        std::cout << "case3 oracle_error=" << oracle_error
                  << " wall_power=" << wall_power
                  << " max_T=" << stage.max_temperature << "\n";
    });

    run_case("production_energy_fail_closed_declarations", [&] {
        const std::size_t m = 6;
        Mesh mesh = chain_mesh(m, dx);
        auto geometry = build_fv_geometry(mesh);
        Field<double, Location::FACE> mass_flux(mesh.n_faces(), "mdot",
                                                "kg/s", 1);
        for (std::size_t f = 0; f <= m; ++f)
            mass_flux(f) = (f == 0 ? -1.0 : 1.0) * rho * u_in;
        // wall faces stay at zero: slip walls carry no mass flux

        const auto base = make_setup(rho, mu, k, cp, "5000",
                                     cfdx::io::BCValueType::WALL_NO_SLIP,
                                     0.0, "", t_in);

        // energy_model gate: accepted declarations vs fail-closed rejects.
        {
            auto s = base;
            EXPECT_TRUE(production_energy_requested(s));
            for (const char* model : {"", "none", "isothermal",
                                      "solid_conduction", "compressible"}) {
                s.energy_model = model;
                EXPECT_TRUE(!production_energy_requested(s));
                validate_production_energy_model(s);
            }
            s.energy_model = "conjugate";
            EXPECT_TRUE(production_energy_requested(s));
            validate_production_energy_model(s);
        }
        {
            auto s = base;
            s.energy_model = "buoyant";
            EXPECT_THROW(validate_production_energy_model(s),
                         std::invalid_argument);
            s.energy_model = "not_a_model";
            EXPECT_THROW(validate_production_energy_model(s),
                         std::invalid_argument);
        }

        // Exactly one fluid material.
        {
            auto s = base;
            s.materials.clear();
            EXPECT_THROW(build_production_energy_case(s),
                         std::invalid_argument);
        }
        {
            auto s = base;
            auto extra = s.materials.front();
            extra.name = "solid";
            s.materials.push_back(extra);
            EXPECT_THROW(build_production_energy_case(s),
                         std::invalid_argument);
        }

        // Finite positive physical properties.
        {
            auto s = base;
            s.materials[0].specific_heat = 0.0;
            EXPECT_THROW(build_production_energy_case(s),
                         std::invalid_argument);
        }
        {
            auto s = base;
            s.materials[0].thermal_conductivity = 0.0;
            EXPECT_THROW(build_production_energy_case(s),
                         std::invalid_argument);
        }
        {
            auto s = base;
            s.materials[0].density = 0.0;
            EXPECT_THROW(build_production_energy_case(s),
                         std::invalid_argument);
        }
        {
            auto s = base;
            s.materials[0].density =
                std::numeric_limits<double>::quiet_NaN();
            EXPECT_THROW(build_production_energy_case(s),
                         std::invalid_argument);
        }
        // The carrier flow solves the same declared fluid: an incomplete
        // material (missing / invalid dynamic viscosity) is rejected whole
        // instead of letting the flow run with silently substituted
        // defaults.
        {
            auto s = base;
            s.materials[0].dynamic_viscosity = 0.0;
            EXPECT_THROW(build_production_energy_case(s),
                         std::invalid_argument);
        }
        {
            auto s = base;
            s.materials[0].dynamic_viscosity =
                std::numeric_limits<double>::quiet_NaN();
            EXPECT_THROW(build_production_energy_case(s),
                         std::invalid_argument);
        }

        // Volumetric power must be numeric and finite.
        {
            auto s = base;
            s.materials[0].extra_properties["volumetric_power"] = "nan";
            EXPECT_THROW(build_production_energy_case(s),
                         std::invalid_argument);
        }
        {
            auto s = base;
            s.materials[0].extra_properties["volumetric_power"] = "abc";
            EXPECT_THROW(build_production_energy_case(s),
                         std::invalid_argument);
        }
        // Strict numeric parsing: std::stod alone accepts a numeric prefix
        // ("500abc" silently becomes 500); the whole text must parse as one
        // finite number, and an explicitly declared empty value is a
        // malformed declaration, not a silent default.
        {
            auto s = base;
            s.materials[0].extra_properties["volumetric_power"] = "500abc";
            EXPECT_THROW(build_production_energy_case(s),
                         std::invalid_argument);
        }
        {
            auto s = base;
            s.materials[0].extra_properties["volumetric_power"] = "";
            EXPECT_THROW(build_production_energy_case(s),
                         std::invalid_argument);
        }
        {
            auto s = base;
            s.materials[0].extra_properties["volumetric_power"] = "inf";
            EXPECT_THROW(build_production_energy_case(s),
                         std::invalid_argument);
        }

        // Wall heat flux must be numeric and finite.
        {
            auto s = base;
            s.boundary_conditions[2].raw_params["heat_flux"] = "abc";
            const auto energy_case = build_production_energy_case(s);
            EXPECT_THROW(build_production_temperature_bcs(mesh, s, energy_case),
                         std::invalid_argument);
        }
        {
            auto s = base;
            s.boundary_conditions[2].raw_params["heat_flux"] = "500abc";
            const auto energy_case = build_production_energy_case(s);
            EXPECT_THROW(build_production_temperature_bcs(mesh, s, energy_case),
                         std::invalid_argument);
        }
        {
            auto s = base;
            s.boundary_conditions[2].raw_params["heat_flux"] =
                "1e308";  // /k overflows despite both operands being finite
            const auto energy_case = build_production_energy_case(s);
            EXPECT_THROW(build_production_temperature_bcs(mesh, s, energy_case),
                         std::invalid_argument);
        }
        // heat_flux is a WALL boundary condition (fiche 14): declaring it
        // on a flow boundary patch is a case-contract violation, rejected
        // instead of silently imposing a temperature gradient on an
        // inlet/outlet.
        {
            auto s = base;
            s.boundary_conditions[0].raw_params["heat_flux"] = "100";
            const auto energy_case = build_production_energy_case(s);
            EXPECT_THROW(build_production_temperature_bcs(mesh, s, energy_case),
                         std::invalid_argument);
        }
        {
            auto s = base;
            s.boundary_conditions[1].raw_params["heat_flux"] = "100";
            const auto energy_case = build_production_energy_case(s);
            EXPECT_THROW(build_production_temperature_bcs(mesh, s, energy_case),
                         std::invalid_argument);
        }
        // A patch declaring both a retained Dirichlet temperature and a
        // heat_flux is contradictory: exactly one thermal boundary
        // condition per patch, no silent priority.
        {
            auto s = base;
            s.boundary_conditions[2].value_type =
                cfdx::io::BCValueType::WALL_THERMAL;
            s.boundary_conditions[2].wall_temperature = 350.0;
            s.boundary_conditions[2].raw_params["heat_flux"] = "100";
            const auto energy_case = build_production_energy_case(s);
            EXPECT_THROW(build_production_temperature_bcs(mesh, s, energy_case),
                         std::invalid_argument);
        }
        {
            auto s = base;
            s.boundary_conditions[0].raw_params["heat_flux"] = "100";
            // inlet temperature t_in is already declared on patch 0
            const auto energy_case = build_production_energy_case(s);
            EXPECT_THROW(build_production_temperature_bcs(mesh, s, energy_case),
                         std::invalid_argument);
        }

        // Declared temperatures are absolute (kelvin): finite and > 0.
        {
            auto s = base;
            s.boundary_conditions[2].value_type =
                cfdx::io::BCValueType::WALL_THERMAL;
            s.boundary_conditions[2].wall_temperature = -5.0;
            const auto energy_case = build_production_energy_case(s);
            EXPECT_THROW(build_production_temperature_bcs(mesh, s, energy_case),
                         std::invalid_argument);
        }
        {
            auto s = base;
            s.boundary_conditions[0].temperature = -1.0;
            const auto energy_case = build_production_energy_case(s);
            EXPECT_THROW(build_production_temperature_bcs(mesh, s, energy_case),
                         std::invalid_argument);
        }

        // At least one Dirichlet patch: all-Neumann is singular.
        {
            auto s = base;
            s.boundary_conditions[0].temperature = 0.0;
            const auto energy_case = build_production_energy_case(s);
            EXPECT_THROW(build_production_temperature_bcs(mesh, s, energy_case),
                         std::invalid_argument);
        }
        {
            // WALL_THERMAL with wall_temperature = 0 means "not declared".
            auto s = base;
            s.boundary_conditions[0].temperature = 0.0;
            s.boundary_conditions[2].value_type =
                cfdx::io::BCValueType::WALL_THERMAL;
            s.boundary_conditions[2].wall_temperature = 0.0;
            const auto energy_case = build_production_energy_case(s);
            EXPECT_THROW(build_production_temperature_bcs(mesh, s, energy_case),
                         std::invalid_argument);
        }

        // The stage itself is fail-closed.
        {
            auto s = base;
            s.energy_model = "none";
            Field<double, Location::CELL> t(m, "temperature", "K", 1);
            t.fill(t_in);
            EXPECT_THROW(run_production_energy_stage(mesh, geometry, s,
                                                     mass_flux, t),
                         std::invalid_argument);
        }
        {
            auto bad_flux = mass_flux;
            bad_flux(0) = std::numeric_limits<double>::quiet_NaN();
            Field<double, Location::CELL> t(m, "temperature", "K", 1);
            t.fill(t_in);
            EXPECT_THROW(run_production_energy_stage(mesh, geometry, base,
                                                     bad_flux, t),
                         std::invalid_argument);
        }
        {
            Field<double, Location::CELL> t(m - 1, "temperature", "K", 1);
            t.fill(t_in);
            EXPECT_THROW(run_production_energy_stage(mesh, geometry, base,
                                                     mass_flux, t),
                         std::invalid_argument);
        }
        {
            Field<double, Location::CELL> t(m, "temperature", "K", 1);
            t.fill(t_in);
            t(0) = std::numeric_limits<double>::quiet_NaN();
            EXPECT_THROW(run_production_energy_stage(mesh, geometry, base,
                                                     mass_flux, t),
                         std::invalid_argument);
        }
        {
            // Absolute temperature contract at the kernel entry point too.
            Field<double, Location::CELL> t(m, "temperature", "K", 1);
            t.fill(t_in);
            t(0) = -1.0;
            EXPECT_THROW(run_production_energy_stage(mesh, geometry, base,
                                                     mass_flux, t),
                         std::invalid_argument);
        }
        {
            // Entry-point geometry contract: the kernel is a public header
            // API and must not trust the caller to hand over a geometry
            // matching the mesh.
            FvGeometry foreign_geometry;
            Field<double, Location::CELL> t(m, "temperature", "K", 1);
            t.fill(t_in);
            EXPECT_THROW(run_production_energy_stage(mesh, foreign_geometry,
                                                     base, mass_flux, t),
                         std::invalid_argument);
        }
        {
            // Non-finite cell centres are rejected, not silently
            // assembled: a degenerate geometry can produce finite but
            // physically wrong results.
            auto bad_geometry = build_fv_geometry(mesh);
            bad_geometry.cell_centres[0] = Vec3{
                std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0};
            Field<double, Location::CELL> t(m, "temperature", "K", 1);
            t.fill(t_in);
            EXPECT_THROW(run_production_energy_stage(mesh, bad_geometry,
                                                     base, mass_flux, t),
                         std::invalid_argument);
        }
        {
            // Non-finite face area vectors are rejected for the same
            // reason: the reconstruction and audits index them.
            auto bad_geometry = build_fv_geometry(mesh);
            bad_geometry.face_area_vectors[0] = Vec3{
                std::numeric_limits<double>::infinity(), 0.0, 0.0};
            Field<double, Location::CELL> t(m, "temperature", "K", 1);
            t.fill(t_in);
            EXPECT_THROW(run_production_energy_stage(mesh, bad_geometry,
                                                     base, mass_flux, t),
                         std::invalid_argument);
        }
        {
            // A boundary face must map to a declared patch: an invalid
            // patch index would otherwise be silently skipped by the
            // per-patch flux audit and its flux would miss from the
            // balance.
            auto bad_geometry = build_fv_geometry(mesh);
            bad_geometry.face_patch[0] =
                mesh.boundary().n_patches() + 7;  // face 0 is the inlet
            Field<double, Location::CELL> t(m, "temperature", "K", 1);
            t.fill(t_in);
            EXPECT_THROW(run_production_energy_stage(mesh, bad_geometry,
                                                     base, mass_flux, t),
                         std::invalid_argument);
        }
    });

    return run_all();
}

