
// Production S2S radiation coupled to the energy stage (capability #885
// stage 2): verification and validation on a frozen flow, one
// delayed-Picard pass.
//
// Verified properties:
//   * The production S2S stage reproduces an independent hand oracle: a
//     two-surface parallel-plate enclosure (bottom / top of a 4x4x1 hexa
//     box, view factors 1 both ways) whose radiosity system is solved by
//     hand with Cramer's rule, and whose net flux is checked against the
//     parallel-plate identity q = sigma (T_b^4 - T_t^4) / (1/eps_b +
//     1/eps_t - 1). Wall temperatures are recomputed INDEPENDENTLY by the
//     test (|Sf|-weighted average of adjacent cell temperatures over the
//     patch faces) from the pass-1 production energy field.
//   * Blocking machine gates: the global enclosure balance closes
//     (balance_passed, balance_normalized <= 1e-10) and the per-surface
//     identity holds (surface_identity_passed).
//   * Sign convention (fiche 11 B3): the hot surface loses (q > 0), the
//     cold surface gains (q < 0), q_b + q_t closes to machine zero.
//   * Re-injection contract: the injected heat fluxes -q_i (positive into
//     the domain, same convention as declared heat_flux) are replayed
//     through build_production_temperature_bcs(mesh, setup, energy_case,
//     extra): the bottom gradient is the SUM of the declared heat_flux and
//     the radiative contribution, the top gradient is the radiative
//     contribution alone; replaying pass 1 (fresh field) then pass 2 (with
//     the oracle extra map) reproduces the stage's final temperature field
//     cell by cell.
//   * The stage's final energy result (pass 2) is converged with all
//     independent blocking audits < 1e-12 (global, gross-scale local and
//     per-cell-normalized balances).
//   * Physical closure: bottom patch flux = -A_b (50 - q_b), top patch
//     flux = A_t q_t, Q_left + Q_right = 50 * A_b (steady global balance,
//     the radiative exchange is internal to the enclosure), adiabatic
//     front/back ~ 0.
//   * max_wall_temperature_lag is a finite, non-negative DIAGNOSTIC (one
//     delayed-Picard pass, no claim on the coupled fixed point).
//   * Fail-closed: S2S without the conjugate energy model, no participant,
//     a single participant, emissivity non-numeric / partially numeric /
//     zero / negative / above one / non-finite / empty, emissivity on a
//     non-WALL spec, a participant with a retained Dirichlet temperature,
//     view factors toward the patch itself, toward a non-participant,
//     declared without emissivity, non-reciprocal / non-closed enclosures,
//     out-of-range or non-numeric view factors, extra wall heat flux on a
//     Dirichlet patch / unknown patch / EMPTY patch / non-finite value,
//     non-finite or non-positive temperature fields, DO / solar
//     declarations, and running the S2S stage without a declared S2S
//     radiation model are all rejected.
#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/physics/production_energy_stage.h"
#include "cfdx/physics/production_radiation_contract.h"
#include "cfdx/physics/production_radiation_s2s_stage.h"
#include "common/test_harness.h"

#include <cmath>
#include <cstddef>
#include <map>
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
// The front patch type is parameterized for the fail-closed EMPTY variant.
// ---------------------------------------------------------------------------
static Mesh box_mesh(std::size_t nx, std::size_t ny, double dx, double dy,
                     PatchType front_type = PatchType::WALL)
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
        for (std::size_t j = 0; j < ny; ++j)
            m.faces().push_face({pid(i, j, 0), pid(i, j, 1),
                                 pid(i, j + 1, 1), pid(i, j + 1, 0)});
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j <= ny; ++j)
            m.faces().push_face({pid(i, j, 0), pid(i + 1, j, 0),
                                 pid(i + 1, j, 1), pid(i, j, 1)});
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t l = 0; l < 2; ++l)
                m.faces().push_face({pid(i, j, l), pid(i + 1, j, l),
                                     pid(i + 1, j + 1, l), pid(i, j + 1, l)});

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
    front.type = front_type;
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
// Case: conjugate energy + S2S radiation. Dirichlet walls left (400 K) and
// right (300 K); radiative participants bottom (eps 0.8, declared heat
// flux 50 W/m^2, view factor 1 to top) and top (eps 0.9, view factor 1 to
// bottom); front / back without spec (adiabatic).
// ---------------------------------------------------------------------------
static CaseSetup make_s2s_setup()
{
    CaseSetup setup;
    setup.energy_model = "conjugate";
    setup.radiation_model = "S2S";

    MaterialSpec water;
    water.name = "water";
    water.density = 1000.0;
    water.dynamic_viscosity = 1e-3;
    water.thermal_conductivity = 1.0;
    water.specific_heat = 1000.0;
    setup.materials.push_back(water);

    BoundarySpec left;
    left.patch_name = "left";
    left.type = BCType::WALL;
    left.value_type = BCValueType::WALL_THERMAL;
    left.wall_temperature = 400.0;
    setup.boundary_conditions.push_back(left);

    BoundarySpec right;
    right.patch_name = "right";
    right.type = BCType::WALL;
    right.value_type = BCValueType::WALL_THERMAL;
    right.wall_temperature = 300.0;
    setup.boundary_conditions.push_back(right);

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
    return setup;
}

// Independent |Sf|-weighted wall temperature of a patch: the only link
// between the resolved thermal field and the radiosity system, recomputed
// here without any production helper.
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

int main()
{
    // -----------------------------------------------------------------
    // Happy path: 4x4x1 box, two-surface parallel-plate enclosure.
    // -----------------------------------------------------------------
    run_case("production_s2s_energy_stage_parallel_plates_oracle", [] {
        const std::size_t nx = 4;
        const std::size_t ny = 4;
        const double dx = 0.25;
        const double dy = 0.25;
        auto mesh = box_mesh(nx, ny, dx, dy);
        const auto geometry = build_fv_geometry(mesh);
        auto setup = make_s2s_setup();

        Field<double, Location::FACE> mass_flux(
            mesh.n_faces(), "mass_flux", "kg/s", 1);
        mass_flux.fill(0.0);  // pure conduction carrier (frozen flow)

        // Independent pass 1 (no radiation) on a fresh 350 K field: the
        // stage's pass 1 is deterministic, so this reproduces T0 exactly.
        Field<double, Location::CELL> t0(mesh.n_cells(), "T0", "K", 1);
        t0.fill(350.0);
        const auto pass1 = run_production_energy_stage(
            mesh, geometry, setup, mass_flux, t0);
        EXPECT_TRUE(pass1.converged);

        // Independent wall temperatures from the pass-1 field.
        const double tb0 =
            patch_wall_temperature(mesh, geometry, "bottom", t0);
        const double tt0 =
            patch_wall_temperature(mesh, geometry, "top", t0);
        EXPECT_TRUE(tb0 > 0.0 && tt0 > 0.0);
        EXPECT_TRUE(tb0 > tt0);  // heated bottom wall radiates to the top

        // Hand radiosity oracle (Cramer) for the two-surface enclosure
        // with F = 1 both ways.
        const double sigma = STEFAN_BOLTZMANN;
        const double eps_b = 0.8;
        const double eps_t = 0.9;
        const double eb = sigma * std::pow(tb0, 4);
        const double et = sigma * std::pow(tt0, 4);
        const double det = 1.0 - (1.0 - eps_b) * (1.0 - eps_t);
        const double jb =
            (eps_b * eb + (1.0 - eps_b) * eps_t * et) / det;
        const double jt =
            (eps_t * et + (1.0 - eps_t) * eps_b * eb) / det;
        const double qb = jb - jt;
        // Parallel-plate identity, computed by a different formula.
        const double q_parallel =
            sigma * (std::pow(tb0, 4) - std::pow(tt0, 4)) /
            (1.0 / eps_b + 1.0 / eps_t - 1.0);
        EXPECT_NEAR(qb, q_parallel, 1e-9);

        // Full stage on a fresh 350 K field (fresh copy, pass-1 field
        // above was mutated in place by run_production_energy_stage).
        Field<double, Location::CELL> temperature(
            mesh.n_cells(), "temperature", "K", 1);
        temperature.fill(350.0);
        const auto stage = run_production_s2s_energy_stage(
            mesh, geometry, setup, mass_flux, temperature);

        EXPECT_TRUE(stage.energy_pass1_converged);
        EXPECT_TRUE(stage.s2s_radiation_converged);
        EXPECT_TRUE(stage.energy_pass2_converged);
        EXPECT_TRUE(stage.converged);

        const auto bi = surface_index(stage.radiation, "bottom");
        const auto ti = surface_index(stage.radiation, "top");
        EXPECT_TRUE(stage.radiation.surface_count == 2);

        // Radiation reproduces the hand oracle.
        EXPECT_NEAR(stage.radiation.wall_temperatures[bi], tb0, 1e-12);
        EXPECT_NEAR(stage.radiation.wall_temperatures[ti], tt0, 1e-12);
        EXPECT_NEAR(stage.radiation.radiosity[bi], jb, 1e-11);
        EXPECT_NEAR(stage.radiation.radiosity[ti], jt, 1e-11);
        EXPECT_NEAR(stage.radiation.net_flux[bi], qb, 1e-11);
        // Irradiation of a two-surface enclosure: G_i = J_j at F = 1.
        EXPECT_NEAR(stage.radiation.irradiation[bi], jt, 1e-11);
        EXPECT_NEAR(stage.radiation.irradiation[ti], jb, 1e-11);

        // Sign convention: the hot bottom surface loses (q > 0), the cold
        // top surface gains (q < 0), and the enclosure closes.
        EXPECT_TRUE(stage.radiation.net_flux[bi] > 0.0);
        EXPECT_TRUE(stage.radiation.net_flux[ti] < 0.0);
        EXPECT_NEAR(stage.radiation.net_flux[ti],
                    -stage.radiation.net_flux[bi], 1e-11);
        EXPECT_TRUE(stage.radiation.balance_passed);
        EXPECT_TRUE(stage.radiation.balance_normalized <= 1e-10);
        EXPECT_TRUE(stage.radiation.surface_identity_passed);
        EXPECT_TRUE(stage.radiation.max_surface_identity_normalized <=
                    1e-10);

        // Injected flux map: -q per patch, positive into the domain.
        EXPECT_NEAR(stage.radiation.injected_heat_flux.at("bottom"),
                    -qb, 1e-11);
        EXPECT_NEAR(stage.radiation.injected_heat_flux.at("top"),
                    -stage.radiation.net_flux[ti], 1e-11);

        // -------------------------------------------------------------
        // Re-injection contract replayed through the BC builder with the
        // ORACLE extra map (never the production values).
        // -------------------------------------------------------------
        const std::map<std::string, double> extra{
            {"bottom", -qb},
            {"top", -stage.radiation.net_flux[ti]}};
        const auto energy_case = build_production_energy_case(setup);
        const auto bcs = build_production_temperature_bcs(
            mesh, setup, energy_case, extra);
        EXPECT_TRUE(bcs.at("bottom").type ==
                    ScalarBoundaryType::FIXED_GRADIENT);
        EXPECT_NEAR(bcs.at("bottom").gradient, (50.0 - qb) / 1.0, 1e-12);
        EXPECT_TRUE(bcs.at("top").type ==
                    ScalarBoundaryType::FIXED_GRADIENT);
        EXPECT_NEAR(bcs.at("top").gradient,
                    -stage.radiation.net_flux[ti] / 1.0, 1e-12);
        EXPECT_TRUE(bcs.at("left").type == ScalarBoundaryType::FIXED_VALUE);
        EXPECT_NEAR(bcs.at("left").value, 400.0, 0.0);
        EXPECT_TRUE(bcs.at("right").type ==
                    ScalarBoundaryType::FIXED_VALUE);
        EXPECT_NEAR(bcs.at("right").value, 300.0, 0.0);

        // Two-pass replay (pass 1 fresh, pass 2 with the oracle extra
        // map) reproduces the stage's final field cell by cell.
        Field<double, Location::CELL> replay(
            mesh.n_cells(), "replay", "K", 1);
        replay.fill(350.0);
        const auto replay_pass1 = run_production_energy_stage(
            mesh, geometry, setup, mass_flux, replay);
        EXPECT_TRUE(replay_pass1.converged);
        const auto replay_pass2 = run_production_energy_stage(
            mesh, geometry, setup, mass_flux, replay, extra);
        EXPECT_TRUE(replay_pass2.converged);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            EXPECT_NEAR(replay(c), temperature(c), 1e-12);

        // -------------------------------------------------------------
        // Final energy verdict (pass 2): independent blocking audits.
        // -------------------------------------------------------------
        EXPECT_TRUE(stage.energy.linear_solver_converged);
        EXPECT_TRUE(stage.energy.global_energy_balance_passed);
        EXPECT_TRUE(stage.energy.local_energy_balance_passed);
        EXPECT_TRUE(stage.energy.global_balance_normalized <= 1e-12);
        EXPECT_TRUE(stage.energy.local_balance_normalized <= 1e-12);
        EXPECT_TRUE(stage.energy.max_cell_residual_normalized <= 1e-12);

        // -------------------------------------------------------------
        // Physical closure: the radiative exchange is internal to the
        // enclosure, so the walls carry the declared 50 W/m^2 in total.
        // Areas: bottom = top = nx * dx * 1 = 1 m^2.
        // -------------------------------------------------------------
        const double area = static_cast<double>(nx) * dx;
        EXPECT_NEAR(stage.radiation.areas[bi], area, 1e-12);
        EXPECT_NEAR(stage.radiation.areas[ti], area, 1e-12);
        EXPECT_NEAR(stage.energy.patch_flux_out.at("bottom"),
                    -area * (50.0 - qb), 1e-6);
        EXPECT_NEAR(stage.energy.patch_flux_out.at("top"),
                    area * stage.radiation.net_flux[ti], 1e-6);
        EXPECT_NEAR(stage.energy.patch_flux_out.at("left") +
                        stage.energy.patch_flux_out.at("right"),
                    50.0 * area, 1e-6);
        EXPECT_NEAR(stage.energy.patch_flux_out.at("front"), 0.0, 1e-6);
        EXPECT_NEAR(stage.energy.patch_flux_out.at("back"), 0.0, 1e-6);

        // Delayed-Picard diagnostic: finite, non-negative, NOT a gate.
        EXPECT_TRUE(std::isfinite(stage.max_wall_temperature_lag));
        EXPECT_TRUE(stage.max_wall_temperature_lag >= 0.0);
    });

    // -----------------------------------------------------------------
    // Fail-closed catalogue on a 2x2x1 box.
    // -----------------------------------------------------------------
    run_case("production_s2s_energy_stage_fail_closed", [] {
        const std::size_t nx = 2;
        const std::size_t ny = 2;
        auto mesh = box_mesh(nx, ny, 0.25, 0.25);
        const auto geometry = build_fv_geometry(mesh);

        auto fresh_field = [&](double value) {
            Field<double, Location::CELL> field(
                mesh.n_cells(), "T", "K", 1);
            field.fill(value);
            return field;
        };

        // S2S without the conjugate energy model: rejected by the
        // selection guard (explicit conjugate requirement)...
        {
            auto setup = make_s2s_setup();
            setup.energy_model = "none";
            bool threw = false;
            try {
                validate_production_radiation_selection(setup);
            } catch (const std::invalid_argument& error) {
                threw = std::string(error.what()).find("conjugate") !=
                        std::string::npos;
            }
            EXPECT_TRUE(threw);
            Field<double, Location::FACE> mass_flux(
                mesh.n_faces(), "mass_flux", "kg/s", 1);
            mass_flux.fill(0.0);
            auto temperature = fresh_field(350.0);
            EXPECT_THROW(
                run_production_s2s_energy_stage(
                    mesh, geometry, setup, mass_flux, temperature),
                std::invalid_argument);
        }

        // ...and DO / solar remain "not coupled" rejections.
        {
            auto setup = make_s2s_setup();
            setup.radiation_model = "DO";
            EXPECT_THROW(validate_production_radiation_selection(setup),
                         std::invalid_argument);
            setup.radiation_model = "solar";
            EXPECT_THROW(validate_production_radiation_selection(setup),
                         std::invalid_argument);
        }

        // Running the S2S stage without a declared S2S model.
        {
            auto setup = make_s2s_setup();
            setup.radiation_model = "none";
            Field<double, Location::FACE> mass_flux(
                mesh.n_faces(), "mass_flux", "kg/s", 1);
            mass_flux.fill(0.0);
            auto temperature = fresh_field(350.0);
            EXPECT_THROW(
                run_production_s2s_energy_stage(
                    mesh, geometry, setup, mass_flux, temperature),
                std::invalid_argument);
            EXPECT_THROW(
                solve_production_s2s_radiation(
                    mesh, geometry, setup, temperature),
                std::invalid_argument);
        }

        auto solve_throws = [&](CaseSetup setup) {
            auto temperature = fresh_field(350.0);
            EXPECT_THROW(
                solve_production_s2s_radiation(
                    mesh, geometry, setup, temperature),
                std::invalid_argument);
        };

        // No participant at all.
        {
            auto setup = make_s2s_setup();
            for (auto& spec : setup.boundary_conditions)
                spec.raw_params.erase("emissivity");
            solve_throws(setup);
        }

        // A single participant.
        {
            auto setup = make_s2s_setup();
            for (auto& spec : setup.boundary_conditions) {
                if (spec.patch_name == "bottom")
                    continue;
                spec.raw_params.erase("emissivity");
                spec.raw_params.erase("s2s_view_factor_to_bottom");
            }
            // A lone participant cannot form an enclosure: rejected.
            solve_throws(setup);
        }

        // Emissivity parsing: strict numeric, (0, 1].
        for (const char* bad : {"abc", "0.8abc", "0", "-0.5", "1.5",
                                "inf", ""}) {
            auto setup = make_s2s_setup();
            for (auto& spec : setup.boundary_conditions)
                if (spec.patch_name == "bottom")
                    spec.raw_params["emissivity"] = bad;
            solve_throws(setup);
        }

        // Emissivity on a non-WALL spec type.
        {
            auto setup = make_s2s_setup();
            for (auto& spec : setup.boundary_conditions)
                if (spec.patch_name == "bottom")
                    spec.type = BCType::SYMMETRY;
            solve_throws(setup);
        }

        // Participant with a retained Dirichlet temperature.
        {
            auto setup = make_s2s_setup();
            for (auto& spec : setup.boundary_conditions)
                if (spec.patch_name == "bottom") {
                    spec.raw_params.clear();
                    spec.value_type = BCValueType::WALL_THERMAL;
                    spec.wall_temperature = 350.0;
                    spec.raw_params["emissivity"] = "0.8";
                    spec.raw_params["s2s_view_factor_to_top"] = "1";
                }
            solve_throws(setup);
        }

        // View factor toward the patch itself.
        {
            auto setup = make_s2s_setup();
            for (auto& spec : setup.boundary_conditions)
                if (spec.patch_name == "bottom") {
                    spec.raw_params["s2s_view_factor_to_bottom"] = "0.5";
                    spec.raw_params.erase("s2s_view_factor_to_top");
                }
            solve_throws(setup);
        }

        // View factor toward a non-participant patch.
        {
            auto setup = make_s2s_setup();
            for (auto& spec : setup.boundary_conditions)
                if (spec.patch_name == "bottom")
                    spec.raw_params["s2s_view_factor_to_left"] = "1";
            solve_throws(setup);
        }

        // View factors declared without an emissivity: the existing left
        // spec (Dirichlet wall, never a participant) gains a view factor
        // key, which a non-participant must not carry.
        {
            auto setup = make_s2s_setup();
            for (auto& spec : setup.boundary_conditions)
                if (spec.patch_name == "left")
                    spec.raw_params["s2s_view_factor_to_bottom"] = "0.5";
            solve_throws(setup);
        }

        // Non-closed enclosure (0.5 / 0.5) and non-reciprocal enclosure
        // (1.0 / 0.9): rejected, never silently patched (fiche 11 B5).
        {
            auto setup = make_s2s_setup();
            for (auto& spec : setup.boundary_conditions) {
                if (spec.patch_name == "bottom")
                    spec.raw_params["s2s_view_factor_to_top"] = "0.5";
                if (spec.patch_name == "top")
                    spec.raw_params["s2s_view_factor_to_bottom"] = "0.5";
            }
            solve_throws(setup);
        }
        {
            auto setup = make_s2s_setup();
            for (auto& spec : setup.boundary_conditions) {
                if (spec.patch_name == "top")
                    spec.raw_params["s2s_view_factor_to_bottom"] = "0.9";
            }
            solve_throws(setup);
        }

        // View factor out of range or non-numeric.
        for (const char* bad : {"-0.1", "1.5", "abc"}) {
            auto setup = make_s2s_setup();
            for (auto& spec : setup.boundary_conditions)
                if (spec.patch_name == "bottom")
                    spec.raw_params["s2s_view_factor_to_top"] = bad;
            solve_throws(setup);
        }

        // Temperature field: non-finite and non-positive are rejected.
        {
            auto setup = make_s2s_setup();
            auto bad = fresh_field(350.0);
            bad(0) = std::nan("");
            EXPECT_THROW(
                solve_production_s2s_radiation(mesh, geometry, setup, bad),
                std::invalid_argument);
            auto negative = fresh_field(350.0);
            negative(1) = -1.0;
            EXPECT_THROW(
                solve_production_s2s_radiation(mesh, geometry, setup,
                                                negative),
                std::invalid_argument);
        }

        // Extra wall heat flux: Dirichlet patch, unknown patch, EMPTY
        // patch and non-finite values are all rejected.
        {
            auto setup = make_s2s_setup();
            const auto energy_case = build_production_energy_case(setup);
            EXPECT_THROW(
                build_production_temperature_bcs(
                    mesh, setup, energy_case, {{"left", 10.0}}),
                std::invalid_argument);
            EXPECT_THROW(
                build_production_temperature_bcs(
                    mesh, setup, energy_case, {{"nope", 10.0}}),
                std::invalid_argument);
            EXPECT_THROW(
                build_production_temperature_bcs(
                    mesh, setup, energy_case,
                    {{"bottom", std::nan("")}}),
                std::invalid_argument);
            auto empty_front_mesh = box_mesh(nx, ny, 0.25, 0.25,
                                             PatchType::EMPTY);
            const auto empty_geometry = build_fv_geometry(empty_front_mesh);
            EXPECT_THROW(
                build_production_temperature_bcs(
                    empty_front_mesh, setup, energy_case,
                    {{"front", 10.0}}),
                std::invalid_argument);
        }
    });

    return run_all();
}
