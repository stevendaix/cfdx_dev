// Production S2S radiation coupled to the energy stage - external
// energy-radiation iteration (capability #885 step 3, issue #885):
// verification on a frozen flow.
//
// Verified properties:
//   * The outer loop converges to the COUPLED FIXED POINT. The final
//     temperature field is the energy solution for the applied radiative
//     flux map, and the radiative response of that same field - recomputed
//     here with an independent hand Cramer oracle on the final wall
//     temperatures - reproduces the applied map: replaying ONE energy pass
//     with the ORACLE map on a fresh field reproduces the loop's final
//     field cell by cell (to the fixed-point tolerance coupled through
//     the wall-flux response).
//   * Separate blocking criteria: the unrelaxed radiative fixed-point
//     residual and the wall-temperature stationarity are both <= 1e-10
//     and reported; every history entry records green inner stages
//     (fiche 14 energy gates, fiche 11 radiosity machine gates).
//   * The loop actually iterates: the declared 50 W/m^2 bottom flux moves
//     the pass-1 wall temperatures away from the uniform initial field,
//     so at least two outer iterations run.
//   * The stage-2 single delayed-Picard pass is NOT the fixed point: the
//     radiative response of its own final field differs from the map it
//     applied by orders of magnitude above the fixed-point tolerance,
//     while the outer loop's does not - the loop adds the missing
//     convergence.
//   * Physical closure: the radiative exchange is internal to the
//     enclosure, so the walls carry the declared 50 W/m^2 in total and
//     the adiabatic front/back patches carry ~0.
//   * Sign convention (fiche 11 B3) on the final response: the hot bottom
//     surface loses, the cold top gains, the pair closes to machine zero.
//   * Fail-closed: running the outer loop without a declared S2S
//     radiation model or without the conjugate energy model is rejected.
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
// Same construction as the stage-2 VV so both tests exercise the identical
// production geometry path.
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
// Case: conjugate energy + S2S radiation. Dirichlet walls left (400 K) and
// right (300 K); radiative participants bottom (eps 0.8, declared heat
// flux 50 W/m^2, view factor 1 to top) and top (eps 0.9, view factor 1 to
// bottom); front / back without spec (adiabatic). Same case as the stage-2
// VV: the outer loop must converge THIS configuration to its fixed point.
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

// Independent boundary-law reconstruction, without the production helper.
static double patch_wall_temperature(
    const Mesh& mesh, const FvGeometry& geometry, const std::string& name,
    const Field<double, Location::CELL>& temperature,
    double heat_flux_into_domain,
    const std::vector<double>* conductivity_values = nullptr)
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
            const auto delta = geometry.face_centres[face] -
                               geometry.cell_centres[owner];
            const double normal_distance = std::abs(
                delta.x * area.x + delta.y * area.y + delta.z * area.z) / weight;
            const double conductivity = conductivity_values
                ? conductivity_values->at(owner) : 1.0;
            weighted += weight * (temperature(owner) +
                heat_flux_into_domain * normal_distance / conductivity);
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
    run_case("production_s2s_outer_loop_isothermal_zero_applied_map", [] {
        auto mesh = box_mesh(2, 2, 0.25, 0.25);
        const auto geometry = build_fv_geometry(mesh);
        auto setup = make_s2s_setup();
        for (auto& spec : setup.boundary_conditions) {
            if (spec.patch_name == "left" || spec.patch_name == "right")
                spec.wall_temperature = 350.0;
            if (spec.patch_name == "bottom") spec.raw_params["heat_flux"] = "0";
        }
        Field<double, Location::FACE> mass_flux(mesh.n_faces(), "mass_flux", "kg/s", 1);
        mass_flux.fill(0.0);
        Field<double, Location::CELL> temperature(mesh.n_cells(), "T", "K", 1);
        temperature.fill(350.0);
        const auto result = run_production_s2s_energy_outer_loop(
            mesh, geometry, setup, mass_flux, temperature);
        EXPECT_TRUE(result.converged);
        EXPECT_TRUE(result.outer_iterations == 1);
        EXPECT_TRUE(result.radiation.applied_heat_flux.size() == 2);
        for (const auto& name : result.radiation.surface_names)
            EXPECT_NEAR(result.radiation.applied_heat_flux.at(name), 0.0, 0.0);
        for (const double flux : result.radiation.net_flux)
            EXPECT_NEAR(flux, 0.0, 1e-10);
    });

    // -----------------------------------------------------------------
    // Happy path: 4x4x1 box, two-surface parallel-plate enclosure, outer
    // energy-radiation iteration to the coupled fixed point.
    // -----------------------------------------------------------------
    run_case("production_s2s_energy_outer_loop_fixed_point", [] {
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

        Field<double, Location::CELL> temperature(
            mesh.n_cells(), "temperature", "K", 1);
        temperature.fill(350.0);
        const auto stage = run_production_s2s_energy_outer_loop(
            mesh, geometry, setup, mass_flux, temperature);

        // Verdict: converged coupled fixed point with BOTH separate
        // criteria green.
        EXPECT_TRUE(stage.converged);
        EXPECT_TRUE(stage.fixed_point_flux_passed);
        EXPECT_TRUE(stage.fixed_point_wall_temperature_passed);
        EXPECT_TRUE(stage.max_flux_fixed_point_normalized <= 1e-10);
        EXPECT_TRUE(stage.max_wall_temperature_change_normalized <=
                    1e-10);
        EXPECT_TRUE(stage.outer_iterations >= 2);
        EXPECT_TRUE(!stage.history.empty());
        for (const auto& entry : stage.history) {
            EXPECT_TRUE(entry.energy_converged);
            EXPECT_TRUE(entry.radiation_converged);
            EXPECT_TRUE(std::isfinite(
                entry.max_flux_fixed_point_normalized));
            EXPECT_TRUE(std::isfinite(
                entry.max_wall_temperature_change_normalized));
        }

        // Inner gates of the final pass (fiche 14 / fiche 11).
        EXPECT_TRUE(stage.energy.converged);
        EXPECT_TRUE(stage.energy.global_balance_normalized <= 1e-12);
        EXPECT_TRUE(stage.energy.local_balance_normalized <= 1e-12);
        EXPECT_TRUE(stage.energy.max_cell_residual_normalized <= 1e-12);
        EXPECT_TRUE(stage.radiation.converged);
        EXPECT_TRUE(stage.radiation.balance_normalized <= 1e-10);
        EXPECT_TRUE(stage.radiation.max_surface_identity_normalized <=
                    1e-10);

        const auto bi = surface_index(stage.radiation, "bottom");
        const auto ti = surface_index(stage.radiation, "top");
        EXPECT_TRUE(stage.radiation.surface_count == 2);

        // Independent wall temperatures of the FINAL field.
        const double tb =
            patch_wall_temperature(mesh, geometry, "bottom", temperature,
                50.0 + stage.radiation.applied_heat_flux.at("bottom"));
        const double tt =
            patch_wall_temperature(mesh, geometry, "top", temperature,
                stage.radiation.applied_heat_flux.at("top"));
        EXPECT_NEAR(stage.radiation.wall_temperatures[bi], tb, 1e-12);
        EXPECT_NEAR(stage.radiation.wall_temperatures[ti], tt, 1e-12);

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
        EXPECT_NEAR(stage.radiation.net_flux[bi], qb, 1e-11);
        EXPECT_NEAR(stage.radiation.net_flux[ti], qt, 1e-11);
        // Fiche 11 B3 sign convention + machine closure of the pair.
        EXPECT_TRUE(stage.radiation.net_flux[bi] > 0.0);
        EXPECT_TRUE(stage.radiation.net_flux[ti] < 0.0);
        EXPECT_NEAR(stage.radiation.net_flux[bi] +
                        stage.radiation.net_flux[ti],
                    0.0, 1e-11);

        // -------------------------------------------------------------
        // Fixed point replayed INDEPENDENTLY: ONE energy pass with the
        // ORACLE map on a fresh field must reproduce the loop's final
        // field. Tolerance note: the applied production map and the
        // oracle response agree to the fixed-point tolerance (1e-10
        // normalized by the gross radiant scale), which couples back
        // through the wall-flux response; 1e-6 K bounds that with
        // orders of margin.
        // -------------------------------------------------------------
        Field<double, Location::CELL> replay(
            mesh.n_cells(), "replay", "K", 1);
        replay.fill(350.0);
        const std::map<std::string, double> oracle_flux{
            {"bottom", -qb}, {"top", -qt}};
        const auto replay_energy = run_production_energy_stage(
            mesh, geometry, setup, mass_flux, replay, oracle_flux);
        EXPECT_TRUE(replay_energy.converged);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            EXPECT_NEAR(replay(c), temperature(c), 1e-6);

        // -------------------------------------------------------------
        // Physical closure (areas: bottom = top = nx * dx = 1 m^2): the
        // radiative exchange is internal to the enclosure, so the walls
        // carry the declared 50 W/m^2 in total.
        // -------------------------------------------------------------
        const double area = static_cast<double>(nx) * dx;
        EXPECT_NEAR(stage.energy.patch_flux_out.at("bottom"),
                    -area * (50.0 - qb), 1e-6);
        EXPECT_NEAR(stage.energy.patch_flux_out.at("top"),
                    area * qt, 1e-6);
        EXPECT_NEAR(stage.energy.patch_flux_out.at("left") +
                        stage.energy.patch_flux_out.at("right"),
                    50.0 * area, 1e-6);
        EXPECT_NEAR(stage.energy.patch_flux_out.at("front"), 0.0, 1e-6);
        EXPECT_NEAR(stage.energy.patch_flux_out.at("back"), 0.0, 1e-6);

        // -------------------------------------------------------------
        // The stage-2 single delayed-Picard pass is NOT the fixed
        // point: the radiative response of its own final field differs
        // from the map it applied by orders of magnitude above the
        // fixed-point tolerance, and its final field differs from the
        // converged outer-loop field.
        // -------------------------------------------------------------
        Field<double, Location::CELL> two_pass_field(
            mesh.n_cells(), "two_pass", "K", 1);
        two_pass_field.fill(350.0);
        const auto two_pass = run_production_s2s_energy_stage(
            mesh, geometry, setup, mass_flux, two_pass_field);
        EXPECT_TRUE(two_pass.converged);  // its own stage-2 gates pass
        const auto response = solve_production_s2s_radiation(
            mesh, geometry, setup, two_pass_field,
            two_pass.radiation.injected_heat_flux);
        double two_pass_residual = 0.0;
        for (std::size_t i = 0; i < response.surface_count; ++i) {
            const std::string& name = response.surface_names[i];
            const double applied =
                two_pass.radiation.injected_heat_flux.at(name);
            const double black_body_power =
                sigma * std::pow(response.wall_temperatures[i], 4);
            const double scale = std::max(
                response.emissivities[i] * black_body_power +
                    response.irradiation[i],
                1.0);
            two_pass_residual = std::max(
                two_pass_residual,
                std::abs(response.injected_heat_flux.at(name) - applied) /
                    scale);
        }
        EXPECT_TRUE(two_pass_residual > 1e-6);
        bool fields_differ = false;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            if (std::abs(two_pass_field(c) - temperature(c)) > 1e-9)
                fields_differ = true;
        EXPECT_TRUE(fields_differ);

        // The last iteration's lag is a finite, non-negative diagnostic.
        EXPECT_TRUE(std::isfinite(stage.max_wall_temperature_lag));
        EXPECT_TRUE(stage.max_wall_temperature_lag >= 0.0);
    });

    run_case("production_s2s_outer_loop_local_conductivity_replay", [] {
        auto mesh = box_mesh(4, 4, 0.25, 0.25);
        const auto geometry = build_fv_geometry(mesh);
        const auto setup = make_s2s_setup();
        Field<double, Location::FACE> mass_flux(mesh.n_faces(), "mass_flux", "kg/s", 1);
        mass_flux.fill(0.0);
        Field<double, Location::CELL> temperature(mesh.n_cells(), "T", "K", 1);
        temperature.fill(350.0);
        std::vector<double> conductivity(mesh.n_cells());
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            conductivity[c] = 2.0 + geometry.cell_centres[c].x;
        const auto result = run_production_s2s_energy_outer_loop(
            mesh, geometry, setup, mass_flux, temperature, &conductivity);
        EXPECT_TRUE(result.converged);
        EXPECT_TRUE(result.fixed_point_flux_passed);
        EXPECT_TRUE(result.fixed_point_wall_temperature_passed);
        EXPECT_TRUE(result.energy.global_balance_normalized <= 1e-12);
        EXPECT_TRUE(result.energy.local_balance_normalized <= 1e-12);
        const double tb = patch_wall_temperature(
            mesh, geometry, "bottom", temperature,
            50.0 + result.radiation.applied_heat_flux.at("bottom"), &conductivity);
        const double tt = patch_wall_temperature(
            mesh, geometry, "top", temperature,
            result.radiation.applied_heat_flux.at("top"), &conductivity);
        const auto bi = surface_index(result.radiation, "bottom");
        const auto ti = surface_index(result.radiation, "top");
        EXPECT_NEAR(result.radiation.wall_temperatures[bi], tb, 1e-12);
        EXPECT_NEAR(result.radiation.wall_temperatures[ti], tt, 1e-12);
        const double oracle = STEFAN_BOLTZMANN *
            (std::pow(tb, 4) - std::pow(tt, 4)) /
            (1.0 / 0.8 + 1.0 / 0.9 - 1.0);
        EXPECT_NEAR(result.radiation.net_flux[bi], oracle, 1e-10);
        EXPECT_NEAR(result.radiation.net_flux[ti], -oracle, 1e-10);
        Field<double, Location::CELL> replay(mesh.n_cells(), "replay", "K", 1);
        replay.fill(350.0);
        const std::map<std::string, double> oracle_flux{
            {"bottom", -oracle}, {"top", oracle}};
        const auto replay_result = run_production_energy_stage(
            mesh, geometry, setup, mass_flux, replay, oracle_flux, &conductivity);
        EXPECT_TRUE(replay_result.converged);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            EXPECT_NEAR(replay(c), temperature(c), 1e-6);
        EXPECT_NEAR(result.energy.patch_flux_out.at("left") +
                    result.energy.patch_flux_out.at("right"), 50.0, 1e-6);
    });

    // -----------------------------------------------------------------
    // Fail-closed catalogue on a 2x2x1 box.
    // -----------------------------------------------------------------
    run_case("production_s2s_energy_outer_loop_fail_closed", [] {
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
        Field<double, Location::FACE> mass_flux(
            mesh.n_faces(), "mass_flux", "kg/s", 1);
        mass_flux.fill(0.0);

        // S2S without the conjugate energy model: rejected by the loop's
        // own guard (and by the selection guard).
        {
            auto setup = make_s2s_setup();
            setup.energy_model = "none";
            auto temperature = fresh_field(350.0);
            EXPECT_THROW(
                run_production_s2s_energy_outer_loop(
                    mesh, geometry, setup, mass_flux, temperature),
                std::invalid_argument);
        }

        // Running the outer loop without a declared S2S model.
        {
            auto setup = make_s2s_setup();
            setup.radiation_model = "none";
            auto temperature = fresh_field(350.0);
            EXPECT_THROW(
                run_production_s2s_energy_outer_loop(
                    mesh, geometry, setup, mass_flux, temperature),
                std::invalid_argument);
        }
    });

    return run_all();
}
