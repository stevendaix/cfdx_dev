// T1 qualification: solid conduction through multi-material zones with
// volumetric power and contact resistance (issues #880, #885).
//
// Verified properties:
//   * Composite wall, no source: the two-point series-resistance
//     discretization reproduces the exact resistance-network solution at
//     cell centres (machine precision gate).
//   * Contact resistance: still exact, R_c simply adds in series at
//     cross-zone interfaces.
//   * Volumetric power: second-order convergence towards the analytic
//     parabola under mesh refinement, exact cumulative-power identity
//     through an internal plane, and energy closure
//     |produced - net boundary flux| / scale < 1e-12.
//   * Fixed-gradient (heat flux) boundary: exact linear solution and closed
//     energy balance.
//   * Fail-closed: invalid conductivity, undeclared zone id, inconsistent
//     double-sided contact resistance, non-finite power, duplicate zone id
//     and a Dirichlet-free (singular) system are all rejected.
#include "cfdx/physics/thermal_zones.h"
#include "common/test_harness.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

// ---------------------------------------------------------------------------
// 1D chain mesh: n cells of width dx along x, unit cross-section. Columns of
// 4 points at each x boundary; one x-face per boundary, four wall faces per
// cell. Patches: "left", "right", "walls" (mirrors test_thermal_vv).
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

// Zone map: cells [0, n/3) -> id a, [n/3, 2n/3) -> id b, rest -> id c.
static std::vector<int> three_zone_map(std::size_t n, int a, int b, int c)
{
    std::vector<int> zone(n, c);
    for (std::size_t i = 0; i < n; ++i) {
        if (i < n / 3) zone[i] = a;
        else if (i < 2 * (n / 3)) zone[i] = b;
    }
    return zone;
}

// Exact cell-centre temperatures of the resistance network for a chain of
// equal cells: piecewise conductivity k[i], Dirichlet ends, optional contact
// resistance r[i] on the internal face between cells i and i+1.
static std::vector<double> exact_chain_temperatures(
    const std::vector<double>& k,
    double dx,
    double t_left,
    double t_right,
    const std::vector<double>& interface_resistance)
{
    const std::size_t n = k.size();
    std::vector<double> r_between(n - 1, 0.0);
    for (std::size_t i = 0; i + 1 < n; ++i)
        r_between[i] = dx / (2.0 * k[i]) +
                       (i < interface_resistance.size() ? interface_resistance[i] : 0.0) +
                       dx / (2.0 * k[i + 1]);
    const double r_wall_left = dx / (2.0 * k.front());
    const double r_wall_right = dx / (2.0 * k.back());
    double r_total = r_wall_left + r_wall_right;
    for (const double r : r_between) r_total += r;
    const double q = (t_left - t_right) / r_total;  // W through the wall

    std::vector<double> t(n, 0.0);
    t[0] = t_left - q * r_wall_left;
    for (std::size_t i = 1; i < n; ++i) t[i] = t[i - 1] - q * r_between[i - 1];
    return t;
}

static ScalarBoundaryConditions chain_bcs(double t_left, double t_right)
{
    ScalarBoundaryConditions bc;
    bc["left"] = {ScalarBoundaryType::FIXED_VALUE, t_left, 0.0};
    bc["right"] = {ScalarBoundaryType::FIXED_VALUE, t_right, 0.0};
    bc["walls"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    return bc;
}

static double max_abs_error(const Vector& temperature,
                            const std::vector<double>& expected)
{
    double error = 0.0;
    for (std::size_t i = 0; i < expected.size(); ++i)
        error = std::max(error, std::abs(temperature(i) - expected[i]));
    return error;
}

int main()
{
    run_case("thermal_zones_composite_wall_exact_network", [] {
        const std::size_t n = 6;
        const double dx = 0.5;
        Mesh m = chain_mesh(n, dx);
        auto g = build_fv_geometry(m);

        // Three zones of two cells each with distinct conductivities.
        std::vector<ThermalZoneSpec> specs = {
            {"mat_a", 10, 1.0, 0.0, 0.0, true},
            {"mat_b", 20, 2.0, 0.0, 0.0, true},
            {"mat_c", 30, 4.0, 0.0, 0.0, true},
        };
        const auto zones = build_thermal_zones(m, three_zone_map(n, 10, 20, 30), specs);
        EXPECT_TRUE(zones.diagnostics.n_zones == 3);
        EXPECT_TRUE(zones.diagnostics.cells_per_zone.at("mat_a") == 2);
        EXPECT_TRUE(zones.diagnostics.n_interface_faces == 2);

        const auto bc = chain_bcs(400.0, 300.0);
        Vector temperature(n, 293.15);
        const auto result =
            solve_thermal_conduction(m, g, zones, bc, temperature,
                                     ScalarSolveControls{5000, 1e-12, 1.0});

        const std::vector<double> k = {1, 1, 2, 2, 4, 4};
        const auto expected = exact_chain_temperatures(k, dx, 400.0, 300.0, {});
        EXPECT_TRUE(result.converged);
        EXPECT_TRUE(max_abs_error(temperature, expected) < 1e-9);
        EXPECT_TRUE(result.energy_imbalance_relative < 1e-12);
        // Closed balance: everything entering on the left leaves on the right.
        EXPECT_NEAR(result.patch_flux_out.at("left"),
                    -result.patch_flux_out.at("right"), 1e-9);
        std::cout << "THERMAL_ZONES_VV: composite_wall"
                  << " max_error=" << max_abs_error(temperature, expected)
                  << " flux_left=" << result.patch_flux_out.at("left")
                  << " imbalance=" << result.energy_imbalance_relative << '\n';
    });

    run_case("thermal_zones_contact_resistance_exact", [] {
        const std::size_t n = 6;
        const double dx = 0.5;
        Mesh m = chain_mesh(n, dx);
        auto g = build_fv_geometry(m);

        // Contact resistance declared on the middle zone applies at both of
        // its interfaces (mat_a|mat_b and mat_b|mat_c).
        std::vector<ThermalZoneSpec> specs = {
            {"mat_a", 10, 1.0, 0.0, 0.0, true},
            {"mat_b", 20, 2.0, 0.0, 0.05, true},
            {"mat_c", 30, 4.0, 0.0, 0.0, true},
        };
        const auto zones = build_thermal_zones(m, three_zone_map(n, 10, 20, 30), specs);
        EXPECT_NEAR(zones.diagnostics.max_contact_resistance, 0.05, 1e-15);

        const auto bc = chain_bcs(400.0, 300.0);
        Vector temperature(n, 293.15);
        const auto result =
            solve_thermal_conduction(m, g, zones, bc, temperature,
                                     ScalarSolveControls{5000, 1e-12, 1.0});

        const std::vector<double> k = {1, 1, 2, 2, 4, 4};
        // Internal faces 1|2 and 3|4 both carry the contact resistance.
        const auto expected = exact_chain_temperatures(
            k, dx, 400.0, 300.0, {0.0, 0.05, 0.0, 0.05, 0.0});
        EXPECT_TRUE(result.converged);
        EXPECT_TRUE(max_abs_error(temperature, expected) < 1e-9);
        EXPECT_TRUE(result.energy_imbalance_relative < 1e-12);
    });

    run_case("thermal_zones_volumetric_power_balance_and_convergence", [] {
        const double power = 1.0e3;  // W/m^3, uniform, k = 1, L = 2 m
        const auto solve_chain = [power](std::size_t n) {
            const double dx = 2.0 / static_cast<double>(n);
            Mesh m = chain_mesh(n, dx);
            auto g = build_fv_geometry(m);
            std::vector<ThermalZoneSpec> specs = {{"mat", 1, 1.0, power, 0.0, true}};
            const std::vector<int> zone(n, 1);
            const auto zones = build_thermal_zones(m, zone, specs);
            const auto bc = chain_bcs(350.0, 300.0);
            Vector temperature(n, 300.0);
            const auto result =
                solve_thermal_conduction(m, g, zones, bc, temperature,
                                         ScalarSolveControls{5000, 1e-12, 1.0});
            // Analytic cell-centre values: T(x) = T_L - (T_L-T_R) x/L +
            // q''' x (L-x) / (2k). The two-point discretization is second
            // order, so this error must shrink ~4x per refinement.
            double error = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                const double x = (static_cast<double>(i) + 0.5) * dx;
                const double expected = 350.0 - 25.0 * x + 500.0 * x * (2.0 - x);
                error = std::max(error, std::abs(temperature(i) - expected));
            }
            // Power audit: the conductive flux through the internal face at
            // x = 1 m must carry exactly half the generated power. This is an
            // exact discrete identity, checked at machine precision.
            double mid_flux_error = -1.0;
            if (n % 2 == 0) {
                const std::size_t mid = n / 2;  // face between cells mid-1, mid
                const double conductance = 1.0 / dx;  // k=1, A=1, uniform cells
                const double flux = conductance * (temperature(mid - 1) - temperature(mid));
                mid_flux_error = std::abs(flux - power * 1.0) / (power * 1.0);
            }
            return std::pair<double, double>(error, mid_flux_error);
        };

        const auto coarse = solve_chain(8);
        const auto fine = solve_chain(32);
        // Discretization error is real and shrinks at second order.
        EXPECT_TRUE(coarse.first > 1.0);
        EXPECT_TRUE(fine.first < 0.25 * coarse.first);
        // Cumulative-power identity through the mid plane: machine precision.
        EXPECT_TRUE(coarse.second >= 0.0 && coarse.second < 1e-9);
        EXPECT_TRUE(fine.second >= 0.0 && fine.second < 1e-9);
        std::cout << "THERMAL_ZONES_VV: volumetric_power"
                  << " error_coarse=" << coarse.first
                  << " error_fine=" << fine.first
                  << " mid_flux_rel_error=" << fine.second << '\n';
    });

    run_case("thermal_zones_energy_closure_with_source", [] {
        const std::size_t n = 16;
        const double dx = 2.0 / static_cast<double>(n);
        Mesh m = chain_mesh(n, dx);
        auto g = build_fv_geometry(m);
        const double power = 1.0e3;  // W/m^3
        std::vector<ThermalZoneSpec> specs = {{"mat", 1, 1.0, power, 0.0, true}};
        const std::vector<int> zone(n, 1);
        const auto zones = build_thermal_zones(m, zone, specs);
        const auto bc = chain_bcs(350.0, 300.0);
        Vector temperature(n, 300.0);
        const auto result =
            solve_thermal_conduction(m, g, zones, bc, temperature,
                                     ScalarSolveControls{5000, 1e-12, 1.0});
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(result.total_source_power, power * 2.0, 1e-9);
        EXPECT_TRUE(result.energy_imbalance_relative < 1e-12);
        // All generated power leaves through the two Dirichlet ends.
        EXPECT_NEAR(result.boundary_flux_out, power * 2.0, 1e-6);
    });

    run_case("thermal_zones_heat_flux_boundary_linear_profile", [] {
        const std::size_t n = 8;
        const double dx = 0.25;  // L = 2 m
        Mesh m = chain_mesh(n, dx);
        auto g = build_fv_geometry(m);

        const double heat_flux = 500.0;  // W/m^2, positive into the domain
        std::vector<ThermalZoneSpec> specs = {{"mat", 1, 2.0, 0.0, 0.0, true}};
        const std::vector<int> zone(n, 1);
        const auto zones = build_thermal_zones(m, zone, specs);

        ScalarBoundaryConditions bc;
        bc["left"] = {ScalarBoundaryType::FIXED_GRADIENT, 0.0, heat_flux / 2.0};
        bc["right"] = {ScalarBoundaryType::FIXED_VALUE, 300.0, 0.0};
        bc["walls"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        Vector temperature(n, 300.0);
        const auto result =
            solve_thermal_conduction(m, g, zones, bc, temperature,
                                     ScalarSolveControls{5000, 1e-12, 1.0});
        // T(x) = T_R + q'' (L - x) / k, exact at cell centres.
        double error = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double x = (static_cast<double>(i) + 0.5) * dx;
            const double expected = 300.0 + heat_flux * (2.0 - x) / 2.0;
            error = std::max(error, std::abs(temperature(i) - expected));
        }
        EXPECT_TRUE(result.converged);
        EXPECT_TRUE(error < 1e-9);
        // No volumetric source: the injected flux must leave through the
        // right patch (audit signs: left negative = entering).
        EXPECT_NEAR(result.patch_flux_out.at("left"), -heat_flux, 1e-9);
        EXPECT_NEAR(result.boundary_flux_out, 0.0, 1e-9);
        EXPECT_TRUE(result.energy_imbalance_relative < 1e-12);
    });

    run_case("thermal_zones_fail_closed_invalid_declarations", [] {
        const std::size_t n = 6;
        Mesh m = chain_mesh(n, 0.5);
        auto g = build_fv_geometry(m);

        bool threw = false;
        try {
            std::vector<ThermalZoneSpec> bad_k = {
                {"mat", 1, -1.0, 0.0, 0.0, true}};
            (void)build_thermal_zones(m, std::vector<int>(n, 1), bad_k);
        } catch (const std::exception&) {
            threw = true;
        }
        EXPECT_TRUE(threw);

        threw = false;
        try {
            std::vector<ThermalZoneSpec> nan_power = {
                {"mat", 1, 1.0, std::numeric_limits<double>::quiet_NaN(), 0.0, true}};
            (void)build_thermal_zones(m, std::vector<int>(n, 1), nan_power);
        } catch (const std::exception&) {
            threw = true;
        }
        EXPECT_TRUE(threw);

        threw = false;
        try {
            // Cell 0 references zone id 7 which is not declared.
            std::vector<int> zone(n, 1);
            zone[0] = 7;
            std::vector<ThermalZoneSpec> specs = {{"mat", 1, 1.0, 0.0, 0.0, true}};
            (void)build_thermal_zones(m, zone, specs);
        } catch (const std::exception&) {
            threw = true;
        }
        EXPECT_TRUE(threw);

        threw = false;
        try {
            // Both sides of the mat_a|mat_b interface declare different
            // contact resistances. All three zone ids are declared so the
            // rejection is specifically the resistance inconsistency.
            std::vector<ThermalZoneSpec> specs = {
                {"mat_a", 10, 1.0, 0.0, 0.01, true},
                {"mat_b", 20, 2.0, 0.0, 0.02, true},
                {"mat_c", 30, 4.0, 0.0, 0.0, true},
            };
            (void)build_thermal_zones(m, three_zone_map(n, 10, 20, 30), specs);
        } catch (const std::exception&) {
            threw = true;
        }
        EXPECT_TRUE(threw);

        threw = false;
        try {
            // No Dirichlet boundary anywhere: the system is singular and the
            // assembly must reject it instead of producing garbage.
            std::vector<ThermalZoneSpec> specs = {{"mat", 1, 1.0, 0.0, 0.0, true}};
            const auto zones = build_thermal_zones(m, std::vector<int>(n, 1), specs);
            ScalarBoundaryConditions adiabatic;
            adiabatic["left"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
            adiabatic["right"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
            adiabatic["walls"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
            (void)assemble_thermal_conduction(m, g, zones, adiabatic);
        } catch (const std::exception&) {
            threw = true;
        }
        EXPECT_TRUE(threw);

        threw = false;
        try {
            // Duplicate zone ids are ambiguous and rejected.
            std::vector<ThermalZoneSpec> dup = {
                {"mat_a", 1, 1.0, 0.0, 0.0, true},
                {"mat_b", 1, 2.0, 0.0, 0.0, true}};
            (void)build_thermal_zones(m, std::vector<int>(n, 1), dup);
        } catch (const std::exception&) {
            threw = true;
        }
        EXPECT_TRUE(threw);
    });

    return run_all();
}
