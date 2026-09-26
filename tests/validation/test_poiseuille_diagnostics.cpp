// Poiseuille Level-A diagnostic validation.
//
// This executable intentionally accumulates all diagnostics and returns a
// non-zero status only after every refinement level and every diagnostic has
// been evaluated. CI therefore reports the complete failure set instead of
// stopping at the first failed assertion.
//
// Scope: scalar finite-volume verification path used by momentum diffusion.
// This is not a pressure-velocity-coupled Navier-Stokes validation.

#include "cfdx/core/solvers/scalar_diffusion.h"
#include "verification_metrics.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::verification;

namespace {

struct MeshCase {
    Mesh mesh;
};

MeshCase make_channel(std::size_t n)
{
    if (n == 0) throw std::invalid_argument("n must be positive");

    Mesh m;
    m.points().resize(8 * n);
    const double dy = 1.0 / static_cast<double>(n);

    for (std::size_t i = 0; i < n; ++i) {
        const double y0 = dy * static_cast<double>(i);
        const double y1 = dy * static_cast<double>(i + 1);
        const std::size_t b = 8 * i;
        const double p[8][3] = {
            {0.0,y0,0.0}, {1.0,y0,0.0}, {1.0,y1,0.0}, {0.0,y1,0.0},
            {0.0,y0,1.0}, {1.0,y0,1.0}, {1.0,y1,1.0}, {0.0,y1,1.0}
        };
        for (std::size_t j = 0; j < 8; ++j)
            m.points().set(b+j, p[j][0], p[j][1], p[j][2]);
    }

    std::vector<std::size_t> bottom, top, x0, x1, z0, z1, internal;
    auto add_face = [&](std::initializer_list<std::size_t> vertices) {
        const std::size_t id = m.faces().n_faces();
        m.faces().push_face(std::vector<FaceIndex>(vertices.begin(), vertices.end()));
        return id;
    };

    bottom.push_back(add_face({0,1,5,4}));
    top.push_back(add_face({8*(n-1)+3,8*(n-1)+7,8*(n-1)+6,8*(n-1)+2}));
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t b = 8*i;
        x0.push_back(add_face({b,b+4,b+7,b+3}));
        x1.push_back(add_face({b+1,b+2,b+6,b+5}));
        z0.push_back(add_face({b,b+3,b+2,b+1}));
        z1.push_back(add_face({b+4,b+5,b+6,b+7}));
    }
    for (std::size_t i = 0; i + 1 < n; ++i) {
        const std::size_t b = 8*i;
        internal.push_back(add_face({b+3,b+7,b+6,b+2}));
    }

    m.ownership().resize(m.n_faces());
    m.ownership().set_owner(bottom[0], 0);
    m.ownership().set_neighbour(bottom[0], FaceOwnership::BOUNDARY);
    m.ownership().set_owner(top[0], n-1);
    m.ownership().set_neighbour(top[0], FaceOwnership::BOUNDARY);

    for (std::size_t i = 0; i < n; ++i) {
        for (const auto f : {x0[i], x1[i], z0[i], z1[i]}) {
            m.ownership().set_owner(f, i);
            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        }
    }
    for (std::size_t i = 0; i + 1 < n; ++i) {
        m.ownership().set_owner(internal[i], i);
        m.ownership().set_neighbour(internal[i], static_cast<std::int64_t>(i+1));
    }

    for (std::size_t i = 0; i < n; ++i) {
        m.cells().push_cell({
            i == 0 ? bottom[0] : internal[i-1],
            i+1 == n ? top[0] : internal[i],
            x0[i], x1[i], z0[i], z1[i]
        });
    }

    auto add_patch = [&](const char* name, const std::vector<std::size_t>& faces) {
        Patch p;
        p.name = name;
        p.type = PatchType::WALL;
        p.face_ids = faces;
        m.boundary().add_patch(p);
    };
    add_patch("bottom", bottom);
    add_patch("top", top);
    add_patch("x0", x0);
    add_patch("x1", x1);
    add_patch("z0", z0);
    add_patch("z1", z1);

    return {std::move(m)};
}

struct Result {
    std::size_t n = 0;
    double l2 = 0.0;
    double linf = 0.0;
    double rel_l2 = 0.0;
    double rel_linf = 0.0;
    double max_u = 0.0;
    double q = 0.0;
    double conservation = 0.0;
    double residual_inf = 0.0;
    double order = std::numeric_limits<double>::quiet_NaN();
    bool converged = false;
};

struct Diagnostics {
    std::vector<std::string> failures;
    std::size_t checks = 0;
    std::size_t passed = 0;

    void check(bool condition, const std::string& name, const std::string& detail)
    {
        ++checks;
        if (condition) {
            ++passed;
            std::cout << "  PASS  " << name;
            if (!detail.empty()) std::cout << " — " << detail;
            std::cout << '\n';
        } else {
            failures.push_back(name + " — " + detail);
            std::cout << "  FAIL  " << name << " — " << detail << '\n';
        }
    }

    void exception(const std::string& name, const std::exception& e)
    {
        ++checks;
        failures.push_back(name + " — exception: " + e.what());
        std::cout << "  FAIL  " << name << " — exception: " << e.what() << '\n';
    }

    void diagnostic(bool condition, const std::string& name, const std::string& detail)
    {
        ++checks;
        if (condition) {
            ++passed;
            std::cout << "  DIAG-PASS  " << name;
        } else {
            std::cout << "  DIAG-FAIL  " << name;
        }
        if (!detail.empty()) std::cout << " — " << detail;
        std::cout << '\n';
    }
};

std::string sci(double x)
{
    std::ostringstream s;
    s << std::scientific << std::setprecision(6) << x;
    return s.str();
}

Result solve_level(std::size_t n, Diagnostics& d)
{
    constexpr double H = 1.0;
    constexpr double mu = 2.0;
    constexpr double G = 4.0;
    const auto mc = make_channel(n);
    const auto geometry = make_geometry_cache(mc.mesh);

    bool geometry_ok = true;
    double volume_sum = 0.0;
    double min_volume = std::numeric_limits<double>::infinity();
    double max_volume = 0.0;
    for (double v : geometry.cell_volumes) {
        geometry_ok = geometry_ok && std::isfinite(v) && v > 0.0;
        volume_sum += v;
        min_volume = std::min(min_volume, v);
        max_volume = std::max(max_volume, v);
    }
    d.check(geometry_ok, "geometry/positive-cell-volumes",
            "min=" + sci(min_volume) + " max=" + sci(max_volume));
    d.check(std::abs(volume_sum - 1.0) < 1e-13,
            "geometry/total-volume", "V=" + sci(volume_sum));

    bool finite_geometry = true;
    for (const auto& c : geometry.cell_centres)
        finite_geometry = finite_geometry && std::isfinite(c.x) &&
                          std::isfinite(c.y) && std::isfinite(c.z);
    for (const auto& s : geometry.face_Sf)
        finite_geometry = finite_geometry && std::isfinite(s.x) &&
                          std::isfinite(s.y) && std::isfinite(s.z);
    d.check(finite_geometry, "geometry/finite-centres-and-face-vectors", "");

    constexpr std::size_t expected_faces_per_cell = 6;
    bool topology_ok = true;
    for (std::size_t c = 0; c < mc.mesh.n_cells(); ++c)
        topology_ok = topology_ok &&
                      mc.mesh.cells().offsets_data()[c+1] -
                      mc.mesh.cells().offsets_data()[c] == expected_faces_per_cell;
    d.check(topology_ok, "topology/six-faces-per-hexa-cell", "");

    std::vector<double> source(n, G);
    DirichletBoundary bc;
    bc.face_values.assign(mc.mesh.n_faces(),
                           std::numeric_limits<double>::quiet_NaN());
    bc.face_values[0] = 0.0;
    bc.face_values[1] = 0.0;
    for (std::size_t f = 2; f < mc.mesh.n_faces(); ++f)
        bc.face_values[f] = std::numeric_limits<double>::quiet_NaN();

    ScalarDiffusionConfig cfg;
    cfg.diffusivity = mu;
    cfg.max_iterations = 5000;
    cfg.tolerance = 1e-13;

    Result r;
    r.n = n;
    const bool accuracy_gate = n >= 32;
    try {
        const auto result = solve_poisson_dirichlet(mc.mesh, bc, source, cfg);
        r.converged = result.linear_result.status == SolverStatus::CONVERGED;
        d.check(r.converged, "solver/converged",
                "status=" + std::to_string(static_cast<int>(result.linear_result.status)) +
                " iter=" + std::to_string(result.linear_result.iterations) +
                " reported_residual=" + sci(result.linear_result.residual));

        bool finite_solution = true;
        std::vector<double> numerical(n);
        std::vector<double> exact(n);
        for (std::size_t c = 0; c < n; ++c) {
            numerical[c] = result.solution(c);
            const double y = geometry.cell_centres[c].y;
            exact[c] = G*y*(H-y)/(2.0*mu);
            finite_solution = finite_solution && std::isfinite(numerical[c]);
        }
        d.check(finite_solution, "solution/all-finite", "");

        double max_u = -std::numeric_limits<double>::infinity();
        double q = 0.0;
        double min_u = std::numeric_limits<double>::infinity();
        double max_profile_symmetry = 0.0;
        for (std::size_t c = 0; c < n; ++c) {
            max_u = std::max(max_u, numerical[c]);
            min_u = std::min(min_u, numerical[c]);
            q += numerical[c] * geometry.cell_volumes[c];
            const std::size_t mirror = n - 1 - c;
            max_profile_symmetry =
                std::max(max_profile_symmetry, std::abs(numerical[c] - numerical[mirror]));
        }
        r.max_u = max_u;
        r.q = q;

        const auto e = error_norms(numerical, exact, geometry.cell_volumes);
        r.l2 = e.l2;
        r.linf = e.linf;
        r.rel_l2 = e.l2_relative;
        r.rel_linf = e.linf_relative;

        const std::string accuracy_tag = accuracy_gate ? "GATE " : "DIAGNOSTIC ";
        d.diagnostic(e.l2_relative < 5e-3, "accuracy/relative-L2",
                     accuracy_tag + "relL2=" + sci(e.l2_relative));
        d.diagnostic(e.linf_relative < 1.0e-2, "accuracy/relative-Linf",
                     accuracy_tag + "relLinf=" + sci(e.linf_relative));
        if (accuracy_gate) {
            d.check(e.l2_relative < 5e-3, "accuracy-gate/relative-L2",
                    "N>=32 relL2=" + sci(e.l2_relative));
            d.check(e.linf_relative < 1.0e-2, "accuracy-gate/relative-Linf",
                    "N>=32 relLinf=" + sci(e.linf_relative));
        }

        const double exact_umax = G*H*H/(8.0*mu);
        const double exact_q = G*H*H*H/(12.0*mu);
        const double rel_umax_error = std::abs(max_u - exact_umax) / exact_umax;
        const double rel_q_error = std::abs(q - exact_q) / exact_q;
        d.diagnostic(rel_umax_error < 5e-3, "QoI/maximum-velocity",
                     accuracy_tag + "relErr=" + sci(rel_umax_error) +
                     " num=" + sci(max_u) + " exact=" + sci(exact_umax));
        d.diagnostic(rel_q_error < 5e-3, "QoI/flow-rate-cell-centered",
                     accuracy_tag + "relErr=" + sci(rel_q_error) +
                     " num=" + sci(q) + " exact=" + sci(exact_q));
        if (accuracy_gate) {
            d.check(rel_umax_error < 5e-3, "accuracy-gate/maximum-velocity",
                    "relErr=" + sci(rel_umax_error));
            d.check(rel_q_error < 5e-3, "accuracy-gate/flow-rate-cell-centered",
                    "relErr=" + sci(rel_q_error));
        }

        const double h = H / static_cast<double>(n);
        const double expected_discrete_q = exact_q + h*h/3.0;
        d.diagnostic(std::abs(q - expected_discrete_q) < 1e-12,
                     "discrete-reference/cell-centered-Q",
                     "Qh=" + sci(q) + " expected=" + sci(expected_discrete_q) +
                     " formula=1/6+h^2/3");

        auto boundary_flux = [&](std::size_t face) {
            const std::size_t owner = mc.mesh.ownership().owner(face);
            const double dpf = (geometry.face_centres[face] - geometry.cell_centres[owner]).mag();
            return mu * (bc.face_values[face] - numerical[owner]) / dpf;
        };
        const double bottom_outward_flux = boundary_flux(0);
        const double top_outward_flux = boundary_flux(1);
        const double bottom_flow = -bottom_outward_flux;
        const double top_flow = -top_outward_flux;
        const double exact_wall_flow = G * H / 2.0;
        const double rel_bottom_flux_error = std::abs(bottom_flow - exact_wall_flow) / exact_wall_flow;
        const double rel_top_flux_error = std::abs(top_flow - exact_wall_flow) / exact_wall_flow;
        d.diagnostic(rel_bottom_flux_error < 5e-3, "QoI/wall-flux-bottom",
                     accuracy_tag + "relErr=" + sci(rel_bottom_flux_error) +
                     " num=" + sci(bottom_flow) + " exact=" + sci(exact_wall_flow));
        d.diagnostic(rel_top_flux_error < 5e-3, "QoI/wall-flux-top",
                     accuracy_tag + "relErr=" + sci(rel_top_flux_error) +
                     " num=" + sci(top_flow) + " exact=" + sci(exact_wall_flow));
        d.diagnostic(std::abs(bottom_flow - top_flow) / exact_wall_flow < 1e-12,
                     "QoI/wall-flux-symmetry",
                     "bottom=" + sci(bottom_flow) + " top=" + sci(top_flow));
        if (accuracy_gate) {
            d.check(rel_bottom_flux_error < 5e-3, "accuracy-gate/wall-flux-bottom",
                    "relErr=" + sci(rel_bottom_flux_error));
            d.check(rel_top_flux_error < 5e-3, "accuracy-gate/wall-flux-top",
                    "relErr=" + sci(rel_top_flux_error));
        }
        d.check(min_u >= -1e-12, "physics/non-negative-profile",
                "minU=" + sci(min_u));
        d.check(max_u > 0.0, "physics/positive-flow-direction",
                "maxU=" + sci(max_u));
        d.check(max_profile_symmetry < 5e-3,
                "profile/symmetry",
                "max|u(y)-u(H-y)|=" + sci(max_profile_symmetry));

        r.residual_inf = 0.0;
        for (std::size_t i = 0; i < result.solution.size(); ++i) {
            double ri = -result.rhs(i);
            const auto begin = result.matrix.row_offsets_data()[i];
            const auto end = result.matrix.row_offsets_data()[i+1];
            for (std::uint32_t k = begin; k < end; ++k)
                ri += result.matrix.values_data()[k] *
                       result.solution(result.matrix.columns_data()[k]);
            r.residual_inf = std::max(r.residual_inf, std::abs(ri));
        }
        d.check(std::isfinite(r.residual_inf) && r.residual_inf < 1e-10,
                "solver/independent-true-residual",
                "||Ax-b||inf=" + sci(r.residual_inf));

        double balance = 0.0;
        try {
            balance = poisson_conservation_balance(
                mc.mesh, bc, source, result.solution, geometry, mu);
            r.conservation = balance;
            d.check(std::isfinite(balance) && std::abs(balance) < 1e-10,
                    "conservation/integral-balance",
                    "balance=" + sci(balance));
        } catch (const std::exception& e) {
            d.exception("conservation/integral-balance", e);
        }

        bool matrix_ok = true;
        double symmetry_error = 0.0;
        double min_diag = std::numeric_limits<double>::infinity();
        double max_offdiag = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < result.matrix.n_rows(); ++i) {
            double diag = 0.0;
            for (std::uint32_t k = result.matrix.row_offsets_data()[i];
                 k < result.matrix.row_offsets_data()[i+1]; ++k) {
                const std::size_t j = result.matrix.columns_data()[k];
                const double aij = result.matrix.values_data()[k];
                if (i == j) diag += aij;
                else max_offdiag = std::max(max_offdiag, aij);
                symmetry_error =
                    std::max(symmetry_error, std::abs(aij - result.matrix(j,i)));
            }
            matrix_ok = matrix_ok && std::isfinite(diag) && diag > 0.0;
            min_diag = std::min(min_diag, diag);
        }
        matrix_ok = matrix_ok && symmetry_error < 1e-12 &&
                    max_offdiag <= 1e-14;
        d.check(matrix_ok, "assembly/SPD-structure",
                "minDiag=" + sci(min_diag) +
                " maxOffDiag=" + sci(max_offdiag) +
                " symmetry=" + sci(symmetry_error));

        return r;
    } catch (const std::exception& e) {
        d.exception("solver/execution", e);
        return r;
    }
}

} // namespace

int main(int argc, char** argv)
{
    const bool quick = argc == 2 && std::string(argv[1]) == "--quick";
    if (argc > 1 && !quick) {
        std::cerr << "usage: test_poiseuille_diagnostics [--quick]\n";
        return 2;
    }

    Diagnostics d;
    const std::vector<std::size_t> levels = quick
        ? std::vector<std::size_t>{32}
        : std::vector<std::size_t>{8, 16, 32, 64};
    std::vector<Result> results;

    std::cout << "=== CFDX Poiseuille Level-A diagnostic validation ===\n";
    std::cout << "Scope: scalar FV diffusion/source verification; not coupled NS.\n";
    std::cout << "All diagnostics run before the final CI decision.\n\n";

    for (const auto n : levels) {
        std::cout << "--- refinement level N=" << n << " ---\n";
        try {
            results.push_back(solve_level(n, d));
        } catch (const std::exception& e) {
            d.exception("level/N=" + std::to_string(n), e);
            results.push_back(Result{});
        }
        std::cout << '\n';
    }

    std::cout << "=== Refinement diagnostics ===\n";
    for (std::size_t i = 1; i < results.size(); ++i) {
        const auto& coarse = results[i-1];
        const auto& fine = results[i];
        const double p_l2 = observed_order(coarse.l2, fine.l2);
        const double qerr_coarse = std::abs(coarse.q - 1.0/6.0);
        const double qerr_fine = std::abs(fine.q - 1.0/6.0);
        const double p_q = observed_order(qerr_coarse, qerr_fine);
        results[i].order = p_l2;
        d.check(coarse.l2 > 0.0 && fine.l2 > 0.0 && std::isfinite(p_l2) && p_l2 >= 1.90,
                "refinement/L2-observed-order-" + std::to_string(coarse.n) +
                    "-to-" + std::to_string(fine.n),
                "p=" + sci(p_l2) + " required>=1.90");
        d.check(qerr_coarse > 0.0 && qerr_fine > 0.0 && std::isfinite(p_q) && p_q >= 1.90,
                "refinement/Q-observed-order-" + std::to_string(coarse.n) +
                    "-to-" + std::to_string(fine.n),
                "p=" + sci(p_q) + " required>=1.90");
    }

    if (results.size() >= 2) {
        const auto& first = results.front();
        const auto& last = results.back();
        d.check(last.l2 < first.l2,
                "refinement/L2-monotonic-decrease",
                "coarse=" + sci(first.l2) + " fine=" + sci(last.l2));
    }

    std::cout << "\n=== Final diagnostic summary ===\n";
    std::cout << "Checks: " << d.checks << "  PASS: " << d.passed
              << "  FAIL: " << d.failures.size() << '\n';

    if (!d.failures.empty()) {
        std::cout << "\nFAILURES (complete list):\n";
        for (const auto& failure : d.failures)
            std::cout << "  - " << failure << '\n';
        std::cout << "\nPOISEUILLE_VALIDATION: FAIL (decision made after all diagnostics)\n";
        return 1;
    }

    std::cout << (quick
        ? "\nPOISEUILLE_QUICK: PASS\n"
        : "\nPOISEUILLE_VALIDATION: PASS\n");
    return 0;
}
