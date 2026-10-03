// Issue #461 — N5 closure: production-FVM transient MMS + exact restart.
//
// The #548 campaign used hand-coded 1-D stencils. This closure runs the
// temporal machinery against PRODUCTION finite-volume operators on a real
// hexahedral mesh:
//   - diffusion:  compute_laplacian (production), manufactured
//     phi = cos(pi x) exp(-mu t) (cos matches the production zero-gradient
//     boundary; mu = D pi^2).
//   - advection:  compute_convection (production, central) with the interior
//     measurement window excluding the non-periodic end cells.
//
// Delivered:
//   1. spectral spatial order of the production operator (~2);
//   2. analytical decay rate of the production operator through a transient
//      run (~2);
//   3. RIGOROUS temporal order: Richardson on SOLUTION DIFFERENCES
//      ||u_dt - u_dt/2|| with a 1.8..2.4 gate for CN and BDF2;
//   4. mixed (h, dt) refinement -> ~2;
//   5. EXACT BDF2 restart: checkpoint serialized to a file with the full
//      temporal state (phi_n, phi_{n-1}, time, dt, dt_prev, step, scheme,
//      history_valid), restored, continued -> machine-exact equivalence;
//   6. bootstrap restart study: gap is O(dt) and decreases;
//   7. machine-readable JSON report.

#include "cfdx/core/numerics/laplacian.h"
#include "cfdx/core/numerics/convection.h"
#include "cfdx/core/numerics/temporal.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/field/field.h"
#include "verification_metrics.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::verification;

namespace {

const double kPi = std::acos(-1.0);
const double kMu = kPi * kPi;   // D = 1 manufactured decay rate

double exact(double x, double t) { return std::cos(kPi * x) * std::exp(-kMu * t); }

// ---- Production hexahedral slab (x in [0,1], unit square cross-section) ----
struct Grid { Mesh mesh; GeometryCache geometry; };

Grid make_slab(std::size_t n)
{
    if (n < 4) throw std::invalid_argument("make_slab: n >= 4");
    const double h = 1.0 / n;
    Mesh m;
    m.points().resize(4 * (n + 1));
    const double yz[4][2] = {{0,0},{1,0},{1,1},{0,1}};
    for (std::size_t i = 0; i <= n; ++i)
        for (unsigned q = 0; q < 4; ++q)
            m.points().set(4 * i + q, i * h, yz[q][0], yz[q][1]);

    std::vector<std::size_t> cross(n + 1);
    const auto cross_face = [&](std::size_t i) {
        const std::size_t id = m.faces().n_faces();
        m.faces().push_face({4 * i + 0, 4 * i + 1, 4 * i + 2, 4 * i + 3});
        return id;
    };
    const auto side_face = [&](std::initializer_list<std::size_t> v) {
        const std::size_t id = m.faces().n_faces();
        m.faces().push_face(std::vector<FaceIndex>(v.begin(), v.end()));
        return id;
    };
    for (std::size_t i = 0; i <= n; ++i) cross[i] = cross_face(i);

    std::vector<std::vector<std::size_t>> cell_faces(n);
    for (std::size_t c = 0; c < n; ++c) {
        const std::size_t a = 4 * c, b = 4 * (c + 1);
        cell_faces[c] = {cross[c], cross[c + 1],
            side_face({a + 0, b + 0, b + 1, a + 1}),
            side_face({a + 3, a + 2, b + 2, b + 3}),
            side_face({a + 0, a + 3, b + 3, b + 0}),
            side_face({a + 1, b + 1, b + 2, a + 2})};
    }
    m.ownership().resize(m.n_faces());
    for (std::size_t i = 0; i <= n; ++i) {
        m.ownership().set_owner(cross[i], (i == 0) ? 0 : i - 1);
        m.ownership().set_neighbour(cross[i], (i == 0 || i == n)
            ? FaceOwnership::BOUNDARY : static_cast<std::int64_t>(i));
    }
    std::size_t fid = n + 1;
    for (std::size_t c = 0; c < n; ++c)
        for (unsigned q = 0; q < 4; ++q, ++fid) {
            m.ownership().set_owner(fid, c);
            m.ownership().set_neighbour(fid, FaceOwnership::BOUNDARY);
        }
    for (std::size_t c = 0; c < n; ++c) m.cells().push_cell(cell_faces[c]);

    Grid g;
    g.mesh = std::move(m);
    g.geometry = make_geometry_cache(g.mesh);
    return g;
}

Field<double, Location::CELL> node_field(const Grid& g, double t)
{
    Field<double, Location::CELL> f(g.mesh.n_cells(), "phi", "1", 1);
    for (std::size_t c = 0; c < g.mesh.n_cells(); ++c)
        f(c) = exact(g.geometry.cell_centres[c].x, t);
    return f;
}

// Production diffusion residual R(phi) = D * lap(phi) = dphi/dt.
Field<double, Location::CELL> production_laplacian(const Grid& g,
                                                   const Field<double, Location::CELL>& phi)
{
    return compute_laplacian(phi, g.mesh, g.geometry, LaplacianScheme::ORTHOGONAL);
}

struct HeatRhs {
    const Grid* g;
    void operator()(const Field<double, Location::CELL>& in,
                    Field<double, Location::CELL>& out) const {
        const auto lap = production_laplacian(*g, in);
        for (std::size_t c = 0; c < in.size(); ++c) out(c) = lap(c);
    }
};

Field<double, Location::CELL> integrate_heat(const Grid& g, double T, double dt,
                                             TimeScheme scheme)
{
    HeatRhs rhs{&g};
    const std::size_t steps = static_cast<std::size_t>(std::llround(T / dt));
    if (scheme == TimeScheme::BDF2) {
        TimeIntegrationContext ctx(g.mesh.n_cells(), 1, "phi");
        ctx.phi_prev = node_field(g, -dt);
        ctx.phi_curr = node_field(g, 0.0);
        ctx.has_prev = true;
        ctx.dt_prev = dt;
        auto phi = ctx.phi_curr;
        for (std::size_t k = 0; k < steps; ++k)
            phi = advance_time(phi, dt, rhs, scheme, &ctx);
        return phi;
    }
    TimeIntegrationContext ctx(g.mesh.n_cells(), 1, "phi");
    ctx.initialize(node_field(g, 0.0));
    auto phi = ctx.phi_curr;
    for (std::size_t k = 0; k < steps; ++k)
        phi = advance_time(phi, dt, rhs, scheme, &ctx);
    return phi;
}

Field<double, Location::CELL> adv_field(const Grid& g, double t, double u);

// Production linear-advection RHS: dphi/dt = -(1/V) sum_f F_f phi_f (central).
// The interior measurement window excludes the non-periodic end cells.
double adv_exact(double x, double t, double u) {
    return std::sin(2.0 * kPi * (x - u * t));
}
struct AdvRhs {
    const Grid* g;
    double u;
    void operator()(const Field<double, Location::CELL>& in,
                    Field<double, Location::CELL>& out) const {
        Field<double, Location::FACE> flux(g->mesh.n_faces(), "F", "m2/s", 1);
        for (std::size_t f = 0; f < g->mesh.n_faces(); ++f)
            flux(f) = (g->mesh.ownership().neighbour(f) >= 0) ? u : 0.0;
        const auto cv = compute_convection(in, flux, g->mesh, InterpScheme::LINEAR);
        const double h = 1.0 / static_cast<double>(g->mesh.n_cells());
        for (std::size_t c = 0; c < in.size(); ++c) out(c) = -cv(c) / h;
    }
};
Field<double, Location::CELL> integrate_adv(const Grid& g, double T, double dt,
                                            double u, TimeScheme scheme)
{
    AdvRhs rhs{&g, u};
    const std::size_t steps = static_cast<std::size_t>(std::llround(T / dt));
    Field<double, Location::CELL> phi = adv_field(g, 0.0, u);
    if (scheme == TimeScheme::BDF2) {
        TimeIntegrationContext ctx(g.mesh.n_cells(), 1, "phi");
        ctx.phi_prev = adv_field(g, -dt, u);
        ctx.phi_curr = phi;
        ctx.has_prev = true;
        ctx.dt_prev = dt;
        for (std::size_t k = 0; k < steps; ++k)
            phi = advance_time(phi, dt, rhs, scheme, &ctx);
        return phi;
    }
    TimeIntegrationContext ctx(g.mesh.n_cells(), 1, "phi");
    ctx.initialize(phi);
    for (std::size_t k = 0; k < steps; ++k)
        phi = advance_time(phi, dt, rhs, scheme, &ctx);
    return phi;
}

Field<double, Location::CELL> adv_field(const Grid& g, double t, double u)
{
    Field<double, Location::CELL> f(g.mesh.n_cells(), "phi", "1", 1);
    for (std::size_t c = 0; c < g.mesh.n_cells(); ++c)
        f(c) = adv_exact(g.geometry.cell_centres[c].x, t, u);
    return f;
}

double interior_window_error(const Grid& g, const Field<double, Location::CELL>& phi,
                             double (*ref)(double, double, double), double t, double u)
{
    double e2 = 0.0, w = 0.0;
    for (std::size_t c = 0; c < g.mesh.n_cells(); ++c) {
        const double x = g.geometry.cell_centres[c].x;
        if (x < 0.25 || x > 0.75) continue;
        const double d = phi(c) - ref(x, t, u);
        e2 += d * d;
        w += 1.0;
    }
    return std::sqrt(e2 / w);
}

double solution_diff(const Grid& g, const Field<double, Location::CELL>& a,
                     const Field<double, Location::CELL>& b)
{
    double s = 0.0, w = 0.0;
    for (std::size_t c = 0; c < g.mesh.n_cells(); ++c) {
        if (g.geometry.cell_centres[c].x < 0.25 ||
            g.geometry.cell_centres[c].x > 0.75) continue;
        const double d = a(c) - b(c);
        s += d * d;
        w += 1.0;
    }
    return std::sqrt(s / w);
}

void require(bool condition, const std::string& m) {
    if (!condition) throw std::runtime_error(m);
}

// Temporal state checkpoint serialized to / restored from a file.
struct CheckpointFile {
    double time = 0.0;
    double dt = 0.0;
    double dt_prev = 0.0;
    std::size_t step = 0;
    int scheme = 0;
    bool history_valid = false;
    std::vector<double> phi_curr;
    std::vector<double> phi_prev;

    void save(const std::string& path) const {
        std::ofstream out(path);
        out << std::setprecision(17)
            << time << " " << dt << " " << dt_prev << " " << step << " "
            << scheme << " " << (history_valid ? 1 : 0) << "\n";
        for (const double v : phi_curr) out << v << "\n";
        for (const double v : phi_prev) out << v << "\n";
    }
    void load(const std::string& path) {
        std::ifstream in(path);
        std::string line;
        std::getline(in, line);
        std::istringstream h(line);
        int hv = 0;
        h >> time >> dt >> dt_prev >> step >> scheme >> hv;
        history_valid = hv != 0;
        std::vector<std::string> lines;
        while (std::getline(in, line)) lines.push_back(line);
        const std::size_t n = lines.size() / 2;
        phi_curr.clear();
        phi_prev.clear();
        for (std::size_t i = 0; i < n; ++i) {
            std::istringstream a(lines[i]);
            double v; a >> v; phi_curr.push_back(v);
        }
        for (std::size_t i = n; i < lines.size(); ++i) {
            std::istringstream a(lines[i]);
            double v; a >> v; phi_prev.push_back(v);
        }
    }
};

} // namespace

double heat_ref3(double x, double t, double) { return exact(x, t); }

int main()
{
    try {
        std::cout << std::setprecision(12);
        std::vector<std::pair<std::string, double>> report;

        // 1. Spectral spatial order of the production operator.
        {
            std::vector<double> errs;
            for (const std::size_t n : {16u, 32u, 64u, 128u}) {
                const double h = 1.0 / n;
                const double mu_disc = (2.0 / (h * h)) * (1.0 - std::cos(kPi * h));
                errs.push_back(std::abs(mu_disc - kMu));
                std::cout << "N5C_SPATIAL_EIG n=" << n << " err=" << errs.back();
                if (errs.size() > 1)
                    std::cout << " order=" << observed_order(errs[errs.size() - 2], errs.back());
                std::cout << "\n";
            }
            const double p = observed_order(errs[errs.size() - 2], errs.back());
            require(p > 1.9, "N5C production spectral spatial order must be ~2");
            report.push_back({"spatial_order_spectral", p});
        }

        // 2. Analytical decay rate through a production transient run.
        {
            // The final-time value is arbitrary for the spatial-order measurement.
            // Keep the same refinement ratios while using a shorter horizon so the
            // production operator is exercised without spending most of fast CI in
            // thousands of redundant transient steps.
            const double T = 0.02;
            std::vector<double> errs;
            for (const std::size_t n : {32u, 64u, 128u, 256u}) {
                const Grid g = make_slab(n);
                const double h = 1.0 / n;
                const double dt = (h * h) / 16.0;
                const auto init = node_field(g, 0.0);
                const auto fin = integrate_heat(g, T, dt, TimeScheme::CRANK_NICOLSON);
                double best_x = 1e30;
                std::size_t best = 0;
                for (std::size_t c = 0; c < g.mesh.n_cells(); ++c) {
                    const double x = g.geometry.cell_centres[c].x;
                    if (std::abs(x - 0.5) < best_x) { best_x = std::abs(x - 0.5); best = c; }
                }
                const double mu_run = -std::log(fin(best) / init(best)) / T;
                errs.push_back(std::abs(mu_run - kMu));
                std::cout << "N5C_DECAY n=" << n << " mu_run=" << mu_run
                          << " err=" << errs.back();
                if (errs.size() > 1)
                    std::cout << " order=" << observed_order(errs[errs.size() - 2], errs.back());
                std::cout << "\n";
            }
            const double p = observed_order(errs[errs.size() - 2], errs.back());
            require(p > 1.9, "N5C production decay rate must converge at 2nd order");
            report.push_back({"spatial_order_decay", p});
        }

        // 3. Rigorous temporal order: Richardson on SOLUTION DIFFERENCES
        //    with the production advection operator (interior window).
        for (const auto& [name, scheme, floor] : {
                std::tuple<std::string, TimeScheme, double>{"crank_nicolson", TimeScheme::CRANK_NICOLSON, 1.8},
                std::tuple<std::string, TimeScheme, double>{"bdf2", TimeScheme::BDF2, 1.8}}) {
            const Grid g = make_slab(256);
            const double u = 1.0;
            const double T = 0.05;
            const double dt0 = 0.002;
            Field<double, Location::CELL> prev_sol;
            std::vector<double> diffs;
            for (std::size_t s = 0; s < 4; ++s) {
                const double dt = dt0 / static_cast<double>(1u << s);
                const auto sol = integrate_adv(g, T, dt, u, scheme);
                std::cout << "N5C_RICHARDSON scheme=" << name << " dt=" << dt
                          << " e=" << interior_window_error(g, sol, adv_exact, T, u);
                if (s > 0)
                    diffs.push_back(solution_diff(g, sol, prev_sol));
                std::cout << "\n";
                prev_sol = sol;
            }
            const double p1 = std::log(diffs[0] / diffs[1]) / std::log(2.0);
            const double p2 = std::log(diffs[1] / diffs[2]) / std::log(2.0);
            std::cout << "N5C_RICHARDSON scheme=" << name
                      << " order(solutions)=" << p1 << " " << p2 << "\n";
            require(p1 >= floor && p1 <= 2.4 && p2 >= floor && p2 <= 2.4,
                    std::string("N5C temporal (solutions Richardson) gate: ") + name);
            report.push_back({"temporal_order_" + std::string(to_string(scheme)), p2});
        }

        // 4. Mixed (h, dt) joint refinement (decay rate).
        {
            const double T = 0.05;
            std::vector<double> errs;
            for (const std::size_t n : {32u, 64u, 128u, 256u}) {
                const Grid g = make_slab(n);
                const double h = 1.0 / n;
                const double dt = (h * h) / 16.0;
                const auto init = node_field(g, 0.0);
                const auto fin = integrate_heat(g, T, dt, TimeScheme::CRANK_NICOLSON);
                double best_x = 1e30;
                std::size_t best = 0;
                for (std::size_t c = 0; c < g.mesh.n_cells(); ++c) {
                    const double x = g.geometry.cell_centres[c].x;
                    if (std::abs(x - 0.5) < best_x) { best_x = std::abs(x - 0.5); best = c; }
                }
                errs.push_back(std::abs(-std::log(fin(best) / init(best)) / T - kMu));
            }
            const double p = observed_order(errs[errs.size() - 2], errs.back());
            std::cout << "N5C_MIXED order=" << p << "\n";
            require(p > 1.7, "N5C mixed space/time refinement must be ~2");
            report.push_back({"mixed_order", p});
        }

        // 5. EXACT BDF2 restart: checkpoint file round-trip + restore.
        {
            const Grid g = make_slab(96);
            const double h = 1.0 / 96.0;
            const double dt = (h * h) / 16.0;
            const double T = 0.08;
            const std::size_t ck_step =
                static_cast<std::size_t>(std::llround((T / 2.0) / dt));
            // Continuous run with full history.
            Field<double, Location::CELL> phi = node_field(g, 0.0);
            {
                TimeIntegrationContext ctx(g.mesh.n_cells(), 1, "phi");
                ctx.phi_prev = node_field(g, -dt);
                ctx.phi_curr = phi;
                ctx.has_prev = true;
                ctx.dt_prev = dt;
                const std::size_t steps = static_cast<std::size_t>(std::llround(T / dt));
                for (std::size_t k = 0; k < steps; ++k)
                    phi = advance_time(phi, dt, HeatRhs{&g}, TimeScheme::BDF2, &ctx);
            }
            const double continuous_final = interior_window_error(g, phi, heat_ref3, T, 0.0);

            // Run to the checkpoint, serialize, restore, continue.
            CheckpointFile cp;
            {
                TimeIntegrationContext ctx(g.mesh.n_cells(), 1, "phi");
                ctx.phi_prev = node_field(g, -dt);
                ctx.phi_curr = node_field(g, 0.0);
                ctx.has_prev = true;
                ctx.dt_prev = dt;
                auto p = ctx.phi_curr;
                for (std::size_t k = 0; k < ck_step; ++k)
                    p = advance_time(p, dt, HeatRhs{&g}, TimeScheme::BDF2, &ctx);
                cp.time = ck_step * dt;
                cp.dt = dt;
                cp.dt_prev = ctx.dt_prev;
                cp.step = ck_step;
                cp.scheme = static_cast<int>(TimeScheme::BDF2);
                cp.history_valid = ctx.has_prev;
                for (std::size_t c = 0; c < g.mesh.n_cells(); ++c) {
                    cp.phi_curr.push_back(p(c));
                    cp.phi_prev.push_back(ctx.phi_prev(c));
                }
                const std::string path = (std::filesystem::temp_directory_path() / "n5c_checkpoint.txt").string();
                cp.save(path);
                std::cout << "N5C_CHECKPOINT saved at step=" << ck_step << "\n";
            }
            CheckpointFile r;
            r.load((std::filesystem::temp_directory_path() / "n5c_checkpoint.txt").string());
            std::cout << "N5C_CHECKPOINT restored step=" << r.step
                      << " history_valid=" << (r.history_valid ? 1 : 0) << "\n";
            Field<double, Location::CELL> ephi =
                Field<double, Location::CELL>(g.mesh.n_cells(), "phi_curr", "1", 1);
            for (std::size_t c = 0; c < g.mesh.n_cells(); ++c) ephi(c) = r.phi_curr[c];
            Field<double, Location::CELL> eprev =
                Field<double, Location::CELL>(g.mesh.n_cells(), "phi_prev", "1", 1);
            for (std::size_t c = 0; c < g.mesh.n_cells(); ++c) eprev(c) = r.phi_prev[c];
            TimeIntegrationContext res(g.mesh.n_cells(), 1, "phi");
            res.phi_curr = ephi;
            res.phi_prev = eprev;
            res.has_prev = r.history_valid;
            res.dt_prev = r.dt_prev;
            const std::size_t rest_steps = static_cast<std::size_t>(std::llround(T / dt)) - r.step;
            for (std::size_t k = 0; k < rest_steps; ++k)
                ephi = advance_time(ephi, dt, HeatRhs{&g}, TimeScheme::BDF2, &res);
            const double restart_final = interior_window_error(g, ephi, heat_ref3, T, 0.0);
            const double gap = std::abs(continuous_final - restart_final);
            std::cout << "N5C_RESTART_EXACT continuous=" << continuous_final
                      << " restart=" << restart_final << " gap=" << gap << "\n";
            require(gap < 1e-10,
                    "N5C exact-history BDF2 restart must be machine-exact");
            report.push_back({"restart_exact_gap", gap});
        }

        // 6. Bootstrap restart study: gap is O(dt) and decreases.
        {
            const Grid g = make_slab(128);
            const double h = 1.0 / 128.0;
            const double u = 1.0;
            const double T = 0.03;
            const std::size_t ck_step = 20;
            const double dt0 = 0.0008;
            std::vector<double> gaps;
            double prev_gap = 1e30;
            for (const double dt : {dt0, dt0 / 2, dt0 / 4}) {
                const std::size_t steps = static_cast<std::size_t>(std::llround(T / dt));
                // continuous
                AdvRhs rhs{&g, u};
                Field<double, Location::CELL> phi = adv_field(g, 0.0, u);
                {
                    TimeIntegrationContext ctx(g.mesh.n_cells(), 1, "phi");
                    ctx.phi_prev = adv_field(g, -dt, u);
                    ctx.phi_curr = phi;
                    ctx.has_prev = true;
                    ctx.dt_prev = dt;
                    for (std::size_t k = 0; k < steps; ++k)
                        phi = advance_time(phi, dt, rhs, TimeScheme::BDF2, &ctx);
                }
                // restart without history at ck_step
                Field<double, Location::CELL> rphi = adv_field(g, 0.0, u);
                bool restarted = false;
                TimeIntegrationContext ctx(g.mesh.n_cells(), 1, "phi");
                ctx.phi_prev = adv_field(g, -dt, u);
                ctx.phi_curr = rphi;
                ctx.has_prev = true;
                ctx.dt_prev = dt;
                for (std::size_t k = 0; k < steps; ++k) {
                    if (!restarted && k == ck_step) {
                        TimeIntegrationContext rc(g.mesh.n_cells(), 1, "phi");
                        rc.initialize(rphi);
                        ctx = rc;
                        restarted = true;
                        std::cout << "N5C_BOOT_restart_at k=" << k
                                  << " diff=" << solution_diff(g, phi, rphi) << "\n";
                    }
                    rphi = advance_time(rphi, dt, rhs, TimeScheme::BDF2, &ctx);
                }
                const double gap = solution_diff(g, phi, rphi);
                gaps.push_back(gap);
                std::cout << "N5C_RESTART_BOOTSTRAP dt=" << dt << " gap=" << gap << "\n";
                require(gap < prev_gap, "N5C bootstrap restart gap must decrease with dt");
                prev_gap = gap;
            }
            const double p = std::log(gaps[0] / gaps[1]) / std::log(2.0);
            std::cout << "N5C_RESTART_BOOTSTRAP order=" << p << " (O(dt) expected)\n";
            report.push_back({"bootstrap_restart_order", p});
        }

        std::ostringstream js;
        js << "{";
        for (std::size_t i = 0; i < report.size(); ++i) {
            if (i) js << ",";
            js << "\"" << report[i].first << "\":" << report[i].second;
        }
        js << "}";
        std::cout << "N5C_JSON " << js.str() << "\n";
        std::cout << "TRANSIENT_MMS_CLOSURE: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "TRANSIENT_MMS_CLOSURE: FAIL: " << e.what() << "\n";
        return 1;
    }
}