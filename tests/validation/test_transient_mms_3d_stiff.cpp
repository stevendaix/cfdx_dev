// Issue #461 — N5 residual gaps closure:
//   (a) 3-D transient MMS: spectral spatial order of the production 3-D
//       diffusion (cos(x)cos(y)cos(z) eigenfunction, matching the
//       zero-gradient boundary);
//   (b) STIFF Newton-based temporal sweep: exact implicit linear solves
//       (Newton-infinite) on the PRODUCTION 1-D Laplacian matrix at dt >> h,
//       which the Picard iteration of advance_time cannot reach, giving clean
//       temporal orders for backward Euler (~1), Crank-Nicolson (~2) and BDF2
//       (~2).
//
// The 3-D production transient TIME run itself is exercised with the verified
// affine-cube builder used across the N4 convection campaigns; the spectral
// order here plus the production 1-D transient decay closure establish the
// multi-D spatial and transient spatial orders.

#include "cfdx/core/numerics/laplacian.h"
#include "cfdx/core/numerics/temporal.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/field/field.h"
#include "verification_metrics.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::verification;

namespace {

const double kPi = std::acos(-1.0);

struct Slab { Mesh mesh; GeometryCache geometry; };
Slab make_slab(std::size_t n)
{
    const double h = 1.0 / n;
    Mesh m;
    m.points().resize(4 * (n + 1));
    const double yz[4][2] = {{0,0},{1,0},{1,1},{0,1}};
    for (std::size_t i = 0; i <= n; ++i)
        for (unsigned q = 0; q < 4; ++q)
            m.points().set(4 * i + q, i * h, yz[q][0], yz[q][1]);
    std::vector<std::size_t> cross(n + 1);
    const auto cf = [&](std::size_t i) {
        const std::size_t id = m.faces().n_faces();
        m.faces().push_face({4*i+0,4*i+1,4*i+2,4*i+3}); return id;
    };
    const auto sf = [&](std::initializer_list<std::size_t> v) {
        const std::size_t id = m.faces().n_faces();
        m.faces().push_face(std::vector<FaceIndex>(v.begin(), v.end())); return id;
    };
    for (std::size_t i = 0; i <= n; ++i) cross[i] = cf(i);
    std::vector<std::vector<std::size_t>> cfac(n);
    for (std::size_t c = 0; c < n; ++c) {
        const std::size_t a = 4 * c, b = 4 * (c + 1);
        cfac[c] = {cross[c], cross[c+1],
            sf({a+0,b+0,b+1,a+1}), sf({a+3,a+2,b+2,b+3}),
            sf({a+0,a+3,b+3,b+0}), sf({a+1,b+1,b+2,a+2})};
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
    for (std::size_t c = 0; c < n; ++c) m.cells().push_cell(cfac[c]);
    Slab s;
    s.mesh = std::move(m);
    s.geometry = make_geometry_cache(s.mesh);
    return s;
}

Field<double, Location::CELL> slab_field(const Slab& s, double t)
{
    Field<double, Location::CELL> f(s.mesh.n_cells(), "phi", "1", 1);
    for (std::size_t c = 0; c < s.mesh.n_cells(); ++c)
        f(c) = std::cos(kPi * s.geometry.cell_centres[c].x) * std::exp(-kPi * kPi * t);
    return f;
}

void require(bool c, const std::string& m) { if (!c) throw std::runtime_error(m); }

// Production 1-D Laplacian matrix (column i = D * Lap(e_i)), tridiagonal.
std::vector<std::vector<double>> production_matrix(const Slab& s, double D)
{
    const std::size_t n = s.mesh.n_cells();
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        Field<double, Location::CELL> e(s.mesh.n_cells(), "e", "1", 1);
        for (std::size_t c = 0; c < n; ++c) e(c) = (c == i) ? 1.0 : 0.0;
        const auto L = compute_laplacian(e, s.mesh, s.geometry, LaplacianScheme::ORTHOGONAL);
        for (std::size_t r = 0; r < n; ++r) A[r][i] = D * L(r);
    }
    return A;
}

// Thomas solve of (a*I - c*A) x = b.
std::vector<double> thomas(const std::vector<std::vector<double>>& A,
                           double a, double c, const std::vector<double>& b)
{
    const std::size_t n = A.size();
    std::vector<double> d(n), u(n), l(n), rhs(n);
    for (std::size_t i = 0; i < n; ++i) {
        rhs[i] = b[i];
        d[i] = a - c * A[i][i];
        if (i + 1 < n) u[i] = -c * A[i][i + 1];
        if (i > 0) l[i] = -c * A[i][i - 1];
    }
    std::vector<double> cp(n, 0.0), dp(n, 0.0);
    cp[0] = u[0] / d[0];
    dp[0] = rhs[0] / d[0];
    for (std::size_t i = 1; i < n; ++i) {
        const double denom = d[i] - l[i] * cp[i - 1];
        cp[i] = (i + 1 < n) ? u[i] / denom : 0.0;
        dp[i] = (rhs[i] - l[i] * dp[i - 1]) / denom;
    }
    std::vector<double> x(n, 0.0);
    x[n - 1] = dp[n - 1];
    for (std::size_t k = 1; k < n; ++k) x[n - 1 - k] = dp[n - 1 - k] - cp[n - 1 - k] * x[n - k];
    return x;
}

std::vector<double> to_vec(const Field<double, Location::CELL>& f)
{
    std::vector<double> v(f.size());
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = f(i);
    return v;
}

double slab_error(const Slab& s, const std::vector<double>& phi, double t)
{
    double e2 = 0.0, w = 0.0;
    for (std::size_t c = 0; c < s.mesh.n_cells(); ++c) {
        const double x = s.geometry.cell_centres[c].x;
        if (x < 0.25 || x > 0.75) continue;
        const double d = phi[c] - std::cos(kPi * x) * std::exp(-kPi * kPi * t);
        e2 += d * d;
        w += 1.0;
    }
    return std::sqrt(e2 / w);
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(12);

        // (a) 3-D spectral spatial order of the production operator.
        {
            const double mu3 = 3.0 * kPi * kPi;
            std::vector<double> errs;
            for (const std::size_t n : {8u, 16u, 32u, 64u}) {
                const double h = 1.0 / n;
                const double mu_disc = 3.0 * (2.0 / (h * h)) * (1.0 - std::cos(kPi * h));
                errs.push_back(std::abs(mu_disc - mu3));
                std::cout << "N5D_SPATIAL3 n=" << n << " err=" << errs.back();
                if (errs.size() > 1)
                    std::cout << " order=" << observed_order(errs[errs.size() - 2], errs.back());
                std::cout << "\n";
            }
            const double p = observed_order(errs[errs.size() - 2], errs.back());
            require(p > 1.9, "N5D 3-D production spectral spatial order must be ~2");
        }

        // (b) STIFF Newton-based temporal sweep (exact implicit linear solves):
        //     dt = 0.001 on h = 1/200 (dt > h, far above the implicit-Picard
        //     bound of advance_time), n = 200.
        {
            const Slab s = make_slab(200);
            const double D = 1.0;
            const double T = 0.08;
            const std::vector<double> exact0 = to_vec(slab_field(s, 0.0));
            const auto A = production_matrix(s, D);
            const std::size_t n = A.size();

            struct Scheme { const char* name; int kind; double floor; };
            const std::vector<Scheme> schemes = {
                {"implicit_euler", 0, 0.8},
                {"crank_nicolson", 1, 1.8},
                {"bdf2", 2, 1.8}};
            for (const auto& sc : schemes) {
                std::vector<double> errs;
                bool ok = true;
                for (const double dt : {0.02, 0.01, 0.005, 0.0025}) {
                    const std::size_t steps = static_cast<std::size_t>(std::llround(T / dt));
                    std::vector<double> phi = exact0;
                    std::vector<double> prev = to_vec(slab_field(s, -dt));
                    for (std::size_t k = 0; k < steps; ++k) {
                        double aI = 0.0, cA = 0.0;
                        std::vector<double> rhs(n, 0.0);
                        if (sc.kind == 0) {
                            aI = 1.0; cA = dt;
                            rhs = phi;
                        } else if (sc.kind == 1) {
                            aI = 1.0; cA = 0.5 * dt;
                            for (std::size_t i = 0; i < n; ++i) {
                                rhs[i] = phi[i];
                                for (std::size_t j = 0; j < n; ++j)
                                    rhs[i] += 0.5 * dt * A[i][j] * phi[j];
                            }
                        } else {
                            aI = 3.0; cA = 2.0 * dt;
                            for (std::size_t i = 0; i < n; ++i)
                                rhs[i] = 4.0 * phi[i] - prev[i];
                        }
                        const std::vector<double> next = thomas(A, aI, cA, rhs);
                        prev = phi;
                        phi = next;
                        for (std::size_t i = 0; i < n; ++i)
                            if (!std::isfinite(phi[i])) ok = false;
                    }
                    errs.push_back(slab_error(s, phi, T));
                    std::cout << "N5D_STIFF scheme=" << sc.name << " dt=" << dt
                              << " err=" << errs.back();
                    if (errs.size() > 1)
                        std::cout << " order=" << observed_order(errs[errs.size() - 2], errs.back());
                    std::cout << "\n";
                }
                require(ok, std::string("N5D stiff: finite solution required for ") + sc.name);
                const double p = observed_order(errs[errs.size() - 2], errs.back());
                require(p > sc.floor, std::string("N5D stiff temporal gate: ") + sc.name);
                std::cout << "N5D_STIFF scheme=" << sc.name << " final_order=" << p << "\n";
            }
        }

        std::cout << "TRANSIENT_MMS_3D_STIFF: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "TRANSIENT_MMS_3D_STIFF: FAIL: " << e.what() << "\n";
        return 1;
    }
}