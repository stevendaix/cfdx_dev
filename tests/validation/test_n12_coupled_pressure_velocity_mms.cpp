// N12 — coupled pressure/velocity manufactured-solution laboratory.
//
// Solver-level oracle using the fully coupled pressure-based path. The exact
// manufactured state is U=(1,0,0), p=x on [0,1]^3, with rho=1, nu=1e-3 and
// a uniform body force f=(1,0,0). The momentum equation is therefore exactly
// balanced: -grad(p)+f=0, while div(U)=0.
//
// Pressure is fixed to the analytical values p=0 at x=0 and p=1 at x=1;
// velocity is fixed to the exact constant state on every boundary patch.
// This deliberately exercises the coupled pressure/velocity solver rather
// than reconstructing an operator-level pressure oracle. The test reports
// field errors, authoritative continuity, momentum-equation residuals and
// finite/non-finite diagnostics. No fallback algorithm is accepted.

#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {
Mesh make_channel(std::size_t n)
{
    if (n == 0) throw std::invalid_argument("n must be positive");
    Mesh m;
    m.points().resize(8 * n);
    const double dy = 1.0 / static_cast<double>(n);
    for (std::size_t j = 0; j < n; ++j) {
        const double y0 = j * dy;
        const double y1 = (j + 1) * dy;
        const std::size_t b = 8 * j;
        const double p[8][3] = {
            {0.0,y0,0.0}, {1.0,y0,0.0}, {1.0,y1,0.0}, {0.0,y1,0.0},
            {0.0,y0,1.0}, {1.0,y0,1.0}, {1.0,y1,1.0}, {0.0,y1,1.0}
        };
        for (std::size_t q = 0; q < 8; ++q)
            m.points().set(b + q, p[q][0], p[q][1], p[q][2]);
    }

    std::vector<std::size_t> bottom, top, x0, x1, z0, z1, internal;
    auto add_face = [&](std::initializer_list<std::size_t> v) {
        const auto id = m.faces().n_faces();
        m.faces().push_face(std::vector<FaceIndex>(v.begin(), v.end()));
        return id;
    };
    bottom.push_back(add_face({0,1,5,4}));
    top.push_back(add_face({8*(n-1)+3,8*(n-1)+7,8*(n-1)+6,8*(n-1)+2}));
    for (std::size_t j = 0; j < n; ++j) {
        const auto b = 8 * j;
        x0.push_back(add_face({b,b+4,b+7,b+3}));
        x1.push_back(add_face({b+1,b+2,b+6,b+5}));
        z0.push_back(add_face({b,b+3,b+2,b+1}));
        z1.push_back(add_face({b+4,b+5,b+6,b+7}));
    }
    for (std::size_t j = 0; j + 1 < n; ++j) {
        const auto b = 8 * j;
        internal.push_back(add_face({b+3,b+7,b+6,b+2}));
    }

    m.ownership().resize(m.n_faces());
    m.ownership().set_owner(bottom[0], 0);
    m.ownership().set_neighbour(bottom[0], FaceOwnership::BOUNDARY);
    m.ownership().set_owner(top[0], n - 1);
    m.ownership().set_neighbour(top[0], FaceOwnership::BOUNDARY);
    for (std::size_t j = 0; j < n; ++j) {
        for (const auto f : {x0[j], x1[j], z0[j], z1[j]}) {
            m.ownership().set_owner(f, j);
            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        }
    }
    for (std::size_t j = 0; j + 1 < n; ++j) {
        m.ownership().set_owner(internal[j], j);
        m.ownership().set_neighbour(internal[j], static_cast<std::int64_t>(j + 1));
    }
    for (std::size_t j = 0; j < n; ++j) {
        m.cells().push_cell({
            j == 0 ? bottom[0] : internal[j-1],
            j + 1 == n ? top[0] : internal[j],
            x0[j], x1[j], z0[j], z1[j]});
    }
    auto patch = [&](const char* name, const std::vector<std::size_t>& faces) {
        Patch p;
        p.name = name;
        p.type = PatchType::WALL;
        p.face_ids = faces;
        m.boundary().add_patch(p);
    };
    patch("bottom", bottom);
    patch("top", top);
    patch("x0", x0);
    patch("x1", x1);
    patch("z0", z0);
    patch("z1", z1);
    return m;
}

struct Metrics {
    double u_linf = 0.0;
    double p_linf = 0.0;
    double continuity_linf = std::numeric_limits<double>::infinity();
    double momentum_rel = std::numeric_limits<double>::infinity();
    double mass_balance = std::numeric_limits<double>::infinity();
    std::size_t iterations = 0;
};

Metrics run_case(std::size_t n)
{
    const Mesh mesh = make_channel(n);
    Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
    Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
    U.fill(0.0);
    p.fill(0.0);

    VelocityBoundaryConditions ubc;
    const Vec3 exact_u{1.0, 0.0, 0.0};
    for (const char* name : {"bottom", "top", "x0", "x1", "z0", "z1"})
        ubc[name] = {VelocityBoundaryCondition::Type::FIXED_VALUE, exact_u};

    ScalarBoundaryConditions pbc;
    pbc["bottom"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    pbc["top"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    pbc["z0"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    pbc["z1"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    pbc["x0"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
    pbc["x1"] = {ScalarBoundaryType::FIXED_VALUE, 1.0, 0.0};

    IncompressibleSolverControls c;
    c.algorithm = PressureVelocityAlgorithm::COUPLED;
    c.coupling.schur_model = CoupledSchurModel::BlockLocal;
    c.coupling.coupled_max_iterations = 500;
    c.coupling.coupled_linear_tolerance = 1e-12;
    c.convergence.max_iterations = 20;
    c.convergence.relative_tolerance = 1e-10;
    c.convergence.continuity_tolerance = 1e-12;
    c.linear_tolerance = 1e-12;
    c.density = 1.0;
    c.kinematic_viscosity = 1.0e-3;
    c.body_force = Vec3{1.0, 0.0, 0.0};
    c.pressure_gauge_policy = PressureGaugePolicy::REFERENCE_CELL;
    c.pressure_reference_cell = 0;
    c.pressure_reference_value = 0.5;
    c.diagnostics.iteration_trace = false;

    const auto result = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
    if (!result.converged)
        throw std::runtime_error("coupled pressure-velocity MMS did not converge: " + result.convergence_reason);
    if (!result.coupled_linear_plan_resolved)
        throw std::runtime_error("coupled pressure-velocity MMS did not resolve a coupled linear plan");
    if (result.resolved_coupled_schur_model != CoupledSchurModel::BlockLocal)
        throw std::runtime_error("coupled pressure-velocity MMS resolved an unexpected Schur model");
    if (result.authoritative_mass_flux.size() != mesh.n_faces())
        throw std::runtime_error("coupled pressure-velocity MMS has no authoritative mass flux evidence");

    Metrics out;
    out.iterations = result.iterations;
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        const auto geometry = build_fv_geometry(mesh);
        const double x = geometry.cell_centres[cell].x;
        out.u_linf = std::max(out.u_linf,
            std::max({std::abs(U(cell,0) - 1.0), std::abs(U(cell,1)), std::abs(U(cell,2))}));
        out.p_linf = std::max(out.p_linf, std::abs(p(cell) - x));
    }
    const auto& h = result.history.back();
    out.continuity_linf = h.corrected_flux_continuity_linf;
    out.momentum_rel = h.momentum_equation_residual_relative;
    out.mass_balance = h.mass_normalized_imbalance;

    std::cout << std::setprecision(12)
              << "N12_COUPLED_PV_MMS n=" << n
              << " U_Linf=" << out.u_linf
              << " p_Linf=" << out.p_linf
              << " continuity_Linf=" << out.continuity_linf
              << " momentum_rel=" << out.momentum_rel
              << " mass_normalized=" << out.mass_balance
              << " iterations=" << out.iterations
              << " schur=block_local\n";

    if (!std::isfinite(out.u_linf) || !std::isfinite(out.p_linf) ||
        !std::isfinite(out.continuity_linf) || !std::isfinite(out.momentum_rel) ||
        !std::isfinite(out.mass_balance))
        throw std::runtime_error("coupled pressure-velocity MMS produced non-finite evidence");
    if (h.mass_nonfinite_faces != 0 || h.boundedness_nonfinite_velocity != 0)
        throw std::runtime_error("coupled pressure-velocity MMS produced non-finite solver state");
    return out;
}
} // namespace

int main()
{
    try {
        const auto m16 = run_case(16);
        const auto m32 = run_case(32);
        const auto m64 = run_case(64);

        const double u_tol = 1e-10;
        const double p_tol = 1e-10;
        const double continuity_tol = 1e-10;
        const double momentum_tol = 1e-10;
        const double mass_tol = 1e-10;
        for (const auto& m : {m16, m32, m64}) {
            if (!(m.u_linf < u_tol && m.p_linf < p_tol &&
                  m.continuity_linf < continuity_tol &&
                  m.momentum_rel < momentum_tol && m.mass_balance < mass_tol))
                throw std::runtime_error("N12 coupled pressure-velocity MMS acceptance criteria failed");
        }

        std::cout << "N12_COUPLED_PV_MMS: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "N12_COUPLED_PV_MMS: FAIL: " << e.what() << "\n";
        return 1;
    }
}
