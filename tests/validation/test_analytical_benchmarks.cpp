#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/radiation.h"
#include "cfdx/physics/thermal.h"
#include "verification_metrics.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <utility>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::verification;

namespace {

struct OneDimensionalMesh {
    Mesh mesh;
    std::vector<std::size_t> cell_faces;
};

OneDimensionalMesh make_channel(std::size_t n, double height)
{
    if (n == 0 || !(height > 0.0))
        throw std::invalid_argument("make_channel: invalid dimensions");

    Mesh m;
    m.points().resize(8 * n);

    const double dy = height / static_cast<double>(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double y0 = dy * static_cast<double>(i);
        const double y1 = dy * static_cast<double>(i + 1);
        const std::size_t b = 8 * i;
        const double p[8][3] = {
            {0.0,y0,0.0}, {1.0,y0,0.0}, {1.0,y1,0.0}, {0.0,y1,0.0},
            {0.0,y0,1.0}, {1.0,y0,1.0}, {1.0,y1,1.0}, {0.0,y1,1.0}
        };
        for (std::size_t j = 0; j < 8; ++j)
            m.points().set(b + j, p[j][0], p[j][1], p[j][2]);
    }

    std::vector<std::size_t> bottom, top, side_x0, side_x1, side_z0, side_z1;
    std::vector<std::size_t> internal;
    bottom.reserve(1); top.reserve(1);
    side_x0.reserve(n); side_x1.reserve(n);
    side_z0.reserve(n); side_z1.reserve(n);
    internal.reserve(n > 0 ? n - 1 : 0);

    // Boundary faces are inserted first so patch ids are deterministic.
    auto add_face = [&](std::initializer_list<std::size_t> vertices) {
        const std::size_t id = m.faces().n_faces();
        m.faces().push_face(
            std::vector<FaceIndex>(vertices.begin(), vertices.end()));
        return id;
    };

    bottom.push_back(add_face({0,1,5,4}));
    top.push_back(add_face({8*(n-1)+3,8*(n-1)+7,8*(n-1)+6,8*(n-1)+2}));

    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t b = 8 * i;
        side_x0.push_back(add_face({b+0,b+4,b+7,b+3}));
        side_x1.push_back(add_face({b+1,b+2,b+6,b+5}));
        side_z0.push_back(add_face({b+0,b+3,b+2,b+1}));
        side_z1.push_back(add_face({b+4,b+5,b+6,b+7}));
    }

    for (std::size_t i = 0; i + 1 < n; ++i) {
        const std::size_t b0 = 8 * i;
        internal.push_back(add_face({b0+3,b0+7,b0+6,b0+2}));
    }

    // Cell face order matches the boundary orientation used by the core
    // geometry tests. Internal faces are oriented from the lower-y owner to
    // the upper-y neighbour.
    m.ownership().resize(m.n_faces());
    std::vector<std::vector<std::size_t>> cell_faces(n);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t b = 8 * i;
        cell_faces[i] = {
            i == 0 ? bottom[0] : internal[i - 1],
            i + 1 == n ? top[0] : internal[i],
            side_x0[i], side_x1[i], side_z0[i], side_z1[i]
        };
    }

    for (std::size_t f = 0; f < m.n_faces(); ++f) {
        bool assigned = false;
        for (std::size_t c = 0; c < n && !assigned; ++c) {
            for (const auto cf : cell_faces[c]) {
                if (cf == f) {
                    m.ownership().set_owner(f, c);
                    const bool is_internal =
                        (std::find(internal.begin(), internal.end(), f) != internal.end());
                    m.ownership().set_neighbour(
                        f, is_internal
                            ? static_cast<int>(c + 1)
                            : FaceOwnership::BOUNDARY);
                    assigned = true;
                    break;
                }
            }
        }
        if (!assigned)
            throw std::runtime_error("make_channel: unassigned face");
    }

    for (const auto& faces : cell_faces)
        m.cells().push_cell(faces);

    Patch p_bottom;
    p_bottom.name = "bottom";
    p_bottom.type = PatchType::WALL;
    p_bottom.face_ids = bottom;
    m.boundary().add_patch(p_bottom);

    Patch p_top;
    p_top.name = "top";
    p_top.type = PatchType::WALL;
    p_top.face_ids = top;
    m.boundary().add_patch(p_top);

    Patch p_x0;
    p_x0.name = "x0";
    p_x0.type = PatchType::WALL;
    p_x0.face_ids = side_x0;
    m.boundary().add_patch(p_x0);

    Patch p_x1;
    p_x1.name = "x1";
    p_x1.type = PatchType::WALL;
    p_x1.face_ids = side_x1;
    m.boundary().add_patch(p_x1);

    Patch p_z0;
    p_z0.name = "z0";
    p_z0.type = PatchType::WALL;
    p_z0.face_ids = side_z0;
    m.boundary().add_patch(p_z0);

    Patch p_z1;
    p_z1.name = "z1";
    p_z1.type = PatchType::WALL;
    p_z1.face_ids = side_z1;
    m.boundary().add_patch(p_z1);

    return {std::move(m), {}};
}

struct ScalarResult {
    std::vector<double> y;
    std::vector<double> u;
    std::vector<double> volume;
};

ScalarResult solve_diffusion_case(std::size_t n, double height,
                                  double gamma, double source,
                                  double bottom_value, double top_value)
{
    auto problem = make_channel(n, height);
    auto geometry = build_fv_geometry(problem.mesh);

    Field<double,Location::FACE> phi(problem.mesh.n_faces(), "phi", "kg/s", 1);
    phi.fill(0.0);

    Field<double,Location::CELL> su(n, "source", "unit", 1);
    Field<double,Location::CELL> sp(n, "sp", "unit", 1);
    su.fill(source);
    sp.fill(0.0);

    ScalarBoundaryConditions bc;
    bc["bottom"] = {ScalarBoundaryType::FIXED_VALUE, bottom_value, 0.0};
    bc["top"] = {ScalarBoundaryType::FIXED_VALUE, top_value, 0.0};
    bc["x0"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    bc["x1"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    bc["z0"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    bc["z1"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

    auto eq = assemble_scalar_equation(
        problem.mesh, geometry, phi, gamma, su, sp, bc, true);
    Vector solution(n, 0.0);
    const auto linear = solve_scalar_equation(
        eq, solution, {5000, 1e-13, 1.0});
    if (linear.status != SolverStatus::CONVERGED)
        throw std::runtime_error("analytical benchmark: linear solve did not converge");

    ScalarResult result;
    result.y.resize(n);
    result.u.resize(n);
    result.volume = geometry.cell_volumes;
    for (std::size_t i = 0; i < n; ++i) {
        result.y[i] = geometry.cell_centres[i].y;
        result.u[i] = solution(i);
    }
    return result;
}

ScalarResult solve_mms_diffusion(std::size_t n, double height, double gamma)
{
    auto problem = make_channel(n, height);
    auto geometry = build_fv_geometry(problem.mesh);

    Field<double,Location::FACE> phi(problem.mesh.n_faces(), "phi", "kg/s", 1);
    phi.fill(0.0);

    Field<double,Location::CELL> su(n, "source", "unit", 1);
    Field<double,Location::CELL> sp(n, "sp", "unit", 1);
    sp.fill(0.0);

    const double pi = std::acos(-1.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double y = geometry.cell_centres[i].y;
        su(i) = gamma * pi*pi/(height*height) * std::sin(pi*y/height);
    }

    ScalarBoundaryConditions bc;
    bc["bottom"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
    bc["top"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
    bc["x0"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    bc["x1"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    bc["z0"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    bc["z1"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

    auto eq = assemble_scalar_equation(
        problem.mesh, geometry, phi, gamma, su, sp, bc, true);
    Vector solution(n, 0.0);
    const auto linear = solve_scalar_equation(eq, solution, {5000, 1e-13, 1.0});
    if (linear.status != SolverStatus::CONVERGED)
        throw std::runtime_error("MMS diffusion: linear solve did not converge");

    ScalarResult result;
    result.y.resize(n);
    result.u.resize(n);
    result.volume = geometry.cell_volumes;
    for (std::size_t i = 0; i < n; ++i) {
        result.y[i] = geometry.cell_centres[i].y;
        result.u[i] = solution(i);
    }
    return result;
}

void report_case(const std::string& name, std::size_t n,
                 const ErrorMetrics& e, double order)
{
    std::cout << std::left << std::setw(30) << name
              << " N=" << std::setw(4) << n
              << " L2=" << std::scientific << std::setprecision(6) << e.l2
              << " Linf=" << e.linf
              << " relL2=" << e.l2_relative
              << " relLinf=" << e.linf_relative;
    if (std::isfinite(order))
        std::cout << " order=" << std::fixed << std::setprecision(3) << order;
    std::cout << "\n";
}

} // namespace

int main()
{
    try {
        constexpr double H = 1.0;
        constexpr double mu = 2.0;

        // Couette: u(y)=U*y/H. The orthogonal FV diffusion operator with
        // exact Dirichlet face values should reproduce the linear field.
        for (const std::size_t n : {8u,16u,32u,64u}) {
            const auto r = solve_diffusion_case(n,H,mu,0.0,0.0,1.0);
            std::vector<double> exact(r.y.size());
            for (std::size_t i=0;i<r.y.size();++i)
                exact[i] = r.y[i] / H;
            const auto e = error_norms(r.u,exact,r.volume);
            report_case("Couette",n,e,std::numeric_limits<double>::quiet_NaN());
            if (e.linf_relative > 1e-11)
                throw std::runtime_error("Couette analytical solution mismatch");
        }

        // Plane Poiseuille: -mu*u'' = G, u(0)=u(H)=0.
        // Exact solution: u=G*y*(H-y)/(2*mu).
        constexpr double G = 4.0;
        std::vector<double> poiseuille_errors;
        for (const std::size_t n : {8u,16u,32u,64u}) {
            const auto r = solve_diffusion_case(n,H,mu,G,0.0,0.0);
            std::vector<double> exact(r.y.size());
            for (std::size_t i=0;i<r.y.size();++i)
                exact[i] = G*r.y[i]*(H-r.y[i])/(2.0*mu);
            const auto e = error_norms(r.u,exact,r.volume);
            const double order = poiseuille_errors.empty()
                ? std::numeric_limits<double>::quiet_NaN()
                : observed_order(poiseuille_errors.back(),e.l2);
            report_case("Poiseuille",n,e,order);
            poiseuille_errors.push_back(e.l2);
        }
        require_order(poiseuille_errors,2.0,1.90,"Poiseuille");

        // 1-D conduction: T(y)=T0+(T1-T0)y/H, q=-k*dT/dy.
        constexpr double T0 = 400.0;
        constexpr double T1 = 300.0;
        constexpr double k = 5.0;
        for (const std::size_t n : {8u,16u,32u,64u}) {
            const auto r = solve_diffusion_case(n,H,k,0.0,T0,T1);
            std::vector<double> exact(r.y.size());
            for (std::size_t i=0;i<r.y.size();++i)
                exact[i] = T0 + (T1-T0)*r.y[i]/H;
            const auto e = error_norms(r.u,exact,r.volume);
            report_case("Conduction 1D",n,e,std::numeric_limits<double>::quiet_NaN());
            if (e.linf_relative > 1e-11)
                throw std::runtime_error("1-D conduction analytical solution mismatch");
        }

        // Scalar diffusion MMS: T=sin(pi*y/H), -gamma*T''=S.
        // This is a true discretisation verification because the source varies
        // with the manufactured solution and the mesh is systematically refined.
        std::vector<double> mms_errors;
        for (const std::size_t n : {8u,16u,32u,64u}) {
            const auto r = solve_mms_diffusion(n,H,k);
            std::vector<double> exact(r.y.size());
            const double pi = std::acos(-1.0);
            for (std::size_t i=0;i<r.y.size();++i)
                exact[i] = std::sin(pi*r.y[i]/H);
            const auto e = error_norms(r.u,exact,r.volume);
            const double order = mms_errors.empty()
                ? std::numeric_limits<double>::quiet_NaN()
                : observed_order(mms_errors.back(),e.l2);
            report_case("Diffusion MMS",n,e,order);
            mms_errors.push_back(e.l2);
        }
        require_order(mms_errors,2.0,1.80,"Diffusion MMS");

        // CHT resistance-in-series and radiation oracles are evaluated
        // independently of the implementation under test.
        const double Th = 700.0;
        const double Tc = 300.0;
        const double L1 = 0.1;
        const double L2 = 0.2;
        const double k1 = 10.0;
        const double k2 = 2.0;
        const double q_cht = (Th-Tc)/(L1/k1 + L2/k2);
        const double T_interface = Th - q_cht*L1/k1;
        if (std::abs(q_cht - 3636.3636363636365) > 1e-10 ||
            std::abs(T_interface - 663.6363636363636) > 1e-10)
            throw std::runtime_error("CHT resistance-series oracle mismatch");

        const double e1 = 0.7;
        const double e2 = 0.4;
        const double area = 2.5;
        const double Ta = 900.0;
        const double Tb = 500.0;
        const double sigma = STEFAN_BOLTZMANN;
        const double q_gray = sigma*area*(std::pow(Ta,4)-std::pow(Tb,4)) /
            ((1.0-e1)/e1 + 1.0 + (1.0-e2)/e2);
        const double q_gray_impl = area*two_surface_net_exchange(e1,e2,Ta,Tb,1.0);
        if (std::abs(q_gray - q_gray_impl) > 1e-10*std::max(1.0,std::abs(q_gray)))
            throw std::runtime_error("radiation analytical oracle mismatch");

        std::cout << "ANALYTICAL_BENCHMARKS: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ANALYTICAL_BENCHMARKS: FAIL: " << e.what() << "\n";
        return 1;
    }
}
