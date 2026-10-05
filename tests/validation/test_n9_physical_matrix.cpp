#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

Mesh make_channel_mesh(std::size_t nx, std::size_t ny, double skew)
{
    Mesh mesh;
    const std::size_t plane = (nx + 1) * (ny + 1);
    mesh.points().resize(2 * plane);
    const auto id = [nx](std::size_t i, std::size_t j, std::size_t k) {
        return (j * (nx + 1) + i) * 2 + k;
    };

    for (std::size_t j = 0; j <= ny; ++j) {
        const double y = static_cast<double>(j) / static_cast<double>(ny);
        const double sx = skew * y * (1.0 - y);
        for (std::size_t i = 0; i <= nx; ++i) {
            const double x = static_cast<double>(i) / static_cast<double>(nx) + sx;
            mesh.points().set(id(i,j,0), x, y, 0.0);
            mesh.points().set(id(i,j,1), x, y, 1.0);
        }
    }

    std::map<std::vector<std::size_t>, std::size_t> face_map;
    std::vector<std::vector<std::size_t>> cell_faces(nx * ny);

    auto add_face = [&](std::initializer_list<std::size_t> vertices,
                        std::size_t cell) {
        std::vector<std::size_t> key(vertices);
        std::sort(key.begin(), key.end());
        const auto it = face_map.find(key);
        if (it != face_map.end()) {
            mesh.ownership().set_neighbour(it->second, static_cast<int>(cell));
            return it->second;
        }
        const std::size_t f = mesh.faces().n_faces();
        mesh.faces().push_face(vertices);
        mesh.ownership().resize(mesh.faces().n_faces());
        mesh.ownership().set_owner(f, cell);
        mesh.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        face_map.emplace(std::move(key), f);
        return f;
    };

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const std::size_t c = j * nx + i;
            const auto a=id(i,j,0), b=id(i+1,j,0), c0=id(i+1,j+1,0), d=id(i,j+1,0);
            const auto e=id(i,j,1), f=id(i+1,j,1), g=id(i+1,j+1,1), h=id(i,j+1,1);
            cell_faces[c] = {
                add_face({a,d,c0,b},c), add_face({e,f,g,h},c),
                add_face({a,b,f,e},c), add_face({d,h,g,c0},c),
                add_face({a,e,h,d},c), add_face({b,c0,g,f},c)};
        }
    }
    for (const auto& faces : cell_faces) mesh.cells().push_cell(faces);

    Patch inlet{"inlet", PatchType::INLET, {}};
    Patch outlet{"outlet", PatchType::OUTLET, {}};
    Patch bottom{"bottom", PatchType::WALL, {}};
    Patch top{"top", PatchType::WALL, {}};
    Patch front{"front", PatchType::EMPTY, {}};
    Patch back{"back", PatchType::EMPTY, {}};

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f) >= 0) continue;
        const auto& v = mesh.faces().vertices();
        const auto begin = v.begin() + static_cast<std::ptrdiff_t>(mesh.faces().face_offset(f));
        const auto end = begin + static_cast<std::ptrdiff_t>(mesh.faces().face_size(f));
        double x=0.0,y=0.0,z=0.0;
        for (auto it=begin; it!=end; ++it) {
            x += mesh.points().x(*it); y += mesh.points().y(*it); z += mesh.points().z(*it);
        }
        const double n = static_cast<double>(mesh.faces().face_size(f));
        x/=n; y/=n; z/=n;
        constexpr double eps=1e-12;
        // The patch assignment is geometric in the unskewed logical coordinates.
        if (std::abs(y) < eps) bottom.face_ids.push_back(f);
        else if (std::abs(y-1.0) < eps) top.face_ids.push_back(f);
        else if (std::abs(z) < eps) front.face_ids.push_back(f);
        else if (std::abs(z-1.0) < eps) back.face_ids.push_back(f);
        else {
            // x is only used here to separate the two end patches. The
            // logical end faces are exactly the first/last i columns.
            bool at_left=true, at_right=true;
            for (auto it=begin; it!=end; ++it) {
                const double xx=mesh.points().x(*it);
                const double yy=mesh.points().y(*it);
                const double logical = xx - skew*yy*(1.0-yy);
                at_left = at_left && std::abs(logical) < eps;
                at_right = at_right && std::abs(logical-1.0) < eps;
            }
            if (at_left) inlet.face_ids.push_back(f);
            else if (at_right) outlet.face_ids.push_back(f);
            else throw std::runtime_error("unclassified channel boundary face");
        }
    }

    mesh.boundary().add_patch(inlet);
    mesh.boundary().add_patch(outlet);
    mesh.boundary().add_patch(bottom);
    mesh.boundary().add_patch(top);
    mesh.boundary().add_patch(front);
    mesh.boundary().add_patch(back);
    return mesh;
}

Mesh make_cavity_mesh(std::size_t nx, std::size_t ny)
{
    return make_channel_mesh(nx, ny, 0.0);
}

struct Run {
    Field<double, Location::CELL> U;
    Field<double, Location::CELL> p;
    IncompressibleSolveResult result;
};

IncompressibleSolverControls controls_for(
    PressureVelocityAlgorithm algorithm,
    double viscosity = 0.1,
    bool transient_projection = false)
{
    IncompressibleSolverControls c;
    c.algorithm = algorithm;
    c.density = 1.0;
    c.kinematic_viscosity = viscosity;
    c.coupling.alpha_u = 0.7;
    c.coupling.alpha_p = 0.3;
    c.coupling.n_pressure_correctors =
        (algorithm == PressureVelocityAlgorithm::PISO ||
         algorithm == PressureVelocityAlgorithm::PIMPLE) ? 2 : 1;
    c.coupling.n_outer_correctors =
        algorithm == PressureVelocityAlgorithm::PIMPLE ? 2 : 1;
    c.coupling.n_fractional_steps =
        algorithm == PressureVelocityAlgorithm::FRACTIONAL_STEP ? 2 : 1;
    c.coupling.coupled_max_iterations = 1000;
    // The block solve must be at least two orders tighter than the nonlinear
    // gate it feeds. solve_gmres stops on a 2-norm relative residual while the
    // acceptance gate is an infinity-norm momentum residual scaled by the RHS
    // magnitude, so a 1e-9 block tolerance leaves an absolute momentum residual
    // of order 1e-8 and COUPLED cannot satisfy a 1e-8 relative gate on
    // Poiseuille no matter how many outer iterations it is given: the state is
    // frozen and the residual floor is set by the Krylov solve, not the
    // nonlinearity. The block solve is a direct solve of the coupled system, so
    // the tolerance is tightened rather than the gate relaxed.
    c.coupling.coupled_linear_tolerance = 1e-12;
    c.convergence.max_iterations = 1500;
    c.convergence.relative_tolerance = 1e-8;
    c.convergence.continuity_tolerance = 1e-8;
    c.linear_max_iterations = 1000;
    c.linear_tolerance = 1e-9;
    c.pressure_reference_cell = 0;
    c.pressure_reference_value = 0.0;
    c.use_bounded_convection = true;
    c.convection_scheme = ConvectionScheme::UPWIND;
    (void)transient_projection;
    return c;
}

VelocityBoundaryConditions channel_velocity_bc(double top_u)
{
    VelocityBoundaryConditions bc;
    bc["inlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    bc["outlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    bc["bottom"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    bc["top"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{top_u,0,0}};
    bc["front"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    bc["back"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    return bc;
}

ScalarBoundaryConditions channel_pressure_bc()
{
    ScalarBoundaryConditions bc;
    for (const char* n : {"inlet","outlet","bottom","top","front","back"})
        bc[n] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
    return bc;
}

Run solve_case(
    Mesh mesh,
    PressureVelocityAlgorithm algorithm,
    const VelocityBoundaryConditions& ubc,
    const ScalarBoundaryConditions& pbc,
    double body_force_x,
    double viscosity = 0.1)
{
    Field<double, Location::CELL> U(mesh.n_cells(),"U","m/s",3);
    Field<double, Location::CELL> p(mesh.n_cells(),"p","Pa",1);
    U.fill(0.0);
    p.fill(17.0);

    auto c = controls_for(algorithm, viscosity);
    c.body_force = {body_force_x,0.0,0.0};
    const auto result = solve_steady_incompressible(mesh,U,p,ubc,pbc,c);
    return {std::move(U),std::move(p),std::move(result)};
}

void require_physical_convergence(const Run& r, const char* name)
{
    if (!r.result.converged || r.result.history.empty())
        throw std::runtime_error(std::string(name)+" did not converge");
    const auto& h=r.result.history.back();
    if (!std::isfinite(h.continuity_linf) ||
        !std::isfinite(h.momentum_residual) ||
        h.continuity_linf > 1e-7 ||
        h.momentum_residual > 1e-7)
        throw std::runtime_error(std::string(name)+" physical residual gate failed");
}

double couette_l2(const Run& r, std::size_t nx, std::size_t ny)
{
    double e2=0.0;
    for (std::size_t c=0;c<r.U.size();++c) {
        const std::size_t i=c%nx, j=c/nx;
        (void)i;
        const double y=(static_cast<double>(j)+0.5)/static_cast<double>(ny);
        const double e=r.U.component_data(0)[c]-y;
        e2 += e*e;
    }
    return std::sqrt(e2/static_cast<double>(r.U.size()));
}

double poiseuille_l2(const Run& r, std::size_t nx, std::size_t ny, double G, double nu)
{
    double e2=0.0;
    for (std::size_t c=0;c<r.U.size();++c) {
        const std::size_t j=c/nx;
        const double y=(static_cast<double>(j)+0.5)/static_cast<double>(ny);
        const double exact=G*y*(1.0-y)/(2.0*nu);
        const double e=r.U.component_data(0)[c]-exact;
        e2 += e*e;
    }
    return std::sqrt(e2/static_cast<double>(r.U.size()));
}

double cavity_profile_linf(const Run& r, std::size_t nx, std::size_t ny)
{
    // Ghia et al. Re=100 centreline reference values. These points are
    // deliberately sparse here because the dedicated Ghia campaign retains
    // the complete 17-point profiles; N9 uses them as an independent physical
    // cross-algorithm gate.
    constexpr double y[] = {0.0625, 0.1719, 0.5000, 0.8516, 0.9609};
    constexpr double u_ref[] = {-0.04192, -0.10150, -0.20581, 0.23151, 0.73722};
    constexpr double x[] = {0.0625, 0.2266, 0.5000, 0.9063, 0.9688};
    constexpr double v_ref[] = {0.09233, 0.17507, 0.05454, -0.16914, -0.05906};
    double max_error = 0.0;
    const auto sample = [&](double qx, double qy, std::size_t component) {
        const double fx = qx*static_cast<double>(nx)-0.5;
        const double fy = qy*static_cast<double>(ny)-0.5;
        const auto clamp_index = [](double q, std::size_t n) {
            return std::clamp(static_cast<long>(std::floor(q)), 0L,
                              static_cast<long>(n)-2L);
        };
        const long ix0=clamp_index(fx,nx), iy0=clamp_index(fy,ny);
        const double tx=std::clamp(fx-static_cast<double>(ix0),0.0,1.0);
        const double ty=std::clamp(fy-static_cast<double>(iy0),0.0,1.0);
        const std::size_t i=static_cast<std::size_t>(ix0), j=static_cast<std::size_t>(iy0);
        const auto at=[&](std::size_t ii,std::size_t jj) { return r.U.component_data(component)[jj*nx+ii]; };
        return (1.0-ty)*((1.0-tx)*at(i,j)+tx*at(i+1,j))
             + ty*((1.0-tx)*at(i,j+1)+tx*at(i+1,j+1));
    };
    for (std::size_t i=0;i<5;++i) {
        max_error=std::max(max_error,std::abs(sample(0.5,y[i],0)-u_ref[i]));
        max_error=std::max(max_error,std::abs(sample(x[i],0.5,1)-v_ref[i]));
    }
    return max_error;
}

double max_difference(const Run& a, const Run& b)
{
    if (a.U.size()!=b.U.size()) throw std::invalid_argument("field size mismatch");
    double e=0.0;
    for (std::size_t c=0;c<a.U.size();++c)
        for (std::size_t d=0;d<3;++d)
            e=std::max(e,std::abs(a.U.component_data(d)[c]-b.U.component_data(d)[c]));
    return e;
}

std::array<PressureVelocityAlgorithm,6> algorithms()
{
    return {PressureVelocityAlgorithm::SIMPLE,PressureVelocityAlgorithm::SIMPLEC,
            PressureVelocityAlgorithm::PISO,PressureVelocityAlgorithm::PIMPLE,
            PressureVelocityAlgorithm::FRACTIONAL_STEP,PressureVelocityAlgorithm::COUPLED};
}

} // namespace

int main()
{
    const auto algs=algorithms();

    // Common physical matrix: Couette, pressure-driven Poiseuille and lid-driven
    // cavity. The same meshes, BCs, viscosity and convergence contract are used
    // for all six pressure-velocity paths.
    std::vector<Run> couette;
    std::vector<Run> poiseuille;
    std::vector<Run> cavity;
    couette.reserve(algs.size());
    poiseuille.reserve(algs.size());
    cavity.reserve(algs.size());

    const auto u_channel=channel_velocity_bc(1.0);
    const auto p_channel=channel_pressure_bc();

    for (std::size_t k=0;k<algs.size();++k) {
        couette.push_back(solve_case(make_channel_mesh(12,16,0.0),algs[k],u_channel,p_channel,0.0));
        require_physical_convergence(couette[k],"Couette");
        if (couette_l2(couette[k],12,16)>2e-3)
            throw std::runtime_error("Couette analytic L2 gate failed");

        poiseuille.push_back(solve_case(make_channel_mesh(12,16,0.0),algs[k],channel_velocity_bc(0.0),p_channel,1.0));
        require_physical_convergence(poiseuille[k],"Poiseuille");
        if (poiseuille_l2(poiseuille[k],12,16,1.0,0.1)>5e-3)
            throw std::runtime_error("Poiseuille analytic L2 gate failed");

        auto cavity_bc=channel_velocity_bc(0.0);
        cavity_bc["inlet"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
        cavity_bc["outlet"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
        cavity_bc["bottom"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
        cavity_bc["top"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{1,0,0}};
        cavity.push_back(solve_case(make_cavity_mesh(32,32),algs[k],cavity_bc,p_channel,0.0,0.01));
        require_physical_convergence(cavity[k],"Ghia Re=100 cavity");
        if (cavity_profile_linf(cavity[k],32,32)>0.20)
            throw std::runtime_error("Ghia Re=100 centreline oracle gate failed");
    }

    // Cross-algorithm physical-equivalence gate on each benchmark.
    for (std::size_t k=1;k<algs.size();++k) {
        if (max_difference(couette[0],couette[k])>2e-4)
            throw std::runtime_error("Couette algorithm equivalence gate failed");
        if (max_difference(poiseuille[0],poiseuille[k])>2e-4)
            throw std::runtime_error("Poiseuille algorithm equivalence gate failed");
        if (max_difference(cavity[0],cavity[k])>2e-2)
            throw std::runtime_error("cavity algorithm equivalence gate failed");
    }

    // Controlled skew/non-orthogonal campaign. The physical solution remains
    // affine Couette, so the analytic gate is independent of algorithm while
    // the mesh exercises non-orthogonal/skew face geometry.
    for (std::size_t k=0;k<algs.size();++k) {
        const auto skew=solve_case(make_channel_mesh(12,16,0.25),algs[k],
                                   u_channel,p_channel,0.0);
        require_physical_convergence(skew,"skew Couette");
        if (couette_l2(skew,12,16)>2e-2)
            throw std::runtime_error("skew Couette analytic L2 gate failed");
        if (max_difference(couette[0],skew)>3e-2)
            throw std::runtime_error("skew Couette physical-equivalence gate failed");
    }

    return 0;
}
