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
    bool transient_projection = false)
{
    IncompressibleSolverControls c;
    c.algorithm = algorithm;
    c.density = 1.0;
    c.kinematic_viscosity = 0.1;
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
    double body_force_x)
{
    Field<double, Location::CELL> U(mesh.n_cells(),"U","m/s",3);
    Field<double, Location::CELL> p(mesh.n_cells(),"p","Pa",1);
    U.fill(0.0);
    p.fill(17.0);

    auto c = controls_for(algorithm);
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

double max_difference(const Run& a, const Run& b)
{
    if (a.U.size()!=b.U.size()) throw std::invalid_argument("field size mismatch");
    double e=0.0;
    for (std::size_t c=0;c<a.U.size();++c)
        for (std::size_t d=0;d<3;++d)
            e=std::max(e,std::abs(a.U.component_data(d)[c]-b.U.component_data(d)[c]));
    return e;
}

// External flow. A uniform stream is admitted at the inlet and leaves through
// an outlet whose static pressure is prescribed. The plate and the far field
// carry no prescribed velocity, so their face flux interpolates the owner value
// and the wall applies no shear; the uniform stream is then an exact discrete
// solution of the momentum and continuity equations with a constant pressure
// that satisfies the prescribed outlet value. The configuration therefore has a
// production-quality oracle, and it is the member of the matrix that pins the
// pressure level instead of leaving it to the gauge, which the closed-channel
// benchmarks never exercise.
VelocityBoundaryConditions external_flow_velocity_bc()
{
    VelocityBoundaryConditions bc;
    bc["inlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{1,0,0}};
    for (const char* n : {"outlet","bottom","top","front","back"})
        bc[n] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    return bc;
}

ScalarBoundaryConditions external_flow_pressure_bc()
{
    ScalarBoundaryConditions bc = channel_pressure_bc();
    bc["outlet"] = {ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
    return bc;
}

// Pressure-driven external flow. The same duct is driven by a prescribed static
// pressure difference across the inlet and the outlet instead of by a body
// force. The exact solution is Poiseuille, so the analytic oracle is unchanged,
// and the pressure level is now determined by two fixed-pressure boundaries
// rather than by the pure-Neumann gauge. This is the branch that separates a
// fixed-pressure boundary from a gauge constraint.
ScalarBoundaryConditions pressure_driven_pressure_bc(double p_in, double p_out)
{
    ScalarBoundaryConditions bc = channel_pressure_bc();
    bc["inlet"] = {ScalarBoundaryType::FIXED_VALUE,p_in,p_in};
    bc["outlet"] = {ScalarBoundaryType::FIXED_VALUE,p_out,p_out};
    return bc;
}

double shear_stream_l2(const Run& r)
{
    double e2=0.0;
    for (std::size_t c=0;c<r.U.size();++c) {
        const double e=r.U.component_data(0)[c]-1.0;
        e2 += e*e;
    }
    return std::sqrt(e2/static_cast<double>(r.U.size()));
}

double pressure_linf(const Run& r, double exact)
{
    double e=0.0;
    for (std::size_t c=0;c<r.p.size();++c)
        e=std::max(e,std::abs(r.p(c)-exact));
    return e;
}

double pressure_driven_l2(
    const Run& r, std::size_t nx, std::size_t ny, double dp, double nu)
{
    double e2=0.0;
    for (std::size_t c=0;c<r.U.size();++c) {
        const std::size_t j=c/nx;
        const double y=(static_cast<double>(j)+0.5)/static_cast<double>(ny);
        const double exact=dp*y*(1.0-y)/(2.0*nu);
        const double e=r.U.component_data(0)[c]-exact;
        e2 += e*e;
    }
    return std::sqrt(e2/static_cast<double>(r.U.size()));
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
        cavity.push_back(solve_case(make_cavity_mesh(16,16),algs[k],cavity_bc,p_channel,0.0));
        require_physical_convergence(cavity[k],"cavity");
    }

    // Cross-algorithm physical-equivalence gate on each benchmark.
    for (std::size_t k=1;k<algs.size();++k) {
        if (max_difference(couette[0],couette[k])>2e-4)
            throw std::runtime_error("Couette algorithm equivalence gate failed");
        if (max_difference(poiseuille[0],poiseuille[k])>2e-4)
            throw std::runtime_error("Poiseuille algorithm equivalence gate failed");
        if (max_difference(cavity[0],cavity[k])>5e-3)
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

    // External flow: uniform stream in, prescribed static pressure out. The
    // oracle is the exact discrete solution, so the gate is set by round-off and
    // not by mesh resolution, and every one of the six algorithms has to land on
    // it independently.
    std::vector<Run> external;
    external.reserve(algs.size());
    const auto u_external=external_flow_velocity_bc();
    const auto p_external=external_flow_pressure_bc();
    for (std::size_t k=0;k<algs.size();++k) {
        external.push_back(solve_case(make_channel_mesh(12,16,0.0),algs[k],
                                      u_external,p_external,0.0));
        require_physical_convergence(external[k],"external flow");
        if (shear_stream_l2(external[k])>1e-6)
            throw std::runtime_error("external flow analytic L2 gate failed");
        if (pressure_linf(external[k],0.0)>1e-6)
            throw std::runtime_error("external flow pressure-level gate failed");
    }
    for (std::size_t k=1;k<algs.size();++k)
        if (max_difference(external[0],external[k])>1e-6)
            throw std::runtime_error("external flow algorithm equivalence gate failed");

    // Pressure-driven external flow: the prescribed pressure difference replaces
    // the body force, so the same analytic profile is reached with the pressure
    // level fixed at both ends and no gauge freedom at all.
    constexpr double dp=0.1;
    std::vector<Run> driven;
    driven.reserve(algs.size());
    const auto u_driven=channel_velocity_bc(0.0);
    const auto p_driven=pressure_driven_pressure_bc(dp,0.0);
    for (std::size_t k=0;k<algs.size();++k) {
        driven.push_back(solve_case(make_channel_mesh(12,16,0.0),algs[k],
                                    u_driven,p_driven,0.0));
        require_physical_convergence(driven[k],"pressure-driven external flow");
        if (pressure_driven_l2(driven[k],12,16,dp,0.1)>5e-3)
            throw std::runtime_error("pressure-driven external flow analytic L2 gate failed");
    }
    for (std::size_t k=1;k<algs.size();++k)
        if (max_difference(driven[0],driven[k])>2e-4)
            throw std::runtime_error("pressure-driven external flow equivalence gate failed");

    return 0;
}
