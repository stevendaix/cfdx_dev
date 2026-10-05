#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {

Mesh make_cavity_mesh(std::size_t nx, std::size_t ny)
{
    Mesh mesh;
    mesh.points().resize((nx + 1) * (ny + 1) * 2);
    const auto id = [nx](std::size_t i, std::size_t j, std::size_t k) {
        return (j * (nx + 1) + i) * 2 + k;
    };
    for (std::size_t j = 0; j <= ny; ++j)
        for (std::size_t i = 0; i <= nx; ++i) {
            const double x = static_cast<double>(i) / static_cast<double>(nx);
            const double y = static_cast<double>(j) / static_cast<double>(ny);
            mesh.points().set(id(i,j,0), x, y, 0.0);
            mesh.points().set(id(i,j,1), x, y, 1.0);
        }

    std::map<std::vector<std::size_t>, std::size_t> faces;
    std::vector<std::vector<std::size_t>> cell_faces(nx * ny);
    auto add_face = [&](std::initializer_list<std::size_t> vertices, std::size_t cell) {
        std::vector<std::size_t> key(vertices);
        std::sort(key.begin(), key.end());
        const auto found = faces.find(key);
        if (found != faces.end()) {
            mesh.ownership().set_neighbour(found->second, static_cast<int>(cell));
            return found->second;
        }
        const std::size_t face = mesh.faces().n_faces();
        mesh.faces().push_face(vertices);
        mesh.ownership().resize(mesh.faces().n_faces());
        mesh.ownership().set_owner(face, cell);
        mesh.ownership().set_neighbour(face, FaceOwnership::BOUNDARY);
        faces.emplace(std::move(key), face);
        return face;
    };

    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i) {
            const std::size_t c = j * nx + i;
            const auto a=id(i,j,0), b=id(i+1,j,0), cc=id(i+1,j+1,0), d=id(i,j+1,0);
            const auto e=id(i,j,1), f=id(i+1,j,1), g=id(i+1,j+1,1), h=id(i,j+1,1);
            cell_faces[c] = {
                add_face({a,d,cc,b},c), add_face({e,f,g,h},c),
                add_face({a,b,f,e},c), add_face({d,h,g,cc},c),
                add_face({a,e,h,d},c), add_face({b,cc,g,f},c)};
        }
    for (const auto& cf : cell_faces)
        mesh.cells().push_cell(cf);

    Patch bottom{"bottom", PatchType::WALL, {}};
    Patch top{"top", PatchType::WALL, {}};
    Patch left{"left", PatchType::WALL, {}};
    Patch right{"right", PatchType::WALL, {}};
    Patch front{"front", PatchType::WALL, {}};
    Patch back{"back", PatchType::WALL, {}};

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f) >= 0)
            continue;
        const auto& vertices = mesh.faces().vertices();
        const auto begin = vertices.begin() + static_cast<std::ptrdiff_t>(mesh.faces().face_offset(f));
        const auto end = begin + static_cast<std::ptrdiff_t>(mesh.faces().face_size(f));
        double x=0.0, y=0.0, z=0.0;
        for (auto p=begin; p!=end; ++p) {
            x += mesh.points().x(*p); y += mesh.points().y(*p); z += mesh.points().z(*p);
        }
        const double n = static_cast<double>(mesh.faces().face_size(f));
        x /= n; y /= n; z /= n;
        constexpr double tol = 1e-12;
        if (std::abs(y) < tol) bottom.face_ids.push_back(f);
        else if (std::abs(y-1.0) < tol) top.face_ids.push_back(f);
        else if (std::abs(x) < tol) left.face_ids.push_back(f);
        else if (std::abs(x-1.0) < tol) right.face_ids.push_back(f);
        else if (std::abs(z) < tol) front.face_ids.push_back(f);
        else if (std::abs(z-1.0) < tol) back.face_ids.push_back(f);
        else throw std::runtime_error("N7 Re=400 cavity: unclassified boundary face");
    }
    mesh.boundary().add_patch(bottom); mesh.boundary().add_patch(top);
    mesh.boundary().add_patch(left); mesh.boundary().add_patch(right);
    mesh.boundary().add_patch(front); mesh.boundary().add_patch(back);
    return mesh;
}

double centre_u(const Field<double,Location::CELL>& U)
{
    constexpr std::size_t n = 16;
    const std::size_t i0=n/2-1, i1=n/2, j0=n/2-1, j1=n/2;
    const double* v=U.component_data(0);
    return 0.25*(v[j0*n+i0]+v[j0*n+i1]+v[j1*n+i0]+v[j1*n+i1]);
}

IncompressibleSolverControls controls()
{
    IncompressibleSolverControls c;
    c.algorithm=PressureVelocityAlgorithm::SIMPLE;
    c.density=1.0;
    c.kinematic_viscosity=0.0025; // Re=400, U=1, L=1.
    c.linear_max_iterations=10000;
    c.linear_tolerance=1e-10;
    c.pressure_reference_cell=128;
    c.pressure_reference_value=0.0;
    c.use_bounded_convection=true;
    c.convection_scheme=ConvectionScheme::UPWIND;
    c.coupling.alpha_u=0.7;
    c.coupling.alpha_p=0.3;
    c.convergence.max_iterations=2500;
    c.convergence.relative_tolerance=1e-8;
    c.convergence.continuity_tolerance=1e-8;
    return c;
}

void set_bcs(VelocityBoundaryConditions& u, ScalarBoundaryConditions& p)
{
    for (const char* name : {"bottom","left","right","front","back"})
        u[name]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    u["top"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{1,0,0}};
    for (const char* name : {"bottom","top","left","right","front","back"})
        p[name]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
}

} // namespace

int main()
{
    return run_case("n7_adaptive_relaxation_re400_16x16", [] {
        const auto mesh=make_cavity_mesh(16,16);
        Field<double,Location::CELL> Ud(mesh.n_cells(),"U","m/s",3);
        Field<double,Location::CELL> pd(mesh.n_cells(),"p","Pa",1);
        Field<double,Location::CELL> Ua(mesh.n_cells(),"U","m/s",3);
        Field<double,Location::CELL> pa(mesh.n_cells(),"p","Pa",1);
        Ud.fill(0.0); pd.fill(0.0); Ua.fill(0.0); pa.fill(0.0);

        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        set_bcs(ubc,pbc);

        auto direct_controls=controls();
        const auto direct=solve_steady_incompressible(mesh,Ud,pd,ubc,pbc,direct_controls);
        EXPECT_TRUE(direct.converged);
        EXPECT_TRUE(direct.convergence_status==ConvergenceStatus::CONVERGED);

        auto adaptive_controls=direct_controls;
        adaptive_controls.adaptive_relaxation.enabled=true;
        adaptive_controls.adaptive_relaxation.min_alpha_u=0.6;
        adaptive_controls.adaptive_relaxation.max_alpha_u=0.8;
        adaptive_controls.adaptive_relaxation.min_alpha_p=0.25;
        adaptive_controls.adaptive_relaxation.max_alpha_p=0.35;

        const auto adaptive=solve_steady_incompressible(
            mesh,Ua,pa,ubc,pbc,adaptive_controls);
        EXPECT_TRUE(adaptive.converged);
        EXPECT_TRUE(adaptive.convergence_status==ConvergenceStatus::CONVERGED);
        EXPECT_TRUE(adaptive.history.size()>=3);

        std::size_t changes=0;
        for (const auto& h : adaptive.history) {
            EXPECT_TRUE(std::isfinite(h.nonlinear_convergence_metric));
            EXPECT_TRUE(h.effective_alpha_u>=0.6 && h.effective_alpha_u<=0.8);
            EXPECT_TRUE(h.effective_alpha_p>=0.25 && h.effective_alpha_p<=0.35);
            if (std::abs(h.effective_alpha_u-0.7)>1e-12 ||
                std::abs(h.effective_alpha_p-0.3)>1e-12)
                ++changes;
        }
        EXPECT_TRUE(changes>0);

        const double scale=std::max(std::abs(centre_u(Ud)),1e-3);
        EXPECT_TRUE(std::abs(centre_u(Ua)-centre_u(Ud)) <= 1e-2*scale);
    });
}
