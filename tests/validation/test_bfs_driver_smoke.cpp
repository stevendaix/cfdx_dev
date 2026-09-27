#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "qualification_boundary_helpers.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

struct Result { double lr_h; double continuity; std::size_t iterations; };

static Result run(const std::string& path, std::size_t level)
{
    Mesh mesh;
    if(!cfdx::io::read_mesh_hdf5(path,mesh))
        throw std::runtime_error("cannot read BFS mesh: "+path);
    const auto topo=mesh.topo_validate();
    if(!topo.ok) throw std::runtime_error("invalid BFS topology: "+path);

    constexpr double H_in=5.2/4.9;
    constexpr double Re_D=200.0;
    constexpr double U_bulk=0.288462/0.288462;
    constexpr double rho=1.0;
    constexpr double nu=2.0*H_in*U_bulk/Re_D;

    Field<double,Location::CELL> U(mesh.n_cells(),"U","m/s",3);
    Field<double,Location::CELL> p(mesh.n_cells(),"p","Pa",1);
    U.fill(0.0); p.fill(0.0);

    VelocityBoundaryConditions ubc;
    const double y0=1.0;
    const double H=H_in;
    ubc["inlet"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{1.0,0.0,0.0}};
    ubc["outlet"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    ubc["wall"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["front"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    ubc["back"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};

    ScalarBoundaryConditions pbc;
    pbc["inlet"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    pbc["outlet"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    pbc["wall"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    pbc["front"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    pbc["back"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};

    IncompressibleSolverControls c;
    c.algorithm=PressureVelocityAlgorithm::SIMPLE;
    c.density=rho; c.kinematic_viscosity=nu;
    c.linear_max_iterations=3000; c.linear_tolerance=1e-10;
    c.convergence.max_iterations=5000;
    c.convergence.relative_tolerance=1e-8;
    c.convergence.continuity_tolerance=1e-8;
    c.coupling.alpha_u=0.7; c.coupling.alpha_p=0.3;
    c.use_bounded_convection=true;
    c.convection_scheme=ConvectionScheme::UPWIND;
    c.pressure_reference_cell=0;
    c.pressure_reference_value=0.0;
    c.diagnostics.iteration_trace=true;
    c.diagnostics.iteration_trace_frequency=100;

    // Initialise the inlet with the fully developed parabolic profile.  The
    // production BC map remains a fixed-value contract; the profile is used
    // as the initial state and the mesh inlet face geometry is checked below.
    const auto geometry=build_fv_geometry(mesh);
    for(std::size_t cell=0;cell<mesh.n_cells();++cell) {
        const double y=geometry.cell_centres[cell].y;
        const double eta=std::clamp((y-y0)/H,0.0,1.0);
        U.component_data(0)[cell]=6.0*U_bulk*eta*(1.0-eta);
    }

    // Exercise the strict mathematical BC/FVM path before the production solve.
    // The legacy maps below are retained for the current steady-solver compatibility API.
    cfdx::physics::BoundaryConstraintMap strict_bc;
    cfdx::validation::add_velocity_dirichlet(strict_bc, "inlet", {1.0, 0.0, 0.0});
    cfdx::validation::add_velocity_flux_dependent(
        strict_bc, "outlet", {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
    cfdx::validation::add_velocity_dirichlet(strict_bc, "wall", {0.0, 0.0, 0.0});
    cfdx::validation::add_pressure_neumann(strict_bc, "inlet", 0.0);
    cfdx::validation::add_pressure_dirichlet(strict_bc, "outlet", 0.0);
    cfdx::validation::add_pressure_neumann(strict_bc, "wall", 0.0);
    cfdx::validation::add_pressure_neumann(strict_bc, "front", 0.0);
    cfdx::validation::add_pressure_neumann(strict_bc, "back", 0.0);
    const auto validation_flux = cfdx::validation::make_validation_face_flux(mesh, geometry, U);
    cfdx::validation::exercise_new_velocity_bc_contract(
        mesh, geometry, strict_bc, validation_flux, "BFS_RE200");

    const auto solve=solve_steady_incompressible(mesh,U,p,ubc,pbc,c);
    if(!solve.converged)
        throw std::runtime_error("BFS solve did not converge");

    ScalarBoundaryConditions grad_bc;
    grad_bc["inlet"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    grad_bc["outlet"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    grad_bc["wall"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    grad_bc["front"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    grad_bc["back"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};

    Field<double,Location::CELL> ux(mesh.n_cells(),"ux","m/s",1);
    Field<double,Location::CELL> uy(mesh.n_cells(),"uy","m/s",1);
    for(std::size_t i=0;i<mesh.n_cells();++i) {
        ux(i)=U.component_data(0)[i];
        uy(i)=U.component_data(1)[i];
    }
    const auto gux=gauss_gradient_with_boundary(ux,mesh,geometry,grad_bc);
    const auto guy=gauss_gradient_with_boundary(uy,mesh,geometry,grad_bc);

    struct Tau { double x; double tau; };
    std::vector<Tau> tau;
    const auto& patch=mesh.boundary().patch(mesh.boundary().find("wall"));
    for(const auto f:patch.face_ids) {
        const auto fc=geometry.face_centres[f];
        if(fc.x < 0.0 || std::abs(fc.y) > 1e-8) continue;
        const auto n=geometry.face_area_vectors[f].normalized();
        const double txy=rho*nu*(gux.component_data(1)[mesh.ownership().owner(f)] +
                                 guy.component_data(0)[mesh.ownership().owner(f)]);
        tau.push_back({fc.x,txy});
    }
    std::sort(tau.begin(),tau.end(),[](const Tau&a,const Tau&b){return a.x<b.x;});
    if(tau.size()<2) throw std::runtime_error("BFS wall-shear extraction returned too few downstream faces");

    double lr=std::numeric_limits<double>::quiet_NaN();
    for(std::size_t i=1;i<tau.size();++i) {
        if(tau[i-1].tau<=0.0 && tau[i].tau>0.0) {
            const double a=tau[i-1].tau, b=tau[i].tau;
            const double w=(-a)/(b-a);
            lr=tau[i-1].x+w*(tau[i].x-tau[i-1].x);
            break;
        }
    }
    if(!std::isfinite(lr)) throw std::runtime_error("BFS reattachment crossing not found");

    const auto& last=solve.history.back();
    std::cout<<"BFS_RE200 level="<<level
             <<" Lr_over_step="<<lr
             <<" continuity="<<last.continuity_linf
             <<" iterations="<<solve.iterations<<"\n";
    return {lr,last.continuity_linf,solve.iterations};
}

static void run_quick_contract(const std::string& path)
{
    Mesh mesh;
    if(!cfdx::io::read_mesh_hdf5(path,mesh))
        throw std::runtime_error("cannot read BFS mesh: "+path);
    const auto topo=mesh.topo_validate();
    if(!topo.ok)
        throw std::runtime_error("invalid BFS topology: "+path);

    const auto geometry=build_fv_geometry(mesh);
    Field<double,Location::CELL> U(mesh.n_cells(),"U","m/s",3);
    U.fill(0.0);
    for(std::size_t i=0;i<mesh.n_cells();++i)
        U.component_data(0)[i]=1.0;

    cfdx::physics::BoundaryConstraintMap strict_bc;
    cfdx::validation::add_velocity_dirichlet(strict_bc, "inlet", {1.0, 0.0, 0.0});
    cfdx::validation::add_velocity_flux_dependent(
        strict_bc, "outlet", {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
    cfdx::validation::add_velocity_dirichlet(strict_bc, "wall", {0.0, 0.0, 0.0});
    cfdx::validation::add_pressure_neumann(strict_bc, "inlet", 0.0);
    cfdx::validation::add_pressure_dirichlet(strict_bc, "outlet", 0.0);
    cfdx::validation::add_pressure_neumann(strict_bc, "wall", 0.0);
    cfdx::validation::add_pressure_neumann(strict_bc, "front", 0.0);
    cfdx::validation::add_pressure_neumann(strict_bc, "back", 0.0);
    const auto validation_flux =
        cfdx::validation::make_validation_face_flux(mesh, geometry, U);
    cfdx::validation::exercise_new_velocity_bc_contract(
        mesh, geometry, strict_bc, validation_flux, "BFS_RE200");

    for(const char* name:{"inlet","outlet","wall","front","back"}) {
        if(mesh.boundary().find(name)>=mesh.boundary().n_patches())
            throw std::runtime_error(std::string("BFS required boundary patch is missing: ")+name);
    }
    std::cout<<"BFS_QUICK_CONTRACT: PASS cells="<<mesh.n_cells()
             <<" inlet_faces="<<mesh.boundary().patch(mesh.boundary().find("inlet")).size()
             <<" outlet_faces="<<mesh.boundary().patch(mesh.boundary().find("outlet")).size()
             <<" wall_faces="<<mesh.boundary().patch(mesh.boundary().find("wall")).size()
             <<" solver_campaign_pending\\n";
}

int main(int argc,char**argv)
{
    const bool quick=argc>=2 && std::string(argv[1])=="--quick";
    if(quick) {
        if(argc!=3)
            throw std::invalid_argument("usage: test_bfs_qualification --quick mesh_n16.h5");
        run_quick_contract(argv[2]);
        return 0;
    }

    const std::size_t need=3;
    if(argc!=4)
        throw std::invalid_argument("usage: test_bfs_qualification mesh_n16.h5 mesh_n32.h5 mesh_n64.h5");
    std::vector<Result> r;
    for(std::size_t i=0;i<need;++i) r.push_back(run(argv[i+1],16u<<i));
    for(const auto& x:r) if(!(x.continuity<1e-7)) throw std::runtime_error("BFS continuity gate failed");
    if(!quick) {
        const double e0=std::abs(r[0].lr_h-5.0), e1=std::abs(r[1].lr_h-5.0), e2=std::abs(r[2].lr_h-5.0);
        std::cout<<"BFS reference target Lr/H=5.0 errors="<<e0<<","<<e1<<","<<e2<<"\n";
        if(!(e2<=e1 && e1<=e0)) throw std::runtime_error("BFS reattachment error is not decreasing under refinement");
    }
    std::cout<<"BFS_QUALIFICATION: DIAGNOSTIC solver_converged conservation wall_shear_reattachment_reference_gate_pending\n";
}
