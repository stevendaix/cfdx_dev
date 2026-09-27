#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

struct Forces { double cd,cl,cm,cdp,cdv,continuity; };

static Forces run(const std::string& path,std::size_t level)
{
    Mesh mesh;
    if(!cfdx::io::read_mesh_hdf5(path,mesh))
        throw std::runtime_error("cannot read NACA0012 mesh: "+path);
    if(!mesh.topo_validate().ok)
        throw std::runtime_error("invalid NACA0012 topology: "+path);

    constexpr double rho=1.0, Uinf=1.0, Re=1000.0, chord=1.0;
    constexpr double nu=Uinf*chord/Re;
    Field<double,Location::CELL> U(mesh.n_cells(),"U","m/s",3);
    Field<double,Location::CELL> p(mesh.n_cells(),"p","Pa",1);
    U.fill(0.0); p.fill(0.0);
    VelocityBoundaryConditions ubc;
    ubc["farfield"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{Uinf,0,0}};
    ubc["airfoil"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["front"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    ubc["back"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    ScalarBoundaryConditions pbc;
    pbc["farfield"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    pbc["airfoil"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    pbc["front"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    pbc["back"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};

    IncompressibleSolverControls c;
    c.algorithm=PressureVelocityAlgorithm::SIMPLE;
    c.density=rho; c.kinematic_viscosity=nu;
    c.linear_max_iterations=4000; c.linear_tolerance=1e-10;
    c.convergence.max_iterations=8000;
    c.convergence.relative_tolerance=1e-8;
    c.convergence.continuity_tolerance=1e-8;
    c.coupling.alpha_u=0.6; c.coupling.alpha_p=0.25;
    c.use_bounded_convection=true;
    c.convection_scheme=ConvectionScheme::UPWIND;
    c.pressure_reference_cell=0;
    c.pressure_reference_value=0.0;
    c.diagnostics.iteration_trace=true;
    c.diagnostics.iteration_trace_frequency=200;

    const auto geometry=build_fv_geometry(mesh);
    for(std::size_t i=0;i<mesh.n_cells();++i)
        U.component_data(0)[i]=Uinf;

    const auto solve=solve_steady_incompressible(mesh,U,p,ubc,pbc,c);
    if(!solve.converged)
        throw std::runtime_error("NACA0012 solver did not converge");

    ScalarBoundaryConditions gx,gy;
    gx["farfield"]={ScalarBoundaryType::FIXED_VALUE,Uinf,0};
    gx["airfoil"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    gx["front"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    gx["back"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    gy["farfield"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    gy["airfoil"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    gy["front"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    gy["back"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    Field<double,Location::CELL> ux(mesh.n_cells(),"ux","m/s",1), uy(mesh.n_cells(),"uy","m/s",1);
    for(std::size_t i=0;i<mesh.n_cells();++i){ux(i)=U.component_data(0)[i];uy(i)=U.component_data(1)[i];}
    const auto gux=gauss_gradient_with_boundary(ux,mesh,geometry,gx);
    const auto guy=gauss_gradient_with_boundary(uy,mesh,geometry,gy);

    const auto& patch=mesh.boundary().patch(mesh.boundary().find("airfoil"));
    double Fx=0,Fy=0,Mz=0,Fxp=0,Fyp=0,Fxv=0,Fyv=0;
    const double q=0.5*rho*Uinf*Uinf;
    for(const auto f:patch.face_ids){
        const std::size_t o=mesh.ownership().owner(f);
        const auto Sf=geometry.face_area_vectors[f];
        const auto n=Sf.normalized();
        const double pp=p(o);
        const double txx=2.0*rho*nu*gux.component_data(0)[o];
        const double tyy=2.0*rho*nu*guy.component_data(1)[o];
        const double txy=rho*nu*(gux.component_data(1)[o]+guy.component_data(0)[o]);
        const double tx=-(pp*n.x+txx*n.x+txy*n.y)*Sf.mag();
        const double ty=-(pp*n.y+txy*n.x+tyy*n.y)*Sf.mag();
        Fx+=tx; Fy+=ty;
        Fxp+=-pp*n.x*Sf.mag(); Fyp+=-pp*n.y*Sf.mag();
        Fxv+=-(txx*n.x+txy*n.y)*Sf.mag();
        Fyv+=-(txy*n.x+tyy*n.y)*Sf.mag();
        const auto r=geometry.face_centres[f]-Vec3{0.25,0,0};
        Mz+=r.x*ty-r.y*tx;
    }
    const double cd=Fx/q, cl=Fy/q, cm=Mz/q;
    const auto& last=solve.history.back();
    std::cout<<"NACA0012_RE1000_A0 level="<<level
             <<" Cd="<<cd<<" Cl="<<cl<<" Cm="<<cm
             <<" Cd_pressure="<<Fxp/q<<" Cd_viscous="<<Fxv/q
             <<" continuity="<<last.continuity_linf
             <<" iterations="<<solve.iterations<<"\n";
    return {cd,cl,cm,Fxp/q,Fxv/q,last.continuity_linf};
}

int main(int argc,char**argv)
{
    const bool quick=argc>=2 && std::string(argv[1])=="--quick";
    const std::size_t need=quick?1:3;
    if(argc!=static_cast<int>(need+(quick?2:1)))
        throw std::invalid_argument("usage: test_naca0012_qualification [--quick] mesh_n64.h5 [mesh_n128.h5 mesh_n256.h5]");
    std::vector<Forces> r;
    for(std::size_t i=0;i<need;++i) r.push_back(run(argv[i+(quick?2:1)],64u<<i));
    for(const auto& x:r) {
        if(!(x.continuity<1e-7)) throw std::runtime_error("NACA0012 continuity gate failed");
        if(!std::isfinite(x.cd) || !std::isfinite(x.cl)) throw std::runtime_error("NACA0012 force integration produced non-finite QoI");
        if(std::abs(x.cl)>2e-2) throw std::runtime_error("NACA0012 alpha=0 symmetry gate failed");
    }
    if(!quick) {
        const double d01=std::abs(r[1].cd-r[0].cd), d12=std::abs(r[2].cd-r[1].cd);
        std::cout<<"NACA0012 Cd mesh differences="<<d01<<","<<d12<<" literature Re=1000 alpha=0 Cd~0.12\n";
        if(!(d12<=d01)) throw std::runtime_error("NACA0012 drag is not stabilizing under refinement");
    }
    std::cout<<"NACA0012_QUALIFICATION: DIAGNOSTIC solver_converged forces_conservation_symmetry_reference_gate_pending\n";
}
