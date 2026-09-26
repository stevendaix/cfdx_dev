#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/io/gmsh/gmsh_importer.h"
#include "cfdx/core/geometry/face_geometry.h"
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
constexpr double D=1.0, RHO=1.0, UINF=1.0, MU=0.01, NU=0.01, RE=100.0;
constexpr double AREA=M_PI*D*D/4.0, CDREF=1.0895;
std::string quote(const std::string&s){std::string q="'";for(char c:s)q+=c=='\''?"'\\''":std::string(1,c);return q+"'";}
void mesh(const std::string&out){
#ifndef CFDX_SOURCE_DIR
 throw std::runtime_error("CFDX_SOURCE_DIR is not defined");
#else
#ifdef CFDX_PYTHON_EXECUTABLE
 const std::string py=CFDX_PYTHON_EXECUTABLE;
#else
 const std::string py="python3";
#endif
 const auto script=std::filesystem::path(CFDX_SOURCE_DIR)/"scripts"/"generate_vmfl036_mesh.py";
 const std::string cmd=quote(py)+" "+quote(script.string())+" "+quote(out);
 std::cout<<"VMFL036_MESH_COMMAND: "<<cmd<<"\n";
 if(std::system(cmd.c_str())!=0)throw std::runtime_error("Gmsh mesh generation failed");
#endif
}
ScalarBoundaryConditions pbc(){ScalarBoundaryConditions b;for(const char*n:{"sphere","inlet","outlet","farfield"})b[n]={ScalarBoundaryType::ZERO_GRADIENT,0,0};return b;}
ScalarBoundaryConditions velocity_scalar_bcs(){ScalarBoundaryConditions b;b["sphere"]={ScalarBoundaryType::FIXED_VALUE,0,0};b["inlet"]={ScalarBoundaryType::FIXED_VALUE,UINF,0};b["outlet"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};b["farfield"]={ScalarBoundaryType::FIXED_VALUE,UINF,0};return b;}
}
int main(){try{
 const auto msh=(std::filesystem::temp_directory_path()/"cfdx_vmfl036_re100.msh").string(); mesh(msh);
 Mesh m;if(!cfdx::io::gmsh::import_gmsh_mesh(msh,m))throw std::runtime_error("Gmsh import failed");
 const auto topo=m.topo_validate();if(!topo.ok)throw std::runtime_error("invalid imported topology");
 for(const char*n:{"sphere","inlet","outlet","farfield"})if(!m.boundary().has_patch(n))throw std::runtime_error(std::string("missing patch ")+n);
 Field<double,Location::CELL> U(m.n_cells(),"U","m/s",3),p(m.n_cells(),"p","Pa",1);U.fill(0);p.fill(0);for(std::size_t c=0;c<m.n_cells();++c)U.component_data(0)[c]=UINF;
 VelocityBoundaryConditions ubc;ubc["sphere"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};ubc["inlet"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{UINF,0,0}};ubc["outlet"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};ubc["farfield"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{UINF,0,0}};
 IncompressibleSolverControls c;c.algorithm=PressureVelocityAlgorithm::SIMPLE;c.density=RHO;c.kinematic_viscosity=NU;c.linear_max_iterations=3000;c.linear_tolerance=1e-9;c.convergence.max_iterations=1500;c.convergence.relative_tolerance=1e-7;c.convergence.continuity_tolerance=1e-7;c.coupling.alpha_u=0.7;c.coupling.alpha_p=0.3;c.use_bounded_convection=true;c.convection_scheme=ConvectionScheme::UPWIND;c.pressure_reference_cell=0;c.pressure_reference_value=0;c.diagnostics.iteration_trace=true;c.diagnostics.iteration_trace_frequency=50;
 std::cout<<std::setprecision(12)<<"VMFL036 CFDX START Re="<<RE<<" cells="<<m.n_cells()<<" faces="<<m.n_faces()<<"\n";
 const auto r=solve_steady_incompressible(m,U,p,ubc,pbc(),c);if(!r.converged)throw std::runtime_error("solver did not converge");
 const auto geom=build_fv_geometry(m);Field<double,Location::CELL> ux(m.n_cells(),"Ux","m/s",1),uy(m.n_cells(),"Uy","m/s",1),uz(m.n_cells(),"Uz","m/s",1);for(std::size_t i=0;i<m.n_cells();++i){ux(i)=U.component_data(0)[i];uy(i)=U.component_data(1)[i];uz(i)=U.component_data(2)[i];}
 const auto gx=gauss_gradient_with_boundary(ux,m,geom,velocity_scalar_bcs());const auto gy=gauss_gradient_with_boundary(uy,m,geom,velocity_scalar_bcs());const auto gz=gauss_gradient_with_boundary(uz,m,geom,velocity_scalar_bcs());
 double fp=0,fv=0;for(const auto f:m.boundary().patch(m.boundary().find("sphere")).face_ids){const auto o=m.ownership().owner(f);const auto&v=m.faces().vertices();const auto*off=m.faces().offsets_data();const auto fg=compute_face_geometry_oriented(m.points().x_data(),m.points().y_data(),m.points().z_data(),v.data(),off[f],off[f+1]-off[f],geom.cell_centres[o]);fp+=p(o)*fg.Sf.x;const double txx=2*MU*gx.component_data(0)[o],txy=MU*(gx.component_data(1)[o]+gy.component_data(0)[o]),txz=MU*(gx.component_data(2)[o]+gz.component_data(0)[o]);fv-=(txx*fg.Sf.x+txy*fg.Sf.y+txz*fg.Sf.z);}
 const double denom=0.5*RHO*UINF*UINF*AREA,ft=fp+fv,cd=ft/denom,err=std::abs(cd-CDREF)/CDREF;const auto&h=r.history.back();
 std::cout<<"VMFL036 RESULT iterations="<<r.iterations<<" continuity="<<h.continuity_linf<<" continuity_norm="<<h.continuity_normalized<<" momentum_eq="<<h.momentum_equation_residual<<" pressure_linear_iterations="<<h.pressure_linear_iterations<<"\n";
 std::cout<<"VMFL036 DRAG F_pressure="<<fp<<" F_viscous="<<fv<<" F_total="<<ft<<" Cd_pressure="<<fp/denom<<" Cd_viscous="<<fv/denom<<" Cd_total="<<cd<<" Cd_reference="<<CDREF<<" Cd_relative_error="<<err<<"\n";
 if(!std::isfinite(cd)||!std::isfinite(ft))throw std::runtime_error("non-finite drag");
 std::cout<<"VMFL036_CFDX_RESULT: PASS (physical run; literature accuracy remains diagnostic)\n";std::error_code ec;std::filesystem::remove(msh,ec);return 0;
 }catch(const std::exception&e){std::cerr<<"VMFL036_CFDX_RESULT: FAIL: "<<e.what()<<"\n";return 1;}}
