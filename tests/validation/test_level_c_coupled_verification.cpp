#include "cfdx/physics/turbulence_solver.h"
#include "cfdx/physics/energy_solver.h"
#include "cfdx/physics/radiation_solver.h"
#include "cfdx/physics/cht_solver.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::core; using namespace cfdx::physics;
namespace {
Mesh cube(const char* name){Mesh m;m.points().resize(8);const double p[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};for(size_t i=0;i<8;++i)m.points().set(i,p[i][0],p[i][1],p[i][2]);m.faces().push_face({0,3,2,1});m.faces().push_face({4,5,6,7});m.faces().push_face({0,1,5,4});m.faces().push_face({3,7,6,2});m.faces().push_face({0,4,7,3});m.faces().push_face({1,2,6,5});m.ownership().resize(6);for(size_t f=0;f<6;++f){m.ownership().set_owner(f,0);m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);}m.cells().push_cell({0,1,2,3,4,5});Patch p;p.name=name;p.type=PatchType::WALL;p.face_ids={0,1,2,3,4,5};m.boundary().add_patch(p);return m;}
std::vector<Direction> dirs(){return{{1,0,0,2*M_PI/3},{-1,0,0,2*M_PI/3},{0,1,0,2*M_PI/3},{0,-1,0,2*M_PI/3},{0,0,1,2*M_PI/3},{0,0,-1,2*M_PI/3}};}
}
int main(){try{
 Mesh m=cube("wall");auto g=build_fv_geometry(m);Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1);phi.fill(0);
 Field<double,Location::CELL> k(1,"k","m2/s2",1),e(1,"e","m2/s3",1),S(1,"S","1/s",1);k(0)=.09;e(0)=.03;S(0)=.02;
 TurbulenceTransportControls tc;tc.model=TurbulenceModel::KEPSILON;tc.density=1;tc.molecular_viscosity=1e-3;ScalarBoundaryConditions bc;bc["wall"]={ScalarBoundaryType::FIXED_VALUE,1e-6,0};
 auto tr=solve_kepsilon_transport(m,g,phi,k,e,S,tc,bc,bc,100,1e-8);if(tr.iterations==0||k(0)<tc.k_min||e(0)<tc.epsilon_min)throw std::runtime_error("turbulence-energy positivity");
 Field<double,Location::CELL> T(1,"T","K",1),src(1,"src","W/m3",1),G(1,"G","W/m2",1);T(0)=800;src(0)=0;G(0)=0;ScalarBoundaryConditions tbc;tbc["wall"]={ScalarBoundaryType::FIXED_VALUE,blackbody_intensity(800),0};
 RadiationEnergyCouplingControls rc;rc.radiation.absorption=.5;rc.radiation.max_iterations=100;rc.radiation.tolerance=1e-10;rc.energy.conductivity=1;rc.energy.density=1;rc.energy.cp=1000;rc.energy.max_iterations=100;rc.energy.tolerance=1e-10;rc.max_outer_iterations=50;rc.tolerance=1e-8;
 auto rr=solve_radiation_energy_coupled(m,g,phi,T,src,G,dirs(),rc,tbc);if(!rr.converged||!std::isfinite(T(0)))throw std::runtime_error("radiation-energy coupling");
 Mesh m1=cube("interface1"),m2=cube("interface2");auto g1=build_fv_geometry(m1),g2=build_fv_geometry(m2);Field<double,Location::FACE> f1(m1.n_faces(),"f1","kg/s",1),f2(m2.n_faces(),"f2","kg/s",1);f1.fill(0);f2.fill(0);Field<double,Location::CELL>T1(1,"T1","K",1),T2(1,"T2","K",1),s1(1,"s1","W/m3",1),s2(1,"s2","W/m3",1);T1(0)=500;T2(0)=300;s1(0)=0;s2(0)=0;EnergySolverControls ec;ec.conductivity=1;ec.max_iterations=100;ec.tolerance=1e-10;ChtInterfaceControls cc;cc.region1_patch="interface1";cc.region2_patch="interface2";cc.conductivity1=1;cc.conductivity2=1;cc.tolerance=1e-8;cc.matching_tolerance=1e-12;auto cr=solve_two_region_cht(m1,g1,m2,g2,f1,f2,T1,T2,s1,s2,ec,ec,cc);if(!cr.converged||cr.interface_imbalance>cc.tolerance)throw std::runtime_error("CHT coupling");
 std::cout<<"LEVEL_C_COUPLED_VERIFICATION: PASS\n";return 0;
}catch(const std::exception&e){std::cerr<<"LEVEL_C_COUPLED_VERIFICATION: FAIL: "<<e.what()<<"\n";return 1;}}
