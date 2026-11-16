// V40: true production FV P1 -> nonconformal CHT energy feedback, with shared
// watts and fixed old physical time level. No external radiation CHT model.
#include "cfdx/physics/cht_p1_coupling_v40.h"
#include "thermal_3d_mesh_v36.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace cfdx::physics;
using namespace cfdx::physics::v37;
using namespace cfdx::physics::v40;
using namespace cfdx::core;
namespace {
void check(bool yes,const char* why){if(!yes)throw std::runtime_error(why);}
struct ChtCase {
    Mesh m1,m2;
    FvGeometry g1,g2;
    Field<double,Location::FACE> phi1,phi2;
    Field<double,Location::CELL> T1,T2,s1,s2,G,Q;
    EnergySolverControls e1,e2;
    ChtInterfaceControls ic;
    ChtP1Controls pc;
    ScalarBoundaryConditions bc1,bc2;
    explicit ChtCase(double warp=0):
        m1(cfdx::testing::thermal_v36::make_mesh(3,4,3,warp,0.,1.)),
        m2(cfdx::testing::thermal_v36::make_mesh(4,3,5,warp,1.,2.)),
        g1(build_fv_geometry(m1)),g2(build_fv_geometry(m2)),
        phi1(m1.n_faces(),"phi1","kg/s",1),phi2(m2.n_faces(),"phi2","kg/s",1),
        T1(m1.n_cells(),"T1","K",1),T2(m2.n_cells(),"T2","K",1),
        s1(m1.n_cells(),"q1","W/m3",1),s2(m2.n_cells(),"q2","W/m3",1),
        G(m1.n_cells(),"G","W/m2",1),Q(m1.n_cells(),"Q","W/m3",1) {
        phi1.fill(0.);phi2.fill(0.);s1.fill(0.);s2.fill(0.);
        T1.fill(380.);T2.fill(310.);G.fill(4*blackbody_emissive_power(380.));Q.fill(0.);
        e1.conductivity=2.;e2.conductivity=.5;
        e1.density=e2.density=1.;e1.cp=e2.cp=100.;
        e1.relaxation=e2.relaxation=1.;
        e1.tolerance=e2.tolerance=1.e-10;
        ic.region1_patch="top";ic.region2_patch="bottom";
        ic.conductivity1=2.;ic.conductivity2=.5;
        ic.tolerance=1.e-10;ic.max_iterations=50;
        ic.area_relative_tolerance=1.e-10;
        bc1["bottom"]={ScalarBoundaryType::FIXED_VALUE,420.,0};
        bc2["top"]={ScalarBoundaryType::FIXED_VALUE,300.,0};
        for(const char* wall:{"west","east","south","north"}){
            bc1[wall]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
            bc2[wall]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
        }
        pc.absorption=.6;pc.scattering=1.3;pc.emissivity=.75;
        pc.max_outer_iterations=250;pc.outer_relaxation=0.2;pc.coupling_tolerance=1.e-10;
        pc.p1.physical_tolerance=1.e-10;
        pc.p1.balance_tolerance=1.e-10;
    }
    ChtP1Report run(){return solve_nonconformal_cht_p1(m1,g1,m2,g2,phi1,phi2,T1,T2,s1,s2,G,Q,e1,e2,ic,pc,bc1,bc2);}
    double total_energy()const {
        long double energy=0.;
        for(std::size_t i=0;i<T1.size();++i)
            energy+=static_cast<long double>(g1.cell_volumes[i])*e1.density*e1.cp*T1(i);
        for(std::size_t i=0;i<T2.size();++i)
            energy+=static_cast<long double>(g2.cell_volumes[i])*e2.density*e2.cp*T2(i);
        return static_cast<double>(energy);
    }
};
void stationary(double warp){
    ChtCase c(warp), control=c;
    const auto base=solve_nonconformal_cht(control.m1,control.g1,control.m2,control.g2,
        control.phi1,control.phi2,control.T1,control.T2,
        control.s1,control.s2,control.e1,control.e2,control.ic,control.bc1,control.bc2);
    check(base.converged,"baseline CHT failed");
    const auto r=c.run();
    std::cout<<std::scientific<<std::setprecision(12)
       <<"CHT_P1_V40 warp="<<warp<<" accepted="<<r.converged<<" it="<<r.outer_iterations
       <<" to_solid="<<r.radiation_power_to_solid_W
       <<" from_fluid="<<r.radiation_power_from_fluid_W
       <<" rad_balance="<<r.radiative_balance
       <<" dT="<<r.max_temperature_relative<<" ds="<<r.max_radiative_source_relative
       <<" cht_cell="<<r.cht_max_cell_residual<<"\n";
    check(r.converged,"CHT P1 outer fixed-point did not converge");
    check(r.interface_overlaps>12,"nonconformal overlaps not exercised");
    check(r.radiative_balance<1.e-10,"radiative power not conserved");
    check(r.cht_max_cell_residual<1.e-10,"CHT residual not closed");
    check(std::abs(r.radiation_power_to_solid_W)>1.e-4,
          "coupled radiative exchange not activated");
    double delta=0;
    for(std::size_t i=0;i<c.T1.size();++i)
        delta=std::max(delta,std::abs(c.T1(i)-control.T1(i)));
    for(std::size_t i=0;i<c.T2.size();++i)
        delta=std::max(delta,std::abs(c.T2(i)-control.T2(i)));
    check(delta>1.e-6,"CHT is insensitive to participating P1 radiation");
}
void rejected_transaction(){
    ChtCase c(.008);
    const auto t1=c.T1,t2=c.T2,g=c.G,q=c.Q;
    c.pc.max_outer_iterations=1;
    auto report=c.run();
    check(!report.converged,"one nonlinear iteration certified incorrectly");
    for(std::size_t i=0;i<c.T1.size();++i)
        check(c.T1(i)==t1(i)&&c.G(i)==g(i)&&c.Q(i)==q(i),
              "rejected P1-CHT changed fluid/radiation");
    for(std::size_t i=0;i<c.T2.size();++i)
        check(c.T2(i)==t2(i),"rejected P1-CHT changed solid");
    std::cout<<"CHT_P1_V40_ROLLBACK PASS\n";
}
void transient(){
    ChtCase c(.008);
    c.e1.dt=c.e2.dt=0.005;
    c.bc1.erase("bottom");c.bc2.erase("top");
    c.T1.fill(390.);c.T2.fill(320.);
    c.G.fill(4.*blackbody_emissive_power(390.));
    const auto energy0=c.total_energy();
    for(int step=0;step<3;++step){
        const auto r=c.run();
        std::cout<<"CHT_P1_V40_TRANSIENT step="<<step+1
          <<" converged="<<r.converged<<" iterations="<<r.outer_iterations
          <<" drift="<<std::abs(c.total_energy()-energy0)/energy0<<"\n";
        check(r.converged,"transient CHT P1 not converged");
        check(std::abs(c.total_energy()-energy0)/energy0<1.e-10,
              "CHT P1 transient total energy not conserved");
    }
}
}
int main(){
    try {stationary(0.);stationary(.008);rejected_transaction();transient();
         std::cout<<"CHT_P1_V40_CERTIFICATE PASS\n";return 0;}
    catch(const std::exception& e){std::cerr<<"CHT_P1_V40_CERTIFICATE FAIL: "<<e.what()<<"\n";return 1;}
}
