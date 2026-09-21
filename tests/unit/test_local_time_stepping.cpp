#include "cfdx/physics/local_time_stepping.h"
#include "common/test_harness.h"
#include <limits>
using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh unit_mesh() {
    Mesh m;
    m.points().resize(8);
    const double p[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for(std::size_t i=0;i<8;++i)m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1});m.faces().push_face({4,5,6,7});m.faces().push_face({0,1,5,4});m.faces().push_face({3,7,6,2});m.faces().push_face({0,4,7,3});m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);for(std::size_t i=0;i<6;++i){m.ownership().set_owner(i,0);m.ownership().set_neighbour(i,FaceOwnership::BOUNDARY);}m.cells().push_cell({0,1,2,3,4,5});return m;
}
int main(){
 run_case("local_dt_matches_mass_flux_cfl",[](){auto m=unit_mesh();Field<double,Location::FACE> mf(6,"mf","kg/s",1);mf.fill(1.0);Field<double,Location::CELL> rho(1,"rho","kg/m3",1);rho(0)=2.0;Field<double,Location::CELL> dt;LocalTimeStepControls c;c.cfl=0.5;compute_local_time_step(m,mf,rho,dt,c);EXPECT_NEAR(dt(0),1.0/6.0,1e-14);});
 run_case("local_dt_rejects_invalid_density",[](){auto m=unit_mesh();Field<double,Location::FACE> mf(6,"mf","kg/s",1);mf.fill(1.0);Field<double,Location::CELL> rho(1,"rho","kg/m3",1);rho(0)=0.0;Field<double,Location::CELL> dt;EXPECT_THROW(compute_local_time_step(m,mf,rho,dt),std::invalid_argument);});
 run_case("local_dt_rejects_nonfinite_flux",[](){auto m=unit_mesh();Field<double,Location::FACE> mf(6,"mf","kg/s",1);mf.fill(1.0);mf(0)=std::numeric_limits<double>::quiet_NaN();Field<double,Location::CELL> rho(1,"rho","kg/m3",1);rho(0)=1.0;Field<double,Location::CELL> dt;EXPECT_THROW(compute_local_time_step(m,mf,rho,dt),std::invalid_argument);});
 return run_all();}
