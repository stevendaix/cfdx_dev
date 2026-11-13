#include "cfdx/physics/bfs_wall_regime_v74.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void need(bool b,const char* msg){if(!b)throw std::runtime_error(msg);}
int main(){
    const std::vector<BfsWallSampleV74> s={
        {2,0,1e-5,0.,1.5e-5,.2,100.,.1},
        {3,0,3e-4,5.,1.5e-5,.2,100.,.1},
        {4,12,1e-3,40.,1.5e-5,.2,100.,.1},
    };
    const auto a=bfs_wall_regime_v74(s);
    need(a.faces.size()==3 && a.corner_owner_count==1,
         "V74 corner owner or wall faces omitted");
    need(a.wall_resolved>=1 && a.max_yplus_velocity>a.min_yplus_velocity,
         "V74 wall regime distribution is inconsistent");
    need(a.wall_resolved+a.buffer_or_transition+a.log_layer+
         a.out_of_model_range==a.faces.size(),
         "V74 yplus regimes are not a partition");
    need(a.faces[0].yplus_velocity==0. && a.faces[0].nut_spalding==0.,
         "V74 stationary wall no-shear reference wrong");
    for(const auto& x:a.faces) {
        need(std::abs(x.viscous_omega_60-10.*x.viscous_omega_6)<1e-9*
             std::max(1.,x.viscous_omega_60),
             "V74 low-Re 6 vs NASA 60 omega incorrect");
        need(x.omega_step>0. && x.omega_blended>0.,
             "V74 wall omega branch invalid");
    }
    need(!a.same_operator_wall_triplet_certified && !a.physical_bfs_validated,
         "V74 wall reference falsely claims OpenFOAM or BFS validation");
    bool rejected=false;
    try{auto bad=s;bad[1].distance=0.;(void)bfs_wall_regime_v74(bad);}
    catch(const std::invalid_argument&){rejected=true;}
    need(rejected,"V74 accepted zero wall normal distance");
    rejected=false;
    try{(void)bfs_wall_regime_v74({});}
    catch(const std::invalid_argument&){rejected=true;}
    need(rejected,"V74 accepted a case without wall faces");
    std::cout<<"V74_WALL_REGIME_PASS faces="<<a.faces.size()
             <<" corners="<<a.corner_owner_count<<"\n";
}
