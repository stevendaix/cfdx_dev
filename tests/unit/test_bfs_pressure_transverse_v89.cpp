@@ -0,0 +1,55 @@
#include "cfdx/physics/bfs_pressure_transverse_trust_v89.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
int main() {
    using namespace cfdx::physics;
    constexpr double u=44.2,rho=1.2;
    const double q=0.5*rho*u*u;
    const auto first=check_bfs_pressure_transverse_pre_v89(0.28*u,22*q,u,rho,true);
    if(!first.may_continue_sst || first.reason!="startup_baseline_not_certified")return 1;
    const auto old2=check_bfs_pressure_transverse_pre_v89(0.07*u,21.6*q,u,rho,false);
    if(old2.may_continue_sst || old2.reason!="gauge_invariant_pressure_increment")return 2;
    const auto old3=check_bfs_pressure_transverse_pre_v89(0.153*u,14*q,u,rho,false);
    if(old3.may_continue_sst || old3.reason!="transverse_momentum_increment")return 3;
    const auto admissible=check_bfs_pressure_transverse_pre_v89(0.10*u,10*q,u,rho,false);
    if(!admissible.may_continue_sst ||
       admissible.reason!="pressure_transverse_preflight_not_certified")return 4;
    // Relative-to-q spans are unchanged by *independent* gauge shifts of
    // old and trial fields. The range of a vector difference, not |p|,
    // is the measured quantity. Numerical verification on explicit vectors:
    const double oldp[]={-5*q,3*q,4*q};
    const double newp[]={-1*q,5*q,8*q};
    double min1=1e50,max1=-1e50,min2=1e50,max2=-1e50;
    for (int i=0;i<3;++i) {
        const double d=newp[i]-oldp[i];
        const double shifted=(newp[i]+17*q)-(oldp[i]-9*q);
        min1=std::min(min1,d);max1=std::max(max1,d);
        min2=std::min(min2,shifted);max2=std::max(max2,shifted);
    }
    if(std::abs((max1-min1)-(max2-min2))>1e-9)return 5;
    if(check_bfs_pressure_transverse_post_v89(6*q,2*q,true,0.9,0.8,true,
        true,true,1e-16,1e-11,u,rho).may_commit==false)return 6;
    if(check_bfs_pressure_transverse_post_v89(6.1*q,2*q,true,0.9,0.8,true,
        true,true,1e-16,1e-11,u,rho).may_commit)return 7;
    if(check_bfs_pressure_transverse_post_v89(3*q,2*q,true,1.6,1.0,true,
        true,true,1e-16,1e-11,u,rho).may_commit)return 8;
    if(check_bfs_pressure_transverse_post_v89(3*q,2*q,true,0.5,1.0,true,
        true,true,1e-8,1e-11,u,rho).may_commit)return 9;
    if(check_bfs_pressure_transverse_post_v89(3*q,2*q,true,0.5,1.0,true,
        false,true,1e-16,1e-11,u,rho).may_commit)return 10;
    if(check_bfs_pressure_transverse_post_v89(3*q,2*q,true,0.5,1.0,true,
        true,false,1e-16,1e-11,u,rho).may_commit)return 11;
    if(check_bfs_pressure_transverse_post_v89(22*q,0,false,0.9,0,false,
        true,true,1e-16,1e-11,u,rho).may_commit==false)return 12;
    for(double v:{-1.,std::numeric_limits<double>::quiet_NaN(),
                 std::numeric_limits<double>::infinity()}) {
        try{(void)check_bfs_pressure_transverse_pre_v89(0.,v,u,rho,false);return 13;}
        catch (const std::invalid_argument&) {}
        try{(void)check_bfs_pressure_transverse_post_v89(v,0,false,0,0,false,
            true,true,0,1e-11,u,rho);return 14;}
        catch (const std::invalid_argument&) {}
    }
    return 0;
}
