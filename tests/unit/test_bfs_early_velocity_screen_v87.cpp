#include "cfdx/physics/bfs_early_velocity_screen_v87.h"
#include <cmath>
#include <limits>
#include <stdexcept>
int main() {
    using cfdx::physics::screen_bfs_velocity_before_sst_v87;
    using cfdx::physics::evaluate_bfs_field_trust_v86;
    constexpr double u=44.2;
    for(double speed: {0.,55.,2*u}) {
        auto r=screen_bfs_velocity_before_sst_v87(speed,u);
        if(!r.may_continue || r.reason!="velocity_preflight_only_not_certified")return 1;
        if(!evaluate_bfs_field_trust_v86(speed,0.,u).accepted)return 2;
    }
    const auto rejected=screen_bfs_velocity_before_sst_v87(148.37,u);
    if(rejected.may_continue || rejected.reason!="velocity_envelope")return 3;
    // Passing early screen NEVER means the updated SST state is valid.
    if(!screen_bfs_velocity_before_sst_v87(60.,u).may_continue)return 4;
    if(evaluate_bfs_field_trust_v86(60.,648.15,u).accepted)return 5;
    for(double x: {-1.,std::numeric_limits<double>::quiet_NaN(),
                   std::numeric_limits<double>::infinity()}) {
        try { (void)screen_bfs_velocity_before_sst_v87(x,u); return 6; }
        catch(const std::invalid_argument&){}
        try { (void)screen_bfs_velocity_before_sst_v87(1.,x); return 7; }
        catch(const std::invalid_argument&){}
    }
    return 0;
}