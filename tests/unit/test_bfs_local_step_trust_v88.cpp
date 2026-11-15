#include "cfdx/physics/bfs_local_step_trust_v88.h"
#include <cmath>
#include <limits>
#include <stdexcept>
int main() {
    using cfdx::physics::check_bfs_local_step_v88;
    const double inlet=44.2, rho=1.2;
    const auto first=check_bfs_local_step_v88(14.3302284171,25313.24005,inlet,rho);
    if (!first.may_continue_sst || first.reason!="delta_velocity_preflight_only_not_certified")return 1;
    if (std::abs(first.max_delta_velocity_ratio-14.3302284171/inlet)>1e-12)return 2;
    const auto third=check_bfs_local_step_v88(84.8479545595,14193.07466,inlet,rho);
    if (third.may_continue_sst || third.reason!="local_velocity_increment_envelope")return 3;
    if (std::abs(third.pressure_increment_span_over_qin-12.10823101)>1e-5)return 4;
    // Case-relative diagnostic must be invariant under arbitrary constant
    // gauge shifts in both old and new pressure fields.
    const auto any_gauge=check_bfs_local_step_v88(15.,1000.,inlet,rho);
    if (!any_gauge.may_continue_sst)return 5;
    if (check_bfs_local_step_v88(0.5*inlet,0.,inlet,rho).may_continue_sst==false)return 6;
    if (check_bfs_local_step_v88(0.50001*inlet,0.,inlet,rho).may_continue_sst)return 7;
    for (double invalid: {-1.,std::numeric_limits<double>::infinity(),
                          std::numeric_limits<double>::quiet_NaN()}) {
        try { (void)check_bfs_local_step_v88(invalid,0.,inlet,rho);return 8; }
        catch(const std::invalid_argument&) {}
        try { (void)check_bfs_local_step_v88(0.,invalid,inlet,rho);return 9; }
        catch(const std::invalid_argument&) {}
        try { (void)check_bfs_local_step_v88(0.,0.,invalid,rho);return 10; }
        catch(const std::invalid_argument&) {}
    }
    return 0;
}
