#include "cfdx/physics/bfs_field_trust_guard_v86.h"
#include <cmath>
#include <limits>
#include <stdexcept>
int main() {
    using cfdx::physics::evaluate_bfs_field_trust_v86;
    const double u=44.2;
    if(!evaluate_bfs_field_trust_v86(55.08,.04911,u).accepted)return 1;
    const auto bad_u=evaluate_bfs_field_trust_v86(148.37,27.43,u);
    if(bad_u.accepted || bad_u.reason!="velocity_envelope")return 2;
    const auto bad_k=evaluate_bfs_field_trust_v86(75.6,648.14,u);
    if(bad_k.accepted || bad_k.reason!="turbulent_energy_envelope")return 3;
    if(!evaluate_bfs_field_trust_v86(2*u,.03*u*u,u).accepted)return 4;
    for(const double x:{-1.,std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::infinity()}) {
       try { (void)evaluate_bfs_field_trust_v86(x,1.,u);return 5; }
       catch(const std::invalid_argument&) {}
       try { (void)evaluate_bfs_field_trust_v86(1.,x,u);return 6; }
       catch(const std::invalid_argument&) {}
    }
    return 0;
}
