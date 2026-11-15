#pragma once
// V87: necessary-only pre-SST velocity screening for the *opt-in* BFS V86
// campaign. If a trial already exceeds V86's speed bound before SST, a
// positive k/omega update cannot repair the velocity in the existing
// split U/p -> SST trial. Reject before the expensive transport. This
// is not a residual certificate and NEVER accepts a trial on its own.
#include "cfdx/physics/bfs_field_trust_guard_v86.h"
#include <cmath>
#include <stdexcept>
#include <string>
namespace cfdx::physics {
struct BfsEarlyVelocityV87 {
    bool may_continue=false;
    std::string reason;
    double speed_ratio=0.0;
};
inline BfsEarlyVelocityV87 screen_bfs_velocity_before_sst_v87(
    double max_speed, double inlet_speed)
{
    // Reuse the SAME velocity envelope as the complete V86 field gate.
    // k=0 merely disables the not-yet-updated SST check: V86 always
    // rechecks the full updated k,omega and global physical residual.
    const auto v86=evaluate_bfs_field_trust_v86(max_speed,0.0,inlet_speed);
    BfsEarlyVelocityV87 result;
    result.may_continue=v86.accepted;
    result.speed_ratio=v86.speed_ratio;
    result.reason=v86.accepted ? "velocity_preflight_only_not_certified"
                               : v86.reason;
    return result;
}
}