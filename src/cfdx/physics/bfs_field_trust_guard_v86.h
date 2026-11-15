#pragma once
// V86 experimental step safeguard for the *specific* TH03 BFS campaign.
// This is a trust-region test, not a universal physical law or convergence
// certificate. It evaluates a fully reassembled private U,k,omega candidate;
// a failure triggers the existing nonlinear transaction rollback/retry.
#include <cmath>
#include <stdexcept>
#include <string>
namespace cfdx::physics {
struct BfsFieldTrustV86 {
    bool accepted=false;
    std::string reason;
    double speed_ratio=0.0;
    double k_fraction_of_uin_squared=0.0;
};
inline BfsFieldTrustV86 evaluate_bfs_field_trust_v86(
    double maximum_velocity_magnitude,
    double maximum_turbulent_kinetic_energy,
    double inlet_speed,
    double speed_ratio_limit=2.0,
    double k_over_u2_limit=0.03)
{
    if (!(std::isfinite(maximum_velocity_magnitude) && maximum_velocity_magnitude>=0.0 &&
          std::isfinite(maximum_turbulent_kinetic_energy) && maximum_turbulent_kinetic_energy>=0.0 &&
          std::isfinite(inlet_speed) && inlet_speed>0.0 &&
          std::isfinite(speed_ratio_limit) && speed_ratio_limit>1.0 &&
          std::isfinite(k_over_u2_limit) && k_over_u2_limit>0.0))
        throw std::invalid_argument("V86 nonfinite/invalid field trust-region data");
    BfsFieldTrustV86 result;
    result.speed_ratio=maximum_velocity_magnitude/inlet_speed;
    result.k_fraction_of_uin_squared=
        (maximum_turbulent_kinetic_energy/inlet_speed)/inlet_speed;
    if (result.speed_ratio>speed_ratio_limit) {
        result.reason="velocity_envelope";
        return result;
    }
    if (result.k_fraction_of_uin_squared>k_over_u2_limit) {
        result.reason="turbulent_energy_envelope";
        return result;
    }
    result.accepted=true;
    result.reason="provisional_step_within_field_trust_region";
    return result;
}
} // namespace cfdx::physics
