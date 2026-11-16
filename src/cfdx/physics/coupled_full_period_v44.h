#pragma once
// V44: strict, bounded full-field initialization from a four-phase orbit.
// The phase-average eddy viscosity is NOT a valid SST closure. The caller
// MUST recompute strain, F2 and nu_t for the trial U/k/omega before assembly.
// The entire trial is then reassembled and subjected to the original physical
// U/p/k/omega/Rhie-Chow gates; this helper does not certify any flow solution.
#include <cmath>
#include <stdexcept>

namespace cfdx::physics {
struct CoupledSSTPhaseCandidateV44 {
    double k = 0.0;
    double omega = 0.0;
};
inline CoupledSSTPhaseCandidateV44 make_coupled_sst_phase_candidate_v44(
    double k_average, double omega_average) {
    if (!std::isfinite(k_average) || k_average < 0.0 ||
        !std::isfinite(omega_average) || !(omega_average > 0.0))
        throw std::invalid_argument("V44 four-phase SST predictor has invalid k/omega");
    return {k_average, omega_average};
}
} // namespace cfdx::physics
