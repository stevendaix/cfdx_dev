#pragma once
// V45: finite, transaction-safe continuation window following an observed
// four-cycle of fully committed COUPLED/SST fields. This is a Picard update
// policy, not a modification of the frozen linear system or an SST model.
// No state changes until the nonlinear transaction commits.
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {

struct CoupledPhaseEscapeTrialV45 {
    bool active = false;
    bool triggered = false;
    std::size_t remaining_including_this = 0;
    double beta = 1.0;
};

inline void validate_coupled_phase_escape_v45(double beta, std::size_t steps) {
    if (!std::isfinite(beta) || beta < 0.0 || beta > 1.0 ||
        steps < 1 || steps > 128)
        throw std::invalid_argument("V45 cycle-escape damping requires beta in [0,1] and 1..128 steps");
}

class CoupledPhaseEscapeLedgerV45 {
public:
    CoupledPhaseEscapeTrialV45 trial(std::size_t iteration,
                                     bool four_cycle_from_committed_fields,
                                     double requested_beta,
                                     std::size_t steps,
                                     double retry_scale) const {
        validate_coupled_phase_escape_v45(requested_beta, steps);
        if (iteration == 0)
            throw std::invalid_argument("V45 cycle-escape iteration must be positive");
        CoupledPhaseEscapeTrialV45 result;
        if (requested_beta == 0.0) return result;
        if (!std::isfinite(retry_scale) || retry_scale <= 0.0 || retry_scale > 1.0)
            throw std::invalid_argument("V45 cycle-escape retry factor invalid");
        if (four_cycle_from_committed_fields &&
            iteration > std::numeric_limits<std::size_t>::max()-steps+1)
            throw std::invalid_argument("V45 cycle-escape iteration overflow");
        const bool continuing = iteration <= active_through_;
        if (!four_cycle_from_committed_fields && !continuing) return result;
        result.active = true;
        result.triggered = four_cycle_from_committed_fields;
        result.remaining_including_this = four_cycle_from_committed_fields ?
            steps : active_through_ - iteration + 1;
        result.beta = requested_beta * retry_scale;
        if (!(result.beta > 0.0 && result.beta <= 1.0) || !std::isfinite(result.beta))
            throw std::invalid_argument("V45 cycle-escape trial beta underflow/overflow");
        return result;
    }

    // Called strictly AFTER successful commit; a rejected candidate must not
    // advance the window. The iteration bound is already tested in trial().
    void commit(const CoupledPhaseEscapeTrialV45& accepted,
                std::size_t iteration, std::size_t steps) noexcept {
        if (accepted.active && accepted.triggered)
            active_through_ = iteration + steps - 1;
    }

    std::size_t active_through() const noexcept { return active_through_; }

private:
    std::size_t active_through_ = 0;
};

} // namespace cfdx::physics
