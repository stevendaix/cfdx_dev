#include "cfdx/physics/pressure_velocity_algorithms.h"
#include "cfdx/physics/timestep_control.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {

struct AcceptedState {
    std::size_t iteration = 0;
    double residual_u = 0.0;
    double residual_p = 0.0;
    double alpha_u = 0.0;
    double alpha_p = 0.0;
};

void append_accepted(std::vector<AcceptedState>& history, AcceptedState state)
{
    if (!history.empty() && state.iteration <= history.back().iteration)
        throw std::invalid_argument("accepted iterations must be strictly increasing");
    history.push_back(state);
}

} // namespace

int main()
{
    run_case("n7_causal_overshoot_retry_preserves_history_and_channel_state", [] {
        AdaptiveRelaxationControls controls;
        controls.enabled = true;
        controls.min_alpha_u = 0.6;
        controls.max_alpha_u = 0.8;
        controls.min_alpha_p = 0.25;
        controls.max_alpha_p = 0.35;
        controls.increase_step = 0.005;
        controls.decrease_step = 0.02;
        controls.severe_degradation_ratio = 2.0;
        controls.recovery_cooldown_windows = 2;
        controls.adaptation_window = 4;
        controls.required_trend_windows = 2;
        validate_adaptive_relaxation_controls(controls);

        std::vector<AcceptedState> history;
        append_accepted(history, {10, 1.0, 1.0, 0.700, 0.300});

        std::size_t improvement_streak = 0;
        std::size_t degradation_streak = 0;
        double alpha_u = 0.700;
        alpha_u = adapt_relaxation_factor_windowed(
            alpha_u, 1.0, 0.9, controls.min_alpha_u, controls.max_alpha_u,
            controls, improvement_streak, degradation_streak);
        EXPECT_NEAR(alpha_u, 0.700, 1e-14);
        alpha_u = adapt_relaxation_factor_windowed(
            alpha_u, 1.0, 0.9, controls.min_alpha_u, controls.max_alpha_u,
            controls, improvement_streak, degradation_streak);
        EXPECT_NEAR(alpha_u, 0.705, 1e-14);

        const double candidate_residual_u = 2.1;
        const bool severe_u = adaptive_relaxation_severe_degradation(
            history.back().residual_u, candidate_residual_u,
            history.back().alpha_u, alpha_u, controls);
        EXPECT_TRUE(severe_u);

        NonlinearRetryControls retry_controls;
        retry_controls.max_retries = 2;
        retry_controls.relaxation_shrink = 0.5;
        retry_controls.minimum_alpha_u = controls.min_alpha_u;
        retry_controls.minimum_alpha_p = controls.min_alpha_p;
        NonlinearRetryController retry(retry_controls);

        alpha_u = std::max(
            controls.min_alpha_u, alpha_u - controls.decrease_step);
        retry.retry_without_relaxation();
        EXPECT_TRUE(retry.retries() == std::size_t{1});
        EXPECT_NEAR(retry.alpha_u(0.705), 0.705, 1e-14);
        EXPECT_NEAR(alpha_u, 0.685, 1e-14);
        EXPECT_TRUE(history.size() == std::size_t{1});

        const double retry_residual_u = 1.2;
        const bool severe_on_retry = adaptive_relaxation_severe_degradation(
            history.back().residual_u, retry_residual_u,
            history.back().alpha_u, alpha_u, controls);
        EXPECT_TRUE(!severe_on_retry);
        append_accepted(history, {11, retry_residual_u, 0.8, alpha_u, 0.300});

        EXPECT_TRUE(history.size() == std::size_t{2});
        EXPECT_TRUE(history.back().iteration == std::size_t{11});
        EXPECT_NEAR(history.back().alpha_u, 0.685, 1e-14);

        const bool severe_p_without_alpha_increase =
            adaptive_relaxation_severe_degradation(
                history.back().residual_p, 2.1,
                history.back().alpha_p, history.back().alpha_p, controls);
        EXPECT_TRUE(!severe_p_without_alpha_increase);

        const bool severe_p_with_increase =
            adaptive_relaxation_severe_degradation(
                history.back().residual_p, 2.1,
                history.back().alpha_p, 0.305, controls);
        EXPECT_TRUE(severe_p_with_increase);
    });

    return run_all();
}
