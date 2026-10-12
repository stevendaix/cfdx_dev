#include "cfdx/physics/timestep_control.h"
#include "common/test_harness.h"

#include <cmath>
#include <limits>

using namespace cfdx::physics;
using namespace cfdx::core;
using namespace cfdx::testing;

static Mesh thin_box()
{
    Mesh m;
    m.points().resize(8);
    m.points().set(0, 0.0, 0.0, 0.0);
    m.points().set(1, 2.0, 0.0, 0.0);
    m.points().set(2, 2.0, 0.1, 0.0);
    m.points().set(3, 0.0, 0.1, 0.0);
    m.points().set(4, 0.0, 0.0, 0.1);
    m.points().set(5, 2.0, 0.0, 0.1);
    m.points().set(6, 2.0, 0.1, 0.1);
    m.points().set(7, 0.0, 0.1, 0.1);
    m.faces().push_face({0, 3, 2, 1});
    m.faces().push_face({4, 5, 6, 7});
    m.faces().push_face({0, 1, 5, 4});
    m.faces().push_face({1, 2, 6, 5});
    m.faces().push_face({2, 3, 7, 6});
    m.faces().push_face({3, 0, 4, 7});
    m.ownership().resize(6);
    for (std::size_t f = 0; f < 6; ++f) {
        m.ownership().set_owner(f, 0);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0, 1, 2, 3, 4, 5});
    return m;
}

int main()
{
    run_case("characteristic_length_is_sensitive_to_thin_cell", [] {
        const auto h = compute_cell_characteristic_lengths(thin_box());
        EXPECT_TRUE(h.size() == 1);
        EXPECT_TRUE(h[0] > 0.0);
        EXPECT_TRUE(h[0] < 0.1);
    });

    run_case("global_and_local_cfl_are_deterministic", [] {
        Mesh m = thin_box();
        Field<double, Location::FACE> flux(6, "phi", "m3/s", 1);
        flux.fill(0.1);
        const std::vector<double> volume{0.02};
        const auto local = compute_local_convective_cfl(m, flux, volume, 0.1);
        const auto summary = summarize_cfl(local);
        EXPECT_NEAR(local[0], 3.0, 1e-14);
        EXPECT_NEAR(summary.global_cfl, 3.0, 1e-14);
        EXPECT_TRUE(summary.limiting_cell == 0);
    });

    run_case("stagnant_cell_has_zero_cfl_without_hidden_fallback", [] {
        Mesh m = thin_box();
        Field<double, Location::FACE> flux(6, "phi", "m3/s", 1);
        flux.fill(0.0);
        const auto local = compute_local_convective_cfl(m, flux, {0.02}, 0.1);
        EXPECT_NEAR(local[0], 0.0, 1e-14);
        const auto summary = summarize_cfl(local);
        EXPECT_NEAR(summary.global_cfl, 0.0, 1e-14);
    });

    run_case("invalid_cell_volume_is_rejected", [] {
        Mesh m = thin_box();
        Field<double, Location::FACE> flux(6, "phi", "m3/s", 1);
        flux.fill(0.1);
        EXPECT_THROW(
            compute_local_convective_cfl(m, flux, {0.0}, 0.1),
            std::invalid_argument);
    });

    run_case("controller_is_deterministic_and_bounded", [] {
        TimeStepControllerControls c;
        c.target_cfl = 1.0;
        c.min_dt = 0.1;
        c.max_dt = 2.0;
        c.growth_limit = 2.0;
        c.shrink_limit = 0.5;
        AuditableTimeStepController controller(c);

        const auto a = controller.propose(1.0, 4.0);
        const auto b = controller.propose(1.0, 4.0);
        EXPECT_NEAR(a.new_dt, 0.5, 1e-14);
        EXPECT_NEAR(a.new_dt, b.new_dt, 1e-14);
        EXPECT_TRUE(a.reason == TimeStepChangeReason::TargetCfl);

        const auto low = controller.propose(1.0, 1.0e-6);
        EXPECT_NEAR(low.new_dt, 2.0, 1e-14);
        const auto retry = controller.rollback(0.4, 1);
        EXPECT_NEAR(retry.new_dt, 0.2, 1e-14);
        EXPECT_TRUE(retry.reason == TimeStepChangeReason::NonlinearFailureRollback);

        EXPECT_TRUE(controller.history().size() == 4);
        EXPECT_TRUE(controller.history().front().reason == TimeStepChangeReason::TargetCfl);
        EXPECT_TRUE(controller.history().back().reason == TimeStepChangeReason::NonlinearFailureRollback);
    });

    run_case("nonlinear_retry_reduces_relaxation_deterministically", [] {
        NonlinearRetryControls c;
        c.max_retries = 2;
        c.relaxation_shrink = 0.5;
        c.minimum_alpha_u = 0.1;
        c.minimum_alpha_p = 0.05;
        NonlinearRetryController retry(c);
        EXPECT_NEAR(retry.alpha_u(0.8), 0.8, 1e-14);
        EXPECT_TRUE(retry.can_retry());
        retry.reject();
        EXPECT_NEAR(retry.alpha_u(0.8), 0.4, 1e-14);
        EXPECT_NEAR(retry.alpha_p(0.3), 0.15, 1e-14);
        retry.reject();
        EXPECT_NEAR(retry.alpha_u(0.8), 0.2, 1e-14);
        EXPECT_TRUE(!retry.can_retry());
    });

    run_case("nonlinear_retry_can_shrink_relaxation_channels_independently", [] {
        NonlinearRetryControls c;
        c.max_retries = 2;
        c.relaxation_shrink = 0.5;
        c.minimum_alpha_u = 0.1;
        c.minimum_alpha_p = 0.05;
        NonlinearRetryController retry(c);

        retry.reject(true, false);
        EXPECT_NEAR(retry.alpha_u(0.8), 0.4, 1e-14);
        EXPECT_NEAR(retry.alpha_p(0.3), 0.3, 1e-14);

        retry.reject(false, true);
        EXPECT_NEAR(retry.alpha_u(0.8), 0.4, 1e-14);
        EXPECT_NEAR(retry.alpha_p(0.3), 0.15, 1e-14);
        EXPECT_TRUE(!retry.can_retry());
    });

    run_case("nonlinear_retry_can_retry_without_generic_relaxation_shrink", [] {
        NonlinearRetryControls c;
        c.max_retries = 1;
        c.relaxation_shrink = 0.5;
        c.minimum_alpha_u = 0.1;
        c.minimum_alpha_p = 0.05;
        NonlinearRetryController retry(c);

        retry.retry_without_relaxation();
        EXPECT_NEAR(retry.alpha_u(0.8), 0.8, 1e-14);
        EXPECT_NEAR(retry.alpha_p(0.3), 0.3, 1e-14);
        EXPECT_TRUE(!retry.can_retry());
    });

    run_case("nonlinear_state_rollback_restores_velocity_and_pressure", [] {
        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(1.0);
        p.fill(2.0);
        NonlinearStateRollback transaction(U, p);
        transaction.begin();
        U(0, 0) = 9.0;
        U(1, 1) = 8.0;
        p(0) = 7.0;
        transaction.reject();
        EXPECT_NEAR(U(0, 0), 1.0, 1e-14);
        EXPECT_NEAR(U(1, 1), 1.0, 1e-14);
        EXPECT_NEAR(p(0), 2.0, 1e-14);
        EXPECT_TRUE(!transaction.active());
    });

    run_case("rollback_restores_state_and_temporal_history", [] {
        Field<double, Location::CELL> phi(2, "phi", "1", 1);
        phi(0) = 1.0;
        phi(1) = 2.0;

        TimeIntegrationContext ctx(2, 1, "phi");
        ctx.initialize(phi);
        ctx.dt_prev = 0.25;

        TimeStepRollback transaction(ctx);
        transaction.begin(phi);

        phi(0) = 7.0;
        phi(1) = 9.0;
        ctx.phi_prev(0) = -3.0;
        ctx.phi_curr(0) = 8.0;
        ctx.has_prev = true;
        ctx.dt_prev = 0.01;

        transaction.reject(phi);
        EXPECT_NEAR(phi(0), 1.0, 1e-14);
        EXPECT_NEAR(phi(1), 2.0, 1e-14);
        EXPECT_NEAR(ctx.dt_prev, 0.25, 1e-14);
        EXPECT_NEAR(ctx.phi_curr(0), 1.0, 1e-14);
        EXPECT_TRUE(!ctx.has_prev);
        EXPECT_TRUE(!transaction.active());
    });

    run_case("rollback_requires_active_transaction", [] {
        Field<double, Location::CELL> phi(1, "phi", "1", 1);
        TimeIntegrationContext ctx(1, 1, "phi");
        TimeStepRollback transaction(ctx);
        EXPECT_THROW(transaction.reject(phi), std::logic_error);
    });

    return run_all();
}
