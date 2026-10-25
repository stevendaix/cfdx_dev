#include "cfdx/core/numerics/temporal.h"
#include "common/test_harness.h"
#include <cmath>
#include <stdexcept>
#include <limits>

using namespace cfdx::core;

int main() {
    using cfdx::testing::run_case;

    run_case("implicit_euler", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        phi(0) = 1.0;
        RhsFunction rhs = [](const ScalarCellField& in, ScalarCellField& out) {
            out(0) = -2.0 * in(0);
        };
        const auto result = advance_time(phi, 0.1, rhs, TimeScheme::EULER_IMPLICIT);
        EXPECT_NEAR(result(0), 1.0 / 1.2, 1e-10);
    });

    run_case("crank_nicolson", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        phi(0) = 1.0;
        RhsFunction rhs = [](const ScalarCellField& in, ScalarCellField& out) {
            out(0) = -2.0 * in(0);
        };
        const auto result = advance_time(phi, 0.1, rhs, TimeScheme::CRANK_NICOLSON);
        EXPECT_NEAR(result(0), 0.9 / 1.1, 1e-10);
    });

    run_case("context_shape_mismatch", [] {
        ScalarCellField phi(2, "phi", "1", 1);
        TimeIntegrationContext ctx(1, 1, "phi");
        ctx.initialize(ScalarCellField(1, "initial", "1", 1));
        RhsFunction rhs = [](const ScalarCellField&, ScalarCellField& out) { out(0) = 0.0; };
        EXPECT_THROW(advance_time(phi, 0.1, rhs, TimeScheme::BDF2, &ctx), std::invalid_argument);
    });

    run_case("nonfinite_rhs_rejected", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        phi(0) = 1.0;
        RhsFunction rhs = [](const ScalarCellField&, ScalarCellField& out) {
            out(0) = std::numeric_limits<double>::quiet_NaN();
        };
        EXPECT_THROW(advance_time(phi, 0.1, rhs, TimeScheme::EULER_EXPLICIT), std::runtime_error);
    });

    run_case("invalid_dt", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        RhsFunction rhs = [](const ScalarCellField&, ScalarCellField& out) {
            out(0) = 0.0;
        };
        EXPECT_THROW(
            advance_time(phi, 0.0, rhs, TimeScheme::EULER_IMPLICIT),
            std::invalid_argument);
    });

    run_case("bdf2_bootstrap_and_second_step", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        phi(0) = 1.0;
        TimeIntegrationContext ctx(1, 1, "phi");
        ctx.initialize(phi);
        RhsFunction rhs = [](const ScalarCellField& in, ScalarCellField& out) {
            out(0) = -2.0 * in(0);
        };

        const auto first = advance_time(phi, 0.1, rhs, TimeScheme::BDF2, &ctx);
        EXPECT_NEAR(first(0), 1.0 / 1.2, 1e-10);

        const auto second = advance_time(first, 0.1, rhs, TimeScheme::BDF2, &ctx);
        const double expected = (4.0 * first(0) - 1.0) / 3.4;
        EXPECT_NEAR(second(0), expected, 1e-10);
        EXPECT_TRUE(std::isfinite(second(0)));
        EXPECT_TRUE(second(0) > 0.0);
    });

    for (const auto scheme : {TimeScheme::EULER_IMPLICIT,
                              TimeScheme::CRANK_NICOLSON, TimeScheme::BDF2}) {
        run_case(std::string("intermediate_nonfinite_rhs_rejected_") + to_string(scheme), [scheme] {
            ScalarCellField phi(1, "phi", "1", 1);
            phi(0) = 1.0;
            TimeIntegrationContext ctx(1, 1, "phi");
            ctx.initialize(phi);
            const RhsFunction rhs = [](const ScalarCellField& in, ScalarCellField& out) {
                out(0) = -std::sqrt(in(0));
            };
            const double dt = scheme == TimeScheme::CRANK_NICOLSON ? 4.0 : 2.0;
            EXPECT_THROW(advance_time(phi, dt, rhs, scheme, &ctx), std::runtime_error);
            EXPECT_NEAR(ctx.phi_curr(0), 1.0, 0.0);
            EXPECT_FALSE(ctx.has_prev);
        });
    }

    run_case("nonfinite_initial_field_rejected_with_finite_rhs", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        phi(0) = std::numeric_limits<double>::quiet_NaN();
        const RhsFunction rhs = [](const ScalarCellField&, ScalarCellField& out) { out(0) = 0.0; };
        EXPECT_THROW(advance_time(phi, 0.1, rhs, TimeScheme::EULER_EXPLICIT), std::runtime_error);
    });

    run_case("explicit_overflow_rejected_before_history_commit", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        phi(0) = std::numeric_limits<double>::max();
        TimeIntegrationContext ctx(1, 1, "phi");
        ctx.initialize(phi);
        const RhsFunction rhs = [](const ScalarCellField&, ScalarCellField& out) {
            out(0) = std::numeric_limits<double>::max();
        };
        EXPECT_THROW(advance_time(phi, 1.0, rhs, TimeScheme::EULER_EXPLICIT, &ctx), std::runtime_error);
        EXPECT_NEAR(ctx.phi_curr(0), phi(0), 0.0);
        EXPECT_FALSE(ctx.has_prev);
    });

    run_case("bdf2_history_path_rejects_nonfinite_rhs_without_committing", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        phi(0) = 1.0;
        TimeIntegrationContext ctx(1, 1, "phi");
        ctx.initialize(phi);
        ctx.shift(phi);
        ctx.dt_prev = 2.0;
        const RhsFunction rhs = [](const ScalarCellField& in, ScalarCellField& out) {
            out(0) = -std::sqrt(in(0));
        };
        EXPECT_THROW(advance_time(phi, 2.0, rhs, TimeScheme::BDF2, &ctx), std::runtime_error);
        EXPECT_NEAR(ctx.phi_curr(0), 1.0, 0.0);
        EXPECT_NEAR(ctx.phi_prev(0), 1.0, 0.0);
        EXPECT_NEAR(ctx.dt_prev, 2.0, 0.0);
        EXPECT_TRUE(ctx.has_prev);
    });

    return cfdx::testing::run_all();
}
