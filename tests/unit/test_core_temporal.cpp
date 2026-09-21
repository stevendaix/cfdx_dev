#include "cfdx/core/numerics/temporal.h"
#include "common/test_harness.h"
#include <cmath>
#include <stdexcept>

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
        EXPECT_TRUE(std::isfinite(second(0)));
        EXPECT_TRUE(second(0) > 0.0);
    });

    return cfdx::testing::run_all();
}
