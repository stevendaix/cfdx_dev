// Issue #461 — N5: BDF2 history lifecycle, restart semantics and
// variable-step consistency.
//
//   - bdf2_coefficients unit checks (uniform ω=1 reduces to 4/3,1/3,2/3);
//   - a single variable-step BDF2 step against the analytic implicit result;
//   - history lifecycle: has_prev transitions and history tracking;
//   - variable-step order: a fixed dt *pattern*, scaled uniformly, converges
//     at ~2nd order;
//   - restart semantics: restarting without history (Euler bootstrap) stays
//     stable and correct, and differs from the continuous run by O(dt).

#include "cfdx/core/numerics/temporal.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

const double lambda = 2.0;
double reaction(double phi) { return -lambda * phi; }

RhsFunction ode_rhs = [](const Field<double, Location::CELL>& x,
                         Field<double, Location::CELL>& rhs) {
    rhs(0) = -lambda * x(0);
};

Field<double, Location::CELL> field1(double u) {
    Field<double, Location::CELL> f(1, "u", "1", 1);
    f(0) = u;
    return f;
}

double exact(double t) { return std::exp(-lambda * t); }

} // namespace

int main() {
    run_case("bdf2_coefficients_uniform", [&]() {
        double a0, a1, denom, b;
        bdf2_coefficients(0.0, 0.1, a0, a1, denom, b);   // unset history -> ω=1
        EXPECT_NEAR(a0 / denom, 4.0 / 3.0, 1e-14);
        EXPECT_NEAR(a1 / denom, 1.0 / 3.0, 1e-14);
        EXPECT_NEAR(b / denom, 2.0 * 0.1 / 3.0, 1e-14);
        bdf2_coefficients(0.1, 0.2, a0, a1, denom, b);   // ω=2
        EXPECT_NEAR(a0, 9.0, 1e-14);
        EXPECT_NEAR(a1, 4.0, 1e-14);
        EXPECT_NEAR(denom, 5.0, 1e-14);
        EXPECT_NEAR(b, 0.6, 1e-14);
    });

    run_case("bdf2_variable_step_analytic", [&]() {
        // Two steps: dt0 then dt1 (dt1 != dt0). Seed exact history.
        const double dt0 = 0.2, dt1 = 0.1;
        TimeIntegrationContext ctx(1, 1, "u");
        ctx.phi_curr = field1(exact(dt0));
        ctx.phi_prev = field1(exact(0.0));
        ctx.has_prev = true;
        ctx.dt_prev = dt0;
        auto result = advance_time(ctx.phi_curr, dt1, ode_rhs, TimeScheme::BDF2, &ctx);
        // Analytic implicit variable-step BDF2 for u' = -λu:
        //   (a0 φn - a1 φprev) / (denom + b λ)
        double a0, a1, denom, b;
        bdf2_coefficients(dt0, dt1, a0, a1, denom, b);
        const double expected =
            (a0 * exact(dt0) - a1 * exact(0.0)) / (denom + b * lambda);
        EXPECT_NEAR(result(0), expected, 1e-12);
        EXPECT_TRUE(std::abs(result(0) - exact(dt0 + dt1)) < 0.02);
    });

    run_case("bdf2_history_lifecycle", [&]() {
        TimeIntegrationContext ctx(1, 1, "u");
        ctx.initialize(field1(1.0));                       // φ^0 = history seed
        EXPECT_TRUE(!ctx.has_prev);
        auto u1 = advance_time(field1(1.0), 0.1, ode_rhs, TimeScheme::BDF2, &ctx);
        EXPECT_TRUE(ctx.has_prev);
        EXPECT_NEAR(ctx.phi_curr(0), u1(0), 1e-15);
        EXPECT_NEAR(ctx.phi_prev(0), 1.0, 1e-15);          // the original φ^0
        auto u2 = advance_time(u1, 0.1, ode_rhs, TimeScheme::BDF2, &ctx);
        EXPECT_NEAR(ctx.phi_curr(0), u2(0), 1e-15);
        EXPECT_NEAR(ctx.phi_prev(0), u1(0), 1e-15);        // φ^1 shifted to prev
    });

    run_case("bdf2_variable_step_order", [&]() {
        // A fixed dt pattern scaled by base; integrate to t=1, expect 2nd order.
        const std::vector<double> pattern = {0.5, 1.0, 1.4, 0.7, 0.9, 1.2, 0.8, 1.1, 0.6, 1.3};
        auto integrate_pattern = [&](double base) {
            // Seed the FIRST step exactly so the start-up does not pollute the
            // order (a first-order bootstrap only would dominate the error).
            TimeIntegrationContext ctx(1, 1, "u");
            const double dt0 = base * pattern[0];
            ctx.phi_prev = field1(exact(0.0));
            ctx.phi_curr = field1(exact(dt0));
            ctx.has_prev = true;
            ctx.dt_prev = dt0;
            auto u = ctx.phi_curr;
            double t = dt0;
            int k = 1;
            while (t < 1.0 - 1e-12) {
                const double dt = std::min(base * pattern[k % pattern.size()], 1.0 - t);
                u = advance_time(u, dt, ode_rhs, TimeScheme::BDF2, &ctx);
                t += dt;
                ++k;
            }
            return u(0);
        };
        double e0 = std::abs(integrate_pattern(0.05) - exact(1.0));
        double e1 = std::abs(integrate_pattern(0.025) - exact(1.0));
        double e2 = std::abs(integrate_pattern(0.0125) - exact(1.0));
        const double p1 = std::log(e0 / e1) / std::log(2.0);
        const double p2 = std::log(e1 / e2) / std::log(2.0);
        EXPECT_TRUE(p1 > 1.5 && p2 > 1.5);
        std::cout << "BDF2_VAR_ORDER p=" << p1 << " " << p2 << "\n";
    });

    run_case("bdf2_restart_semantics", [&]() {
        // Continuous run vs restart-at-midpoint (fresh context -> Euler boot).
        auto run_continuous = []() {
            TimeIntegrationContext ctx(1, 1, "u");
            ctx.initialize(field1(1.0));
            auto u = field1(1.0);
            double t = 0.0;
            const double dt = 0.01;
            while (t < 1.0 - 1e-12) { u = advance_time(u, dt, ode_rhs, TimeScheme::BDF2, &ctx); t += dt; }
            return u(0);
        };
        auto run_restarted = []() {
            TimeIntegrationContext ctx(1, 1, "u");
            ctx.initialize(field1(1.0));
            auto u = field1(1.0);
            double t = 0.0;
            const double dt = 0.01;
            while (t < 0.5 - 1e-12) { u = advance_time(u, dt, ode_rhs, TimeScheme::BDF2, &ctx); t += dt; }
            // Restart WITHOUT stored history: seed the restart context with the
            // current field; the first step after restart degrades to an
            // implicit-Euler bootstrap (first order), then BDF2 resumes.
            TimeIntegrationContext ctx_restart(1, 1, "u");
            ctx_restart.initialize(u);
            while (t < 1.0 - 1e-12) { u = advance_time(u, dt, ode_rhs, TimeScheme::BDF2, &ctx_restart); t += dt; }
            return u(0);
        };
        const double c = run_continuous();
        const double r = run_restarted();
        EXPECT_TRUE(std::isfinite(r) && r > 0.0);
        EXPECT_NEAR(r, exact(1.0), 0.01);                       // correct to ~1% at dt=0.01
        const double diff = std::abs(r - c);
        EXPECT_TRUE(diff < 0.02);                               // restart gap is bounded (= one Euler step at restart)
        std::cout << "BDF2_RESTART continuous=" << c << " restarted=" << r
                  << " diff=" << diff << "\n";
    });

    return run_all();
}