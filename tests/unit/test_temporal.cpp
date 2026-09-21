// M0.13-T02 — Temporal Discretization Tests (Physics module)
// Tests for explicit_euler_step, implicit_euler_step, crank_nicolson_step, bdf2_step

#include "cfdx/physics/temporal.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::physics;
using namespace cfdx::testing;

const double lambda = 2.0;
const double phi0 = 1.0;

double reaction(double phi) { return -lambda * phi; }

int main() {
    run_case("explicit_euler_decay", [&]() {
        double phi_n = 1.0;
        double dt = 0.1;
        double result = explicit_euler_step(phi_n, dt, reaction);
        double expected = phi_n + dt * reaction(phi_n);
        EXPECT_NEAR(result, expected, 1e-12);
    });

    run_case("implicit_euler_decay", [&]() {
        double phi_n = 1.0;
        double dt = 0.1;
        double result = implicit_euler_step(phi_n, dt, reaction);
        EXPECT_NEAR(result, 1.0 / 1.2, 1e-10);
    });

    run_case("crank_nicolson_decay", [&]() {
        double phi_n = 1.0;
        double dt = 0.1;
        double result = crank_nicolson_step(phi_n, dt, reaction);
        double expected = (1.0 - 0.1) / (1.0 + 0.1);
        EXPECT_NEAR(result, expected, 1e-10);
    });

    run_case("bdf2_decay", [&]() {
        double phi_n = 1.0;
        double phi_prev = 1.0;
        double dt = 0.1;
        double result = bdf2_step(phi_n, dt, phi_prev, reaction);
        double expected = 3.0 / 3.4;
        EXPECT_NEAR(result, expected, 1e-10);
    });

    run_case("derivative_finite_diff", [&]() {
        double h = 1e-6;
        double phi = 1.0;
        double deriv = derivative(reaction, phi, h);
        double expected = -lambda;
        EXPECT_NEAR(deriv, expected, 1e-3);
    });

    return run_all();
}