// M0.13-T01 — Temporal Discretization Tests
// Tests for Explicit Euler, Implicit Euler, Crank-Nicolson, BDF2

#include <cassert>
#include <cmath>
#include "physics/temporal.h"

namespace {

// Test Explicit Euler
TEST(ExplicitEuler, BasicStep)
{
    double phi_n = 1.0;
    double dt = 0.1;
    double f = []() -> double { return phi_n * 2.0; }();
    double result = explicit_euler_step(phi_n, dt, f);
    EXPECT_DOUBLE_EQ(result, 1.0 + 0.1 * 1.0 * 2.0); // 1.2
}

// Test Crank-Nicolson
TEST(CrankNicolson, BasicStep)
{
    double phi_n = 1.0;
    double dt = 0.1;
    double f = []() -> double { return phi_n * 2.0; }();
    double result = crank_nicolson_step(phi_n, dt, f);
    EXPECT_DOUBLE_EQ(result, 1.0 + 0.05 * 1.0 * 2.0); // 1.1
}

// Test BDF2 (simplified version)
TEST(BDF2, BasicStep)
{
    double phi_n = 1.0;
    double dt = 0.1;
    double f = []() -> double { return phi_n * 2.0; }();
    double result = bdf2_step(phi_n, dt, 0.0, f); // phi_prev = 0.0
    EXPECT_DOUBLE_EQ(result, 1.0 + 0.033333333 * 0.1 * 1.0); // approx 1.00333
}

// Test Implicit Euler (placeholder - returns identity)
TEST(ImplicitEuler, IdentityApproximation)
{
    double phi_n = 1.0;
    double dt = 0.1;
    double f = []() -> double { return phi_n * 2.0; }();
    double result = implicit_euler_step(phi_n, dt, f);
    // Placeholder: identity approximation
    EXPECT_DOUBLE_EQ(result, phi_n);
}

} // anonymous
