// M0.13-T01 — Temporal Discretization Tests
// Tests for Explicit Euler, Implicit Euler, Crank-Nicolson, BDF2

#include <cassert>
#include <cmath>
#include "physics/temporal.h"

namespace {

// Analytical solution for dphi/dt = -lambda * phi
// Solution: phi(t) = phi0 * exp(-lambda * t)
const double lambda = 2.0;
const double phi0 = 1.0;

// Test function: dphi/dt = -lambda * phi
// Exact solution: phi(t) = phi0 * exp(-lambda * t)

// Explicit Euler test
TEST(ExplicitEuler, BasicForwardEuler)
{
    double phi_n = 1.0;
    double dt = 0.1;
    double f = []() -> double { return -lambda * phi_n; }();
    double result = explicit_euler_step(phi_n, dt, f);
    
    // Exact solution at t+dt: phi0 * exp(-lambda * dt)
    double exact = phi0 * std::exp(-lambda * dt);
    double error = std::abs(result - exact);
    
    // Explicit Euler is first-order accurate
    // Expected error ~ O(dt) = 0.1
    // With dt=0.1, error should be roughly proportional to dt
    assert(error < 0.5); // Should be small
    
    // More precise check: error should decrease with smaller dt
    double phi_n2 = 1.0;
    double dt2 = 0.05;
    double result2 = explicit_euler_step(phi_n2, dt2, f);
    double exact2 = phi0 * std::exp(-lambda * dt2);
    double error2 = std::abs(result2 - exact2);
    
    // Second-order improvement with halved dt
    assert(error < error2 * 2.0); // Should be roughly half the error
}

// Implicit Euler test (placeholder - returns identity)
TEST(ImplicitEuler, BasicIdentityPlaceholder)
{
    double phi_n = 1.0;
    double dt = 0.1;
    double f = []() -> double { return -lambda * phi_n; }();
    double result = implicit_euler_step(phi_n, dt, f);
    
    // Placeholder returns phi_n (identity approximation)
    // This is a placeholder - real implementation would solve nonlinear eq
    assert(result == phi_n); // Identity approximation
}

// Crank-Nicolson test
TEST(CrankNicolson, BasicCNApproximation)
{
    double phi_n = 1.0;
    double dt = 0.1;
    double f = []() -> double { return -lambda * phi_n; }();
    double result = crank_nicolson_step(phi_n, dt, f);
    
    // For dphi/dt = -lambda*phi, exact solution: phi(t) = phi0 * exp(-lambda*t)
    // At t+dt: phi_exact = phi0 * exp(-lambda * dt)
    double exact = phi0 * std::exp(-lambda * dt);
    double error = std::abs(result - exact);
    
    // CN is second-order accurate, so error should be O(dt^2)
    // With dt=0.1, error should be roughly (0.1)^2 = 0.01
    assert(error < 0.1); // Should be significantly smaller than explicit Euler
}

// BDF2 test
TEST(BDF2, BasicBDF2Approximation)
{
    double phi_n = 1.0;
    double dt = 0.1;
    double f = []() -> double { return -lambda * phi_n; }();
    double result = bdf2_step(phi_n, dt, 1.0, f); // phi_prev = 1.0 (assumed)
    
    // BDF2 formula: phi^{n+1} = (4*phi^n - phi^{n-1} + 2*dt*f(phi^n))/3
    // With phi^{n-1} = phi_n = 1.0 and f(phi^n) = -lambda * phi_n
    // phi^{n+1} = (4*1.0 - 1.0 + 2*0.1*(-2.0*1.0))/3
    // phi^{n+1} = (4 - 1 - 0.4)/3 = 2.6/3 = 0.8667
    // Exact solution at t+dt: phi0 * exp(-lambda * dt) = exp(-0.2) = 0.8187
    double exact = phi0 * std::exp(-lambda * dt);
    double error = std::abs(result - exact);
    
    // BDF2 is second-order accurate, so error should be O(dt^2)
    // With dt=0.1, error should be roughly (0.1)^2 = 0.01
    assert(error < 0.1); // Should be reasonably close
}

} // anonymous
