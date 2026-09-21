// M0.13-T01 — Temporal Discretization
// Explicit Euler, Implicit Euler, Crank-Nicolson, BDF2
#pragma once

#include "cfdx/core/field/field.h"
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace cfdx {
namespace physics {

// Compute derivative of f at phi for Newton-Raphson.
// A centered finite difference gives second-order truncation error.
template<class Function>
inline double derivative(Function&& f, double phi, double h) {
    if (!(h > 0.0) || !std::isfinite(h)) {
        throw std::invalid_argument("derivative: step must be positive and finite");
    }
    return (f(phi + h) - f(phi - h)) / (2.0 * h);
}

// ============================================================================
// Explicit Euler
// phi^{n+1} = phi^n + dt * f(phi^n)
// ============================================================================

template<class Function>
inline double explicit_euler_step(double phi_n, double dt, Function&& f) {
    return phi_n + dt * f(phi_n);
}

// ============================================================================
// Implicit Euler
// phi^{n+1} = phi^n + dt * f(phi^{n+1})
// ============================================================================

template<class Function>
inline double implicit_euler_step(double phi_n, double dt, Function&& f,
                                  std::size_t max_iter = 50,
                                  double tolerance = 1e-12) {
    if (dt < 0.0 || max_iter == 0 || tolerance <= 0.0) {
        throw std::invalid_argument("implicit_euler_step: invalid parameters");
    }

    double x = phi_n + dt * f(phi_n);
    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        const double g = x - phi_n - dt * f(x);
        if (std::abs(g) <= tolerance * std::max(1.0, std::abs(phi_n))) return x;

        const double h = std::sqrt(std::numeric_limits<double>::epsilon()) *
                         std::max(1.0, std::abs(x));
        const double dg = 1.0 - dt * derivative(f, x, h);
        if (!std::isfinite(dg) || std::abs(dg) <= 1e-14) {
            throw std::runtime_error(
                "implicit_euler_step: singular Newton derivative");
        }

        const double dx = g / dg;
        x -= dx;
        if (std::abs(dx) <= tolerance * std::max(1.0, std::abs(x))) return x;
    }

    throw std::runtime_error(
        "implicit_euler_step: Newton iteration did not converge");
}

// ============================================================================
// Crank-Nicolson
// (phi^{n+1} - phi^n)/dt = 0.5*(f(phi^{n+1}) + f(phi^n))
// ============================================================================

template<class Function>
inline double crank_nicolson_step(double phi_n, double dt, Function&& f,
                                  std::size_t max_iter = 50,
                                  double tolerance = 1e-12) {
    if (dt < 0.0 || max_iter == 0 || tolerance <= 0.0) {
        throw std::invalid_argument("crank_nicolson_step: invalid parameters");
    }

    const double fn = f(phi_n);
    double x = phi_n + dt * fn;
    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        const double g = x - phi_n - 0.5 * dt * (f(x) + fn);
        if (std::abs(g) <= tolerance * std::max(1.0, std::abs(phi_n))) return x;

        const double h = std::sqrt(std::numeric_limits<double>::epsilon()) *
                         std::max(1.0, std::abs(x));
        const double dg = 1.0 - 0.5 * dt * derivative(f, x, h);
        if (!std::isfinite(dg) || std::abs(dg) <= 1e-14) {
            throw std::runtime_error(
                "crank_nicolson_step: singular Newton derivative");
        }

        const double dx = g / dg;
        x -= dx;
        if (std::abs(dx) <= tolerance * std::max(1.0, std::abs(x))) return x;
    }

    throw std::runtime_error(
        "crank_nicolson_step: Newton iteration did not converge");
}

// ============================================================================
// BDF2
// (3*phi^{n+1} - 4*phi^n + phi^{n-1})/(2*dt) = f(phi^{n+1})
// ============================================================================

template<class Function>
inline double bdf2_step(double phi_n, double dt, double phi_prev, Function&& f,
                        std::size_t max_iter = 50,
                        double tolerance = 1e-12) {
    if (dt <= 0.0 || max_iter == 0 || tolerance <= 0.0) {
        throw std::invalid_argument("bdf2_step: invalid parameters");
    }

    double x = phi_n + dt * f(phi_n);
    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        const double g =
            3.0 * x - 4.0 * phi_n + phi_prev - 2.0 * dt * f(x);
        if (std::abs(g) <= tolerance * std::max(1.0, std::abs(phi_n))) return x;

        const double h = std::sqrt(std::numeric_limits<double>::epsilon()) *
                         std::max(1.0, std::abs(x));
        const double dg = 3.0 - 2.0 * dt * derivative(f, x, h);
        if (!std::isfinite(dg) || std::abs(dg) <= 1e-14) {
            throw std::runtime_error(
                "bdf2_step: singular Newton derivative");
        }

        const double dx = g / dg;
        x -= dx;
        if (std::abs(dx) <= tolerance * std::max(1.0, std::abs(x))) return x;
    }

    throw std::runtime_error(
        "bdf2_step: Newton iteration did not converge");
}

} // namespace physics
} // namespace cfdx
