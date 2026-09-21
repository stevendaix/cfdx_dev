// M0.13-T01 — Temporal Discretization
// Scalar reference implementations for time-integration schemes.
// Production field integration is provided by cfdx/core/numerics/temporal.h.
#pragma once

#include "cfdx/core/field/field.h"
#include <cmath>
#include <stdexcept>
#include <vector>

namespace cfdx {
namespace physics {

namespace detail {

inline void validate_dt(double dt) {
    if (!(dt > 0.0) || !std::isfinite(dt)) {
        throw std::invalid_argument("time step must be finite and strictly positive");
    }
}

template<class Function>
inline double fixed_point(
    double initial,
    Function&& f,
    double dt,
    double rhs_scale,
    double constant,
    int max_iterations = 100,
    double tolerance = 1e-12)
{
    double x = initial;
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        const double x_new = constant + rhs_scale * f(x);
        if (!std::isfinite(x_new)) {
            throw std::runtime_error("temporal fixed-point iteration produced a non-finite value");
        }
        if (std::abs(x_new - x) <= tolerance * std::max(1.0, std::abs(x_new))) {
            return x_new;
        }
        x = x_new;
    }
    throw std::runtime_error("temporal fixed-point iteration did not converge");
}

} // namespace detail

template<class Function>
inline double explicit_euler_step(double phi_n, double dt, Function&& f) {
    detail::validate_dt(dt);
    return phi_n + dt * f(phi_n);
}

// phi^{n+1} = phi^n + dt f(phi^{n+1})
template<class Function>
inline double implicit_euler_step(double phi_n, double dt, Function&& f) {
    detail::validate_dt(dt);
    return detail::fixed_point(phi_n, std::forward<Function>(f), dt, dt, phi_n);
}

// phi^{n+1} = phi^n + dt/2 [f(phi^n) + f(phi^{n+1})]
template<class Function>
inline double crank_nicolson_step(double phi_n, double dt, Function&& f) {
    detail::validate_dt(dt);
    const double rhs_n = f(phi_n);
    return detail::fixed_point(
        phi_n,
        std::forward<Function>(f),
        dt,
        0.5 * dt,
        phi_n + 0.5 * dt * rhs_n);
}

// (3 phi^{n+1} - 4 phi^n + phi^{n-1})/(2 dt) = f(phi^{n+1})
template<class Function>
inline double bdf2_step(
    double phi_n,
    double dt,
    double phi_prev,
    Function&& f)
{
    detail::validate_dt(dt);
    const double constant = (4.0 * phi_n - phi_prev) / 3.0;
    const double rhs_scale = 2.0 * dt / 3.0;
    return detail::fixed_point(
        phi_n,
        std::forward<Function>(f),
        dt,
        rhs_scale,
        constant);
}

template<class Function>
inline double derivative(Function&& f, double phi, double h) {
    if (!(h > 0.0) || !std::isfinite(h)) {
        throw std::invalid_argument("finite-difference step must be finite and strictly positive");
    }
    return (f(phi + h) - f(phi - h)) / (2.0 * h);
}

} // namespace physics
} // namespace cfdx
