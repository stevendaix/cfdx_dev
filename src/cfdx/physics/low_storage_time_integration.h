#pragma once
#include "cfdx/core/field/field.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::physics {

template<class RHS>
inline void ssprk3_step(
    cfdx::core::Field<double,cfdx::core::Location::CELL>& u,
    double dt, RHS&& rhs) {
    if (!(dt > 0.0) || !std::isfinite(dt))
        throw std::invalid_argument("invalid dt");
    const std::size_t n = u.size();
    cfdx::core::Field<double,cfdx::core::Location::CELL> k(n, "rhs");
    cfdx::core::Field<double,cfdx::core::Location::CELL> u0 = u;

    rhs(u, k);
    for (std::size_t i = 0; i < n; ++i) u(i) = u0(i) + dt * k(i);

    rhs(u, k);
    for (std::size_t i = 0; i < n; ++i)
        u(i) = 0.75 * u0(i) + 0.25 * (u(i) + dt * k(i));

    rhs(u, k);
    for (std::size_t i = 0; i < n; ++i)
        u(i) = (1.0 / 3.0) * u0(i) + (2.0 / 3.0) * (u(i) + dt * k(i));
}

// Williamson 3-stage third-order low-storage RK. This is genuinely 2-register:
// the solution and one residual/register field are retained.
template<class RHS>
inline void low_storage_rk3_step(
    cfdx::core::Field<double,cfdx::core::Location::CELL>& u,
    double dt, RHS&& rhs) {
    if (!(dt > 0.0) || !std::isfinite(dt))
        throw std::invalid_argument("invalid dt");
    const std::size_t n = u.size();
    cfdx::core::Field<double,cfdx::core::Location::CELL> r(n, "rk3_register");
    r.fill(0.0);

    constexpr double alpha[3] = {0.0, -5.0 / 9.0, -153.0 / 128.0};
    constexpr double beta[3]  = {1.0 / 3.0, 15.0 / 16.0, 8.0 / 15.0};

    cfdx::core::Field<double,cfdx::core::Location::CELL> k(n, "rhs");
    for (int s = 0; s < 3; ++s) {
        rhs(u, k);
        for (std::size_t i = 0; i < n; ++i) {
            r(i) = alpha[s] * r(i) + dt * k(i);
            u(i) += beta[s] * r(i);
        }
    }
}

template<class RHS>
inline void low_storage_rk2_step(
    cfdx::core::Field<double,cfdx::core::Location::CELL>& u,
    double dt, RHS&& rhs) {
    if (!(dt > 0.0) || !std::isfinite(dt))
        throw std::invalid_argument("invalid dt");
    const std::size_t n = u.size();
    cfdx::core::Field<double,cfdx::core::Location::CELL> k(n, "rhs");
    rhs(u, k);
    for (std::size_t i = 0; i < n; ++i) u(i) += 0.5 * dt * k(i);
    rhs(u, k);
    for (std::size_t i = 0; i < n; ++i) u(i) += 0.5 * dt * k(i);
}

} // namespace cfdx::physics
