#pragma once
#include "cfdx/core/field/field.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::physics {

inline void validate_rk_rhs(
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& u,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k) {
    if(u.size()!=k.size() || u.dimension()!=k.dimension())
        throw std::invalid_argument("RK RHS shape mismatch");
    for(std::size_t d=0;d<u.dimension();++d)
        for(std::size_t i=0;i<u.size();++i)
            if(!std::isfinite(u(i,d)) || !std::isfinite(k(i,d)))
                throw std::runtime_error("RK non-finite input or RHS");
}

template<class RHS>
inline void ssprk3_step(
    cfdx::core::Field<double,cfdx::core::Location::CELL>& u,
    double dt, RHS&& rhs) {
    if (!(dt > 0.0) || !std::isfinite(dt))
        throw std::invalid_argument("invalid dt");
    const std::size_t n = u.size();
    cfdx::core::Field<double,cfdx::core::Location::CELL> k(n, "rhs", "", u.dimension());
    cfdx::core::Field<double,cfdx::core::Location::CELL> u0 = u;

    rhs(u, k);
    validate_rk_rhs(u,k);
    for (std::size_t d = 0; d < u.dimension(); ++d)
        for (std::size_t i = 0; i < n; ++i) u(i,d) = u0(i,d) + dt * k(i,d);

    rhs(u, k);
    validate_rk_rhs(u,k);
    for (std::size_t d = 0; d < u.dimension(); ++d)
        for (std::size_t i = 0; i < n; ++i)
        u(i,d) = 0.75 * u0(i,d) + 0.25 * (u(i,d) + dt * k(i,d));

    rhs(u, k);
    validate_rk_rhs(u,k);
    for (std::size_t d = 0; d < u.dimension(); ++d)
        for (std::size_t i = 0; i < n; ++i)
        u(i,d) = (1.0 / 3.0) * u0(i,d) + (2.0 / 3.0) * (u(i,d) + dt * k(i,d));
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
    cfdx::core::Field<double,cfdx::core::Location::CELL> r(n, "rk3_register", "", u.dimension());
    r.fill(0.0);

    constexpr double alpha[3] = {0.0, -5.0 / 9.0, -153.0 / 128.0};
    constexpr double beta[3]  = {1.0 / 3.0, 15.0 / 16.0, 8.0 / 15.0};

    cfdx::core::Field<double,cfdx::core::Location::CELL> k(n, "rhs", "", u.dimension());
    for (int s = 0; s < 3; ++s) {
        rhs(u, k);
    validate_rk_rhs(u,k);
        for (std::size_t d = 0; d < u.dimension(); ++d)
        for (std::size_t i = 0; i < n; ++i) {
            r(i,d) = alpha[s] * r(i,d) + dt * k(i,d);
            u(i,d) += beta[s] * r(i,d);
        }
    }
}

// Explicit RK2 (midpoint rule). Two stages with a half-step predictor:
//   k1 = f(u0);  u_star = u0 + dt/2 k1;  k2 = f(u_star);  u1 = u0 + dt k2
// The b-vector is (0, 1); the ORDER-CONDITION b.c = 1/2 holds, so the method
// is genuinely second order. (A previous version applied dt/2 to BOTH stages,
// which violates b.c = 1/2 and is only first order — fixed, see #461 N5.)
template<class RHS>
inline void low_storage_rk2_step(
    cfdx::core::Field<double,cfdx::core::Location::CELL>& u,
    double dt, RHS&& rhs) {
    if (!(dt > 0.0) || !std::isfinite(dt))
        throw std::invalid_argument("invalid dt");
    const std::size_t n = u.size();
    cfdx::core::Field<double,cfdx::core::Location::CELL> u0 = u;
    cfdx::core::Field<double,cfdx::core::Location::CELL> k(n, "rhs", "", u.dimension());
    rhs(u, k);
    validate_rk_rhs(u,k);
    for (std::size_t d = 0; d < u.dimension(); ++d)
        for (std::size_t i = 0; i < n; ++i) u(i,d) = u0(i,d) + 0.5 * dt * k(i,d);
    rhs(u, k);
    validate_rk_rhs(u,k);
    for (std::size_t d = 0; d < u.dimension(); ++d)
        for (std::size_t i = 0; i < n; ++i) u(i,d) = u0(i,d) + dt * k(i,d);
}

} // namespace cfdx::physics
