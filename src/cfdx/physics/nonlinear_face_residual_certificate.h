#pragma once

#include "cfdx/core/numerics/conservation.h"
#include "cfdx/core/numerics/temporal.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {

// The matrix equation A(U_old) U - b(U_old) is a Picard operator. This
// certificate instead consumes a face-wise re-evaluation of transport on the
// accepted U,p state. In particular, it does NOT use the global signed balance
// (which can be zero when two cells have equal and opposite defects).
struct NonlinearFaceResidualCertificate {
    double local_linf = std::numeric_limits<double>::infinity();
    double local_relative = std::numeric_limits<double>::infinity();
    double local_scale = 0.0;
    std::size_t worst_cell = 0;
    bool finite = false;
};

inline NonlinearFaceResidualCertificate certify_nonlinear_face_residual(
    const cfdx::core::IntegratedBalanceDiagnostics& balance,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& flux,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& source,
    const std::vector<double>& volumes,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* accumulation = nullptr)
{
    if (source.dimension() != 1 || flux.dimension() != 1 ||
        source.size() != volumes.size() ||
        (accumulation && (accumulation->size() != volumes.size() ||
                          accumulation->dimension() != 1)))
        throw std::invalid_argument("nonlinear face residual: shape mismatch");
    NonlinearFaceResidualCertificate out;
    if (!balance.finite() || !std::isfinite(balance.max_cell_residual) ||
        balance.max_cell_residual < 0.0)
        return out;
    double scale = 0.0;
    for (std::size_t face = 0; face < flux.size(); ++face) {
        const double value = flux(face);
        if (!std::isfinite(value)) return out;
        scale = std::max(scale, std::abs(value));
    }
    for (std::size_t cell = 0; cell < source.size(); ++cell) {
        const double volume = volumes[cell];
        const double value = source(cell);
        if (!std::isfinite(value) || !std::isfinite(volume) || volume <= 0.0)
            return out;
        scale = std::max(scale, std::abs(value * volume));
        if (accumulation) {
            const double temporal = (*accumulation)(cell);
            if (!std::isfinite(temporal)) return out;
            scale = std::max(scale, std::abs(temporal));
        }
    }
    out.local_linf = balance.max_cell_residual;
    out.local_scale = scale;
    out.worst_cell = balance.worst_cell;
    // Exact zero-transport equations have zero residual, not 0/0. A
    // nonzero residual with a zero physical scale must never be certified.
    out.local_relative = scale > 0.0
        ? out.local_linf / scale
        : (out.local_linf == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
    out.finite = std::isfinite(out.local_relative);
    return out;
}

// A nonlinear iterate is not certified by a Picard matrix residual alone.
// Use the same UNCHANGED relative tolerance as the carrier monitor. Every
// momentum component must pass, including a nominally zero transverse one.
inline bool independent_momentum_face_residual_gate(
    const std::array<double, 3>& local_relative, double tolerance)
{
    if (!std::isfinite(tolerance) || tolerance <= 0.0)
        throw std::invalid_argument("independent face residual: invalid tolerance");
    return std::all_of(local_relative.begin(), local_relative.end(),
        [tolerance](double value) {
            return std::isfinite(value) && value >= 0.0 && value <= tolerance;
        });
}

// audit_integrated_balance uses div(flux) + source - accumulation = 0.
// In momentum, d(rho U)/dt is on the LEFT, so the argument to subtract is
// the *negative* discrete time derivative. The exact same BE/variable-step
// BDF2 coefficients as the production predictor are used; no pseudo-time
// momentum diagonal belongs to this physical balance.
inline double momentum_audit_negative_accumulation(
    double u_new, double u_previous, double u_previous_previous,
    double density, double volume, double dt, double dt_previous, bool bdf2)
{
    if (!(density > 0.) || !(volume > 0.) || !(dt > 0.) ||
        !std::isfinite(density) || !std::isfinite(volume) ||
        !std::isfinite(dt) ||
        !std::isfinite(u_new) || !std::isfinite(u_previous) ||
        !std::isfinite(u_previous_previous) ||
        (bdf2 && (!(dt_previous > 0.) || !std::isfinite(dt_previous))))
        throw std::invalid_argument("nonlinear face residual: invalid transient state");
    double a0 = 1.0, a1 = 0.0, denom = 1.0, b = dt;
    if (bdf2)
        cfdx::core::bdf2_coefficients(dt_previous, dt, a0, a1, denom, b);
    const double value = -density * volume *
        (denom * u_new - a0 * u_previous + a1 * u_previous_previous) / b;
    if (!std::isfinite(value))
        throw std::invalid_argument("nonlinear face residual: overflow in time derivative");
    return value;
}

} // namespace cfdx::physics
