#pragma once

#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/vector.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {

// Opt-in, diagnostic-only verification of an ASSEMBLED saddle-point Jacobian.
// Caller must supply the actual discrete residual, including its boundary
// equations, Rhie-Chow terms and gauge row, at the same linearization state.
// The checker deliberately does not infer that callback from a matrix: doing
// so would create a circular oracle that could not detect assembly defects.
struct PressureVelocityJacobianProbeControls {
    double relative_tolerance = 1e-6;
    double absolute_tolerance = 1e-8;
    // 0 checks every column. A positive bound samples a deterministic set
    // containing first and last columns; it is not exhaustive qualification.
    std::size_t max_columns = 0;
};

struct PressureVelocityJacobianCertificate {
    bool passed = false;
    std::size_t columns_checked = 0;
    double max_absolute_error = 0.0;
    double max_tolerance_ratio = 0.0;
    std::size_t worst_row = 0;
    std::size_t worst_column = 0;
};

using PressureVelocityDiscreteResidual =
    std::function<bool(const cfdx::core::Vector&, cfdx::core::Vector&)>;

inline PressureVelocityJacobianCertificate certify_pressure_velocity_jacobian(
    const cfdx::core::BlockOperator& assembled,
    const cfdx::core::Vector& linearization_state,
    const PressureVelocityDiscreteResidual& evaluate_discrete_residual,
    const PressureVelocityJacobianProbeControls& settings = {})
{
    using cfdx::core::Vector;
    assembled.validate();
    const std::size_t nu = assembled.velocity_size();
    const std::size_t np = assembled.pressure_size();
    const std::size_t n = nu + np;
    if (n == 0 || linearization_state.size() != n ||
        !linearization_state.is_valid() || !evaluate_discrete_residual ||
        !std::isfinite(settings.absolute_tolerance) ||
        !std::isfinite(settings.relative_tolerance) ||
        settings.absolute_tolerance < 0.0 ||
        settings.relative_tolerance < 0.0 ||
        (settings.absolute_tolerance == 0.0 &&
         settings.relative_tolerance == 0.0))
        throw std::invalid_argument("pressure-velocity Jacobian certificate: invalid inputs");

    PressureVelocityJacobianCertificate result;
    Vector baseline(n, 0.0);
    if (!evaluate_discrete_residual(linearization_state, baseline) ||
        baseline.size() != n || !baseline.is_valid())
        return result;

    const std::size_t count = settings.max_columns == 0
        ? n : std::min(n, settings.max_columns);
    const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
    for (std::size_t sample = 0; sample < count; ++sample) {
        const std::size_t column = count == 1 ? 0 :
            sample * (n - 1) / (count - 1);
        const double step = eps * (1.0 + std::abs(linearization_state(column)));
        Vector plus = linearization_state;
        Vector minus = linearization_state;
        plus(column) += step;
        minus(column) -= step;
        Vector plus_residual(n, 0.0);
        Vector minus_residual(n, 0.0);
        if (!evaluate_discrete_residual(plus, plus_residual) ||
            !evaluate_discrete_residual(minus, minus_residual) ||
            plus_residual.size() != n || minus_residual.size() != n ||
            !plus_residual.is_valid() || !minus_residual.is_valid())
            return result; // Fail-closed, no partial certificate is promoted.

        for (std::size_t row = 0; row < n; ++row) {
            const double analytic = row < nu
                ? (column < nu
                   ? assembled.Auu()(row, column)
                   : assembled.G()(row, column - nu))
                : (column < nu
                   ? assembled.D()(row - nu, column)
                   : assembled.C()(row - nu, column - nu));
            const double numeric =
                (plus_residual(row) - minus_residual(row)) / (2.0 * step);
            const double error = std::abs(analytic - numeric);
            const double scale = settings.absolute_tolerance +
                settings.relative_tolerance *
                    std::max(std::abs(analytic), std::abs(numeric));
            const double ratio = scale > 0.0 ? error / scale :
                (error == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
            if (!std::isfinite(ratio))
                return result;
            result.max_absolute_error = std::max(result.max_absolute_error, error);
            if (ratio > result.max_tolerance_ratio) {
                result.max_tolerance_ratio = ratio;
                result.worst_row = row;
                result.worst_column = column;
            }
        }
        ++result.columns_checked;
    }
    result.passed = result.columns_checked == count &&
        result.max_tolerance_ratio <= 1.0;
    return result;
}


// Matrix-free directional test inspired by PETSc's SNES Jacobian checking.
// Cost: 2 * directions evaluations of the *independent* FV residual rather
// than 2 * unknowns for a full column scan. This stochastic-looking but
// deterministic probing is a SCREENING test, not exhaustive certification.
// It finds wrong block signs and omitted Rhie-Chow pressure derivatives
// without constructing a dense finite-difference Jacobian.
struct PressureVelocityDirectionalControls {
    std::size_t directions = 8;
    double relative_tolerance = 1e-6;
    double absolute_tolerance = 1e-8;
    std::uint64_t seed = 0x9e3779b97f4a7c15ULL;
};

struct PressureVelocityDirectionalCertificate {
    bool passed = false;
    std::size_t directions_checked = 0;
    std::size_t residual_evaluations = 0;
    double max_absolute_error = 0.0;
    double max_tolerance_ratio = 0.0;
    std::size_t worst_row = 0;
    std::size_t worst_direction = 0;
};

inline PressureVelocityDirectionalCertificate certify_pressure_velocity_jacobian_directions(
    const cfdx::core::BlockOperator& assembled,
    const cfdx::core::Vector& state,
    const PressureVelocityDiscreteResidual& residual,
    const PressureVelocityDirectionalControls& controls = {})
{
    using cfdx::core::Vector;
    assembled.validate();
    const std::size_t nu = assembled.velocity_size();
    const std::size_t np = assembled.pressure_size();
    const std::size_t n = nu + np;
    if (n == 0 || nu == 0 || np == 0 || state.size() != n ||
        !state.is_valid() || !residual || controls.directions == 0 ||
        !std::isfinite(controls.relative_tolerance) ||
        !std::isfinite(controls.absolute_tolerance) ||
        controls.relative_tolerance < 0.0 || controls.absolute_tolerance < 0.0 ||
        (controls.relative_tolerance == 0.0 && controls.absolute_tolerance == 0.0))
        throw std::invalid_argument("directional pressure-velocity Jacobian: invalid inputs");

    PressureVelocityDirectionalCertificate result;
    // Scale h to the state norm; direction is normalized independently.
    double state_norm = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        state_norm = std::hypot(state_norm, state(i));
    const double step = std::cbrt(std::numeric_limits<double>::epsilon()) *
        std::max(1.0, state_norm);
    if (!std::isfinite(step) || !(step > 0.0)) return result;

    for (std::size_t probe = 0; probe < controls.directions; ++probe) {
        Vector direction(n, 0.0);
        // Splitmix64: immutable seed and indexes; no global RNG state and
        // identical probes in CI, DebugSanitizers and field captures.
        for (std::size_t i = 0; i < n; ++i) {
            std::uint64_t v = controls.seed ^
                (static_cast<std::uint64_t>(probe + 1) * 0xbf58476d1ce4e5b9ULL) ^
                (static_cast<std::uint64_t>(i + 1) * 0x94d049bb133111ebULL);
            v += 0x9e3779b97f4a7c15ULL;
            v = (v ^ (v >> 30)) * 0xbf58476d1ce4e5b9ULL;
            v = (v ^ (v >> 27)) * 0x94d049bb133111ebULL;
            v ^= v >> 31;
            direction(i) = (v & 1ULL) ? 1.0 : -1.0;
        }
        const double inv_norm = 1.0 / std::sqrt(static_cast<double>(n));
        for (std::size_t i = 0; i < n; ++i) direction(i) *= inv_norm;

        Vector plus = state, minus = state;
        for (std::size_t i = 0; i < n; ++i) {
            plus(i) += step * direction(i);
            minus(i) -= step * direction(i);
            if (plus(i) == minus(i)) return result; // perturbation lost
        }
        Vector r_plus(n, 0.0), r_minus(n, 0.0);
        if (!residual(plus, r_plus) || !residual(minus, r_minus) ||
            r_plus.size() != n || r_minus.size() != n ||
            !r_plus.is_valid() || !r_minus.is_valid())
            return result; // partial evidence never promoted
        result.residual_evaluations += 2;

        Vector du(nu, 0.0), dp(np, 0.0);
        for (std::size_t i = 0; i < nu; ++i) du(i) = direction(i);
        for (std::size_t i = 0; i < np; ++i) dp(i) = direction(nu + i);
        Vector Mu(nu, 0.0), Gp(nu, 0.0), Du(np, 0.0), Cp(np, 0.0);
        assembled.Auu().matvec_into(du, Mu);
        assembled.G().matvec_into(dp, Gp);
        assembled.D().matvec_into(du, Du);
        assembled.C().matvec_into(dp, Cp);
        if (!Mu.is_valid() || !Gp.is_valid() ||
            !Du.is_valid() || !Cp.is_valid()) return result;
        for (std::size_t row = 0; row < n; ++row) {
            const double analytic = row < nu
                ? Mu(row) + Gp(row) : Du(row - nu) + Cp(row - nu);
            const double numeric = (r_plus(row) - r_minus(row)) / (2.0 * step);
            const double error = std::abs(analytic - numeric);
            const double tolerance = controls.absolute_tolerance +
                controls.relative_tolerance *
                    std::max(std::abs(analytic), std::abs(numeric));
            const double ratio = tolerance > 0.0 ? error / tolerance :
                (error == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
            if (!std::isfinite(ratio)) return result;
            result.max_absolute_error = std::max(result.max_absolute_error, error);
            if (ratio > result.max_tolerance_ratio) {
                result.max_tolerance_ratio = ratio;
                result.worst_row = row;
                result.worst_direction = probe;
            }
        }
        ++result.directions_checked;
    }
    result.passed = result.directions_checked == controls.directions &&
        result.max_tolerance_ratio <= 1.0;
    return result;
}

} // namespace cfdx::physics
