#pragma once

#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/vector.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
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

} // namespace cfdx::physics
