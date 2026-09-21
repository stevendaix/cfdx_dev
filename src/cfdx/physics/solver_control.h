#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

struct ConvergenceCriteria {
    double absolute_tolerance = 1e-10;
    double relative_tolerance = 1e-6;
    double continuity_tolerance = 1e-10;
    double energy_tolerance = 1e-10;
    std::size_t max_iterations = 1000;
    bool require_turbulence = false;
    bool require_energy = false;
};

struct IterationMetrics {
    double momentum_residual = std::numeric_limits<double>::infinity();
    double pressure_residual = std::numeric_limits<double>::infinity();
    double turbulence_residual = std::numeric_limits<double>::infinity();
    double energy_residual = std::numeric_limits<double>::infinity();
    double continuity_imbalance = std::numeric_limits<double>::infinity();
    double energy_imbalance = std::numeric_limits<double>::infinity();
    std::size_t iteration = 0;
};

inline void validate_convergence_criteria(const ConvergenceCriteria& c)
{
    if (!(c.absolute_tolerance > 0.0) ||
        !(c.relative_tolerance > 0.0) ||
        !(c.continuity_tolerance > 0.0) ||
        !(c.energy_tolerance > 0.0) ||
        c.max_iterations == 0)
        throw std::invalid_argument("invalid convergence criteria");
}

inline bool residual_converged(double initial_residual, double residual,
                               const ConvergenceCriteria& c)
{
    if (!std::isfinite(initial_residual) || !std::isfinite(residual))
        return false;
    const double scale = std::max(std::abs(initial_residual), 1.0);
    return std::abs(residual) <= c.absolute_tolerance ||
           std::abs(residual) <= c.relative_tolerance * scale;
}

inline bool conservation_converged(const IterationMetrics& m,
                                   const ConvergenceCriteria& c)
{
    return std::isfinite(m.continuity_imbalance) &&
           std::isfinite(m.energy_imbalance) &&
           std::abs(m.continuity_imbalance) <= c.continuity_tolerance &&
           std::abs(m.energy_imbalance) <= c.energy_tolerance;
}

inline bool converged(const IterationMetrics& m,
                      const IterationMetrics& initial,
                      const ConvergenceCriteria& c)
{
    const bool turbulence_ok =
        !c.require_turbulence ||
        residual_converged(initial.turbulence_residual, m.turbulence_residual, c);
    const bool energy_ok =
        !c.require_energy ||
        residual_converged(initial.energy_residual, m.energy_residual, c);
    const bool conservation_ok =
        std::isfinite(m.continuity_imbalance) &&
        std::abs(m.continuity_imbalance) <= c.continuity_tolerance &&
        (!c.require_energy ||
         (std::isfinite(m.energy_imbalance) &&
          std::abs(m.energy_imbalance) <= c.energy_tolerance));
    return residual_converged(initial.momentum_residual, m.momentum_residual, c) &&
           residual_converged(initial.pressure_residual, m.pressure_residual, c) &&
           turbulence_ok && energy_ok && conservation_ok;
}

} // namespace cfdx::physics
