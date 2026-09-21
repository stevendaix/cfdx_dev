#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct FsiControls {
    double relaxation = 0.5;
    double aitken_min = 0.05;
    double aitken_max = 1.0;
    std::size_t max_sub_iterations = 50;
    double displacement_tolerance = 1e-8;
    double force_tolerance = 1e-8;
};

struct FsiInterfaceState {
    std::vector<double> displacement;
    std::vector<double> force;
};

inline double interface_work(
    const std::vector<double>& force,
    const std::vector<double>& displacement_increment)
{
    if (force.size() != displacement_increment.size())
        throw std::invalid_argument("FSI interface vector size mismatch");
    double work = 0.0;
    for (std::size_t i = 0; i < force.size(); ++i)
        work += force[i] * displacement_increment[i];
    return work;
}

inline double displacement_residual(
    const std::vector<double>& current,
    const std::vector<double>& previous)
{
    if (current.size() != previous.size())
        throw std::invalid_argument("FSI displacement size mismatch");
    double norm = 0.0;
    for (std::size_t i = 0; i < current.size(); ++i)
        norm = std::max(norm, std::abs(current[i] - previous[i]));
    return norm;
}

inline double aitken_relaxation(
    double previous_relaxation,
    const std::vector<double>& residual,
    const std::vector<double>& previous_residual,
    double minimum,
    double maximum)
{
    if (residual.size() != previous_residual.size() || residual.empty())
        throw std::invalid_argument("invalid Aitken residual vectors");
    if (!(minimum > 0.0) || !(maximum >= minimum))
        throw std::invalid_argument("invalid Aitken bounds");

    double denominator = 0.0;
    double numerator = 0.0;
    for (std::size_t i = 0; i < residual.size(); ++i) {
        const double delta = residual[i] - previous_residual[i];
        denominator += delta * delta;
        numerator += previous_residual[i] * delta;
    }
    if (!(denominator > 0.0) || !std::isfinite(denominator))
        return std::clamp(previous_relaxation, minimum, maximum);

    const double candidate = -previous_relaxation * numerator / denominator;
    if (!std::isfinite(candidate))
        return std::clamp(previous_relaxation, minimum, maximum);
    return std::clamp(candidate, minimum, maximum);
}

inline bool fsi_converged(
    double displacement_residual_value,
    double force_residual_value,
    const FsiControls& controls)
{
    return std::isfinite(displacement_residual_value) &&
           std::isfinite(force_residual_value) &&
           displacement_residual_value <= controls.displacement_tolerance &&
           force_residual_value <= controls.force_tolerance;
}

} // namespace cfdx::physics
