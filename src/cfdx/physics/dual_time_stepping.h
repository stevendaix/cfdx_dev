#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>\n#include <functional>
#include <string>
#include <vector>

namespace cfdx::physics {

enum class DualTimePhysicalScheme { BACKWARD_EULER, BDF2 };

struct DualTimeStepControls {
    DualTimePhysicalScheme physical_scheme = DualTimePhysicalScheme::BDF2;
    std::size_t max_pseudo_iterations = 100;
    double pseudo_dt_initial = 1.0;
    double pseudo_dt_min = 1.0e-8;
    double pseudo_dt_max = 1.0e6;
    double pseudo_dt_growth = 1.5;
    double pseudo_dt_shrink = 0.5;
    double absolute_tolerance = 1.0e-10;
    double relative_tolerance = 1.0e-8;
};

struct DualTimeIteration {
    std::size_t iteration = 0;
    double residual_norm = 0.0;
    double relative_residual = 0.0;
    double pseudo_dt = 0.0;
};

struct DualTimeStepReport {
    bool converged = false;
    std::size_t pseudo_iterations = 0;
    double initial_residual = 0.0;
    double final_residual = 0.0;
    std::vector<DualTimeIteration> history;
};

class DualTimeConvergenceFailure : public std::runtime_error {
public:
    explicit DualTimeConvergenceFailure(const std::string& message)
        : std::runtime_error(message) {}
};

inline void validate_dual_time_controls(const DualTimeStepControls& c) {
    if (c.max_pseudo_iterations == 0) throw std::invalid_argument("dual-time iteration budget must be positive");
    if (!(c.pseudo_dt_initial > 0.0) || !std::isfinite(c.pseudo_dt_initial)) throw std::invalid_argument("invalid initial pseudo-time step");
    if (!(c.pseudo_dt_min > 0.0) || !std::isfinite(c.pseudo_dt_min)) throw std::invalid_argument("invalid minimum pseudo-time step");
    if (!(c.pseudo_dt_max >= c.pseudo_dt_min) || !std::isfinite(c.pseudo_dt_max)) throw std::invalid_argument("invalid maximum pseudo-time step");
    if (c.pseudo_dt_initial < c.pseudo_dt_min || c.pseudo_dt_initial > c.pseudo_dt_max) throw std::invalid_argument("initial pseudo-time step outside bounds");
    if (!(c.pseudo_dt_growth > 1.0) || !std::isfinite(c.pseudo_dt_growth)) throw std::invalid_argument("invalid pseudo-time growth factor");
    if (!(c.pseudo_dt_shrink > 0.0 && c.pseudo_dt_shrink < 1.0) || !std::isfinite(c.pseudo_dt_shrink)) throw std::invalid_argument("invalid pseudo-time shrink factor");
    if (!(c.absolute_tolerance >= 0.0) || !std::isfinite(c.absolute_tolerance)) throw std::invalid_argument("invalid absolute tolerance");
    if (!(c.relative_tolerance >= 0.0) || !std::isfinite(c.relative_tolerance)) throw std::invalid_argument("invalid relative tolerance");
    if (c.absolute_tolerance == 0.0 && c.relative_tolerance == 0.0) throw std::invalid_argument("dual-time requires a convergence tolerance");
}

inline double dual_time_norm(const std::vector<double>& r) {
    long double sum = 0.0L;
    for (double value : r) {
        if (!std::isfinite(value)) throw std::runtime_error("dual-time residual is non-finite");
        sum += static_cast<long double>(value) * static_cast<long double>(value);
    }
    return std::sqrt(static_cast<double>(sum));
}

} // namespace cfdx::physics


inline std::vector<double> dual_time_physical_residual(
    const std::vector<double>& state,
    const std::vector<double>& state_n,
    const std::vector<double>& state_nm1,
    double dt,
    DualTimePhysicalScheme scheme,
    const std::function<std::vector<double>(const std::vector<double>&)>& rhs)
{
    if (!(dt > 0.0) || !std::isfinite(dt))
        throw std::invalid_argument("dual-time physical dt must be finite and positive");
    if (state.size() != state_n.size())
        throw std::invalid_argument("dual-time state/history size mismatch");
    if (scheme == DualTimePhysicalScheme::BDF2 && state.size() != state_nm1.size())
        throw std::invalid_argument("dual-time BDF2 history size mismatch");
    const auto f = rhs(state);
    if (f.size() != state.size())
        throw std::invalid_argument("dual-time RHS size mismatch");

    std::vector<double> residual(state.size(), 0.0);
    if (scheme == DualTimePhysicalScheme::BACKWARD_EULER) {
        for (std::size_t i = 0; i < state.size(); ++i)
            residual[i] = state[i] - state_n[i] - dt * f[i];
    } else {
        for (std::size_t i = 0; i < state.size(); ++i)
            residual[i] = 1.5 * state[i] - 2.0 * state_n[i]
                        + 0.5 * state_nm1[i] - dt * f[i];
    }
    return residual;
}

} // namespace cfdx::physics
