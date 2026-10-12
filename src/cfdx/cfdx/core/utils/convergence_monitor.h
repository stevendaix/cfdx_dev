#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::core {

enum class ConvergenceStatus { CONTINUE, CONVERGED, STAGNATED, DIVERGED, MAX_ITERATIONS };

inline const char* to_string(ConvergenceStatus status) {
    switch (status) {
    case ConvergenceStatus::CONTINUE: return "CONTINUE";
    case ConvergenceStatus::CONVERGED: return "CONVERGED";
    case ConvergenceStatus::STAGNATED: return "STAGNATED";
    case ConvergenceStatus::DIVERGED: return "DIVERGED";
    case ConvergenceStatus::MAX_ITERATIONS: return "MAX_ITERATIONS";
    }
    return "UNKNOWN";
}

struct ConvergenceGate { double absolute_tolerance = 1e-10; double relative_tolerance = 1e-6; };

struct QoIGate {
    std::string name;
    double absolute_tolerance = 1e-10;
    double relative_tolerance = 1e-6;
    double target = std::numeric_limits<double>::quiet_NaN();
    bool require_target = false;
};

struct ConvergenceMonitorControls {
    ConvergenceGate residual;
    double conservation_tolerance = 1e-10;
    std::size_t minimum_iterations = 1;
    std::size_t max_iterations = 1000;
    std::size_t stagnation_window = 10;
    double stagnation_relative_improvement = 1e-3;
    std::size_t divergence_window = 3;
    double divergence_growth_factor = 10.0;
    std::vector<QoIGate> qoi_gates;
};

struct ConvergenceSample {
    std::size_t iteration = 0;
    double residual = std::numeric_limits<double>::infinity();
    double conservation = std::numeric_limits<double>::infinity();
    std::map<std::string, double> qoi_values;
};

struct ConvergenceReport {
    ConvergenceStatus status = ConvergenceStatus::CONTINUE;
    std::size_t iteration = 0;
    double residual = std::numeric_limits<double>::infinity();
    double normalized_residual = std::numeric_limits<double>::infinity();
    double conservation = std::numeric_limits<double>::infinity();
    std::string reason;
    std::size_t stagnation_count = 0;
    std::size_t divergence_count = 0;
    bool residual_gate = false;
    bool conservation_gate = false;
    bool qoi_gate = false;
};

inline void validate_convergence_monitor_controls(const ConvergenceMonitorControls& c) {
    if (!(c.residual.absolute_tolerance > 0.0) ||
        !(c.residual.relative_tolerance > 0.0) ||
        !(c.conservation_tolerance > 0.0) ||
        c.max_iterations == 0 || c.minimum_iterations > c.max_iterations ||
        c.stagnation_window == 0 || c.stagnation_relative_improvement < 0.0 ||
        c.divergence_window == 0 || c.divergence_growth_factor <= 1.0 ||
        !std::isfinite(c.residual.absolute_tolerance) ||
        !std::isfinite(c.residual.relative_tolerance) ||
        !std::isfinite(c.conservation_tolerance)) {
        throw std::invalid_argument("invalid convergence monitor controls");
    }
    for (const auto& qoi : c.qoi_gates) {
        if (qoi.name.empty() || qoi.absolute_tolerance <= 0.0 ||
            qoi.relative_tolerance <= 0.0 || !std::isfinite(qoi.absolute_tolerance) ||
            !std::isfinite(qoi.relative_tolerance) ||
            (qoi.require_target && !std::isfinite(qoi.target))) {
            throw std::invalid_argument("invalid QoI convergence gate");
        }
    }
}

class ConvergenceMonitor {
public:
    explicit ConvergenceMonitor(ConvergenceMonitorControls controls = {})
        : controls_(std::move(controls)) {
        validate_convergence_monitor_controls(controls_);
    }

    void reset() {
        history_.clear();
        initial_residual_ = std::numeric_limits<double>::quiet_NaN();
        status_ = ConvergenceStatus::CONTINUE;
    }

    ConvergenceReport update(const ConvergenceSample& sample) {
        validate_sample(sample);
        if (!history_.empty() && sample.iteration <= history_.back().iteration)
            throw std::invalid_argument("convergence iterations must be strictly increasing");
        if (history_.empty())
            initial_residual_ = std::max(sample.residual, std::numeric_limits<double>::min());

        history_.push_back(sample);
        const double normalized =
            sample.residual / std::max(initial_residual_, std::numeric_limits<double>::min());
        const bool residual_ok =
            sample.residual <= controls_.residual.absolute_tolerance ||
            normalized <= controls_.residual.relative_tolerance;
        const bool conservation_ok = sample.conservation <= controls_.conservation_tolerance;
        const bool qoi_ok = qois_converged(sample);
        const std::size_t stagnation_count = count_stagnation();
        const std::size_t divergence_count = count_divergence();

        ConvergenceReport report;
        report.iteration = sample.iteration;
        report.residual = sample.residual;
        report.normalized_residual = normalized;
        report.conservation = sample.conservation;
        report.residual_gate = residual_ok;
        report.conservation_gate = conservation_ok;
        report.qoi_gate = qoi_ok;
        report.stagnation_count = stagnation_count;
        report.divergence_count = divergence_count;

        if (sample.iteration >= controls_.minimum_iterations &&
            residual_ok && conservation_ok && qoi_ok) {
            status_ = ConvergenceStatus::CONVERGED;
            report.reason = "all configured convergence gates satisfied";
        } else if (divergence_count >= controls_.divergence_window) {
            status_ = ConvergenceStatus::DIVERGED;
            report.reason = "residual growth exceeded the configured divergence criterion";
        } else if (stagnation_count >= controls_.stagnation_window) {
            status_ = ConvergenceStatus::STAGNATED;
            report.reason = "insufficient residual improvement over the configured stagnation window";
        } else if (sample.iteration >= controls_.max_iterations) {
            status_ = ConvergenceStatus::MAX_ITERATIONS;
            report.reason = "maximum nonlinear iteration count reached without satisfying all gates";
        } else {
            status_ = ConvergenceStatus::CONTINUE;
            report.reason = "convergence gates not yet satisfied";
        }
        report.status = status_;
        return report;
    }

    const std::vector<ConvergenceSample>& history() const { return history_; }
    ConvergenceStatus status() const { return status_; }

    std::string deterministic_report() const {
        std::ostringstream out;
        out << std::setprecision(17) << "status=" << to_string(status_);
        if (!history_.empty()) {
            const auto& s = history_.back();
            out << "\niteration=" << s.iteration
                << "\nresidual=" << s.residual
                << "\nconservation=" << s.conservation;
            for (const auto& [name, value] : s.qoi_values)
                out << "\nqoi." << name << "=" << value;
        }
        return out.str();
    }

private:
    static void validate_sample(const ConvergenceSample& sample) {
        if (!std::isfinite(sample.residual) || sample.residual < 0.0 ||
            !std::isfinite(sample.conservation) || sample.conservation < 0.0)
            throw std::invalid_argument("non-finite or negative convergence metric");
        for (const auto& [name, value] : sample.qoi_values)
            if (name.empty() || !std::isfinite(value))
                throw std::invalid_argument("invalid QoI convergence sample");
    }

    bool qois_converged(const ConvergenceSample& sample) const {
        for (const auto& gate : controls_.qoi_gates) {
            const auto it = sample.qoi_values.find(gate.name);
            if (it == sample.qoi_values.end())
                return false;
            const double value = it->second;
            if (gate.require_target) {
                const double scale = std::max(std::abs(gate.target), 1.0);
                if (std::abs(value - gate.target) >
                    gate.absolute_tolerance + gate.relative_tolerance * scale)
                    return false;
            }
            if (history_.size() < 2)
                return false;
            const auto& previous = history_[history_.size() - 2];
            const auto pit = previous.qoi_values.find(gate.name);
            if (pit == previous.qoi_values.end())
                return false;
            const double scale = std::max(std::abs(value), 1.0);
            if (std::abs(value - pit->second) >
                gate.absolute_tolerance + gate.relative_tolerance * scale)
                return false;
        }
        return true;
    }

    std::size_t count_stagnation() const {
        if (history_.size() < controls_.stagnation_window + 1)
            return 0;
        std::size_t count = 0;
        const std::size_t begin = history_.size() - controls_.stagnation_window;
        for (std::size_t i = begin; i < history_.size(); ++i) {
            const double previous = history_[i - 1].residual;
            const double current = history_[i].residual;
            const double improvement =
                (previous - current) / std::max(previous, std::numeric_limits<double>::min());
            if (improvement < controls_.stagnation_relative_improvement)
                ++count;
        }
        return count;
    }

    std::size_t count_divergence() const {
        if (history_.size() < 2)
            return 0;
        std::size_t count = 0;
        const std::size_t begin =
            history_.size() > controls_.divergence_window
                ? history_.size() - controls_.divergence_window : 1;
        for (std::size_t i = begin; i < history_.size(); ++i)
            if (history_[i].residual >
                controls_.divergence_growth_factor * history_[i - 1].residual)
                ++count;
        return count;
    }

    ConvergenceMonitorControls controls_;
    std::vector<ConvergenceSample> history_;
    double initial_residual_ = std::numeric_limits<double>::quiet_NaN();
    ConvergenceStatus status_ = ConvergenceStatus::CONTINUE;
};

} // namespace cfdx::core
