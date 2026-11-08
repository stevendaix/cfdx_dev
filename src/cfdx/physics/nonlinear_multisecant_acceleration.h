#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <deque>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics {

// Experimental preconditioned nonlinear multisecant *predictor* for SIMPLE.
// It never replaces an accepted solution; the ordinary SIMPLE pressure
// projection and the unfiltered TVD physical-residual gates still run.
// The states must contain the complete (cell velocity, cell pressure,
// conservative face mass flux) triple in a consistent flattening order.
enum class NonlinearMultisecantMethod { OFF, ANDERSON, NGMRES };

struct NonlinearMultisecantControls {
    NonlinearMultisecantMethod method = NonlinearMultisecantMethod::OFF;
    std::size_t depth = 4;
    std::size_t start_iteration = 15;
    double regularization = 1e-7;
    double mixing = 0.7;
    double maximum_coefficient = 2.5;
    // Guard against unlimited extrapolation of the last accepted state.
    // A proposal's infinity-norm change cannot exceed this factor times
    // the largest recent accepted SIMPLE state displacement (in scaled units).
    double maximum_step_multiplier = 2.0;
    // A trial is not committed if either of these physical channels grows
    // beyond this multiple of its last accepted value (above gate noise).
    double maximum_physical_growth = 1.5;
};

inline void validate_nonlinear_multisecant(
    const NonlinearMultisecantControls& c)
{
    if (c.depth < 2 || c.depth > 10 || c.start_iteration < 3 ||
        !std::isfinite(c.regularization) || c.regularization <= 0.0 ||
        c.regularization > 1e-1 ||
        !std::isfinite(c.mixing) || c.mixing <= 0.0 || c.mixing > 1.0 ||
        !std::isfinite(c.maximum_coefficient) || c.maximum_coefficient < 1.0 ||
        c.maximum_coefficient > 20.0 ||
        !std::isfinite(c.maximum_step_multiplier) ||
        c.maximum_step_multiplier < 0.1 || c.maximum_step_multiplier > 20.0 ||
        !std::isfinite(c.maximum_physical_growth) ||
        c.maximum_physical_growth < 1.0 ||
        c.maximum_physical_growth > 10.0)
        throw std::invalid_argument("invalid nonlinear multisecant controls");
}

struct MultisecantSample {
    std::vector<double> input;               // X_k (after any predictor)
    std::vector<double> output;              // G(X_k), accepted U/p/phi
    std::vector<double> physical_residual;   // F(G(X_k)), SIGNED, full TVD
};

struct MultisecantProposal {
    bool available = false;
    std::vector<double> state;
    std::vector<double> coefficients;
    double relative_step = 0.0;
    double reference_step = 0.0;
    double estimated_objective = 0.0;
    std::string reason;
};

// A tiny pivoted linear solver; never invoke a large dense numerical package
// for a history of at most 10 samples. Ill-conditioned fits fail closed.
inline bool multisecant_small_solve(
    std::vector<std::vector<double>> matrix,
    std::vector<double> rhs,
    std::vector<double>& solution)
{
    const std::size_t n = rhs.size();
    if (matrix.size() != n) return false;
    for (const auto& row : matrix) if (row.size() != n) return false;
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t r = k + 1; r < n; ++r)
            if (std::abs(matrix[r][k]) > std::abs(matrix[pivot][k]))
                pivot = r;
        if (!std::isfinite(matrix[pivot][k]) ||
            std::abs(matrix[pivot][k]) < 1e-14) return false;
        std::swap(matrix[pivot], matrix[k]);
        std::swap(rhs[pivot], rhs[k]);
        const double inv = 1.0 / matrix[k][k];
        for (std::size_t j = k; j < n; ++j) matrix[k][j] *= inv;
        rhs[k] *= inv;
        for (std::size_t i = 0; i < n; ++i) {
            if (i == k) continue;
            const double v = matrix[i][k];
            for (std::size_t j = k; j < n; ++j)
                matrix[i][j] -= v * matrix[k][j];
            rhs[i] -= v * rhs[k];
        }
    }
    for (const double x : rhs) if (!std::isfinite(x)) return false;
    solution = std::move(rhs);
    return true;
}

class NonlinearMultisecantHistory {
public:
    explicit NonlinearMultisecantHistory(NonlinearMultisecantControls controls)
        : controls_(controls)
    {
        validate_nonlinear_multisecant(controls_);
    }

    void reset() { samples_.clear(); scales_.clear(); }
    std::size_t size() const { return samples_.size(); }

    // "scales" is a fixed per-entry physical scale, not the varying trial
    // residual or RHS norm. This avoids winning a fit by inflating a scale.
    void append(MultisecantSample sample,
                const std::vector<double>& scales)
    {
        if (sample.input.empty() ||
            sample.input.size() != sample.output.size() ||
            sample.physical_residual.empty() ||
            scales.size() != sample.input.size())
            throw std::invalid_argument("multisecant: inconsistent sample dimensions");
        if (!samples_.empty() &&
            (samples_.back().input.size() != sample.input.size() ||
             samples_.back().physical_residual.size() !=
                 sample.physical_residual.size()))
            throw std::invalid_argument("multisecant: changed problem dimension");
        for (std::size_t j = 0; j < scales.size(); ++j)
            if (!(scales[j] > 0.0) || !std::isfinite(scales[j]))
                throw std::invalid_argument("multisecant: nonpositive state scale");
        for (const auto* v : {&sample.input, &sample.output,
                              &sample.physical_residual})
            for (const double x : *v)
                if (!std::isfinite(x))
                    throw std::invalid_argument("multisecant: nonfinite sample");
        if (scales_.empty()) scales_ = scales;
        // In a single run the scaling is frozen from the first accepted state.
        else if (scales_.size() != scales.size())
            throw std::invalid_argument("multisecant: changed state scaling");
        samples_.push_back(std::move(sample));
        while (samples_.size() > controls_.depth) samples_.pop_front();
    }

    MultisecantProposal propose() const
    {
        MultisecantProposal proposal;
        if (controls_.method == NonlinearMultisecantMethod::OFF ||
            samples_.size() < 2) {
            proposal.reason = "insufficient accepted history";
            return proposal;
        }
        const std::size_t m = samples_.size();
        const std::size_t n = samples_.back().output.size();
        const std::size_t residual_n =
            controls_.method == NonlinearMultisecantMethod::ANDERSON
                ? n : samples_.back().physical_residual.size();
        std::vector<std::vector<double>> residuals(m);
        for (std::size_t i = 0; i < m; ++i) {
            residuals[i].resize(residual_n);
            for (std::size_t j = 0; j < residual_n; ++j) {
                residuals[i][j] =
                    controls_.method == NonlinearMultisecantMethod::ANDERSON
                    ? (samples_[i].output[j] - samples_[i].input[j]) / scales_[j]
                    : samples_[i].physical_residual[j];
            }
        }
        const double normalizer = static_cast<double>(residual_n);
        std::vector<std::vector<double>> system(m + 1,
                                                 std::vector<double>(m + 1, 0.0));
        double maxdiag = 0.0;
        for (std::size_t i = 0; i < m; ++i) {
            for (std::size_t j = 0; j < m; ++j) {
                double dot = 0.0;
                for (std::size_t k = 0; k < residual_n; ++k)
                    dot += residuals[i][k] * residuals[j][k];
                system[i][j] = dot / normalizer;
            }
            maxdiag = std::max(maxdiag, system[i][i]);
            system[i][m] = 1.0;
            system[m][i] = 1.0;
        }
        if (!(maxdiag > 1e-32) || !std::isfinite(maxdiag)) {
            proposal.reason = "vanishing or nonfinite history residual";
            return proposal;
        }
        for (std::size_t i = 0; i < m; ++i)
            system[i][i] += controls_.regularization * maxdiag;
        std::vector<double> rhs(m + 1, 0.0), coeff;
        rhs[m] = 1.0;
        if (!multisecant_small_solve(system, rhs, coeff)) {
            proposal.reason = "ill-conditioned least squares";
            return proposal;
        }
        coeff.resize(m);
        for (const double a : coeff)
            if (!std::isfinite(a) || std::abs(a) > controls_.maximum_coefficient) {
                proposal.reason = "unbounded multisecant coefficients";
                return proposal;
            }
        const auto& last = samples_.back();
        double old_step = 0.0;
        for (std::size_t k = 0; k < n; ++k)
            old_step = std::max(old_step,
                                std::abs(last.output[k] - last.input[k]) / scales_[k]);
        if (!(old_step > 1e-14) || !std::isfinite(old_step)) {
            proposal.reason = "fixed point or unresolvable step";
            return proposal;
        }
        proposal.state = last.output;
        double step = 0.0;
        for (std::size_t k = 0; k < n; ++k) {
            double combination = 0.0;
            for (std::size_t i = 0; i < m; ++i)
                combination += coeff[i] * samples_[i].output[k];
            proposal.state[k] += controls_.mixing *
                                  (combination - last.output[k]);
            if (!std::isfinite(proposal.state[k])) {
                proposal.state.clear();
                proposal.reason = "nonfinite proposed state";
                return proposal;
            }
            step = std::max(step,
                std::abs(proposal.state[k] - last.output[k]) / scales_[k]);
        }
        // Trust radius protects both physical fields and face flux. A cap is
        // applied to the ENTIRE state vector, never to individual face fluxes;
        // therefore linear continuity conservation remains intact.
        const double cap = controls_.maximum_step_multiplier * old_step;
        if (step > cap) {
            const double theta = cap / step;
            for (std::size_t k = 0; k < n; ++k)
                proposal.state[k] = last.output[k] +
                                    theta * (proposal.state[k] - last.output[k]);
            step = cap;
        }
        if (!(step > 1e-14)) {
            proposal.state.clear();
            proposal.reason = "vanishing proposal";
            return proposal;
        }
        double objective = 0.0;
        for (std::size_t k = 0; k < residual_n; ++k) {
            double r = 0.0;
            for (std::size_t i = 0; i < m; ++i)
                r += coeff[i] * residuals[i][k];
            objective += r * r;
        }
        proposal.available = true;
        proposal.coefficients = std::move(coeff);
        proposal.relative_step = step;
        proposal.reference_step = old_step;
        proposal.estimated_objective = std::sqrt(objective / normalizer);
        proposal.reason = "candidate predictor; requires SIMPLE projection and physical gate";
        return proposal;
    }

private:
    NonlinearMultisecantControls controls_;
    std::deque<MultisecantSample> samples_;
    std::vector<double> scales_;
};

// A strongly nonlinear but low-overhead manufactured residual can be passed
// to this Jv oracle. The caller must re-evaluate its COMPLETE target TVD
// discretization (including limiter changes), not the frozen UPWIND matrix.
struct DirectionalResidualControls {
    double relative_step = 1e-5;
    double minimum_step = 1e-9;
    double maximum_step = 1e-3;
};

inline void validate_directional_residual(const DirectionalResidualControls& c)
{
    if (!std::isfinite(c.relative_step) || c.relative_step <= 0.0 ||
        !std::isfinite(c.minimum_step) || c.minimum_step <= 0.0 ||
        !std::isfinite(c.maximum_step) || c.maximum_step < c.minimum_step)
        throw std::invalid_argument("invalid directional residual step");
}

template<class Evaluator>
std::vector<double> nonlinear_directional_jvp(
    Evaluator&& residual, const std::vector<double>& state,
    const std::vector<double>& direction,
    const DirectionalResidualControls& controls = {})
{
    validate_directional_residual(controls);
    if (state.empty() || state.size() != direction.size())
        throw std::invalid_argument("Jv: state and direction dimensions differ");
    double norm_x = 0.0, norm_v = 0.0;
    for (std::size_t i=0;i<state.size();++i) {
        if (!std::isfinite(state[i]) || !std::isfinite(direction[i]))
            throw std::invalid_argument("Jv: nonfinite state or direction");
        norm_x = std::max(norm_x, std::abs(state[i]));
        norm_v = std::max(norm_v, std::abs(direction[i]));
    }
    if (!(norm_v > 0.0))
        throw std::invalid_argument("Jv: zero direction");
    const double step = std::clamp(
        controls.relative_step * (1.0 + norm_x) / norm_v,
        controls.minimum_step / norm_v, controls.maximum_step / norm_v);
    std::vector<double> plus = state, minus = state;
    for (std::size_t i=0;i<state.size();++i) {
        plus[i] += step * direction[i];
        minus[i] -= step * direction[i];
    }
    const auto fp = residual(plus);
    const auto fm = residual(minus);
    if (fp.empty() || fp.size() != fm.size())
        throw std::invalid_argument("Jv: residual size changed between trials");
    std::vector<double> result(fp.size());
    for (std::size_t i=0;i<fp.size();++i) {
        result[i] = (fp[i] - fm[i]) / (2.0 * step);
        if (!std::isfinite(result[i]))
            throw std::runtime_error("Jv: nonfinite directional product");
    }
    return result;
}

} // namespace cfdx::physics
