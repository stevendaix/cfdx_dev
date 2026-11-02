#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace cfdx_n256 {

enum class Phase { robust_start, full_nonorth, second_order };

struct Config {
    double cfl_initial = 0.5;
    double cfl_min = 0.05;
    double cfl_max = 50.0;
    double cfl_growth = 1.15;
    double cfl_shrink = 0.5;
    double good_ratio = 0.85;
    double bad_ratio = 1.05;
    unsigned good_needed = 8;
    unsigned full_nonorth_needed = 12;
    unsigned final_needed = 15;
};

struct Controls {
    double cfl;
    double nonorth_factor;
    double second_order_weight;
    Phase phase;
    bool qualification_ready;
};

class Continuation {
    Config cfg_;
    double cfl_;
    double prev_ = -1.0;
    unsigned good_ = 0;
    unsigned stable_full_ = 0;
    Phase phase_ = Phase::robust_start;
    bool ready_ = false;
    double nonorth_ = 0.0;
    double order2_ = 0.0;

public:
    explicit Continuation(Config c = {}) : cfg_(c), cfl_(c.cfl_initial) {
        if (!(std::isfinite(c.cfl_min) && std::isfinite(c.cfl_initial) &&
              std::isfinite(c.cfl_max) && std::isfinite(c.cfl_growth) &&
              std::isfinite(c.cfl_shrink) && std::isfinite(c.good_ratio) &&
              std::isfinite(c.bad_ratio) &&
              c.cfl_min > 0.0 && c.cfl_min <= c.cfl_initial &&
              c.cfl_initial <= c.cfl_max && c.cfl_growth > 1.0 &&
              c.cfl_shrink > 0.0 && c.cfl_shrink < 1.0 &&
              c.good_ratio > 0.0 && c.good_ratio < 1.0 &&
              c.bad_ratio > 1.0 && c.good_needed > 0 &&
              c.full_nonorth_needed > 0 && c.final_needed > 0)) {
            throw std::invalid_argument("invalid continuation config");
        }
    }

    Controls controls() const {
        return {cfl_, nonorth_, order2_, phase_, ready_};
    }

    // Supply the true, physically normalized, complete-equation nonlinear
    // residual. accepted must be false for a rejected/rolled-back candidate.
    Controls observe(double true_residual, bool accepted,
                    bool physically_converged = false) {
        if (!std::isfinite(true_residual) || true_residual < 0.0) {
            throw std::invalid_argument("nonfinite/negative residual");
        }

        if (!accepted) {
            cfl_ = std::max(cfg_.cfl_min, cfl_ * cfg_.cfl_shrink);
            good_ = stable_full_ = 0;
            ready_ = false;
            return controls();
        }

        double ratio = 1.0;
        if (prev_ == 0.0) {
            ratio = true_residual == 0.0
                ? 0.0
                : std::numeric_limits<double>::infinity();
        } else if (prev_ > 0.0) {
            ratio = true_residual / prev_;
        }
        prev_ = true_residual;

        if (ratio > cfg_.bad_ratio) {
            cfl_ = std::max(cfg_.cfl_min, cfl_ * cfg_.cfl_shrink);
            good_ = stable_full_ = 0;
            ready_ = false;
        } else if (ratio < cfg_.good_ratio) {
            cfl_ = std::min(cfg_.cfl_max, cfl_ * cfg_.cfl_growth);
            ++good_;
        } else {
            good_ = 0;
        }

        if (phase_ == Phase::robust_start && good_ >= cfg_.good_needed) {
            nonorth_ = std::min(1.0, nonorth_ + 0.25);
            good_ = 0;
            if (nonorth_ == 1.0) {
                phase_ = Phase::full_nonorth;
            }
        } else if (phase_ == Phase::full_nonorth &&
                   good_ >= cfg_.full_nonorth_needed) {
            order2_ = std::min(1.0, order2_ + 0.25);
            good_ = 0;
            if (order2_ == 1.0) {
                phase_ = Phase::second_order;
                stable_full_ = 0;
            }
        }

        if (phase_ == Phase::second_order) {
            if (ratio <= cfg_.bad_ratio) {
                ++stable_full_;
            } else {
                stable_full_ = 0;
            }
            // Readiness describes the current accepted state, not a permanent
            // latch: losing the physical gate revokes qualification readiness.
            ready_ = stable_full_ >= cfg_.final_needed && physically_converged;
        } else {
            ready_ = false;
        }
        return controls();
    }
};

// Ratio of smaller/larger adjacent positive cell volumes; low means an
// abrupt jump between neighbouring cells.
inline double neighbor_volume_ratio(double a, double b) {
    if (!std::isfinite(a) || !std::isfinite(b) || a <= 0.0 || b <= 0.0) {
        throw std::invalid_argument("invalid cell volume");
    }
    return std::min(a, b) / std::max(a, b);
}

// For a wall-normal positive-spacing sequence, return the maximum successive
// growth ratio in either direction.
template<class It>
double max_adjacent_growth(It begin, It end) {
    if (begin == end) {
        throw std::invalid_argument("empty spacing");
    }
    double old = *begin++;
    if (!std::isfinite(old) || old <= 0.0) {
        throw std::invalid_argument("invalid spacing");
    }
    double worst = 1.0;
    for (; begin != end; ++begin) {
        const double next = *begin;
        if (!std::isfinite(next) || next <= 0.0) {
            throw std::invalid_argument("invalid spacing");
        }
        worst = std::max(worst, std::max(next / old, old / next));
        old = next;
    }
    return worst;
}

} // namespace cfdx_n256
