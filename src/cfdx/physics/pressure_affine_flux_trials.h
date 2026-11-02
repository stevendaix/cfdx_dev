#pragma once

// Optional acceleration of the Rhie--Chow line-search, valid only while the
// predictor, momentum response and boundary operators are frozen.  The exact
// authoritative face flux is ALWAYS reconstructed after accepting a trial and
// checked against this affine prediction before the result may be used.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

class PressureAffineFluxTrials {
public:
    PressureAffineFluxTrials(const std::vector<double>& flux_at_zero,
                             const std::vector<double>& flux_at_nominal,
                             const std::vector<std::size_t>& owner,
                             const std::vector<std::int64_t>& neighbour,
                             std::size_t n_cells, double nominal_alpha)
        : phi_zero_(flux_at_zero), phi_delta_(flux_at_nominal.size()),
          balance_zero_(n_cells, 0.0), balance_delta_(n_cells, 0.0),
          nominal_alpha_(nominal_alpha)
    {
        if (!(nominal_alpha_ > 0.0) || !std::isfinite(nominal_alpha_) ||
            n_cells == 0 || flux_at_zero.empty() ||
            flux_at_zero.size() != flux_at_nominal.size() ||
            flux_at_zero.size() != owner.size() ||
            owner.size() != neighbour.size())
            throw std::invalid_argument("affine Rhie-Chow: invalid dimensions/alpha");
        for (std::size_t f = 0; f < owner.size(); ++f) {
            const double a = flux_at_zero[f], b = flux_at_nominal[f];
            if (!std::isfinite(a) || !std::isfinite(b) || owner[f] >= n_cells ||
                neighbour[f] < -1 ||
                (neighbour[f] >= 0 && static_cast<std::size_t>(neighbour[f]) >= n_cells))
                throw std::invalid_argument("affine Rhie-Chow: invalid face");
            const double d = b - a;
            if (!std::isfinite(d))
                throw std::invalid_argument("affine Rhie-Chow: nonfinite delta");
            phi_delta_[f] = d;
            balance_zero_[owner[f]] += a;
            balance_delta_[owner[f]] += d;
            if (neighbour[f] >= 0) {
                const auto n = static_cast<std::size_t>(neighbour[f]);
                balance_zero_[n] -= a;
                balance_delta_[n] -= d;
            }
        }
    }

    double continuity_norm(double alpha) const {
        const double t = ratio(alpha);
        double sum = 0.0;
        for (std::size_t c = 0; c < balance_zero_.size(); ++c) {
            const double b = std::fma(t, balance_delta_[c], balance_zero_[c]);
            if (!std::isfinite(b)) return std::numeric_limits<double>::infinity();
            sum += std::abs(b);
        }
        return std::isfinite(sum) ? sum : std::numeric_limits<double>::infinity();
    }

    // Exact flux comparison, including ALL boundary faces.  Normalisation is
    // relative to the scale of the two endpoint fluxes, not the small trial
    // flux (which may be close to zero due to cancellation).
    bool agrees_with_flux(double alpha, const std::vector<double>& actual,
                          double relative_tolerance = 2e-12) const {
        const double t = ratio(alpha);
        if (!(relative_tolerance >= 0.0) || !std::isfinite(relative_tolerance) ||
            actual.size() != phi_zero_.size()) return false;
        for (std::size_t f = 0; f < actual.size(); ++f) {
            const double predicted = std::fma(t, phi_delta_[f], phi_zero_[f]);
            const double scale = std::max({1.0, std::abs(phi_zero_[f]),
                                          std::abs(phi_delta_[f]), std::abs(actual[f])});
            if (!std::isfinite(predicted) || !std::isfinite(actual[f]) ||
                std::abs(predicted - actual[f]) > relative_tolerance*scale)
                return false;
        }
        return true;
    }
private:
    double ratio(double alpha) const {
        if (!std::isfinite(alpha) || alpha < 0.0 || alpha > nominal_alpha_)
            throw std::invalid_argument("affine Rhie-Chow: trial outside interval");
        return alpha / nominal_alpha_;
    }
    std::vector<double> phi_zero_, phi_delta_;
    std::vector<double> balance_zero_, balance_delta_;
    double nominal_alpha_;
};

}  // namespace cfdx::physics
