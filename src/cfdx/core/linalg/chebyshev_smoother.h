#pragma once

#include "cfdx/core/linalg/linear_operator.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::core {

class ChebyshevSmoother {
public:
    struct Controls {
        std::size_t degree = 2;
        double lambda_min = 0.05;
        double lambda_max = 0.0; // <= 0: estimate from operator
        double damping = 1.0;
        std::size_t spectral_iterations = 8;
    };

    ChebyshevSmoother() = default;

    explicit ChebyshevSmoother(const Controls& controls) : c_(controls)
    {
        validate_controls();
    }

    void apply(const LinearOperatorBase& A,
               const Vector& rhs,
               Vector& x) const
    {
        if (rhs.size() != A.rows() || A.rows() != A.cols())
            throw std::invalid_argument("Chebyshev dimension mismatch");
        if (x.size() != rhs.size()) {
            x.resize(rhs.size());
            x.fill(0.0);
        }

        const double lmax =
            c_.lambda_max > 0.0 ? c_.lambda_max : estimate_lambda_max(A);
        validate_spectrum(lmax);

        const double theta = 0.5 * (lmax + c_.lambda_min);
        const double delta = 0.5 * (lmax - c_.lambda_min);
        Vector r(rhs.size());
        Vector Ar(rhs.size());
        Vector xprev(x);

        double alpha = 1.0 / theta;
        double beta = 0.0;
        for (std::size_t k = 0; k < c_.degree; ++k) {
            A.apply(x, Ar);
            for (std::size_t i = 0; i < rhs.size(); ++i)
                r(i) = rhs(i) - Ar(i);

            if (k == 0) {
                alpha = 1.0 / theta;
                beta = 0.0;
            } else {
                beta = std::pow(0.5 * delta * alpha, 2.0);
                alpha = 1.0 / (theta - beta * delta);
            }

            for (std::size_t i = 0; i < rhs.size(); ++i) {
                const double old = x(i);
                x(i) += c_.damping * alpha * r(i) +
                        beta * (x(i) - xprev(i));
                xprev(i) = old;
            }
        }
    }

    double lambda_max(const LinearOperatorBase& A) const
    {
        return c_.lambda_max > 0.0 ? c_.lambda_max : estimate_lambda_max(A);
    }

private:
    void validate_controls() const
    {
        if (c_.degree == 0 || !(c_.lambda_min > 0.0) ||
            !(c_.lambda_max >= 0.0) || !(c_.damping > 0.0) ||
            c_.spectral_iterations == 0)
            throw std::invalid_argument("invalid Chebyshev controls");
    }

    void validate_spectrum(double lmax) const
    {
        if (!(lmax >= c_.lambda_min) || !std::isfinite(lmax))
            throw std::invalid_argument("invalid Chebyshev spectrum");
    }

    double estimate_lambda_max(const LinearOperatorBase& A) const
    {
        Vector q(A.cols(), 1.0);
        Vector Aq(A.rows());
        const double inv_sqrt_n =
            1.0 / std::sqrt(static_cast<double>(q.size()));
        for (std::size_t i = 0; i < q.size(); ++i)
            q(i) *= inv_sqrt_n;

        double estimate = c_.lambda_min;
        for (std::size_t iter = 0; iter < c_.spectral_iterations; ++iter) {
            A.apply(q, Aq);
            const double norm = Aq.norm2();
            if (!(norm > 0.0) || !std::isfinite(norm))
                return c_.lambda_min;
            estimate = std::max(estimate, std::abs(q.dot(Aq)));
            for (std::size_t i = 0; i < q.size(); ++i)
                q(i) = Aq(i) / norm;
        }
        return std::max(estimate, c_.lambda_min);
    }

    Controls c_{};
};

} // namespace cfdx::core
