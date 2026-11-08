#pragma once

// Optional extended-precision *assembly*, not an extended-precision pressure
// solution. This only postpones rounding of face contributions to the final
// diagonal and right-hand side. No acceptance gate, pressure field or
// Rhie-Chow face flux is changed by this accumulator.
//
// Each cell accumulator is independent, enabling later thread-local or
// segmented GPU reductions without changing the finite-volume equations.

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

class PressureFaceAssemblyAccumulator {
public:
    explicit PressureFaceAssemblyAccumulator(std::size_t n)
        : internal_diagonal_(n, 0.0L),
          boundary_diagonal_(n, 0.0L),
          continuity_(n, 0.0L), boundary_rhs_(n, 0.0L) {
        if (n == 0 ||
            std::numeric_limits<long double>::digits <=
                std::numeric_limits<double>::digits)
            throw std::invalid_argument(
                "compensated pressure assembly requires a nonempty mesh "
                "and wider-than-double accumulation");
    }

    void add_continuity(std::size_t cell, double flux) {
        continuity_.at(cell) += static_cast<long double>(flux);
    }

    void add_internal_pair(std::size_t owner, std::size_t neighbour,
                           double coefficient) {
        if (owner == neighbour || !(coefficient > 0.0) ||
            !std::isfinite(coefficient))
            throw std::invalid_argument(
                "compensated pressure assembly: invalid face coefficient");
        internal_diagonal_.at(owner) += static_cast<long double>(coefficient);
        internal_diagonal_.at(neighbour) += static_cast<long double>(coefficient);
    }

    void add_pressure_boundary(std::size_t owner, double coefficient,
                               double rhs_source) {
        if (!(coefficient > 0.0) || !std::isfinite(coefficient) ||
            !std::isfinite(rhs_source))
            throw std::invalid_argument(
                "compensated pressure assembly: invalid boundary coefficient/source");
        boundary_diagonal_.at(owner) += static_cast<long double>(coefficient);
        boundary_rhs_.at(owner) += static_cast<long double>(rhs_source);
    }

    double diagonal(std::size_t cell) const {
        return narrow(internal_diagonal_.at(cell) + boundary_diagonal_.at(cell));
    }

    double rhs(std::size_t cell) const {
        return narrow(-continuity_.at(cell) + boundary_rhs_.at(cell));
    }

private:
    static double narrow(long double value) {
        if (!std::isfinite(value) ||
            std::abs(value) > std::numeric_limits<double>::max())
            throw std::overflow_error(
                "compensated pressure assembly: diagonal or RHS overflow");
        return static_cast<double>(value);
    }
    std::vector<long double> internal_diagonal_;
    std::vector<long double> boundary_diagonal_;
    std::vector<long double> continuity_;
    std::vector<long double> boundary_rhs_;
};

} // namespace cfdx::physics
