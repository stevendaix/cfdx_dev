#pragma once

// Optional double-double compensated *assembly*, not an extended-precision pressure
// solution. This postpones rounding of face contributions to the final
// diagonal and right-hand side. No acceptance gate, pressure field or
// Rhie-Chow face flux is changed by this accumulator.
//
// Each cell accumulator is independent, enabling later thread-local or
// segmented GPU reductions without changing the finite-volume equations.

#include "cfdx/core/numerics/portable_double_double.h"
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

class PressureFaceAssemblyAccumulator {
public:
    explicit PressureFaceAssemblyAccumulator(std::size_t n)
        : internal_diagonal_(n),
          boundary_diagonal_(n),
          continuity_(n), boundary_rhs_(n) {
        if (n == 0)
            throw std::invalid_argument(
                "compensated pressure assembly requires a nonempty mesh");
    }

    void add_continuity(std::size_t cell, double flux) {
        continuity_.at(cell).add(flux);
    }

    void add_internal_pair(std::size_t owner, std::size_t neighbour,
                           double coefficient) {
        if (owner == neighbour || !(coefficient > 0.0) ||
            !std::isfinite(coefficient))
            throw std::invalid_argument(
                "compensated pressure assembly: invalid face coefficient");
        internal_diagonal_.at(owner).add(coefficient);
        internal_diagonal_.at(neighbour).add(coefficient);
    }

    void add_pressure_boundary(std::size_t owner, double coefficient,
                               double rhs_source) {
        if (!(coefficient > 0.0) || !std::isfinite(coefficient) ||
            !std::isfinite(rhs_source))
            throw std::invalid_argument(
                "compensated pressure assembly: invalid boundary coefficient/source");
        boundary_diagonal_.at(owner).add(coefficient);
        boundary_rhs_.at(owner).add(rhs_source);
    }

    double diagonal(std::size_t cell) const {
        auto sum = internal_diagonal_.at(cell);
        sum.add(boundary_diagonal_.at(cell));
        return sum.narrow_to_double();
    }

    double rhs(std::size_t cell) const {
        auto sum = -continuity_.at(cell);
        sum.add(boundary_rhs_.at(cell));
        return sum.narrow_to_double();
    }

private:
    using Accumulator = cfdx::core::numerics::PortableDoubleDouble;
    std::vector<Accumulator> internal_diagonal_;
    std::vector<Accumulator> boundary_diagonal_;
    std::vector<Accumulator> continuity_;
    std::vector<Accumulator> boundary_rhs_;
};

} // namespace cfdx::physics
