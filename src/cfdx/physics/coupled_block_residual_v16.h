#pragma once

#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {

// Read-only diagnostic for the *original* coupled system A x = b.  Velocity
// and pressure/continuity row norms are reported separately, but acceptance
// always uses the full-system true residual.  A gauge row remains included in
// this algebraic residual; it is NOT mistaken for a physical continuity row.
struct CoupledBlockResidualV16 {
    struct Block {
        long double absolute_l2 = 0.0L;
        long double rhs_l2 = 0.0L;
        long double linf = 0.0L;
        std::size_t worst_row = 0;

        long double relative_to_block_rhs() const noexcept {
            return absolute_l2 / (rhs_l2 > 0.0L ? rhs_l2 : 1.0L);
        }
    };
    Block momentum;
    Block continuity_and_gauge;
    long double global_l2 = 0.0L;
    long double global_rhs_l2 = 0.0L;

    long double relative_to_global_rhs() const noexcept {
        return global_l2 / (global_rhs_l2 > 0.0L ? global_rhs_l2 : 1.0L);
    }
    long double momentum_fraction_of_squared_defect() const noexcept {
        return global_l2 > 0.0L
            ? (momentum.absolute_l2 / global_l2) * (momentum.absolute_l2 / global_l2)
            : 0.0L;
    }
};

inline CoupledBlockResidualV16 independent_coupled_block_residual_v16(
    const cfdx::core::SparseMatrix& A,
    const cfdx::core::Vector& b,
    const cfdx::core::Vector& x,
    std::size_t velocity_row_count)
{
    if (!A.is_consistent() || A.n_rows() != A.n_cols() ||
        A.n_rows() != b.size() || A.n_cols() != x.size() ||
        velocity_row_count == 0 || velocity_row_count >= A.n_rows())
        throw std::invalid_argument("COUPLED V16 block certificate: invalid matrix or dimensions");

    CoupledBlockResidualV16 result;
    for (std::size_t row = 0; row < A.n_rows(); ++row) {
        auto& block = row < velocity_row_count
            ? result.momentum : result.continuity_and_gauge;
        long double ax = 0.0L;
        for (std::size_t k = A.row_offsets_data()[row];
             k < A.row_offsets_data()[row + 1]; ++k) {
            const double value = A.values_data()[k];
            const double xv = x(A.columns_data()[k]);
            if (!std::isfinite(value) || !std::isfinite(xv))
                throw std::runtime_error("COUPLED V16 block certificate: nonfinite matrix or solution");
            ax += static_cast<long double>(value) * xv;
        }
        if (!std::isfinite(b(row)))
            throw std::runtime_error("COUPLED V16 block certificate: nonfinite right-hand side");
        const long double defect = static_cast<long double>(b(row)) - ax;
        const long double magnitude = std::abs(defect);
        if (magnitude > block.linf) {
            block.linf = magnitude;
            block.worst_row = row;
        }
        block.absolute_l2 = std::hypot(block.absolute_l2, defect);
        block.rhs_l2 = std::hypot(block.rhs_l2, static_cast<long double>(b(row)));
    }
    result.global_l2 = std::hypot(result.momentum.absolute_l2,
                                  result.continuity_and_gauge.absolute_l2);
    result.global_rhs_l2 = std::hypot(result.momentum.rhs_l2,
                                     result.continuity_and_gauge.rhs_l2);
    if (!std::isfinite(result.global_l2) ||
        !std::isfinite(result.global_rhs_l2) ||
        !std::isfinite(static_cast<double>(result.relative_to_global_rhs())))
        throw std::runtime_error("COUPLED V16 block certificate: nonfinite residual norm");
    return result;
}

} // namespace cfdx::physics
