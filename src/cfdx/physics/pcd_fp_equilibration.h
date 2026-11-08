#pragma once

// Numerical experiment for the nested PCD pressure convection-diffusion solve.
// Left row equilibration changes conditioning, not the unknown or exact solution.
// A quality certificate must always be recomputed against the ORIGINAL Fp.

#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/linalg/preconditioner.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct PcdFpEquilibration {
    cfdx::core::SparseMatrix matrix;
    std::vector<double> row_multipliers;
    double smallest_row_max = std::numeric_limits<double>::infinity();
    double largest_row_max = 0.0;
    double smallest_diagonal_ratio = std::numeric_limits<double>::infinity();
    double largest_diagonal_ratio = 0.0;
};

inline PcdFpEquilibration equilibrate_pcd_fp_rows(
    const cfdx::core::SparseMatrix& original)
{
    const std::size_t n = original.n_rows();
    if (n == 0 || original.n_cols() != n)
        throw std::invalid_argument("PCD Fp row equilibration requires square nonempty matrix");
    PcdFpEquilibration result;
    result.row_multipliers.reserve(n);
    const auto* offsets = original.row_offsets_data();
    const auto* columns = original.columns_data();
    const auto* values = original.values_data();
    std::vector<cfdx::core::SparseMatrix::Index> ro(offsets, offsets + n + 1);
    std::vector<cfdx::core::SparseMatrix::Index> co(columns, columns + original.nnz());
    std::vector<double> va(values, values + original.nnz());
    for (std::size_t row = 0; row < n; ++row) {
        double max_entry = 0.0;
        double diagonal = 0.0;
        for (std::size_t k = offsets[row]; k < offsets[row + 1]; ++k) {
            max_entry = std::max(max_entry, std::abs(values[k]));
            if (columns[k] == row) diagonal = std::abs(values[k]);
        }
        if (!std::isfinite(max_entry) || !(max_entry > 1e-290))
            throw std::invalid_argument("PCD Fp has zero, nonfinite or underflow-scale row");
        const double scale = 1.0 / max_entry;
        result.row_multipliers.push_back(scale);
        result.smallest_row_max = std::min(result.smallest_row_max, max_entry);
        result.largest_row_max = std::max(result.largest_row_max, max_entry);
        result.smallest_diagonal_ratio = std::min(result.smallest_diagonal_ratio, diagonal / max_entry);
        result.largest_diagonal_ratio = std::max(result.largest_diagonal_ratio, diagonal / max_entry);
        for (std::size_t k = offsets[row]; k < offsets[row + 1]; ++k)
            va[k] *= scale;
    }
    result.matrix = cfdx::core::SparseMatrix::from_sorted_csr(n, n,
        std::move(ro), std::move(co), std::move(va));
    return result;
}

inline cfdx::core::Vector pcd_fp_scaled_rhs(
    const cfdx::core::Vector& rhs, const std::vector<double>& multipliers)
{
    if (rhs.size() != multipliers.size())
        throw std::invalid_argument("PCD Fp scaled RHS has inconsistent dimensions");
    cfdx::core::Vector scaled(rhs.size(), 0.0);
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        scaled(i) = rhs(i) * multipliers[i];
        if (!std::isfinite(scaled(i)))
            throw std::invalid_argument("PCD Fp scaled RHS is nonfinite");
    }
    return scaled;
}

// This gate diagnoses a preconditioner that catastrophically amplifies the
// image A*M^-1*r. It does not claim to predict a Krylov convergence rate.
// The threshold is intentionally broad: the original NACA N64 AMG action
// amplified the image by ~3e14, whereas ILU0 was O(1).
struct PcdFpActionProbe {
    bool applied = false;
    bool finite = false;
    double correction_l2 = 0.0;
    double image_l2 = 0.0;
    double residual_l2 = 0.0;
    double amplification = std::numeric_limits<double>::infinity();
};

inline PcdFpActionProbe pcd_fp_action_probe(
    const cfdx::core::SparseMatrix& matrix,
    cfdx::core::Preconditioner& preconditioner,
    const cfdx::core::Vector& residual) {
    if (matrix.n_rows() != matrix.n_cols() || residual.size() != matrix.n_rows())
        throw std::invalid_argument("PCD Fp action probe dimension mismatch");
    PcdFpActionProbe result;
    result.residual_l2 = residual.norm2();
    if (!std::isfinite(result.residual_l2)) return result;
    if (result.residual_l2 == 0.0) {
        result.applied = result.finite = true;
        result.amplification = 0.0;
        return result;
    }
    cfdx::core::Vector correction(residual.size(), 0.0);
    result.applied = preconditioner.apply(residual, correction);
    if (!result.applied || !correction.is_valid()) return result;
    result.correction_l2 = correction.norm2();
    double image_norm = 0.0;
    const auto* ro = matrix.row_offsets_data();
    const auto* co = matrix.columns_data();
    const auto* va = matrix.values_data();
    for (std::size_t row = 0; row < matrix.n_rows(); ++row) {
        long double image = 0.0L;
        for (std::size_t k = ro[row]; k < ro[row + 1]; ++k)
            image += static_cast<long double>(va[k]) * correction(co[k]);
        image_norm = std::hypot(image_norm, static_cast<double>(image));
    }
    result.image_l2 = image_norm;
    result.amplification = image_norm / result.residual_l2;
    result.finite = std::isfinite(result.correction_l2) &&
        std::isfinite(result.image_l2) && std::isfinite(result.amplification);
    return result;
}

struct PcdFpTrueResidual {
    double initial_l2 = 0.0;
    double final_l2 = 0.0;
    double rhs_l2 = 0.0;
    double final_linf = 0.0;
    double relative_initial = 0.0;
    double relative_rhs = 0.0;
};

inline PcdFpTrueResidual pcd_fp_true_residual(
    const cfdx::core::SparseMatrix& matrix,
    const cfdx::core::Vector& rhs,
    const cfdx::core::Vector& initial,
    const cfdx::core::Vector& solution)
{
    if (matrix.n_rows() != matrix.n_cols() ||
        rhs.size() != matrix.n_rows() ||
        initial.size() != rhs.size() || solution.size() != rhs.size())
        throw std::invalid_argument("PCD Fp residual has inconsistent dimensions");
    PcdFpTrueResidual r;
    const auto* offsets = matrix.row_offsets_data();
    const auto* columns = matrix.columns_data();
    const auto* values = matrix.values_data();
    for (std::size_t row = 0; row < rhs.size(); ++row) {
        long double ax_initial = 0.0L;
        long double ax = 0.0L;
        for (std::size_t k = offsets[row]; k < offsets[row + 1]; ++k) {
            ax_initial += static_cast<long double>(values[k]) * initial(columns[k]);
            ax += static_cast<long double>(values[k]) * solution(columns[k]);
        }
        const double e0 = static_cast<double>(static_cast<long double>(rhs(row)) - ax_initial);
        const double e = static_cast<double>(static_cast<long double>(rhs(row)) - ax);
        r.initial_l2 = std::hypot(r.initial_l2, e0);
        r.final_l2 = std::hypot(r.final_l2, e);
        r.rhs_l2 = std::hypot(r.rhs_l2, rhs(row));
        r.final_linf = std::max(r.final_linf, std::abs(e));
    }
    r.relative_initial = r.final_l2 / (r.initial_l2 > 0.0 ? r.initial_l2 : 1.0);
    r.relative_rhs = r.final_l2 / (r.rhs_l2 > 0.0 ? r.rhs_l2 : 1.0);
    return r;
}

} // namespace cfdx::physics
