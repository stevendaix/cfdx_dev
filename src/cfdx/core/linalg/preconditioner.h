#pragma once

#include "sparse_matrix.h"
#include "vector.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace cfdx {
namespace core {

/** Jacobi, Gauss-Seidel and ILU(0) preconditioners. */
class JacobiPreconditioner {
public:
    bool setup(const SparseMatrix& A) {
        if (A.n_rows() != A.n_cols()) return false;
        diag_.assign(A.n_rows(), 0.0);
        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            diag_[i] = A(i, i);
            if (diag_[i] == 0.0 || !std::isfinite(diag_[i])) return false;
        }
        return true;
    }

    bool apply(const Vector& r, Vector& z) const {
        if (r.size() != diag_.size()) return false;
        z.resize(r.size());
        for (std::size_t i = 0; i < r.size(); ++i) z(i) = r(i) / diag_[i];
        return true;
    }

private:
    std::vector<double> diag_;
};

class GaussSeidelPreconditioner {
public:
    bool setup(const SparseMatrix& A) {
        A_ = &A;
        return A.n_rows() == A.n_cols();
    }

    bool apply(const Vector& r, Vector& z) const {
        if (!A_ || r.size() != A_->n_rows()) return false;
        z.resize(r.size());

        const auto* v = A_->values_data();
        const auto* c = A_->columns_data();
        const auto* ro = A_->row_offsets_data();

        for (std::size_t i = 0; i < A_->n_rows(); ++i) {
            double rhs = r(i);
            double diag = 0.0;
            for (std::size_t k = ro[i]; k < ro[i + 1]; ++k) {
                if (c[k] == i) diag = v[k];
                else if (c[k] < i) rhs -= v[k] * z(c[k]);
            }
            if (diag == 0.0 || !std::isfinite(diag)) return false;
            z(i) = rhs / diag;
        }
        return true;
    }

private:
    const SparseMatrix* A_ = nullptr;
};

/**
 * ILU(0) preconditioner.
 *
 * The factorization preserves exactly the CSR sparsity pattern of A.
 * L has an implicit unit diagonal; U contains the diagonal and upper
 * triangular entries. Factor values are stored contiguously in one array
 * and row/column metadata are contiguous as well.
 */
class ILU0Preconditioner {
public:
    bool setup(const SparseMatrix& A) {
        if (A.n_rows() != A.n_cols() || !A.is_consistent()) return false;

        n_ = A.n_rows();
        row_offsets_.assign(
            A.row_offsets_data(), A.row_offsets_data() + n_ + 1);
        columns_.assign(
            A.columns_data(), A.columns_data() + A.nnz());
        lu_values_.assign(
            A.values_data(), A.values_data() + A.nnz());
        diagonal_.assign(n_, invalid_index());

        // Locate all diagonal entries once. CSR columns are sorted by contract.
        for (std::size_t i = 0; i < n_; ++i) {
            const auto begin = row_offsets_[i];
            const auto end = row_offsets_[i + 1];
            const auto it = std::lower_bound(
                columns_.begin() + begin,
                columns_.begin() + end,
                static_cast<SparseMatrix::Index>(i));
            if (it == columns_.begin() + end ||
                *it != static_cast<SparseMatrix::Index>(i)) {
                return false;
            }
            const std::size_t p = static_cast<std::size_t>(
                std::distance(columns_.begin(), it));
            diagonal_[i] = p;
            if (!std::isfinite(lu_values_[p]) ||
                std::abs(lu_values_[p]) <= pivot_tolerance_) {
                return false;
            }
        }

        // In-place ILU(0). Only entries already present in A may be updated.
        for (std::size_t i = 0; i < n_; ++i) {
            const auto row_begin = row_offsets_[i];
            const auto row_end = row_offsets_[i + 1];

            for (std::size_t p = row_begin; p < row_end; ++p) {
                const std::size_t j = columns_[p];
                if (j >= i) break;  // L part is strictly lower triangular.

                const std::size_t dj = diagonal_[j];
                const double ujj = lu_values_[dj];
                if (!std::isfinite(ujj) || std::abs(ujj) <= pivot_tolerance_) {
                    return false;
                }

                const double lij = lu_values_[p] / ujj;
                lu_values_[p] = lij;

                // A_ik <- A_ik - L_ij U_jk for existing (i,k).
                for (std::size_t q = dj + 1; q < row_offsets_[j + 1]; ++q) {
                    const std::size_t k = columns_[q];
                    const auto target = std::lower_bound(
                        columns_.begin() + row_begin,
                        columns_.begin() + row_end,
                        static_cast<SparseMatrix::Index>(k));
                    if (target != columns_.begin() + row_end &&
                        *target == static_cast<SparseMatrix::Index>(k)) {
                        const std::size_t t = static_cast<std::size_t>(
                            std::distance(columns_.begin(), target));
                        lu_values_[t] -= lij * lu_values_[q];
                    }
                }
            }

            const double diag = lu_values_[diagonal_[i]];
            if (!std::isfinite(diag) || std::abs(diag) <= pivot_tolerance_) {
                return false;
            }
        }

        return true;
    }

    bool apply(const Vector& r, Vector& z) const {
        if (r.size() != n_ || diagonal_.size() != n_) return false;

        z.resize(n_);

        // Forward solve: L y = r, with unit diagonal.
        for (std::size_t i = 0; i < n_; ++i) {
            double value = r(i);
            for (std::size_t p = row_offsets_[i];
                 p < diagonal_[i]; ++p) {
                value -= lu_values_[p] * z(columns_[p]);
            }
            z(i) = value;
        }

        // Backward solve: U z = y.
        for (std::size_t ii = n_; ii-- > 0;) {
            const std::size_t i = ii;
            double value = z(i);
            for (std::size_t p = diagonal_[i] + 1;
                 p < row_offsets_[i + 1]; ++p) {
                value -= lu_values_[p] * z(columns_[p]);
            }
            const double diag = lu_values_[diagonal_[i]];
            if (!std::isfinite(diag) || std::abs(diag) <= pivot_tolerance_) {
                return false;
            }
            z(i) = value / diag;
        }

        return true;
    }

private:
    static constexpr double pivot_tolerance_ = 1e-30;

    static constexpr std::uint32_t invalid_index() noexcept {
        return std::numeric_limits<std::uint32_t>::max();
    }

    std::size_t n_ = 0;
    std::vector<std::uint32_t> row_offsets_;
    std::vector<SparseMatrix::Index> columns_;
    std::vector<double> lu_values_;
    std::vector<std::uint32_t> diagonal_;
};

} // namespace core
} // namespace cfdx
