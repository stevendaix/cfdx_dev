#pragma once

#include "cfdx/core/linalg/preconditioner.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace cfdx::core {

// Production-quality stationary Gauss-Seidel preconditioner.
class GaussSeidelPreconditioner final : public Preconditioner {
public:
    bool setup(const SparseMatrix& A) override {
        if (A.n_rows() != A.n_cols() || !A.is_consistent()) return false;
        n_ = A.n_rows();
        diag_.assign(n_, 0.0);
        lower_rows_.assign(n_, {});

        const auto* ro = A.row_offsets_data();
        const auto* ci = A.columns_data();
        const auto* av = A.values_data();

        for (std::size_t i = 0; i < n_; ++i) {
            bool found = false;
            std::uint32_t previous = 0;
            for (std::size_t k = ro[i]; k < ro[i + 1]; ++k) {
                const std::uint32_t j = ci[k];
                if (!std::isfinite(av[k])) {
                    clear();
                    return false;
                }
                if (k > ro[i] && j <= previous) {
                    clear();
                    return false;
                }
                previous = j;
                if (j == i) {
                    if (found) {
                        clear();
                        return false;
                    }
                    diag_[i] = av[k];
                    found = true;
                } else if (j < i) {
                    lower_rows_[i].push_back({j, av[k]});
                }
            }
            if (!found || std::abs(diag_[i]) <= std::numeric_limits<double>::epsilon()) {
                clear();
                return false;
            }
        }
        return true;
    }

    bool apply(const Vector& r, Vector& z) const override {
        if (r.size() != n_ || z.size() != n_) return false;
        for (std::size_t i = 0; i < n_; ++i) {
            double rhs = r(i);
            if (!std::isfinite(rhs)) return false;
            for (const auto& [j, aij] : lower_rows_[i]) rhs -= aij * z(j);
            z(i) = rhs / diag_[i];
            if (!std::isfinite(z(i))) return false;
        }
        return true;
    }

    const char* name() const override { return "Gauss-Seidel"; }

private:
    void clear() {
        n_ = 0;
        diag_.clear();
        lower_rows_.clear();
    }

    std::size_t n_ = 0;
    std::vector<double> diag_;
    std::vector<std::vector<std::pair<std::size_t, double>>> lower_rows_;
};

// ILU(0) in the existing CSR sparsity pattern. The factorization is stored
// in-place as L (strict lower) and U (diagonal + upper). No fill is created.
class ILU0Preconditioner final : public Preconditioner {
public:
    bool setup(const SparseMatrix& A) override {
        if (A.n_rows() != A.n_cols() || !A.is_consistent()) {
            clear();
            return false;
        }

        n_ = A.n_rows();
        row_offsets_.assign(A.row_offsets_data(), A.row_offsets_data() + n_ + 1);
        columns_.assign(A.columns_data(), A.columns_data() + A.nnz());
        lu_.assign(A.values_data(), A.values_data() + A.nnz());
        diag_pos_.assign(n_, invalid_);

        for (std::size_t i = 0; i < n_; ++i) {
            std::uint32_t previous = 0;
            bool found = false;
            for (std::size_t k = row_offsets_[i]; k < row_offsets_[i + 1]; ++k) {
                const auto col = columns_[k];
                if (!std::isfinite(lu_[k])) {
                    clear();
                    return false;
                }
                if (k > row_offsets_[i] && col <= previous) {
                    clear();
                    return false;
                }
                previous = col;
                if (col == i) {
                    if (found) {
                        clear();
                        return false;
                    }
                    diag_pos_[i] = k;
                    found = true;
                }
            }
            if (!found) {
                clear();
                return false;
            }
        }

        for (std::size_t i = 0; i < n_; ++i) {
            for (std::size_t kk = row_offsets_[i]; kk < row_offsets_[i + 1]; ++kk) {
                const std::size_t kcol = columns_[kk];
                if (kcol >= i) break;

                const std::size_t dk = diag_pos_[kcol];
                if (!std::isfinite(lu_[dk]) ||
                    std::abs(lu_[dk]) <= std::numeric_limits<double>::epsilon()) {
                    clear();
                    return false;
                }

                lu_[kk] /= lu_[dk];
                if (!std::isfinite(lu_[kk])) {
                    clear();
                    return false;
                }

                const double lik = lu_[kk];
                for (std::size_t kj = row_offsets_[kcol]; kj < row_offsets_[kcol + 1]; ++kj) {
                    const std::size_t j = columns_[kj];
                    if (j <= kcol) continue;
                    const auto pos = find_position(i, j);
                    if (pos != invalid_) {
                        lu_[pos] -= lik * lu_[kj];
                        if (!std::isfinite(lu_[pos])) {
                            clear();
                            return false;
                        }
                    }
                }
            }

            const double diag = lu_[diag_pos_[i]];
            if (!std::isfinite(diag) ||
                std::abs(diag) <= std::numeric_limits<double>::epsilon()) {
                clear();
                return false;
            }
        }
        return true;
    }

    bool apply(const Vector& r, Vector& z) const override {
        if (r.size() != n_ || z.size() != n_) return false;

        std::vector<double> y(n_, 0.0);
        for (std::size_t i = 0; i < n_; ++i) {
            double v = r(i);
            if (!std::isfinite(v)) return false;
            for (std::size_t k = row_offsets_[i]; k < row_offsets_[i + 1]; ++k) {
                const std::size_t j = columns_[k];
                if (j >= i) break;
                v -= lu_[k] * y[j];
            }
            if (!std::isfinite(v)) return false;
            y[i] = v;
        }

        for (std::size_t ii = n_; ii-- > 0;) {
            const std::size_t i = ii;
            double v = y[i];
            for (std::size_t k = row_offsets_[i]; k < row_offsets_[i + 1]; ++k) {
                const std::size_t j = columns_[k];
                if (j <= i) continue;
                v -= lu_[k] * z(j);
            }
            z(i) = v / lu_[diag_pos_[i]];
            if (!std::isfinite(z(i))) return false;
        }
        return true;
    }

    const char* name() const override { return "ILU(0)"; }

private:
    static constexpr std::size_t invalid_ = std::numeric_limits<std::size_t>::max();

    std::size_t find_position(std::size_t row, std::size_t col) const {
        for (std::size_t k = row_offsets_[row]; k < row_offsets_[row + 1]; ++k) {
            if (columns_[k] == col) return k;
            if (columns_[k] > col) break;
        }
        return invalid_;
    }

    void clear() {
        n_ = 0;
        row_offsets_.clear();
        columns_.clear();
        lu_.clear();
        diag_pos_.clear();
    }

    std::size_t n_ = 0;
    std::vector<std::size_t> row_offsets_;
    std::vector<std::uint32_t> columns_;
    std::vector<double> lu_;
    std::vector<std::size_t> diag_pos_;
};

} // namespace cfdx::core
