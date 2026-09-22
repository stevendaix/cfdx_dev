#pragma once
#include "cfdx/core/linalg/linear_operator.h"
#include "cfdx/core/linalg/preconditioner.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <vector>

namespace cfdx::core {

// Multilevel agglomerated AMG using piecewise-constant aggregation and
// Galerkin coarse operators. The finest operator remains matrix-free during
// apply(); SparseMatrix is used during setup to build the hierarchy.
// This is a lightweight AMG implementation: aggregation is pairwise and
// coarsest solves use damped Jacobi rather than a direct sparse factorization.
class MatrixFreeVcyclePreconditioner final : public Preconditioner {
public:
    MatrixFreeVcyclePreconditioner(const LinearOperatorBase& op,
                                   double omega = 0.7,
                                   std::size_t pre = 2,
                                   std::size_t post = 2)
        : op_(op), omega_(omega), pre_(pre), post_(post)
    {
        if (op.rows() != op.cols() || !(omega_ > 0.0 && omega_ < 2.0)) {
            throw std::invalid_argument("invalid AMG operator");
        }
    }

    bool setup(const SparseMatrix& A) override
    {
        if (A.n_rows() != op_.rows() || A.n_cols() != op_.cols() ||
            A.n_rows() == 0 || !matrix_is_valid(A)) {
            levels_.clear();
            first_aggregate_.clear();
            return false;
        }

        levels_.clear();
        first_aggregate_.clear();

        Level fine;
        fine.A = A;
        if (!build_diagonal(fine)) {
            levels_.clear();
            return false;
        }
        levels_.push_back(std::move(fine));

        // Build a genuine multilevel hierarchy until the coarse problem is
        // small. Keep the first aggregation public for diagnostics/tests.
        while (levels_.back().A.n_rows() > 2) {
            std::vector<std::size_t> aggregate;
            std::size_t coarse_n = 0;
            build_pair_aggregation(levels_.back(), aggregate, coarse_n);
            if (coarse_n >= levels_.back().A.n_rows() || coarse_n == 0) {
                break;
            }

            SparseMatrix coarse = galerkin_coarse(levels_.back().A, aggregate, coarse_n);
            if (!matrix_is_valid(coarse)) {
                levels_.clear();
                first_aggregate_.clear();
                return false;
            }

            if (levels_.size() == 1) {
                first_aggregate_ = aggregate;
            }

            Level next;
            next.A = std::move(coarse);
            if (!build_diagonal(next)) {
                levels_.clear();
                first_aggregate_.clear();
                return false;
            }
            levels_.back().aggregate = std::move(aggregate);
            levels_.push_back(std::move(next));
        }

        if (levels_.empty() || levels_.front().A.n_rows() == 0) {
            levels_.clear();
            return false;
        }
        return true;
    }

    std::size_t coarse_size() const noexcept
    {
        if (levels_.size() < 2) return levels_.empty() ? 0 : levels_.front().A.n_rows();
        return levels_[1].A.n_rows();
    }

    std::size_t aggregate_of(std::size_t fine_cell) const
    {
        if (fine_cell >= first_aggregate_.size()) {
            throw std::out_of_range("aggregate_of: fine-cell index out of range");
        }
        return first_aggregate_[fine_cell];
    }

    bool apply(const Vector& r, Vector& z) const override
    {
        if (levels_.empty() || r.size() != op_.rows() || !r.is_valid()) {
            return false;
        }
        if (z.size() != r.size()) z.resize(r.size());
        z.fill(0.0);
        return vcycle(0, r, z);
    }

    const char* name() const override { return "galerkin-agglomerated-amg"; }

private:
    struct Level {
        SparseMatrix A;
        std::vector<double> inv_diag;
        std::vector<std::size_t> aggregate;
    };

    static bool matrix_is_valid(const SparseMatrix& A)
    {
        if (A.n_rows() != A.n_cols() || A.n_rows() == 0 || !A.is_consistent()) {
            return false;
        }
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            if (row[i] > row[i + 1]) return false;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                if (col[k] >= A.n_cols() || !std::isfinite(val[k])) return false;
                if (k > row[i] && col[k] < col[k - 1]) return false;
            }
        }
        return true;
    }

    static bool build_diagonal(Level& level)
    {
        const std::size_t n = level.A.n_rows();
        const auto* row = level.A.row_offsets_data();
        const auto* col = level.A.columns_data();
        const auto* val = level.A.values_data();
        level.inv_diag.assign(n, 0.0);

        for (std::size_t i = 0; i < n; ++i) {
            double diagonal = 0.0;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                if (col[k] == i) diagonal += val[k];
            }
            if (!std::isfinite(diagonal) || std::abs(diagonal) <= 1e-30) {
                return false;
            }
            level.inv_diag[i] = 1.0 / diagonal;
            if (!std::isfinite(level.inv_diag[i])) return false;
        }
        return true;
    }

    static void build_pair_aggregation(const Level& level,
                                       std::vector<std::size_t>& aggregate,
                                       std::size_t& coarse_n)
    {
        const std::size_t n = level.A.n_rows();
        const auto* row = level.A.row_offsets_data();
        const auto* col = level.A.columns_data();
        const auto* val = level.A.values_data();

        aggregate.assign(n, 0);
        std::vector<bool> matched(n, false);
        coarse_n = 0;

        for (std::size_t i = 0; i < n; ++i) {
            if (matched[i]) continue;

            std::size_t best = n;
            double best_strength = -1.0;
            const double di = 1.0 / level.inv_diag[i];

            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                const std::size_t j = col[k];
                if (j == i || j >= n || matched[j]) continue;
                const double dj = 1.0 / level.inv_diag[j];
                const double denom = std::sqrt(std::max(std::abs(di * dj), 1e-60));
                const double strength = std::abs(val[k]) / denom;
                if (std::isfinite(strength) && strength > best_strength) {
                    best_strength = strength;
                    best = j;
                }
            }

            aggregate[i] = coarse_n;
            matched[i] = true;
            if (best < n) {
                aggregate[best] = coarse_n;
                matched[best] = true;
            }
            ++coarse_n;
        }
    }

    static SparseMatrix galerkin_coarse(const SparseMatrix& A,
                                        const std::vector<std::size_t>& aggregate,
                                        std::size_t coarse_n)
    {
        std::vector<std::map<std::size_t, double>> rows(coarse_n);
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();

        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            const std::size_t ci = aggregate[i];
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                const std::size_t cj = aggregate[col[k]];
                rows[ci][cj] += val[k];
            }
        }

        SparseMatrix coarse(coarse_n, coarse_n);
        for (std::size_t i = 0; i < coarse_n; ++i) {
            for (const auto& [j, value] : rows[i]) {
                if (std::isfinite(value) && std::abs(value) > 1e-30) {
                    coarse.push_back(i, j, value);
                }
            }
        }
        coarse.finalize();
        return coarse;
    }

    bool apply_operator(std::size_t level, const Vector& x, Vector& y) const
    {
        if (level == 0) {
            op_.apply(x, y);
            return y.is_valid();
        }
        const auto& A = levels_[level].A;
        const auto result = A.matvec(x);
        if (y.size() != result.size()) y.resize(result.size());
        for (std::size_t i = 0; i < result.size(); ++i) y(i) = result[i];
        return y.is_valid();
    }

    bool smooth(std::size_t level, const Vector& r, Vector& x, std::size_t sweeps) const
    {
        Vector Ax(r.size());
        const auto& inv_diag = levels_[level].inv_diag;
        for (std::size_t s = 0; s < sweeps; ++s) {
            if (!apply_operator(level, x, Ax)) return false;
            for (std::size_t i = 0; i < r.size(); ++i) {
                x(i) += omega_ * inv_diag[i] * (r(i) - Ax(i));
                if (!std::isfinite(x(i))) return false;
            }
        }
        return true;
    }

    bool vcycle(std::size_t level, const Vector& r, Vector& x) const
    {
        const auto& current = levels_[level];
        if (level + 1 == levels_.size()) {
            return smooth_coarsest(level, r, x);
        }

        if (!smooth(level, r, x, pre_)) return false;

        Vector Ax(r.size());
        if (!apply_operator(level, x, Ax)) return false;

        const auto& aggregate = current.aggregate;
        const std::size_t nc = levels_[level + 1].A.n_rows();
        Vector coarse_r(nc, 0.0);
        for (std::size_t i = 0; i < r.size(); ++i) {
            coarse_r(aggregate[i]) += r(i) - Ax(i);
        }

        Vector coarse_x(nc, 0.0);
        if (!vcycle(level + 1, coarse_r, coarse_x)) return false;

        for (std::size_t i = 0; i < r.size(); ++i) {
            x(i) += coarse_x(aggregate[i]);
            if (!std::isfinite(x(i))) return false;
        }

        return smooth(level, r, x, post_);
    }

    bool smooth_coarsest(std::size_t level, const Vector& r, Vector& x) const
    {
        // The hierarchy is deliberately kept sparse. On the coarsest level,
        // use more damped Jacobi sweeps to provide a stable small-system solve.
        return smooth(level, r, x, 24);
    }

    const LinearOperatorBase& op_;
    double omega_;
    std::size_t pre_, post_;
    std::vector<Level> levels_;
    std::vector<std::size_t> first_aggregate_;
};

using AgglomeratedAMGPreconditioner = MatrixFreeVcyclePreconditioner;
} // namespace cfdx::core
