#pragma once
#include "cfdx/core/linalg/linear_operator.h"
#include "cfdx/core/linalg/preconditioner.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::core {

// Multilevel AMG with a serial C/F splitting, normalized direct interpolation
// and Galerkin coarse operators. The setup follows the transferable parts of
// BoomerAMG's setup pipeline while deliberately omitting its distributed HMIS
// and extended+i machinery. The finest operator remains matrix-free during
// apply(); SparseMatrix is used during setup to build the hierarchy.
class MatrixFreeVcyclePreconditioner final : public Preconditioner {
public:
    MatrixFreeVcyclePreconditioner(const LinearOperatorBase& op,
                                   double omega = 0.7,
                                   std::size_t pre = 4,
                                   std::size_t post = 4,
                                   double strength_threshold = 0.25,
                                   std::size_t max_levels = 25)
        : op_(op), omega_(omega), pre_(pre), post_(post),
          strength_threshold_(strength_threshold), max_levels_(max_levels)
    {
        if (op.rows() != op.cols() || !(omega_ > 0.0 && omega_ < 2.0) ||
            !std::isfinite(strength_threshold_) || strength_threshold_ < 0.0 ||
            strength_threshold_ > 1.0 || max_levels_ == 0) {
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
        while (levels_.back().A.n_rows() > 2 && levels_.size() < max_levels_) {
            std::vector<TransferRow> prolongation;
            std::vector<std::size_t> aggregate;
            std::size_t coarse_n = 0;
            build_cf_interpolation(levels_.back(), strength_threshold_,
                                   prolongation, aggregate, coarse_n);
            if (coarse_n >= levels_.back().A.n_rows() || coarse_n == 0) {
                break;
            }

            SparseMatrix coarse = galerkin_coarse(
                levels_.back().A, prolongation, coarse_n);
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
            levels_.back().prolongation = std::move(prolongation);
            levels_.push_back(std::move(next));
        }

        if (levels_.empty() || levels_.front().A.n_rows() == 0) {
            levels_.clear();
            return false;
        }
        return true;
    }

    bool update_values(const SparseMatrix& A) override
    {
        if (levels_.empty() || A.n_rows() != op_.rows() ||
            A.n_cols() != op_.cols() || !same_pattern(levels_.front().A, A) ||
            !matrix_is_valid(A)) {
            return false;
        }

        std::vector<Level> refreshed = levels_;
        refreshed.front().A = A;
        if (!build_diagonal(refreshed.front())) return false;

        for (std::size_t level = 0; level + 1 < refreshed.size(); ++level) {
            const std::size_t coarse_n = refreshed[level + 1].A.n_rows();
            SparseMatrix coarse = galerkin_coarse(
                refreshed[level].A, refreshed[level].prolongation, coarse_n);
            if (!matrix_is_valid(coarse)) return false;
            refreshed[level + 1].A = std::move(coarse);
            if (!build_diagonal(refreshed[level + 1])) return false;
        }

        levels_ = std::move(refreshed);
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
    using TransferRow = std::vector<std::pair<std::size_t, double>>;

    struct Level {
        SparseMatrix A;
        std::vector<double> inv_diag;
        std::vector<std::size_t> aggregate;
        std::vector<TransferRow> prolongation;
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

    static bool same_pattern(const SparseMatrix& a, const SparseMatrix& b)
    {
        if (a.n_rows() != b.n_rows() || a.n_cols() != b.n_cols() ||
            a.nnz() != b.nnz() || !a.is_consistent() || !b.is_consistent()) {
            return false;
        }
        return std::equal(a.row_offsets_data(),
                          a.row_offsets_data() + a.n_rows() + 1,
                          b.row_offsets_data()) &&
               std::equal(a.columns_data(), a.columns_data() + a.nnz(),
                          b.columns_data());
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

    static std::vector<std::vector<std::size_t>> build_strength_graph(
        const Level& level, double strength_threshold)
    {
        const std::size_t n = level.A.n_rows();
        const auto* row = level.A.row_offsets_data();
        const auto* col = level.A.columns_data();
        const auto* val = level.A.values_data();
        std::vector<std::vector<std::size_t>> strong(n);

        // HYPRE's standard strength graph compares entries having the sign
        // opposite to the diagonal against the largest such entry in the row.
        // Its default max-row-sum guard (0.9) suppresses coarsening for rows
        // that are already strongly diagonally dominant.
        constexpr double max_row_sum = 0.9;
        for (std::size_t i = 0; i < n; ++i) {
            const double di = 1.0 / level.inv_diag[i];
            double row_sum = di;
            double opposite_max = 0.0;
            double absolute_max = 0.0;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                if (col[k] == i) continue;
                row_sum += val[k];
                absolute_max = std::max(absolute_max, std::abs(val[k]));
                if (di * val[k] < 0.0)
                    opposite_max = std::max(opposite_max, std::abs(val[k]));
            }
            if (std::abs(row_sum) > std::abs(di) * max_row_sum) continue;

            const bool use_opposite_sign = opposite_max > 0.0;
            const double scale = use_opposite_sign ? opposite_max : absolute_max;
            if (!(scale > 0.0) || !std::isfinite(scale)) continue;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                const std::size_t j = col[k];
                if (j == i || j >= n) continue;
                if (use_opposite_sign && di * val[k] >= 0.0) continue;
                if (std::abs(val[k]) + std::numeric_limits<double>::epsilon() <
                    strength_threshold * scale) continue;
                strong[i].push_back(j);
            }
            std::sort(strong[i].begin(), strong[i].end());
            strong[i].erase(std::unique(strong[i].begin(), strong[i].end()),
                            strong[i].end());
        }
        return strong;
    }

    static void build_cf_interpolation(const Level& level,
                                       double strength_threshold,
                                       std::vector<TransferRow>& prolongation,
                                       std::vector<std::size_t>& aggregate,
                                       std::size_t& coarse_n)
    {
        const std::size_t n = level.A.n_rows();
        const auto strong = build_strength_graph(level, strength_threshold);
        std::vector<std::vector<std::size_t>> neighborhood = strong;
        for (std::size_t i = 0; i < n; ++i) {
            for (const std::size_t j : strong[i]) neighborhood[j].push_back(i);
        }
        for (auto& neighbors : neighborhood) {
            std::sort(neighbors.begin(), neighbors.end());
            neighbors.erase(std::unique(neighbors.begin(), neighbors.end()),
                            neighbors.end());
        }

        // Deterministic serial maximal-independent-set splitting. This is the
        // local analogue of the independent-set completion used by HMIS: pick
        // the most connected undecided point as C and mark its graph neighbors F.
        enum class Point : unsigned char { Undecided, Fine, Coarse };
        std::vector<Point> point(n, Point::Undecided);
        std::size_t undecided = n;
        while (undecided > 0) {
            std::size_t best = n;
            std::size_t best_measure = 0;
            for (std::size_t i = 0; i < n; ++i) {
                if (point[i] != Point::Undecided) continue;
                std::size_t measure = 0;
                for (const std::size_t j : neighborhood[i])
                    if (point[j] == Point::Undecided) ++measure;
                if (best == n || measure > best_measure ||
                    (measure == best_measure && i < best)) {
                    best = i;
                    best_measure = measure;
                }
            }
            point[best] = Point::Coarse;
            --undecided;
            for (const std::size_t j : neighborhood[best]) {
                if (point[j] == Point::Undecided) {
                    point[j] = Point::Fine;
                    --undecided;
                }
            }
        }

        // Every F point needs a direct C interpolatory neighbor. Promote any
        // exceptional isolated point instead of manufacturing a hidden fallback.
        for (std::size_t i = 0; i < n; ++i) {
            if (point[i] != Point::Fine) continue;
            const bool has_coarse = std::any_of(
                neighborhood[i].begin(), neighborhood[i].end(),
                [&](std::size_t j) { return point[j] == Point::Coarse; });
            if (!has_coarse) point[i] = Point::Coarse;
        }

        std::vector<std::size_t> coarse_index(n, n);
        coarse_n = 0;
        for (std::size_t i = 0; i < n; ++i)
            if (point[i] == Point::Coarse) coarse_index[i] = coarse_n++;

        const auto* row = level.A.row_offsets_data();
        const auto* col = level.A.columns_data();
        const auto* val = level.A.values_data();
        prolongation.assign(n, {});
        aggregate.assign(n, 0);
        for (std::size_t i = 0; i < n; ++i) {
            if (point[i] == Point::Coarse) {
                prolongation[i].push_back({coarse_index[i], 1.0});
                aggregate[i] = coarse_index[i];
                continue;
            }

            std::vector<std::pair<std::size_t, double>> raw;
            double raw_sum = 0.0;
            double magnitude_sum = 0.0;
            const double diagonal = 1.0 / level.inv_diag[i];
            for (const std::size_t j : neighborhood[i]) {
                if (point[j] != Point::Coarse) continue;
                double aij = 0.0;
                for (std::size_t k = row[i]; k < row[i + 1]; ++k)
                    if (col[k] == j) aij += val[k];
                const double weight = -aij / diagonal;
                raw.push_back({coarse_index[j], weight});
                raw_sum += weight;
                magnitude_sum += std::abs(weight);
            }

            if (std::abs(raw_sum) <= 1e-30 || !std::isfinite(raw_sum)) {
                raw_sum = magnitude_sum;
                for (auto& entry : raw) entry.second = std::abs(entry.second);
            }
            if (!(std::abs(raw_sum) > 1e-30) || !std::isfinite(raw_sum)) {
                const double uniform = 1.0 / static_cast<double>(raw.size());
                for (auto& entry : raw) entry.second = uniform;
            } else {
                for (auto& entry : raw) entry.second /= raw_sum;
            }

            double largest = -1.0;
            for (const auto& entry : raw) {
                prolongation[i].push_back(entry);
                if (std::abs(entry.second) >= largest) {
                    largest = std::abs(entry.second);
                    aggregate[i] = entry.first;
                }
            }
        }
    }

    static SparseMatrix galerkin_coarse(const SparseMatrix& A,
                                        const std::vector<TransferRow>& prolongation,
                                        std::size_t coarse_n)
    {
        std::vector<std::map<std::size_t, double>> rows(coarse_n);
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();

        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                const std::size_t j = col[k];
                for (const auto& [ci, wi] : prolongation[i])
                    for (const auto& [cj, wj] : prolongation[j])
                        rows[ci][cj] += wi * val[k] * wj;
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

        const auto& prolongation = current.prolongation;
        const std::size_t nc = levels_[level + 1].A.n_rows();
        Vector coarse_r(nc, 0.0);
        for (std::size_t i = 0; i < r.size(); ++i) {
            const double residual = r(i) - Ax(i);
            for (const auto& [coarse, weight] : prolongation[i])
                coarse_r(coarse) += weight * residual;
        }

        Vector coarse_x(nc, 0.0);
        if (!vcycle(level + 1, coarse_r, coarse_x)) return false;

        for (std::size_t i = 0; i < r.size(); ++i) {
            for (const auto& [coarse, weight] : prolongation[i])
                x(i) += weight * coarse_x(coarse);
            if (!std::isfinite(x(i))) return false;
        }

        return smooth(level, r, x, post_);
    }

    bool smooth_coarsest(std::size_t level, const Vector& r, Vector& x) const
    {
        // Solve the tiny coarsest problem directly. Using a fixed number of
        // Jacobi sweeps leaves low-frequency error on small Poisson systems
        // and makes the V-cycle quality depend on the arbitrary sweep count.
        // The coarsest level is intentionally small, so a dense pivoted solve
        // is both robust and negligible compared with the fine-grid work.
        const auto& A = levels_[level].A;
        const std::size_t n = A.n_rows();
        if (r.size() != n || n == 0) return false;

        std::vector<double> m(n * n, 0.0);
        std::vector<double> b(n, 0.0);
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        for (std::size_t i = 0; i < n; ++i) {
            b[i] = r(i);
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                m[i * n + col[k]] += val[k];
            }
        }

        constexpr double pivot_tol = 1e-14;
        for (std::size_t k = 0; k < n; ++k) {
            std::size_t pivot = k;
            double pivot_abs = std::abs(m[k * n + k]);
            for (std::size_t i = k + 1; i < n; ++i) {
                const double candidate = std::abs(m[i * n + k]);
                if (candidate > pivot_abs) {
                    pivot_abs = candidate;
                    pivot = i;
                }
            }
            if (!std::isfinite(pivot_abs) || pivot_abs <= pivot_tol) return false;

            if (pivot != k) {
                for (std::size_t j = k; j < n; ++j) {
                    std::swap(m[k * n + j], m[pivot * n + j]);
                }
                std::swap(b[k], b[pivot]);
            }

            const double diagonal = m[k * n + k];
            for (std::size_t i = k + 1; i < n; ++i) {
                const double factor = m[i * n + k] / diagonal;
                if (!std::isfinite(factor)) return false;
                m[i * n + k] = 0.0;
                for (std::size_t j = k + 1; j < n; ++j) {
                    m[i * n + j] -= factor * m[k * n + j];
                }
                b[i] -= factor * b[k];
            }
        }

        if (x.size() != n) x.resize(n);
        for (std::size_t ii = n; ii-- > 0;) {
            double sum = b[ii];
            for (std::size_t j = ii + 1; j < n; ++j) {
                sum -= m[ii * n + j] * x(j);
            }
            const double diagonal = m[ii * n + ii];
            if (!std::isfinite(diagonal) || std::abs(diagonal) <= pivot_tol) return false;
            x(ii) = sum / diagonal;
            if (!std::isfinite(x(ii))) return false;
        }
        return x.is_valid();
    }

    const LinearOperatorBase& op_;
    double omega_;
    std::size_t pre_, post_;
    double strength_threshold_;
    std::size_t max_levels_;
    std::vector<Level> levels_;
    std::vector<std::size_t> first_aggregate_;
};

using AgglomeratedAMGPreconditioner = MatrixFreeVcyclePreconditioner;
} // namespace cfdx::core
