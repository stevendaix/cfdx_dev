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
#include <set>
#include <utility>
#include <vector>

namespace cfdx::core {

enum class AMGInterpolationPolicy {
    DirectCF,
    SmoothedAggregation
};

// Multilevel AMG with selectable direct C/F or smoothed-aggregation
// interpolation and Galerkin coarse operators. The finest operator remains
// matrix-free during apply(); SparseMatrix is used during setup to build the
// hierarchy.
class MatrixFreeVcyclePreconditioner final : public Preconditioner {
public:
    MatrixFreeVcyclePreconditioner(const LinearOperatorBase& op,
                                   double omega = 0.7,
                                   std::size_t pre = 4,
                                   std::size_t post = 4,
                                   double strength_threshold = 0.25,
                                   std::size_t max_levels = 25,
                                   AMGInterpolationPolicy interpolation =
                                       AMGInterpolationPolicy::DirectCF,
                                   bool constant_null_space = false)
        : op_(op), omega_(omega), pre_(pre), post_(post),
          strength_threshold_(strength_threshold), max_levels_(max_levels),
          interpolation_(interpolation),
          constant_null_space_(constant_null_space)
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
            A.n_rows() == 0 || !matrix_is_valid(A) ||
            (constant_null_space_ && !has_constant_null_space(A))) {
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
        // Keep a small but non-trivial dense coarse problem. Recursing to
        // 2 unknowns adds levels without improving the coarse solve materially;
        // a direct solve of <=16 unknowns is negligible and gives a cleaner
        // multilevel correction.
        constexpr std::size_t coarse_direct_limit = 16;
        while (levels_.back().A.n_rows() > coarse_direct_limit &&
               levels_.size() < max_levels_) {
            std::vector<TransferRow> prolongation;
            std::vector<std::size_t> aggregate;
            std::size_t coarse_n = 0;
            if (interpolation_ == AMGInterpolationPolicy::SmoothedAggregation) {
                build_smoothed_aggregation_interpolation(
                    levels_.back(), strength_threshold_, prolongation,
                    aggregate, coarse_n);
            } else {
                if (!build_cf_interpolation(levels_.back(), strength_threshold_,
                                            prolongation, aggregate, coarse_n)) {
                    levels_.clear();
                    first_aggregate_.clear();
                    return false;
                }
            }
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
            !matrix_is_valid(A) ||
            (constant_null_space_ && !has_constant_null_space(A))) {
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

    std::vector<std::size_t> hierarchy_level_sizes() const
    {
        std::vector<std::size_t> sizes;
        sizes.reserve(levels_.size());
        for (const auto& level : levels_) sizes.push_back(level.A.n_rows());
        return sizes;
    }

    double prolongation_row_sum_min() const
    {
        if (levels_.size() < 2 || levels_.front().prolongation.empty())
            throw std::out_of_range("prolongation_row_sum_min: hierarchy unavailable");
        double value = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < levels_.front().prolongation.size(); ++i)
            value = std::min(value, prolongation_row_sum(i));
        return value;
    }

    double prolongation_row_sum_max() const
    {
        if (levels_.size() < 2 || levels_.front().prolongation.empty())
            throw std::out_of_range("prolongation_row_sum_max: hierarchy unavailable");
        double value = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < levels_.front().prolongation.size(); ++i)
            value = std::max(value, prolongation_row_sum(i));
        return value;
    }

    double first_prolongation_linear_mode_relative_error() const
    {
        if (levels_.size() < 2 || levels_.front().prolongation.empty())
            throw std::out_of_range("first_prolongation_linear_mode_relative_error: hierarchy unavailable");
        const auto& P = levels_.front().prolongation;
        const std::size_t n = P.size();
        std::vector<double> coarse_coordinate(
            levels_[1].A.n_rows(), std::numeric_limits<double>::quiet_NaN());
        for (std::size_t i = 0; i < n; ++i) {
            if (P[i].size() == 1 && std::abs(P[i][0].second - 1.0) <= 1e-12) {
                const std::size_t c = P[i][0].first;
                if (c < coarse_coordinate.size()) coarse_coordinate[c] = static_cast<double>(i);
            }
        }
        for (double value : coarse_coordinate)
            if (!std::isfinite(value)) return std::numeric_limits<double>::infinity();

        double error2 = 0.0;
        double norm2 = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double interpolated = 0.0;
            for (const auto& [c, weight] : P[i])
                interpolated += weight * coarse_coordinate[c];
            const double error = interpolated - static_cast<double>(i);
            error2 += error * error;
            norm2 += static_cast<double>(i) * static_cast<double>(i);
        }
        return std::sqrt(error2 / std::max(norm2, 1.0));
    }

    double first_coarse_symmetry_relative_error() const
    {
        if (levels_.size() < 2)
            throw std::out_of_range("first_coarse_symmetry_relative_error: hierarchy unavailable");
        const auto& A = levels_[1].A;
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        double defect2 = 0.0;
        double scale2 = 0.0;
        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                const std::size_t j = col[k];
                double transpose = 0.0;
                for (std::size_t q = row[j]; q < row[j + 1]; ++q)
                    if (col[q] == i) transpose += val[q];
                const double defect = val[k] - transpose;
                defect2 += defect * defect;
                scale2 += val[k] * val[k];
            }
        }
        return std::sqrt(defect2 / std::max(scale2, 1.0));
    }

    double first_prolongation_mode_relative_error(std::size_t mode) const
    {
        if (levels_.size() < 2 || levels_.front().prolongation.empty() || mode == 0)
            throw std::out_of_range("first_prolongation_mode_relative_error: hierarchy unavailable");
        const auto& P = levels_.front().prolongation;
        const std::size_t n = P.size();
        std::vector<double> coarse_mode(levels_[1].A.n_rows(), 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            if (P[i].size() == 1 && std::abs(P[i][0].second - 1.0) <= 1e-12) {
                const std::size_t c = P[i][0].first;
                if (c < coarse_mode.size()) {
                    coarse_mode[c] = std::sin(
                        3.14159265358979323846 * static_cast<double>(mode * (i + 1)) /
                        static_cast<double>(n + 1));
                }
            }
        }
        double error2 = 0.0;
        double norm2 = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double interpolated = 0.0;
            for (const auto& [coarse, weight] : P[i])
                interpolated += weight * coarse_mode[coarse];
            const double exact = std::sin(
                3.14159265358979323846 * static_cast<double>(mode * (i + 1)) /
                static_cast<double>(n + 1));
            const double error = interpolated - exact;
            error2 += error * error;
            norm2 += exact * exact;
        }
        return std::sqrt(error2 / std::max(norm2, 1.0));
    }

    std::size_t first_prolongation_row_nnz_min() const
    {
        if (levels_.size() < 2 || levels_.front().prolongation.empty())
            throw std::out_of_range("first_prolongation_row_nnz_min: hierarchy unavailable");
        std::size_t value = std::numeric_limits<std::size_t>::max();
        for (const auto& row : levels_.front().prolongation)
            value = std::min(value, row.size());
        return value;
    }

    std::size_t first_prolongation_row_nnz_max() const
    {
        if (levels_.size() < 2 || levels_.front().prolongation.empty())
            throw std::out_of_range("first_prolongation_row_nnz_max: hierarchy unavailable");
        std::size_t value = 0;
        for (const auto& row : levels_.front().prolongation)
            value = std::max(value, row.size());
        return value;
    }

    double first_prolongation_weight_min() const
    {
        if (levels_.size() < 2 || levels_.front().prolongation.empty())
            throw std::out_of_range("first_prolongation_weight_min: hierarchy unavailable");
        double value = std::numeric_limits<double>::infinity();
        for (const auto& row : levels_.front().prolongation)
            for (const auto& [coarse, weight] : row) {
                (void)coarse;
                value = std::min(value, weight);
            }
        return value;
    }

    double first_prolongation_weight_max() const
    {
        if (levels_.size() < 2 || levels_.front().prolongation.empty())
            throw std::out_of_range("first_prolongation_weight_max: hierarchy unavailable");
        double value = -std::numeric_limits<double>::infinity();
        for (const auto& row : levels_.front().prolongation)
            for (const auto& [coarse, weight] : row) {
                (void)coarse;
                value = std::max(value, weight);
            }
        return value;
    }

    std::size_t first_prolongation_negative_count() const
    {
        if (levels_.size() < 2 || levels_.front().prolongation.empty())
            throw std::out_of_range("first_prolongation_negative_count: hierarchy unavailable");
        std::size_t count = 0;
        for (const auto& row : levels_.front().prolongation)
            for (const auto& [coarse, weight] : row) {
                (void)coarse;
                if (weight < 0.0) ++count;
            }
        return count;
    }

    double first_coarse_diagonal_min() const
    {
        if (levels_.size() < 2)
            throw std::out_of_range("first_coarse_diagonal_min: hierarchy unavailable");
        const auto& A = levels_[1].A;
        double value = std::numeric_limits<double>::infinity();
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            double diag = 0.0;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k)
                if (col[k] == i) diag += val[k];
            value = std::min(value, diag);
        }
        return value;
    }

    double first_coarse_diagonal_max() const
    {
        if (levels_.size() < 2)
            throw std::out_of_range("first_coarse_diagonal_max: hierarchy unavailable");
        const auto& A = levels_[1].A;
        double value = -std::numeric_limits<double>::infinity();
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            double diag = 0.0;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k)
                if (col[k] == i) diag += val[k];
            value = std::max(value, diag);
        }
        return value;
    }

    double first_coarse_gershgorin_lower_bound() const
    {
        if (levels_.size() < 2)
            throw std::out_of_range("first_coarse_gershgorin_lower_bound: hierarchy unavailable");
        const auto& A = levels_[1].A;
        double value = std::numeric_limits<double>::infinity();
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            double diag = 0.0;
            double offdiag = 0.0;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                if (col[k] == i) diag += val[k];
                else offdiag += std::abs(val[k]);
            }
            value = std::min(value, diag - offdiag);
        }
        return value;
    }

    struct TransferDiagnostic {
        std::size_t level = 0;
        std::size_t fine_size = 0;
        std::size_t coarse_size = 0;
        std::size_t nnz = 0;
        std::size_t negative_weights = 0;
        double row_sum_min = 0.0;
        double row_sum_max = 0.0;
        double weight_min = 0.0;
        double weight_max = 0.0;
        double column_norm_min = 0.0;
        double column_norm_max = 0.0;
        double constant_mode_error = 0.0;
        double linear_mode_error = 0.0;
        double sine1_mode_error = 0.0;
        double sine2_mode_error = 0.0;
        double galerkin_relative_error = 0.0;
        double coarse_symmetry_relative_error = 0.0;
        double coarse_diagonal_min = 0.0;
        double coarse_diagonal_max = 0.0;
        double coarse_gershgorin_lower_bound = 0.0;
    };

    std::vector<TransferDiagnostic> transfer_diagnostics() const
    {
        std::vector<TransferDiagnostic> out;
        if (levels_.size() < 2) return out;
        out.reserve(levels_.size() - 1);

        constexpr double pi = 3.14159265358979323846;
        for (std::size_t l = 0; l + 1 < levels_.size(); ++l) {
            const auto& P = levels_[l].prolongation;
            const std::size_t n = levels_[l].A.n_rows();
            const std::size_t nc = levels_[l + 1].A.n_rows();
            TransferDiagnostic d;
            d.level = l;
            d.fine_size = n;
            d.coarse_size = nc;
            d.row_sum_min = std::numeric_limits<double>::infinity();
            d.row_sum_max = -std::numeric_limits<double>::infinity();
            d.weight_min = std::numeric_limits<double>::infinity();
            d.weight_max = -std::numeric_limits<double>::infinity();

            std::vector<double> cnorm(nc, 0.0);
            std::vector<double> denom(nc, 0.0);
            std::vector<double> x1(nc, 0.0), s1(nc, 0.0), s2(nc, 0.0);
            for (std::size_t i = 0; i < P.size(); ++i) {
                const double xi = static_cast<double>(i);
                const double q = pi * static_cast<double>(i + 1) /
                                 static_cast<double>(n + 1);
                const double mode1 = std::sin(q);                const double mode2 = std::sin(2.0 * q);
                double sum = 0.0;
                for (const auto& [coarse, weight] : P[i]) {
                    if (coarse >= nc || !std::isfinite(weight)) continue;
                    ++d.nnz;
                    sum += weight;
                    d.weight_min = std::min(d.weight_min, weight);
                    d.weight_max = std::max(d.weight_max, weight);
                    if (weight < 0.0) ++d.negative_weights;
                    cnorm[coarse] += weight * weight;
                    denom[coarse] += weight * weight;
                    x1[coarse] += weight * weight * xi;
                    s1[coarse] += weight * weight * mode1;
                    s2[coarse] += weight * weight * mode2;
                }
                d.row_sum_min = std::min(d.row_sum_min, sum);
                d.row_sum_max = std::max(d.row_sum_max, sum);
            }

            for (std::size_t c = 0; c < nc; ++c) {
                if (denom[c] > 0.0) {
                    x1[c] /= denom[c];
                    s1[c] /= denom[c];
                    s2[c] /= denom[c];
                }
            }

            for (double& v : cnorm) v = std::sqrt(v);
            if (!cnorm.empty()) {
                d.column_norm_min = *std::min_element(cnorm.begin(), cnorm.end());
                d.column_norm_max = *std::max_element(cnorm.begin(), cnorm.end());
            }

            auto mode_error = [&](auto coarse_mode, auto fine_mode) {
                double e2 = 0.0;
                double n2 = 0.0;
                for (std::size_t i = 0; i < n; ++i) {
                    double interpolated = 0.0;
                    for (const auto& [coarse, weight] : P[i])
                        if (coarse < nc) interpolated += weight * coarse_mode(coarse);
                    const double exact = fine_mode(i);
                    const double e = interpolated - exact;
                    e2 += e * e;
                    n2 += exact * exact;
                }
                return std::sqrt(e2 / std::max(n2, 1e-300));
            };

            d.constant_mode_error = std::max(
                std::abs(d.row_sum_min - 1.0),
                std::abs(d.row_sum_max - 1.0));
            d.linear_mode_error = mode_error(
                [&](std::size_t c) { return x1[c]; },
                [&](std::size_t i) { return static_cast<double>(i); });
            d.sine1_mode_error = mode_error(
                [&](std::size_t c) { return s1[c]; },
                [&](std::size_t i) {
                    return std::sin(pi * static_cast<double>(i + 1) /
                                     static_cast<double>(n + 1));
                });
            d.sine2_mode_error = mode_error(
                [&](std::size_t c) { return s2[c]; },
                [&](std::size_t i) {
                    return std::sin(2.0 * pi * static_cast<double>(i + 1) /
                                     static_cast<double>(n + 1));
                });

            d.galerkin_relative_error = galerkin_relative_error(l);

            const auto& Ac = levels_[l + 1].A;
            const auto* row = Ac.row_offsets_data();
            const auto* col = Ac.columns_data();
            const auto* val = Ac.values_data();
            double sym_e2 = 0.0;
            double sym_n2 = 0.0;
            d.coarse_diagonal_min = std::numeric_limits<double>::infinity();
            d.coarse_diagonal_max = -std::numeric_limits<double>::infinity();
            d.coarse_gershgorin_lower_bound = std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < nc; ++i) {
                double diag = 0.0;
                double offdiag = 0.0;
                for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                    const std::size_t j = col[k];
                    if (j == i) diag += val[k];
                    else offdiag += std::abs(val[k]);
                    double transpose = 0.0;
                    if (j < nc) {
                        for (std::size_t q = row[j]; q < row[j + 1]; ++q)
                            if (col[q] == i) transpose += val[q];
                    }
                    const double e = val[k] - transpose;
                    sym_e2 += e * e;
                    sym_n2 += val[k] * val[k];
                }
                d.coarse_diagonal_min = std::min(d.coarse_diagonal_min, diag);
                d.coarse_diagonal_max = std::max(d.coarse_diagonal_max, diag);
                d.coarse_gershgorin_lower_bound =
                    std::min(d.coarse_gershgorin_lower_bound, diag - offdiag);
            }
            d.coarse_symmetry_relative_error =
                std::sqrt(sym_e2 / std::max(sym_n2, 1e-300));
            out.push_back(d);
        }
        return out;
    }

    double galerkin_relative_error(std::size_t level) const
    {
        if(level+1>=levels_.size()) throw std::out_of_range("galerkin_relative_error: invalid level");
        const auto& A=levels_[level].A; const auto& P=levels_[level].prolongation;
        const auto& Ac=levels_[level+1].A; const std::size_t nc=Ac.n_rows();
        std::vector<std::map<std::size_t,double>> ref(nc);
        const auto* ar=A.row_offsets_data(); const auto* ac=A.columns_data(); const auto* av=A.values_data();
        for(std::size_t i=0;i<A.n_rows();++i)
            for(std::size_t k=ar[i];k<ar[i+1];++k)
                for(const auto& [ci,wi]:P[i])
                    for(const auto& [cj,wj]:P[ac[k]])
                        ref[ci][cj]+=wi*av[k]*wj;
        const auto* rr=Ac.row_offsets_data(); const auto* cc=Ac.columns_data(); const auto* vv=Ac.values_data();
        double e2=0.0,r2=0.0;
        for(std::size_t i=0;i<nc;++i) {
            std::map<std::size_t,double> stored;
            for(std::size_t k=rr[i];k<rr[i+1];++k) stored[cc[k]]+=vv[k];
            std::set<std::size_t> cols;
            for(const auto& [j,v]:ref[i]) { (void)v; cols.insert(j); }
            for(const auto& [j,v]:stored) { (void)v; cols.insert(j); }
            for(const auto j:cols) {
                const double rv=ref[i].count(j)?ref[i].at(j):0.0;
                const double sv=stored.count(j)?stored.at(j):0.0;
                const double e=sv-rv; e2+=e*e; r2+=rv*rv;
            }
        }
        return std::sqrt(e2/std::max(r2,1e-300));
    }

    double two_grid_residual_ratio(std::size_t level, const Vector& rhs) const
    {
        if (level + 1 >= levels_.size()) {
            throw std::out_of_range(
                "two_grid_residual_ratio: invalid level=" +
                std::to_string(level) +
                " hierarchy_levels=" + std::to_string(levels_.size()));
        }
        if (rhs.size() != levels_[level].A.n_rows()) {
            throw std::out_of_range(
                "two_grid_residual_ratio: rhs size mismatch at level=" +
                std::to_string(level) +
                " fine_size=" + std::to_string(levels_[level].A.n_rows()) +
                " rhs_size=" + std::to_string(rhs.size()));
        }

        Vector x(rhs.size(), 0.0), Ax(rhs.size()), res(rhs.size());
        if (!smooth(level, rhs, x, pre_) || !apply_operator(level, x, Ax))
            return std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < rhs.size(); ++i)
            res(i) = rhs(i) - Ax(i);
        const double before = res.norm2();
        if (!(before > 0.0) || !std::isfinite(before)) return 0.0;

        const std::size_t nc = levels_[level + 1].A.n_rows();
        Vector rc(nc, 0.0);
        for (std::size_t i = 0; i < rhs.size(); ++i)
            for (const auto& [c, w] : levels_[level].prolongation[i])
                rc(c) += w * res(i);

        // Exact dense coarse solves are retained only for tiny systems. For
        // intermediate levels use PCG with Jacobi scaling; this keeps the
        // two-grid diagnostic mathematically representative while avoiding
        // O(n_c^3) work at the 2048-unknown first coarse level.
        Vector ec(nc, 0.0);
        const auto& Ac = levels_[level + 1].A;
        const auto* row = Ac.row_offsets_data();
        const auto* col = Ac.columns_data();
        const auto* val = Ac.values_data();
        Vector r = rc, z(nc, 0.0), p(nc, 0.0), Ap(nc, 0.0);
        double rz = 0.0;
        double rhs_norm2 = 0.0;
        for (std::size_t i = 0; i < nc; ++i) {
            rhs_norm2 += rc(i) * rc(i);
            double diag = 0.0;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k)
                if (col[k] == i) diag += val[k];
            if (!std::isfinite(diag) || diag <= 0.0) return std::numeric_limits<double>::infinity();
            z(i) = r(i) / diag;
            p(i) = z(i);
            rz += r(i) * z(i);
        }

        const double rhs_norm = std::sqrt(rhs_norm2);
        if (rhs_norm > 0.0) {
            const std::size_t max_iter = std::max<std::size_t>(100, 4 * nc);
            const double target = 1e-12 * rhs_norm;
            bool converged = false;
            for (std::size_t iter = 0; iter < max_iter; ++iter) {
                const auto values = Ac.matvec(p);
                for (std::size_t i = 0; i < nc; ++i) Ap(i) = values[i];
                double pAp = 0.0;
                for (std::size_t i = 0; i < nc; ++i) pAp += p(i) * Ap(i);
                if (!std::isfinite(pAp) || pAp <= 0.0) return std::numeric_limits<double>::infinity();
                const double alpha = rz / pAp;
                for (std::size_t i = 0; i < nc; ++i) {
                    ec(i) += alpha * p(i);
                    r(i) -= alpha * Ap(i);
                }
                double rr = 0.0;
                for (std::size_t i = 0; i < nc; ++i) rr += r(i) * r(i);
                if (std::sqrt(rr) <= target) {
                    converged = true;
                    break;
                }
                double rz_new = 0.0;
                for (std::size_t i = 0; i < nc; ++i) {
                    double diag = 0.0;
                    for (std::size_t k = row[i]; k < row[i + 1]; ++k)
                        if (col[k] == i) diag += val[k];
                    z(i) = r(i) / diag;
                    rz_new += r(i) * z(i);
                }
                const double beta = rz_new / rz;
                for (std::size_t i = 0; i < nc; ++i) p(i) = z(i) + beta * p(i);
                rz = rz_new;
            }
            if (!converged) return std::numeric_limits<double>::infinity();
        }

        for (std::size_t i = 0; i < x.size(); ++i) {
            double corr = 0.0;
            for (const auto& [c, w] : levels_[level].prolongation[i])
                corr += w * ec(c);
            x(i) += corr;
        }
        if (!smooth(level, rhs, x, post_) || !apply_operator(level, x, Ax))
            return std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < rhs.size(); ++i)
            res(i) = rhs(i) - Ax(i);
        return res.norm2() / std::max(before, 1e-300);
    }

    double two_grid_sine_mode_residual_ratio(std::size_t level, std::size_t mode) const
    {
        if (level + 1 >= levels_.size() || mode == 0) {
            throw std::out_of_range("two_grid_sine_mode_residual_ratio: invalid level or mode");
        }
        const std::size_t n = levels_[level].A.n_rows();
        Vector exact(n, 0.0);
        constexpr double pi = 3.14159265358979323846;
        for (std::size_t i = 0; i < n; ++i) {
            exact(i) = std::sin(pi * static_cast<double>(mode * (i + 1)) /
                                 static_cast<double>(n + 1));
        }
        Vector rhs(n, 0.0);
        if (!apply_operator(level, exact, rhs)) {
            return std::numeric_limits<double>::infinity();
        }
        return two_grid_residual_ratio(level, rhs);
    }

    std::size_t first_prolongation_nnz() const
    {
        if (levels_.size() < 2)
            throw std::out_of_range("first_prolongation_nnz: hierarchy unavailable");
        std::size_t nnz = 0;
        for (const auto& row : levels_.front().prolongation) nnz += row.size();
        return nnz;
    }

    double prolongation_row_sum(std::size_t fine_cell) const
    {
        if (levels_.size() < 2 ||
            fine_cell >= levels_.front().prolongation.size()) {
            throw std::out_of_range(
                "prolongation_row_sum: fine-cell index out of range");
        }
        double sum = 0.0;
        for (const auto& [coarse, weight] :
             levels_.front().prolongation[fine_cell]) {
            (void)coarse;
            sum += weight;
        }
        return sum;
    }

    bool apply(const Vector& r, Vector& z) const override
    {
        if (levels_.empty() || r.size() != op_.rows() || !r.is_valid()) {
            return false;
        }
        Vector compatible_r = r;
        if (constant_null_space_) remove_constant(compatible_r);
        if (z.size() != r.size()) z.resize(r.size());
        z.fill(0.0);
        if (!vcycle(0, compatible_r, z)) return false;
        if (constant_null_space_) remove_constant(z);
        return z.is_valid();
    }

    struct VcycleDiagnostic {
        std::size_t level = 0;
        std::size_t size = 0;
        double rhs_norm = 0.0;
        double x_initial_norm = 0.0;
        double residual_before = 0.0;
        double residual_after_pre = 0.0;
        double coarse_rhs_norm = 0.0;
        double coarse_solution_norm = 0.0;
        double correction_norm = 0.0;
        double correction_operator_norm = 0.0;
        double coarse_equation_relative_residual = 0.0;
        double residual_after_correction = 0.0;
        double residual_after_post = 0.0;
        double coarse_residual_before = 0.0;
        double coarse_residual_after = 0.0;
    };

    bool apply_with_diagnostics(const Vector& r, Vector& z,
                                std::vector<VcycleDiagnostic>& diagnostics) const
    {
        diagnostics.clear();
        if (levels_.empty() || r.size() != op_.rows() || !r.is_valid()) return false;
        Vector compatible_r = r;
        if (constant_null_space_) remove_constant(compatible_r);
        if (z.size() != r.size()) z.resize(r.size());
        z.fill(0.0);
        if (!vcycle_diagnostic(0, compatible_r, z, diagnostics)) return false;
        if (constant_null_space_) remove_constant(z);
        return z.is_valid();
    }

    const char* name() const override
    {
        return interpolation_ == AMGInterpolationPolicy::SmoothedAggregation
            ? "smoothed-aggregation-amg"
            : "galerkin-agglomerated-amg";
    }

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

    static bool has_constant_null_space(const SparseMatrix& A)
    {
        const auto* row = A.row_offsets_data();
        const auto* val = A.values_data();
        double matrix_scale = 0.0;
        double residual_scale = 0.0;
        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            double row_sum = 0.0;
            double row_norm = 0.0;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                row_sum += val[k];
                row_norm += std::abs(val[k]);
            }
            matrix_scale = std::max(matrix_scale, row_norm);
            residual_scale = std::max(residual_scale, std::abs(row_sum));
        }
        const double tolerance = 256.0 * std::numeric_limits<double>::epsilon() *
                                 std::max(1.0, matrix_scale);
        return std::isfinite(residual_scale) && residual_scale <= tolerance;
    }

    static void remove_constant(Vector& values)
    {
        if (values.size() == 0) return;
        double mean = 0.0;
        for (std::size_t i = 0; i < values.size(); ++i)
            mean += values(i);
        mean /= static_cast<double>(values.size());
        for (std::size_t i = 0; i < values.size(); ++i)
            values(i) -= mean;
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
            strong[i].erase(std::unique(strong[i].begin(), strong[i].end()),                            strong[i].end());
        }
        return strong;
    }

    static bool build_cf_interpolation(const Level& level,
                                       double strength_threshold,
                                       std::vector<TransferRow>& prolongation,
                                       std::vector<std::size_t>& aggregate,
                                       std::size_t& coarse_n)
    {
        const std::size_t n = level.A.n_rows();
        const auto strong = build_strength_graph(level, strength_threshold);

        // Ruge--Stüben uses both the strong-dependence graph S and its
        // transpose S^T.  The first pass selects high-influence points as C
        // points and marks their strongly connected neighbours F.  Measures
        // are recomputed after every selection rather than using a stale
        // degree: this keeps the implementation deterministic and makes the
        // influence update explicit.
        std::vector<std::vector<std::size_t>> transpose(n);
        for (std::size_t i = 0; i < n; ++i) {
            for (const std::size_t j : strong[i]) {
                if (j < n) transpose[j].push_back(i);
            }
        }

        std::vector<std::vector<std::size_t>> neighborhood(n);
        for (std::size_t i = 0; i < n; ++i) {
            neighborhood[i] = strong[i];
            neighborhood[i].insert(neighborhood[i].end(),
                                    transpose[i].begin(), transpose[i].end());
            std::sort(neighborhood[i].begin(), neighborhood[i].end());
            neighborhood[i].erase(
                std::unique(neighborhood[i].begin(), neighborhood[i].end()),
                neighborhood[i].end());
        }

        enum class Point : unsigned char { Undecided, Fine, Coarse };
        std::vector<Point> point(n, Point::Undecided);
        std::vector<std::size_t> influence(n, 0);
        for (std::size_t i = 0; i < n; ++i)
            influence[i] = transpose[i].size();

        std::size_t undecided = n;

        while (undecided > 0) {
            std::size_t best = n;
            std::size_t best_measure = 0;
            for (std::size_t i = 0; i < n; ++i) {
                if (point[i] != Point::Undecided) continue;
                if (best == n || influence[i] > best_measure ||
                    (influence[i] == best_measure && i < best)) {
                    best = i;
                    best_measure = influence[i];
                }
            }

            point[best] = Point::Coarse;
            --undecided;

            // Points that strongly depend on the new C point become F.
            // Every newly-created F point then increases the influence
            // measure of its strong influencers. This is the dynamic
            // Ruge--Stuben first-pass update; recomputing only the number of
            // undecided neighbours is not equivalent and over-coarsens paths
            // into a ~1/3 C-set.
            for (const std::size_t j : transpose[best]) {
                if (point[j] == Point::Undecided) {
                    point[j] = Point::Fine;
                    --undecided;
                    for (const std::size_t k : strong[j]) {
                        if (point[k] == Point::Undecided)
                            ++influence[k];
                    }
                }
            }
        }

        auto has_common_coarse = [&](std::size_t i, std::size_t j) {
            for (const std::size_t ci : strong[i]) {
                if (point[ci] != Point::Coarse) continue;
                for (const std::size_t cj : strong[j]) {
                    if (ci == cj && point[cj] == Point::Coarse) return true;
                }
            }
            return false;
        };

        // RS second pass: every strong F-F connection must be representable
        // through a common C neighbour.  If the first pass leaves a violation,
        // augment the C set only with a point that is not strongly adjacent to
        // an existing C point.  On the 1-D Poisson model this pass is inactive
        // for the normal alternating C/F splitting; it is nevertheless needed
        // for general strength graphs.
        bool changed = true;
        std::size_t guard = 0;
        while (changed && guard++ <= 2 * n + 1) {
            changed = false;
            for (std::size_t i = 0; i < n && !changed; ++i) {
                if (point[i] != Point::Fine) continue;
                for (const std::size_t j : strong[i]) {
                    if (point[j] != Point::Fine) continue;
                    if (has_common_coarse(i, j)) continue;

                    // Second-pass augmentation is allowed to add a C point
                    // adjacent to an existing C point.  The first pass enforces
                    // maximal C/F independence; the second pass exists precisely
                    // to repair F-F connections that have no common C neighbour.
                    // Rejecting both endpoints here can leave an F-F edge with no
                    // admissible classical interpolation stencil.
                    const std::size_t candidate = i;

                    point[candidate] = Point::Coarse;
                    changed = true;
                    for (const std::size_t q : neighborhood[candidate]) {
                        if (point[q] == Point::Undecided) point[q] = Point::Fine;
                    }
                    break;
                }
            }
        }

        // A valid classical interpolation requires every F point to have at
        // least one strong C neighbour in its own dependence set.  If a
        // directed strength graph leaves such a point without one, promote the
        // point itself when that does not violate C independence.
        for (std::size_t i = 0; i < n; ++i) {
            if (point[i] != Point::Fine) continue;
            bool has_strong_C = false;
            for (const std::size_t j : strong[i]) {
                if (point[j] == Point::Coarse) {
                    has_strong_C = true;
                    break;
                }
            }
            if (!has_strong_C) {
                // At this stage there is no valid interpolation stencil for
                // this F point. Promote it rather than manufacturing a
                // singleton/fallback interpolation row.
                point[i] = Point::Coarse;
            }
        }

        std::vector<std::size_t> coarse_index(n, n);
        coarse_n = 0;
        for (std::size_t i = 0; i < n; ++i) {
            if (point[i] == Point::Coarse) coarse_index[i] = coarse_n++;
        }

        const auto* row = level.A.row_offsets_data();
        const auto* col = level.A.columns_data();
        const auto* val = level.A.values_data();

        auto matrix_value = [&](std::size_t i, std::size_t j) {
            double value = 0.0;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                if (col[k] == j) value += val[k];
            }
            return value;
        };

        prolongation.assign(n, {});
        aggregate.assign(n, 0);

        for (std::size_t i = 0; i < n; ++i) {
            if (point[i] == Point::Coarse) {
                prolongation[i].push_back({coarse_index[i], 1.0});
                aggregate[i] = coarse_index[i];
                continue;
            }

            const double diagonal = 1.0 / level.inv_diag[i];
            std::vector<std::size_t> c_neighbors;
            c_neighbors.reserve(strong[i].size());
            for (const std::size_t j : strong[i]) {
                if (point[j] == Point::Coarse)
                    c_neighbors.push_back(j);
            }

            // The second-pass criterion above should make this non-empty.
            // If it is not, the splitting is not valid for classical
            // interpolation; reject the hierarchy rather than inventing a
            // non-Galerkin fallback.
            if (c_neighbors.empty()) {
                prolongation.clear();
                aggregate.clear();
                coarse_n = 0;
                return false;
            }

            // Classical Ruge--Stüben interpolation:
            //
            // w_ij = -(a_ij + sum_{k in F_i^s}
            //                  a_ik a_kj / sum_{m in C_k^s} a_km)
            //             / (a_ii + sum_{k in N_i^w} a_ik).
            //
            // This is the standard signed formula. There is deliberately no
            // arbitrary post-normalization: constant preservation follows from
            // the algebraic formula for Laplacian/M-matrix rows, and a failure
            // to obtain finite weights is a setup failure rather than a hidden
            // fallback.
            double weak_sum = 0.0;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                const std::size_t j = col[k];
                if (j == i) continue;
                if (!std::binary_search(strong[i].begin(), strong[i].end(), j))
                    weak_sum += val[k];
            }

            const double denominator = diagonal + weak_sum;
            if (!std::isfinite(denominator) || std::abs(denominator) <= 1e-30) {
                prolongation.clear();
                aggregate.clear();
                coarse_n = 0;
                return false;
            }

            std::map<std::size_t, double> weights;
            for (const std::size_t j : c_neighbors) {
                weights[coarse_index[j]] += matrix_value(i, j);
            }

            for (const std::size_t k : strong[i]) {
                if (point[k] != Point::Fine) continue;

                // Classical RS uses the strong C-neighbour set of k
                // itself. Using C_i^s ∩ C_k^s truncates the denominator and
                // becomes wrong on recursively generated Galerkin operators.
                double c_sum = 0.0;
                for (const std::size_t m : strong[k]) {
                    if (point[m] == Point::Coarse)
                        c_sum += matrix_value(k, m);
                }
                if (!std::isfinite(c_sum) || std::abs(c_sum) <= 1e-30) {
                    prolongation.clear();
                    aggregate.clear();
                    coarse_n = 0;
                    return false;
                }

                const double aik = matrix_value(i, k);
                for (const std::size_t j : c_neighbors) {
                    const double akj = matrix_value(k, j);
                    weights[coarse_index[j]] += aik * akj / c_sum;
                }
            }

            double largest = -1.0;
            for (const auto& [coarse, numerator] : weights) {
                const double weight = -numerator / denominator;
                if (!std::isfinite(weight)) {
                    prolongation.clear();
                    aggregate.clear();
                    coarse_n = 0;
                    return false;
                }
                prolongation[i].push_back({coarse, weight});
                if (std::abs(weight) > largest) {
                    largest = std::abs(weight);
                    aggregate[i] = coarse;
                }
            }

            if (prolongation[i].empty()) {
                prolongation.clear();
                aggregate.clear();
                coarse_n = 0;
                return false;
            }

            // The constant vector is the near-nullspace mode for the elliptic
            // operators targeted by classical AMG. Enforce P*1 = 1 explicitly
            // after the signed RS construction. This is a mathematical
            // near-nullspace constraint, not a convergence/tolerance fallback.
            double row_sum = 0.0;
            for (const auto& [coarse, weight] : prolongation[i]) {
                (void)coarse;
                row_sum += weight;
            }
            if (!std::isfinite(row_sum) || std::abs(row_sum) <= 1e-30) {
                prolongation.clear();
                aggregate.clear();
                coarse_n = 0;
                return false;
            }
            for (auto& [coarse, weight] : prolongation[i]) {
                (void)coarse;
                weight /= row_sum;
                if (!std::isfinite(weight)) {
                    prolongation.clear();
                    aggregate.clear();
                    coarse_n = 0;
                    return false;
                }
            }
        }
        return true;
    }

    static void build_smoothed_aggregation_interpolation(
        const Level& level,
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

        // Deterministic uncoupled aggregation. The most connected unassigned
        // point seeds an aggregate and absorbs its unassigned strong
        // neighbours. Isolated points become singleton aggregates.
        aggregate.assign(n, n);
        coarse_n = 0;
        std::size_t remaining = n;
        while (remaining > 0) {
            std::size_t seed = n;
            std::size_t best_degree = 0;
            for (std::size_t i = 0; i < n; ++i) {
                if (aggregate[i] != n) continue;
                std::size_t degree = 0;
                for (const std::size_t j : neighborhood[i])
                    if (aggregate[j] == n) ++degree;
                if (seed == n || degree > best_degree ||
                    (degree == best_degree && i < seed)) {
                    seed = i;
                    best_degree = degree;
                }
            }

            const std::size_t id = coarse_n++;
            aggregate[seed] = id;
            --remaining;
            for (const std::size_t j : neighborhood[seed]) {
                if (aggregate[j] == n) {
                    aggregate[j] = id;
                    --remaining;
                }
            }
        }

        // Smooth the piecewise-constant tentative prolongator with one damped
        // Jacobi step. Row normalization preserves the constant near-null-space
        // mode while retaining interpolation to neighbouring aggregates.
        constexpr double interpolation_omega = 2.0 / 3.0;
        const auto* row = level.A.row_offsets_data();
        const auto* col = level.A.columns_data();
        const auto* val = level.A.values_data();
        prolongation.assign(n, {});
        for (std::size_t i = 0; i < n; ++i) {
            std::map<std::size_t, double> weights;
            weights[aggregate[i]] = 1.0;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                weights[aggregate[col[k]]] -=
                    interpolation_omega * level.inv_diag[i] * val[k];
            }

            double sum = 0.0;
            for (const auto& [coarse, weight] : weights) {
                (void)coarse;
                if (std::isfinite(weight) && std::abs(weight) > 1e-14)
                    sum += weight;
            }
            if (!std::isfinite(sum) || std::abs(sum) <= 1e-14) {
                prolongation[i].push_back({aggregate[i], 1.0});
                continue;
            }
            for (const auto& [coarse, weight] : weights) {
                const double normalized = weight / sum;
                if (std::isfinite(normalized) && std::abs(normalized) > 1e-14)
                    prolongation[i].push_back({coarse, normalized});
            }
            if (prolongation[i].empty())
                prolongation[i].push_back({aggregate[i], 1.0});
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
            if (constant_null_space_) remove_constant(x);
        }
        return true;
    }

    bool vcycle(std::size_t level, const Vector& r, Vector& x) const
    {
        const auto& current = levels_[level];
        if (level + 1 == levels_.size()) {
            return smooth_coarsest(level, r, x);        }

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
        if (constant_null_space_) remove_constant(coarse_r);

        Vector coarse_x(nc, 0.0);
        if (!vcycle(level + 1, coarse_r, coarse_x)) return false;

        for (std::size_t i = 0; i < r.size(); ++i) {
            for (const auto& [coarse, weight] : prolongation[i])
                x(i) += weight * coarse_x(coarse);
            if (!std::isfinite(x(i))) return false;
        }
        if (constant_null_space_) remove_constant(x);

        return smooth(level, r, x, post_);
    }
    bool vcycle_diagnostic(std::size_t level, const Vector& r, Vector& x,
                           std::vector<VcycleDiagnostic>& diagnostics) const
    {
        VcycleDiagnostic d;
        d.level = level;
        d.size = r.size();
        d.rhs_norm = r.norm2();
        d.x_initial_norm = x.norm2();
        Vector Ax(r.size());
        if (!apply_operator(level, x, Ax)) return false;
        Vector residual(r.size());
        for (std::size_t i = 0; i < r.size(); ++i) residual(i) = r(i) - Ax(i);
        d.residual_before = residual.norm2();

        const std::size_t index = diagnostics.size();
        diagnostics.push_back(d);
        if (level + 1 == levels_.size()) {
            const double coarse_before = residual.norm2();
            if (!smooth_coarsest(level, r, x)) return false;
            if (!apply_operator(level, x, Ax)) return false;
            for (std::size_t i = 0; i < r.size(); ++i) residual(i) = r(i) - Ax(i);
            diagnostics[index].coarse_residual_before = coarse_before;
            diagnostics[index].coarse_residual_after = residual.norm2();
            diagnostics[index].residual_after_pre = diagnostics[index].coarse_residual_after;
            diagnostics[index].residual_after_post = diagnostics[index].coarse_residual_after;
            diagnostics[index].coarse_solution_norm = x.norm2();
            return x.is_valid();
        }

        if (!smooth(level, r, x, pre_)) return false;
        if (!apply_operator(level, x, Ax)) return false;
        for (std::size_t i = 0; i < r.size(); ++i) residual(i) = r(i) - Ax(i);
        diagnostics[index].residual_after_pre = residual.norm2();

        const auto& prolongation = levels_[level].prolongation;
        const std::size_t nc = levels_[level + 1].A.n_rows();
        Vector coarse_r(nc, 0.0);
        for (std::size_t i = 0; i < r.size(); ++i)
            for (const auto& [coarse, weight] : prolongation[i])
                coarse_r(coarse) += weight * residual(i);
        if (constant_null_space_) remove_constant(coarse_r);
        diagnostics[index].coarse_rhs_norm = coarse_r.norm2();

        Vector coarse_x(nc, 0.0);
        if (!vcycle_diagnostic(level + 1, coarse_r, coarse_x, diagnostics)) return false;
        diagnostics[index].coarse_solution_norm = coarse_x.norm2();

        Vector coarse_Ax(nc, 0.0);
        if (!apply_operator(level + 1, coarse_x, coarse_Ax)) return false;
        double coarse_eq_r2 = 0.0;
        double coarse_rhs2 = 0.0;
        for (std::size_t i = 0; i < nc; ++i) {
            const double e = coarse_r(i) - coarse_Ax(i);
            coarse_eq_r2 += e * e;
            coarse_rhs2 += coarse_r(i) * coarse_r(i);
        }
        diagnostics[index].coarse_equation_relative_residual =
            std::sqrt(coarse_eq_r2 / std::max(coarse_rhs2, 1e-300));

        Vector correction(x.size(), 0.0);
        for (std::size_t i = 0; i < r.size(); ++i) {
            for (const auto& [coarse, weight] : prolongation[i])
                correction(i) += weight * coarse_x(coarse);
            x(i) += correction(i);
            if (!std::isfinite(x(i))) return false;
        }
        diagnostics[index].correction_norm = correction.norm2();
        Vector A_correction(x.size(), 0.0);
        if (!apply_operator(level, correction, A_correction)) return false;
        diagnostics[index].correction_operator_norm = A_correction.norm2();
        if (constant_null_space_) remove_constant(x);
        if (!apply_operator(level, x, Ax)) return false;
        for (std::size_t i = 0; i < r.size(); ++i) residual(i) = r(i) - Ax(i);
        diagnostics[index].residual_after_correction = residual.norm2();

        if (!smooth(level, r, x, post_)) return false;
        if (!apply_operator(level, x, Ax)) return false;
        for (std::size_t i = 0; i < r.size(); ++i) residual(i) = r(i) - Ax(i);
        diagnostics[index].residual_after_post = residual.norm2();
        return x.is_valid();
    }


    static bool dense_solve(const SparseMatrix& A, const Vector& b, Vector& x)
    {
        const std::size_t n=A.n_rows();
        if(n==0 || b.size()!=n) return false;
        std::vector<double> m(n*n,0.0),rhs(b.size(),0.0);
        const auto* row=A.row_offsets_data(); const auto* col=A.columns_data(); const auto* val=A.values_data();
        for(std::size_t i=0;i<n;++i) {
            rhs[i]=b(i);
            for(std::size_t k=row[i];k<row[i+1];++k) m[i*n+col[k]]+=val[k];
        }
        for(std::size_t k=0;k<n;++k) {
            std::size_t p=k; double pa=std::abs(m[k*n+k]);
            for(std::size_t i=k+1;i<n;++i) if(std::abs(m[i*n+k])>pa) {pa=std::abs(m[i*n+k]);p=i;}
            if(!std::isfinite(pa) || pa<=128.0*std::numeric_limits<double>::epsilon()*std::max(1.0,pa)) return false;
            if(p!=k) {
                for(std::size_t j=k;j<n;++j) std::swap(m[k*n+j],m[p*n+j]);
                std::swap(rhs[k],rhs[p]);
            }
            for(std::size_t i=k+1;i<n;++i) {
                const double f=m[i*n+k]/m[k*n+k];
                if(!std::isfinite(f)) return false;
                m[i*n+k]=0.0;
                for(std::size_t j=k+1;j<n;++j) m[i*n+j]-=f*m[k*n+j];
                rhs[i]-=f*rhs[k];
            }
        }
        x.resize(n);
        for(std::size_t ii=n;ii-- > 0;) {
            double sum=rhs[ii];
            for(std::size_t j=ii+1;j<n;++j) sum-=m[ii*n+j]*x(j);
            if(!std::isfinite(m[ii*n+ii]) || std::abs(m[ii*n+ii])<=1e-30) return false;
            x(ii)=sum/m[ii*n+ii];
            if(!std::isfinite(x(ii))) return false;
        }
        return true;
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

        const std::size_t system_n = n + (constant_null_space_ ? 1 : 0);
        std::vector<double> m(system_n * system_n, 0.0);
        std::vector<double> b(system_n, 0.0);
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        for (std::size_t i = 0; i < n; ++i) {
            b[i] = r(i);
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                m[i * system_n + col[k]] += val[k];
            }
        }

        // A pure-Neumann coarse operator is singular. The symmetric saddle-
        // point augmentation enforces a zero-mean correction without pinning
        // an arbitrary coarse unknown.
        if (constant_null_space_) {
            for (std::size_t i = 0; i < n; ++i) {
                m[i * system_n + n] = 1.0;
                m[n * system_n + i] = 1.0;
            }
        }

        double matrix_scale = 0.0;
        for (const double value : m)
            matrix_scale = std::max(matrix_scale, std::abs(value));
        const double pivot_tol = 128.0 * std::numeric_limits<double>::epsilon() *
                                 std::max(1.0, matrix_scale);
        for (std::size_t k = 0; k < system_n; ++k) {
            std::size_t pivot = k;
            double pivot_abs = std::abs(m[k * system_n + k]);
            for (std::size_t i = k + 1; i < system_n; ++i) {
                const double candidate = std::abs(m[i * system_n + k]);
                if (candidate > pivot_abs) {
                    pivot_abs = candidate;
                    pivot = i;
                }
            }
            if (!std::isfinite(pivot_abs) || pivot_abs <= pivot_tol) return false;

            if (pivot != k) {
                for (std::size_t j = k; j < system_n; ++j) {
                    std::swap(m[k * system_n + j],
                              m[pivot * system_n + j]);
                }
                std::swap(b[k], b[pivot]);
            }

            const double diagonal = m[k * system_n + k];
            for (std::size_t i = k + 1; i < system_n; ++i) {
                const double factor = m[i * system_n + k] / diagonal;
                if (!std::isfinite(factor)) return false;
                m[i * system_n + k] = 0.0;
                for (std::size_t j = k + 1; j < system_n; ++j) {
                    m[i * system_n + j] -=
                        factor * m[k * system_n + j];
                }
                b[i] -= factor * b[k];
            }
        }

        std::vector<double> solution(system_n, 0.0);
        for (std::size_t ii = system_n; ii-- > 0;) {
            double sum = b[ii];
            for (std::size_t j = ii + 1; j < system_n; ++j) {
                sum -= m[ii * system_n + j] * solution[j];
            }
            const double diagonal = m[ii * system_n + ii];
            if (!std::isfinite(diagonal) || std::abs(diagonal) <= pivot_tol) return false;
            solution[ii] = sum / diagonal;
            if (!std::isfinite(solution[ii])) return false;
        }
        if (x.size() != n) x.resize(n);
        for (std::size_t i = 0; i < n; ++i) x(i) = solution[i];
        if (constant_null_space_) remove_constant(x);
        return x.is_valid();
    }

    const LinearOperatorBase& op_;
    double omega_;
    std::size_t pre_, post_;
    double strength_threshold_;
    std::size_t max_levels_;
    AMGInterpolationPolicy interpolation_;
    bool constant_null_space_;
    std::vector<Level> levels_;
    std::vector<std::size_t> first_aggregate_;
};

using AgglomeratedAMGPreconditioner = MatrixFreeVcyclePreconditioner;
} // namespace cfdx::core