#pragma once
#include "cfdx/core/linalg/amg_setup_footprint.h"
#include "cfdx/core/linalg/galerkin_csr_product.h"
#include "cfdx/core/linalg/algebraic_line_preconditioner.h"
#include "cfdx/core/linalg/linear_operator.h"
#include "cfdx/core/linalg/preconditioner.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <queue>
#include <stdexcept>
#include <string>
#include <set>
#include <utility>
#include <vector>

namespace cfdx::core {

enum class AMGInterpolationPolicy {
    DirectCF,
    SmoothedAggregation
};

enum class AMGSmootherPolicy {
    DampedJacobi,
    AlgebraicLines
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
                                   bool constant_null_space = false,
                                   bool require_smoother_certificate = true,
                                   AMGSmootherPolicy smoother = AMGSmootherPolicy::DampedJacobi,
                                   double prolongator_filter = 0.0,
                                   std::size_t max_cf_interpolation_entries = 0)
        : op_(op), configured_omega_(omega), omega_(omega), pre_(pre), post_(post),
          strength_threshold_(strength_threshold), max_levels_(max_levels),
          interpolation_(interpolation),
          constant_null_space_(constant_null_space),
          require_smoother_certificate_(require_smoother_certificate),
          smoother_policy_(smoother), prolongator_filter_(prolongator_filter),
          max_cf_interpolation_entries_(max_cf_interpolation_entries)
    {
        if (op.rows() != op.cols() || !(omega_ > 0.0 && omega_ < 2.0) ||
            !std::isfinite(strength_threshold_) || strength_threshold_ < 0.0 ||
            strength_threshold_ > 1.0 || max_levels_ == 0 ||
            !std::isfinite(prolongator_filter_) || prolongator_filter_ < 0.0 ||
            prolongator_filter_ >= 1.0 ||
            max_cf_interpolation_entries_ == 1) {
            throw std::invalid_argument("invalid AMG operator");
        }
    }

    bool setup(const SparseMatrix& A) override
    {
        last_setup_error_.clear();
        last_footprint_ = AMGSetupFootprint{};
        line_factor_fallback_levels_ = 0;
        numeric_line_fallbacks_ = 0;
        if (A.n_rows() != op_.rows() || A.n_cols() != op_.cols() ||
            A.n_rows() == 0 || !matrix_is_valid(A) ||
            (constant_null_space_ && !has_constant_null_space(A))) {
            last_setup_error_ = "invalid finest matrix, dimensions, CSR, or null-space contract";
            levels_.clear();
            first_aggregate_.clear();
            return false;
        }

        levels_.clear();
        first_aggregate_.clear();

        Level fine;
        fine.A = A;
        if (!build_diagonal(fine)) {
            last_setup_error_ = "invalid or zero diagonal on finest level";
            levels_.clear();
            return false;
        }
        if (!build_level_lines(fine)) {
            last_setup_error_ = "invalid finest algebraic-line factorization";
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
                    aggregate, coarse_n,
                    &last_footprint_.fallback_interpolation_rows, prolongator_filter_);
                if (levels_.size() == 1) {
                    std::vector<std::size_t> members(coarse_n, 0);
                    for (const std::size_t id : aggregate) {
                        if (id < coarse_n) ++members[id];
                    }
                    for (const std::size_t member_count : members) {
                        if (member_count == 1)
                            ++last_footprint_.singleton_aggregates;
                    }
                }
            } else {
                if (!build_cf_interpolation(levels_.back(), strength_threshold_,
                                            prolongation, aggregate, coarse_n,
                                            max_cf_interpolation_entries_)) {
                    last_setup_error_ = "Direct-CF interpolation construction failed at level " +
                                        std::to_string(levels_.size() - 1);
                    levels_.clear();
                    first_aggregate_.clear();
                    return false;
                }
            }
            if (coarse_n >= levels_.back().A.n_rows() || coarse_n == 0) {
                break;
            }

            GalerkinCSRProduct product;
            SparseMatrix coarse;
            if (!product.symbolic(levels_.back().A, prolongation, coarse_n) ||
                !product.numeric(levels_.back().A, prolongation, coarse,
                                 &last_footprint_.galerkin_dropped_entries) ||
                !matrix_is_valid(coarse)) {
                last_setup_error_ = "invalid Galerkin coarse matrix from level " +
                                    std::to_string(levels_.size() - 1);
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
                last_setup_error_ = "invalid or zero diagonal on coarse level " +
                                    std::to_string(levels_.size());
                levels_.clear();
                first_aggregate_.clear();
                return false;
            }
            if (!build_level_lines(next)) {
                last_setup_error_ = "invalid coarse algebraic-line factorization";
                levels_.clear();
                return false;
            }
            levels_.back().aggregate = std::move(aggregate);
            levels_.back().prolongation = std::move(prolongation);
            levels_.back().coarse_product = std::move(product);
            levels_.push_back(std::move(next));
        }

        if (levels_.empty() || levels_.front().A.n_rows() == 0) {
            last_setup_error_ = "empty hierarchy after coarsening";
            levels_.clear();
            return false;
        }

        // A fresh hierarchy needs the same damping certification as a
        // numerically refreshed hierarchy. Recompute from the constructor
        // value so a previous difficult matrix cannot permanently poison a
        // later full rebuild. The certificate is a CG-safety contract: it
        // stays fail-closed on the SPD/CG path, while the general
        // nonsymmetric preconditioner path (outer GMRES/FGMRES/BiCGStab,
        // e.g. the nested LSC/BFBt pressure solve on H = D Q^-1 G - C)
        // keeps the pre-certificate behaviour -- configured damping, the
        // outer Krylov true-residual checks remaining authoritative. Every
        // structural check above is fail-closed on both paths; only the
        // damping certificate is policy-dependent.
        double safe_omega = 0.0;
        if (compute_stable_smoother_omega(levels_, safe_omega, &last_setup_error_)) {
            omega_ = safe_omega;
        } else if (require_smoother_certificate_) {
            levels_.clear();
            first_aggregate_.clear();
            return false;
        } else {
            omega_ = configured_omega_;
            last_setup_error_.clear();
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

        // Transactionally stage only numeric state.  Copying `levels_` here
        // used to clone every aggregate and every prolongation row on each
        // pressure corrector even though interpolation is unchanged.  This
        // is particularly expensive for NACA N256's anisotropic hierarchy.
        // Keep the original transfer operators in levels_ until the staged
        // Galerkin matrices, diagonals, line factorizations and damping have
        // all passed their checks. Failure leaves the old hierarchy intact.
        std::vector<Level> refreshed(levels_.size());
        std::size_t refreshed_line_fallbacks = 0;
        for (std::size_t level = 0; level < refreshed.size(); ++level) {
            refreshed[level].use_lines = levels_[level].use_lines;
            // The line topology and its existing factors are used only when
            // a line smoother is active.  Its numeric refresh is itself
            // transactional, so preserve an independent copy until commit.
            if (refreshed[level].use_lines)
                refreshed[level].line_pc = levels_[level].line_pc;
        }
        refreshed.front().A = A;
        if (!build_diagonal(refreshed.front()) ||
            !refresh_level_lines(refreshed.front(), true, refreshed_line_fallbacks))
            return false;

        std::size_t dropped_entries = 0;
        for (std::size_t level = 0; level + 1 < refreshed.size(); ++level) {
            const std::size_t coarse_n = levels_[level + 1].A.n_rows();
            // Transactional fixed-graph Galerkin refresh: duplicate the old
            // coarse CSR once into staging, then recompute only its values.
            // Unlike a full numeric() build this does not reallocate or
            // re-copy the row offsets and columns after every multiplication.
            refreshed[level + 1].A = levels_[level + 1].A;
            if (!levels_[level].coarse_product.numeric_into_existing(
                    refreshed[level].A, levels_[level].prolongation,
                    refreshed[level + 1].A, &dropped_entries) ||
                !matrix_is_valid(refreshed[level + 1].A) ||
                refreshed[level + 1].A.n_rows() != coarse_n) return false;
            if (!build_diagonal(refreshed[level + 1]) ||
                !refresh_level_lines(refreshed[level + 1], false, refreshed_line_fallbacks))
                return false;
        }

        double refreshed_omega = 0.0;
        if (!compute_stable_smoother_omega(refreshed, refreshed_omega)) {
            if (require_smoother_certificate_) return false;
            refreshed_omega = omega_;
        }

        // Numeric-only commit: aggregate and interpolation storage is never
        // moved or copied.  Move-assignments are noexcept for these standard
        // containers, so no half-updated hierarchy can be observed.
        for (std::size_t level = 0; level < levels_.size(); ++level) {
            levels_[level].A = std::move(refreshed[level].A);
            levels_[level].inv_diag = std::move(refreshed[level].inv_diag);
            levels_[level].use_lines = refreshed[level].use_lines;
            if (levels_[level].use_lines)
                levels_[level].line_pc = std::move(refreshed[level].line_pc);
        }
        numeric_line_fallbacks_ = refreshed_line_fallbacks;
        omega_ = refreshed_omega;
        last_footprint_.galerkin_dropped_entries += dropped_entries;
        return true;
    }

    double smoother_omega() const noexcept { return omega_; }
    std::size_t line_factor_fallback_levels() const noexcept {
        return line_factor_fallback_levels_;
    }
    std::size_t numeric_line_fallbacks() const noexcept {
        return numeric_line_fallbacks_;
    }
    std::size_t finest_line_coverage() const noexcept {
        return levels_.empty() || !levels_.front().use_lines ? 0 :
               levels_.front().line_pc.covered_rows();
    }
    std::size_t finest_longest_line() const noexcept {
        return levels_.empty() || !levels_.front().use_lines ? 0 :
               levels_.front().line_pc.longest_line();
    }

    const std::string& last_setup_error() const noexcept { return last_setup_error_; }

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

    // Complexity diagnostics for AMG on stretched meshes.  High operator
    // complexity (sum(nnz(A_l))/nnz(A_0)) often exposes overly dense coarse
    // Galerkin stencils even if the fine-grid Krylov count looks good.
    std::vector<std::size_t> hierarchy_level_nnz() const
    {
        std::vector<std::size_t> sizes;
        sizes.reserve(levels_.size());
        for (const auto& level : levels_) sizes.push_back(level.A.nnz());
        return sizes;
    }

    double operator_complexity() const noexcept
    {
        if (levels_.empty() || levels_.front().A.nnz() == 0) return 0.0;
        double total = 0.0;
        for (const auto& level : levels_) total += static_cast<double>(level.A.nnz());
        return total / static_cast<double>(levels_.front().A.nnz());
    }

    double grid_complexity() const noexcept
    {
        if (levels_.empty() || levels_.front().A.n_rows() == 0) return 0.0;
        double total = 0.0;
        for (const auto& level : levels_) total += static_cast<double>(level.A.n_rows());
        return total / static_cast<double>(levels_.front().A.n_rows());
    }

    const AMGSetupFootprint& setup_footprint() const noexcept
    {
        // Reporting-only counters of the AMG setup decisions that can degrade
        // the coarse-grid quality: singleton smoothed-aggregation aggregates,
        // fallback interpolation rows and Galerkin entries dropped as
        // numerically zero. The counters are reset at the start of setup()
        // and are deliberately kept when a setup attempt fails, so the
        // footprint of a failed attempt stays visible for diagnostics.
        // update_values() accumulates Galerkin drops only after a numeric
        // hierarchy refresh succeeds; failed refreshes do not commit counts. No tolerance or acceptance criterion reads them.
        return last_footprint_;
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

    struct ProlongationRowDiagnostic {
        std::size_t fine_index = 0;
        std::vector<std::pair<std::size_t, double>> entries;
    };

    // Diagnostic-only access to the first/last interpolation rows. This exposes
    // the actual stored P weights without changing the hierarchy construction.
    std::vector<ProlongationRowDiagnostic> prolongation_boundary_rows(
        std::size_t level, std::size_t count) const
    {
        std::vector<ProlongationRowDiagnostic> out;
        if (level + 1 >= levels_.size() || count == 0) return out;
        const auto& P = levels_[level].prolongation;
        const std::size_t n = P.size();
        const std::size_t take = std::min(count, n);
        out.reserve(2 * take);
        for (std::size_t i = 0; i < take; ++i)
            out.push_back({i, P[i]});
        const std::size_t begin = n - take;
        for (std::size_t i = begin; i < n; ++i)
            if (i >= take)
                out.push_back({i, P[i]});
        return out;
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

    struct TransferColumnCoverage {
        std::size_t zero_columns = 0;
        std::size_t min_nnz = 0;
        std::size_t max_nnz = 0;
        std::size_t anchored_columns = 0;
    };

    TransferColumnCoverage transfer_column_coverage(std::size_t level) const
    {
        if (level + 1 >= levels_.size())
            throw std::out_of_range("transfer_column_coverage: invalid level");
        const auto& P = levels_[level].prolongation;
        const std::size_t nc = levels_[level + 1].A.n_rows();
        std::vector<std::size_t> nnz(nc, 0), anchors(nc, 0);
        for (const auto& row : P) {
            for (const auto& [c, w] : row) {
                if (c >= nc) continue;
                ++nnz[c];
                if (row.size() == 1 && std::abs(w - 1.0) <= 1e-14) ++anchors[c];
            }
        }
        TransferColumnCoverage out;
        out.min_nnz = nnz.empty() ? 0 : *std::min_element(nnz.begin(), nnz.end());
        out.max_nnz = nnz.empty() ? 0 : *std::max_element(nnz.begin(), nnz.end());
        out.zero_columns = static_cast<std::size_t>(std::count(nnz.begin(), nnz.end(), 0));
        out.anchored_columns = static_cast<std::size_t>(std::count_if(anchors.begin(), anchors.end(),
            [](std::size_t n) { return n > 0; }));
        return out;
    }

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

    double sine_mode_vcycle_energy_ratio(std::size_t level,
                                            std::size_t mode) const
    {
        if (level >= levels_.size() || mode == 0)
            throw std::out_of_range("sine_mode_vcycle_energy_ratio: invalid level or mode");
        const std::size_t n = levels_[level].A.n_rows();
        Vector exact(n, 0.0);
        constexpr double pi = 3.14159265358979323846;
        for (std::size_t i = 0; i < n; ++i)
            exact(i) = std::sin(pi * static_cast<double>(mode * (i + 1)) /
                                 static_cast<double>(n + 1));
        Vector rhs(n, 0.0);
        if (!apply_operator(level, exact, rhs))
            return std::numeric_limits<double>::infinity();
        Vector x(n, 0.0);
        if (!vcycle(level, rhs, x))
            return std::numeric_limits<double>::infinity();

        Vector error(n, 0.0), Aerror(n, 0.0), Aexact(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) error(i) = exact(i) - x(i);
        if (!apply_operator(level, error, Aerror) ||
            !apply_operator(level, exact, Aexact))
            return std::numeric_limits<double>::infinity();

        double eAe = 0.0;
        double exactAe = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            eAe += error(i) * Aerror(i);
            exactAe += exact(i) * Aexact(i);
        }
        if (!std::isfinite(eAe) || !std::isfinite(exactAe) || exactAe <= 0.0)
            return std::numeric_limits<double>::infinity();
        return std::sqrt(std::max(0.0, eAe) / exactAe);
    }

    double two_grid_sine_mode_energy_ratio(std::size_t level, std::size_t mode) const
    {
        if (level + 1 >= levels_.size() || mode == 0)
            throw std::out_of_range("two_grid_sine_mode_energy_ratio: invalid level or mode");
        const std::size_t n = levels_[level].A.n_rows();
        Vector exact(n, 0.0);
        constexpr double pi = 3.14159265358979323846;
        for (std::size_t i = 0; i < n; ++i)
            exact(i) = std::sin(pi * static_cast<double>(mode * (i + 1)) /
                                 static_cast<double>(n + 1));
        Vector rhs(n, 0.0);
        if (!apply_operator(level, exact, rhs))
            return std::numeric_limits<double>::infinity();

        Vector x(n, 0.0), Ax(n), res(n);
        if (!smooth(level, rhs, x, pre_) || !apply_operator(level, x, Ax))
            return std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < n; ++i) res(i) = rhs(i) - Ax(i);

        const std::size_t nc = levels_[level + 1].A.n_rows();
        Vector rc(nc, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (const auto& [col, weight] : levels_[level].prolongation[i])
                rc(col) += weight * res(i);

        Vector ec(nc, 0.0), r = rc, z(nc, 0.0), p(nc, 0.0), Ap(nc, 0.0);
        const auto& Ac = levels_[level + 1].A;
        const auto* row = Ac.row_offsets_data();
        const auto* col = Ac.columns_data();
        const auto* val = Ac.values_data();
        double rz = 0.0, rhs_norm2 = 0.0;
        for (std::size_t i = 0; i < nc; ++i) {
            rhs_norm2 += rc(i) * rc(i);
            double diag = 0.0;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k)
                if (col[k] == i) diag += val[k];
            if (!std::isfinite(diag) || diag <= 0.0)
                return std::numeric_limits<double>::infinity();
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
                if (!std::isfinite(pAp) || pAp <= 0.0)
                    return std::numeric_limits<double>::infinity();
                const double alpha = rz / pAp;
                for (std::size_t i = 0; i < nc; ++i) {
                    ec(i) += alpha * p(i);
                    r(i) -= alpha * Ap(i);
                }
                double rr = 0.0;
                for (std::size_t i = 0; i < nc; ++i) rr += r(i) * r(i);
                if (std::sqrt(rr) <= target) { converged = true; break; }
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
        for (std::size_t i = 0; i < n; ++i) {
            double correction = 0.0;
            for (const auto& [coarse, weight] : levels_[level].prolongation[i])
                correction += weight * ec(coarse);
            x(i) += correction;
        }
        if (!smooth(level, rhs, x, post_))
            return std::numeric_limits<double>::infinity();

        Vector error(n, 0.0), Aerror(n), Aexact(n);
        for (std::size_t i = 0; i < n; ++i) error(i) = exact(i) - x(i);
        if (!apply_operator(level, error, Aerror) ||
            !apply_operator(level, exact, Aexact))
            return std::numeric_limits<double>::infinity();
        double eAe = 0.0, exactAe = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            eAe += error(i) * Aerror(i);
            exactAe += exact(i) * Aexact(i);
        }
        if (!std::isfinite(eAe) || !std::isfinite(exactAe) || exactAe <= 0.0)
            return std::numeric_limits<double>::infinity();
        return std::sqrt(std::max(0.0, eAe) / exactAe);
    }

    double sine_mode_smoother_residual_ratio(std::size_t level,
                                             std::size_t mode,
                                             std::size_t sweeps) const
    {
        if (level >= levels_.size() || mode == 0) {
            throw std::out_of_range("sine_mode_smoother_residual_ratio: invalid level or mode");
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
        const double before = rhs.norm2();
        if (!(before > 0.0) || !std::isfinite(before)) return 0.0;
        Vector x(n, 0.0);
        if (!smooth(level, rhs, x, sweeps)) {
            return std::numeric_limits<double>::infinity();
        }
        Vector Ax(n, 0.0);
        if (!apply_operator(level, x, Ax)) {
            return std::numeric_limits<double>::infinity();
        }
        Vector residual(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) residual(i) = rhs(i) - Ax(i);
        return residual.norm2() / std::max(before, 1e-300);
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
        AlgebraicLinePreconditioner line_pc;
        bool use_lines = false;
        std::vector<std::size_t> aggregate;
        std::vector<TransferRow> prolongation;
        GalerkinCSRProduct coarse_product;
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

    // A symmetric Jacobi smoother is stable when omega * lambda_max(D^-1/2 A D^-1/2) < 2.
    // The maximum absolute normalized row sum is a conservative Gershgorin upper bound,
    // evaluated on the refreshed Galerkin levels (not the stale pre-update operators).
    static double normalized_gershgorin_upper_bound(const Level& level)
    {
        const SparseMatrix& A = level.A;
        const std::size_t n = A.n_rows();
        if (level.inv_diag.size() != n || n == 0)
            return std::numeric_limits<double>::infinity();
        std::vector<double> diagonal(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            if (!(level.inv_diag[i] > 0.0) || !std::isfinite(level.inv_diag[i]))
                return std::numeric_limits<double>::infinity();
            diagonal[i] = 1.0 / level.inv_diag[i];
            if (!(diagonal[i] > 0.0) || !std::isfinite(diagonal[i]))
                return std::numeric_limits<double>::infinity();
        }
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        double upper_bound = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double row_sum = 0.0;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                const std::size_t j = col[k];
                // Separate square roots avoid overflow/underflow in
                // diagonal[i] * diagonal[j] when the normalized scale itself
                // is representable.
                const double scale =
                    std::sqrt(diagonal[i]) * std::sqrt(diagonal[j]);
                if (!(scale > 0.0) || !std::isfinite(scale))
                    return std::numeric_limits<double>::infinity();
                row_sum += std::abs(val[k]) / scale;
            }
            if (!std::isfinite(row_sum))
                return std::numeric_limits<double>::infinity();
            upper_bound = std::max(upper_bound, row_sum);
        }
        return upper_bound;
    }

    bool compute_stable_smoother_omega(
        const std::vector<Level>& hierarchy, double& safe_omega,
        std::string* failure_reason = nullptr) const
    {
        const auto fail = [failure_reason](const std::string& reason) {
            if (failure_reason) *failure_reason = reason;
            return false;
        };
        if (hierarchy.empty()) return fail("empty AMG hierarchy");
        double maximum_bound = 0.0;
        for (std::size_t level_index = 0; level_index < hierarchy.size(); ++level_index) {
            const auto& level = hierarchy[level_index];
            if (level.inv_diag.size() != level.A.n_rows())
                return fail("diagonal size mismatch at level " + std::to_string(level_index));
            for (std::size_t row = 0; row < level.inv_diag.size(); ++row) {
                if (!(level.inv_diag[row] > 0.0) || !std::isfinite(level.inv_diag[row]))
                    return fail("Jacobi stability certificate requires positive finite diagonal; level=" +
                                std::to_string(level_index) + " row=" + std::to_string(row) +
                                " inverse_diagonal=" + std::to_string(level.inv_diag[row]));
            }
            const double bound = normalized_gershgorin_upper_bound(level);
            if (!std::isfinite(bound) || !(bound > 0.0))
                return fail("normalized Gershgorin bound is invalid at level " +
                            std::to_string(level_index) + " bound=" + std::to_string(bound));
            maximum_bound = std::max(maximum_bound, bound);
        }
        if (!(maximum_bound > 0.0) || !std::isfinite(maximum_bound))
            return fail("maximum normalized Gershgorin bound is invalid");

        // Re-evaluate against the configured damping on every setup/refresh.
        // This permits recovery when a later matrix is less restrictive,
        // while retaining the conservative 10% margin below the Jacobi
        // stability limit for the current hierarchy.
        safe_omega = std::min(configured_omega_, 1.8 / maximum_bound);
        if (!std::isfinite(safe_omega) || !(safe_omega > 0.0 && safe_omega < 2.0))
            return fail("computed Jacobi damping is outside (0,2): " + std::to_string(safe_omega));
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

    bool build_level_lines(Level& level) {
        if (smoother_policy_ != AMGSmootherPolicy::AlgebraicLines) return true;
        // A Galerkin coarse operator can be SPD yet exhibit locally
        // ill-conditioned or non-M-matrix blocks. In that case the optional
        // line smoother MUST NOT be applied: fall back to certified Jacobi
        // while retaining the full hierarchy and physical-residual gate.
        if (!level.line_pc.setup(level.A)) {
            level.use_lines = false;
            ++line_factor_fallback_levels_;
            return true;
        }
        // On a well-conditioned/isotropic grid retain the original Jacobi
        // smoother: only activate the block solver for coherent long paths.
        level.use_lines = level.line_pc.longest_line() >= 4 &&
                          level.line_pc.covered_rows() * 4 >= level.A.n_rows();
        return true;
    }

    bool refresh_level_lines(Level& level, bool finest,
                             std::size_t& fallback_count) const {
        if (smoother_policy_ != AMGSmootherPolicy::AlgebraicLines ||
            !level.use_lines) return true;
        // On the physical (finest) matrix, a changed directional strength
        // graph invalidates numeric-only reuse: require full setup(). On a
        // Galerkin level, whose graph is itself only an approximation, a
        // rejected line factor can safely fall back to Jacobi. This avoids
        // poisoning the active hierarchy while retaining the outer true-
        // residual contract. Record the reversion for diagnostics.
        if (level.line_pc.update_values(level.A)) return true;
        if (finest) return false;
        level.use_lines = false;
        ++fallback_count;
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
