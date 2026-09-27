#pragma once

#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/advanced_preconditioners.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::core {

// Native reduction preconditioner inspired by the reduction/Schur hierarchy
// used by HYPRE MGR. The fine variables are relaxed with a local F operator,
// while the retained C variables are solved with an AMG-preconditioned
// approximate coarse operator
//
//   A = [AFF AFC; ACF ACC]
//   S~ = ACC - ACF diag(AFF)^-1 AFC.
//
// The implementation deliberately keeps the reduction explicit and
// dependency-free. It is not a wrapper around HYPRE MGR and does not claim
// feature or performance parity with it.
class NativeMGRPreconditioner final : public Preconditioner {
public:
    enum class FineRelaxation {
        Diagonal,
        ILU0
    };

    using ReductionPartition =
        std::pair<std::vector<std::size_t>, std::vector<std::size_t>>;

    NativeMGRPreconditioner(
        std::vector<std::size_t> fine_variables,
        std::vector<std::size_t> coarse_variables,
        FineRelaxation fine_relaxation = FineRelaxation::Diagonal,
        std::vector<ReductionPartition> additional_levels = {})
        : fine_(std::move(fine_variables)),
          coarse_(std::move(coarse_variables)),
          fine_relaxation_(fine_relaxation),
          additional_levels_(std::move(additional_levels)) {}

    bool setup(const SparseMatrix& A) override {
        clear();
        n_ = A.n_rows();
        if (!validate_partition(A)) {
            last_error_ = "invalid MGR fine/coarse partition";
            return false;
        }
        if (!build_reduction(A)) return false;
        if (fine_relaxation_ == FineRelaxation::ILU0) {
            fine_ilu_ = std::make_unique<ILU0Preconditioner>();
            if (!fine_ilu_->setup(fine_operator_)) {
                last_error_ = "ILU(0) setup failed on MGR F operator";
                clear_state_only();
                return false;
            }
        }
        if (!setup_coarse_solver()) {
            last_error_ = "coarse solver setup failed on MGR reduced operator";
            clear_state_only();
            return false;
        }
        ready_ = true;
        ++setup_count_;
        last_error_.clear();
        return true;
    }

    bool update_values(const SparseMatrix& A) override {
        if (!ready_) return setup(A);
        if (!validate_partition(A)) {
            last_error_ = "invalid MGR fine/coarse partition";
            return false;
        }
        SparseMatrix next;
        SparseMatrix next_fine_operator;
        std::vector<double> next_fine_diag;
        if (!assemble_reduced(A, next, next_fine_diag, next_fine_operator)) return false;

        const bool fine_pattern_ok =
            fine_relaxation_ == FineRelaxation::Diagonal ||
            same_pattern(fine_operator_, next_fine_operator);
        if (same_pattern(coarse_operator_, next) && fine_pattern_ok) {
            coarse_operator_ = std::move(next);
            fine_operator_ = std::move(next_fine_operator);
            fine_inv_diag_ = std::move(next_fine_diag);
            if (fine_relaxation_ == FineRelaxation::ILU0) {
                if (!fine_ilu_ || !fine_ilu_->setup(fine_operator_)) {
                    last_error_ = "ILU(0) numeric refresh failed on MGR F operator";
                    return false;
                }
            }
            if (!update_coarse_solver()) {
                last_error_ = "coarse solver numeric refresh failed on MGR reduced operator";
                return false;
            }
            ++numeric_updates_;
            last_error_.clear();
            return true;
        }

        last_error_ = fine_pattern_ok
            ? "MGR reduced CSR pattern changed; explicit setup() required"
            : "MGR F-operator CSR pattern changed; explicit setup() required";
        return false;
    }

    bool apply(const Vector& residual, Vector& correction) const override {
        if (!ready_ || residual.size() != n_ || correction.size() != n_)
            return false;

        correction.fill(0.0);
        Vector fine_rhs(fine_.size(), 0.0);
        Vector fine_correction(fine_.size(), 0.0);
        Vector coarse_rhs(coarse_.size(), 0.0);

        for (std::size_t i = 0; i < fine_.size(); ++i)
            fine_rhs(i) = residual(fine_[i]);

        // F-relaxation. Diagonal Jacobi is the baseline; ILU(0) is the
        // stronger local option used to test whether weak F relaxation is the
        // limiting factor on coupled velocity blocks.
        if (fine_relaxation_ == FineRelaxation::ILU0) {
            if (!fine_ilu_ || !fine_ilu_->apply(fine_rhs, fine_correction))
                return false;
        } else {
            for (std::size_t i = 0; i < fine_.size(); ++i)
                fine_correction(i) = fine_inv_diag_[i] * fine_rhs(i);
        }
        for (std::size_t i = 0; i < fine_.size(); ++i)
            correction(fine_[i]) = fine_correction(i);

        // Coarse residual: r_C - A_CF z_F.
        for (std::size_t i = 0; i < coarse_.size(); ++i) {
            double value = residual(coarse_[i]);
            for (const auto& [j, aij] : coarse_to_fine_[i])
                value -= aij * correction(fine_[j]);
            coarse_rhs(i) = value;
        }

        Vector coarse_correction(coarse_.size(), 0.0);
        if (!apply_coarse_solver(coarse_rhs, coarse_correction))
            return false;

        // Prolongate and apply the approximate block elimination correction:
        // z_F <- z_F - diag(A_FF)^-1 A_FC z_C.
        for (std::size_t i = 0; i < coarse_.size(); ++i)
            correction(coarse_[i]) += coarse_correction(i);

        Vector back_substitution_rhs(fine_.size(), 0.0);
        for (std::size_t i = 0; i < fine_.size(); ++i) {
            double coupling = 0.0;
            for (const auto& [j, aij] : fine_to_coarse_[i])
                coupling += aij * coarse_correction(j);
            back_substitution_rhs(i) = fine_rhs(i) - coupling;
        }
        if (fine_relaxation_ == FineRelaxation::ILU0) {
            if (!fine_ilu_ || !fine_ilu_->apply(
                    back_substitution_rhs, fine_correction))
                return false;
        } else {
            for (std::size_t i = 0; i < fine_.size(); ++i)
                fine_correction(i) =
                    fine_inv_diag_[i] * back_substitution_rhs(i);
        }
        for (std::size_t i = 0; i < fine_.size(); ++i)
            correction(fine_[i]) = fine_correction(i);

        return correction.is_valid();
    }

    const char* name() const override { return "NativeMGR-Reduction"; }
    const std::string& last_error() const noexcept { return last_error_; }
    std::size_t setup_count() const noexcept { return setup_count_; }
    std::size_t numeric_updates() const noexcept { return numeric_updates_; }
    std::size_t fine_size() const noexcept { return fine_.size(); }
    std::size_t coarse_size() const noexcept { return coarse_.size(); }
    FineRelaxation fine_relaxation() const noexcept { return fine_relaxation_; }
    std::size_t reduction_levels() const noexcept {
        return 1 + additional_levels_.size();
    }

private:
    bool validate_partition(const SparseMatrix& A) const {
        if (A.n_rows() != A.n_cols() || !A.is_consistent() ||
            fine_.empty() || coarse_.empty())
            return false;
        std::vector<bool> seen(A.n_rows(), false);
        for (const auto i : fine_) {
            if (i >= A.n_rows() || seen[i]) return false;
            seen[i] = true;
        }
        for (const auto i : coarse_) {
            if (i >= A.n_rows() || seen[i]) return false;
            seen[i] = true;
        }
        return std::all_of(seen.begin(), seen.end(), [](bool v) { return v; });
    }

    bool build_reduction(const SparseMatrix& A) {
        return assemble_reduced(
            A, coarse_operator_, fine_inv_diag_, fine_operator_);
    }

    bool setup_coarse_solver() {
        if (additional_levels_.empty()) {
            coarse_mgr_.reset();
            return coarse_amg_.setup(coarse_operator_);
        }
        const auto& level = additional_levels_.front();
        auto nested = std::make_unique<NativeMGRPreconditioner>(
            level.first, level.second, fine_relaxation_,
            std::vector<ReductionPartition>(
                additional_levels_.begin() + 1, additional_levels_.end()));
        if (!nested->setup(coarse_operator_))
            return false;
        coarse_mgr_ = std::move(nested);
        return true;
    }

    bool update_coarse_solver() {
        if (coarse_mgr_) return coarse_mgr_->update_values(coarse_operator_);
        return coarse_amg_.update_values(coarse_operator_);
    }

    bool apply_coarse_solver(
        const Vector& rhs, Vector& correction) const {
        if (coarse_mgr_) return coarse_mgr_->apply(rhs, correction);
        return coarse_amg_.apply(rhs, correction);
    }

    bool assemble_reduced(const SparseMatrix& A, SparseMatrix& reduced,
                          std::vector<double>& inverse_diag,
                          SparseMatrix& fine_operator) {
        const std::size_t nf = fine_.size();
        const std::size_t nc = coarse_.size();
        std::vector<std::size_t> fine_pos(A.n_rows(), nf);
        std::vector<std::size_t> coarse_pos(A.n_rows(), nc);
        for (std::size_t i = 0; i < nf; ++i) fine_pos[fine_[i]] = i;
        for (std::size_t i = 0; i < nc; ++i) coarse_pos[coarse_[i]] = i;

        inverse_diag.assign(nf, 0.0);
        coarse_to_fine_.assign(nc, {});
        fine_to_coarse_.assign(nf, {});
        std::vector<std::vector<std::pair<std::size_t,double>>> acc(nc);

        const auto* ro = A.row_offsets_data();
        const auto* ci = A.columns_data();
        const auto* av = A.values_data();

        for (std::size_t i = 0; i < nf; ++i) {
            const auto row = fine_[i];
            bool found = false;
            for (std::size_t k = ro[row]; k < ro[row + 1]; ++k) {
                const auto col = ci[k];
                if (col == row) {
                    inverse_diag[i] = av[k];
                    found = true;
                } else if (coarse_pos[col] < nc) {
                    fine_to_coarse_[i].push_back({coarse_pos[col], av[k]});
                } else if (fine_pos[col] < nf) {
                    fine_rows[i].push_back({fine_pos[col], av[k]});
                }
            }
            if (!found || !std::isfinite(inverse_diag[i]) ||
                std::abs(inverse_diag[i]) <= std::numeric_limits<double>::epsilon()) {
                last_error_ = "singular fine diagonal in MGR reduction";
                return false;
            }
            inverse_diag[i] = 1.0 / inverse_diag[i];
        }

        // ACC and ACF are assembled directly; the diagonal F approximation
        // contributes -ACF diag(AFF)^-1 AFC.
        std::vector<std::vector<std::pair<std::size_t,double>>> fine_rows(nf);
        std::vector<std::vector<std::pair<std::size_t,double>>> rows(nc);
        for (std::size_t i = 0; i < nc; ++i) {
            const auto row = coarse_[i];
            for (std::size_t k = ro[row]; k < ro[row + 1]; ++k) {
                const auto col = ci[k];
                if (coarse_pos[col] < nc)
                    rows[i].push_back({coarse_pos[col], av[k]});
                else if (fine_pos[col] < nf)
                    coarse_to_fine_[i].push_back({fine_pos[col], av[k]});
            }
        }

        for (std::size_t i = 0; i < nc; ++i) {
            for (const auto& [f, acf] : coarse_to_fine_[i]) {
                for (const auto& [c, afc] : fine_to_coarse_[f]) {
                    rows[i].push_back({c, -acf * inverse_diag[f] * afc});
                }
            }
        }

        fine_operator = SparseMatrix(nf, nf);
        for (std::size_t i = 0; i < nf; ++i) {
            auto& row = fine_rows[i];
            std::sort(row.begin(), row.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });
            for (std::size_t k = 0; k < row.size();) {
                const std::size_t column = row[k].first;
                double value = 0.0;
                while (k < row.size() && row[k].first == column)
                    value += row[k++].second;
                if (!std::isfinite(value)) {
                    last_error_ = "non-finite coefficient in MGR F operator";
                    return false;
                }
                if (value != 0.0)
                    fine_operator.push_back(i, column, value);
            }
        }
        fine_operator.finalize();
        if (!fine_operator.is_consistent()) {
            last_error_ = "invalid MGR F operator";
            return false;
        }

        reduced = SparseMatrix(nc, nc);
        for (std::size_t i = 0; i < nc; ++i) {
            auto& row = rows[i];
            std::sort(row.begin(), row.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });
            for (std::size_t k = 0; k < row.size();) {
                const std::size_t column = row[k].first;
                double value = 0.0;
                while (k < row.size() && row[k].first == column)
                    value += row[k++].second;
                if (!std::isfinite(value)) {
                    last_error_ = "non-finite coefficient in MGR reduced operator";
                    return false;
                }
                if (value != 0.0)
                    reduced.push_back(i, column, value);
            }
        }
        reduced.finalize();
        return reduced.is_consistent();
    }

    static bool same_pattern(const SparseMatrix& a, const SparseMatrix& b) {
        if (a.n_rows() != b.n_rows() || a.n_cols() != b.n_cols() ||
            a.nnz() != b.nnz())
            return false;
        for (std::size_t i = 0; i <= a.n_rows(); ++i)
            if (a.row_offsets_data()[i] != b.row_offsets_data()[i]) return false;
        for (std::size_t k = 0; k < a.nnz(); ++k)
            if (a.columns_data()[k] != b.columns_data()[k]) return false;
        return true;
    }

    void clear_state_only() {
        ready_ = false;
        coarse_mgr_.reset();
        coarse_operator_ = SparseMatrix();
        fine_operator_ = SparseMatrix();
        fine_inv_diag_.clear();
        fine_ilu_.reset();
        coarse_to_fine_.clear();
        fine_to_coarse_.clear();
    }

    void clear() {
        clear_state_only();
        setup_count_ = 0;
        numeric_updates_ = 0;
    }

    std::vector<std::size_t> fine_;
    std::vector<std::size_t> coarse_;
    std::size_t n_{0};
    SparseMatrix coarse_operator_;
    SparseMatrix fine_operator_;
    std::vector<double> fine_inv_diag_;
    FineRelaxation fine_relaxation_{FineRelaxation::Diagonal};
    std::vector<ReductionPartition> additional_levels_;
    std::unique_ptr<ILU0Preconditioner> fine_ilu_;
    std::unique_ptr<NativeMGRPreconditioner> coarse_mgr_;
    std::vector<std::vector<std::pair<std::size_t,double>>> coarse_to_fine_;
    std::vector<std::vector<std::pair<std::size_t,double>>> fine_to_coarse_;
    NativeBoomerAMGPreconditioner coarse_amg_;
    bool ready_{false};
    std::size_t setup_count_{0};
    std::size_t numeric_updates_{0};
    mutable std::string last_error_;
};

} // namespace cfdx::core
