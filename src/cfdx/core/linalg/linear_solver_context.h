#pragma once

#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/gmres_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace cfdx::core {

struct LinearSolverContextStats {
    std::size_t full_setups = 0;
    std::size_t numeric_updates = 0;
    std::size_t unchanged_reuses = 0;
    std::size_t solves = 0;
};

// Tracks symbolic and numeric preconditioner state across repeated matrices.
// SparseMatrix permits direct coefficient access, so explicit snapshots are
// currently more reliable than mutation epochs.
class ReusablePreconditionerContext {
public:
    explicit ReusablePreconditionerContext(Preconditioner& preconditioner)
        : preconditioner_(preconditioner) {}

    bool prepare(const SparseMatrix& matrix) {
        if (matrix.n_rows() == 0 || matrix.n_rows() != matrix.n_cols() ||
            !matrix.is_consistent()) {
            return false;
        }

        const bool same_structure = structure_matches(matrix);
        const bool same_values = same_structure && values_match(matrix);
        bool ready = false;
        if (!same_structure) {
            ready = preconditioner_.setup(matrix);
            if (ready) ++stats_.full_setups;
        } else if (!same_values) {
            ready = preconditioner_.update_values(matrix);
            if (ready) {
                ++stats_.numeric_updates;
            } else {
                // This context owns lifecycle decisions. If a numeric refresh
                // is rejected (for example, a preconditioner can no longer
                // certify its retained hierarchy), rebuild from the current
                // matrix before failing the linear solve. This is explicit
                // full setup, not a silent solver/preconditioner substitution.
                ready = preconditioner_.setup(matrix);
                if (ready) ++stats_.full_setups;
            }
        } else {
            ready = true;
            ++stats_.unchanged_reuses;
        }
        if (!ready) {
            invalidate();
            return false;
        }
        snapshot(matrix);
        return true;
    }

    void invalidate() {
        rows_.clear();
        columns_.clear();
        values_.clear();
    }

    void record_solve() noexcept { ++stats_.solves; }
    Preconditioner& preconditioner() noexcept { return preconditioner_; }
    const LinearSolverContextStats& stats() const noexcept { return stats_; }

private:
    bool structure_matches(const SparseMatrix& matrix) const {
        return rows_.size() == matrix.n_rows() + 1 &&
               columns_.size() == matrix.nnz() &&
               std::equal(rows_.begin(), rows_.end(), matrix.row_offsets_data()) &&
               std::equal(columns_.begin(), columns_.end(), matrix.columns_data());
    }

    bool values_match(const SparseMatrix& matrix) const {
        return values_.size() == matrix.nnz() &&
               std::equal(values_.begin(), values_.end(), matrix.values_data());
    }

    void snapshot(const SparseMatrix& matrix) {
        rows_.assign(matrix.row_offsets_data(),
                     matrix.row_offsets_data() + matrix.n_rows() + 1);
        columns_.assign(matrix.columns_data(),
                        matrix.columns_data() + matrix.nnz());
        values_.assign(matrix.values_data(), matrix.values_data() + matrix.nnz());
    }

    Preconditioner& preconditioner_;
    LinearSolverContextStats stats_;
    std::vector<std::uint32_t> rows_;
    std::vector<std::uint32_t> columns_;
    std::vector<double> values_;
};

// Owns reusable state around repeated GMRES solves.
class ReusableGmresContext {
public:
    explicit ReusableGmresContext(Preconditioner& preconditioner)
        : state_(preconditioner) {}

    bool prepare(const SparseMatrix& matrix) { return state_.prepare(matrix); }

    SolverResult solve(const SparseMatrix& matrix,
                       const Vector& rhs,
                       Vector& solution,
                       int restart = 30,
                       std::size_t max_iterations = 1000,
                       double tolerance = 1e-12,
                       KrylovControls controls = {},
                       PrecisionPolicy precision = {}) {
        if (!state_.prepare(matrix)) {
            SolverResult result;
            result.status = SolverStatus::NOT_APPLICABLE;
            return result;
        }

        LinearOperator op;
        op.size = matrix.n_rows();
        op.apply = [&matrix, precision](const Vector& input, Vector& output) {
            if (output.size() != matrix.n_rows()) output.resize(matrix.n_rows());
            if (precision.enabled &&
                precision.operator_precision == SolverPrecision::FP32) {
                mixed_precision_matvec(matrix, input, output, SolverPrecision::FP32);
                return;
            }
            const auto* row = matrix.row_offsets_data();
            const auto* col = matrix.columns_data();
            const auto* val = matrix.values_data();
            for (std::size_t i = 0; i < matrix.n_rows(); ++i) {
                double sum = 0.0;
                for (std::size_t k = row[i]; k < row[i + 1]; ++k)
                    sum += val[k] * input(col[k]);
                output(i) = sum;
            }
        };

        state_.record_solve();
        return solve_gmres(
            op, rhs, solution, restart, max_iterations, tolerance,
            &state_.preconditioner(), controls, &workspace_);
    }

    void invalidate() { state_.invalidate(); }
    const LinearSolverContextStats& stats() const noexcept {
        return state_.stats();
    }
    std::size_t workspace_reallocations() const noexcept {
        return workspace_.reallocations;
    }

private:
    ReusablePreconditionerContext state_;
    GmresWorkspace workspace_;
};

// Reuses preconditioner symbolic state for repeated SPD pressure/diffusion
// systems while CG retains its true-residual convergence checks.
class ReusableCgContext {
public:
    explicit ReusableCgContext(Preconditioner& preconditioner)
        : state_(preconditioner) {}

    bool prepare(const SparseMatrix& matrix) { return state_.prepare(matrix); }

    SolverResult solve(
        const SparseMatrix& matrix,
        const Vector& rhs,
        Vector& solution,
        std::size_t max_iterations = 1000,
        double tolerance = 1e-12,
        PrecisionPolicy precision = {},
        KrylovReductionPolicy reduction = {},
        const NullSpaceProjector* null_space = nullptr) {
        // Match solve_cg_impl's zero-iteration contract before asking the
        // reusable context to build/update a preconditioner. A converged
        // initial guess needs no Krylov direction and must not fail because
        // an otherwise-unused preconditioner cannot be prepared.
        if (matrix.n_rows() == matrix.n_cols() &&
            rhs.size() == matrix.n_rows() &&
            solution.size() == matrix.n_cols() && matrix.is_consistent()) {
            bool finite_input = true;
            for (std::size_t k = 0; k < matrix.nnz(); ++k)
                finite_input = finite_input &&
                    std::isfinite(matrix.values_data()[k]);
            for (std::size_t i = 0; i < rhs.size(); ++i)
                finite_input = finite_input && std::isfinite(rhs(i)) &&
                    std::isfinite(solution(i));
            const bool compatible_null_space = !null_space ||
                (null_space->dimension() == matrix.n_rows() &&
                 null_space->is_compatible(rhs) &&
                 null_space->is_null_space(matrix));
            if (finite_input && compatible_null_space) {
                if (null_space) null_space->remove(solution);
                Vector initial_residual(rhs.size());
                mixed_precision_true_residual(
                    matrix, rhs, solution, initial_residual);
                if (null_space) null_space->remove(initial_residual);
                const double rhs_norm = krylov_norm2(
                    rhs, SolverPrecision::FP64, reduction);
                const double residual_norm = krylov_norm2(
                    initial_residual, SolverPrecision::FP64, reduction);
                const double tolerance_absolute =
                    tolerance * std::max(rhs_norm, 1e-300);
                if (residual_norm <= tolerance_absolute) {
                    SolverResult result;
                    result.status = SolverStatus::CONVERGED;
                    result.iterations = 0;
                    result.residual = residual_norm;
                    result.residual_relative = rhs_norm > 0.0
                        ? residual_norm / rhs_norm : 0.0;
                    return result;
                }
            }
        }

        if (!state_.prepare(matrix)) {
            SolverResult result;
            result.status = SolverStatus::NOT_APPLICABLE;
            result.iterations = 0;
            if (matrix.n_rows() == matrix.n_cols() &&
                rhs.size() == matrix.n_rows() &&
                solution.size() == matrix.n_cols() && matrix.is_consistent()) {
                Vector failure_residual(rhs.size());
                mixed_precision_true_residual(
                    matrix, rhs, solution, failure_residual);
                result.residual = krylov_norm2(
                    failure_residual, SolverPrecision::FP64, reduction);
                const double rhs_norm = krylov_norm2(
                    rhs, SolverPrecision::FP64, reduction);
                result.residual_relative = rhs_norm > 0.0
                    ? result.residual / rhs_norm
                    : (result.residual == 0.0
                        ? 0.0 : std::numeric_limits<double>::infinity());
            }
            return result;
        }
        state_.record_solve();
        return solve_cg_prepared(
            matrix, rhs, solution, state_.preconditioner(),
            max_iterations, tolerance, precision, reduction, null_space);
    }

    void invalidate() { state_.invalidate(); }
    const LinearSolverContextStats& stats() const noexcept {
        return state_.stats();
    }

private:
    ReusablePreconditionerContext state_;
};

} // namespace cfdx::core
