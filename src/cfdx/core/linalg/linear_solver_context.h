#pragma once

#include "cfdx/core/linalg/gmres_solver.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace cfdx::core {

struct LinearSolverContextStats {
    std::size_t full_setups = 0;
    std::size_t numeric_updates = 0;
    std::size_t unchanged_reuses = 0;
    std::size_t solves = 0;
};

// Owns the reusable state around repeated GMRES solves. Matrix structure and
// values are compared explicitly because SparseMatrix currently permits direct
// coefficient access and therefore cannot provide reliable mutation epochs.
class ReusableGmresContext {
public:
    explicit ReusableGmresContext(Preconditioner& preconditioner)
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
            if (ready) ++stats_.numeric_updates;
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

    SolverResult solve(const SparseMatrix& matrix,
                       const Vector& rhs,
                       Vector& solution,
                       int restart = 30,
                       std::size_t max_iterations = 1000,
                       double tolerance = 1e-12,
                       KrylovControls controls = {},
                       PrecisionPolicy precision = {}) {
        if (!prepare(matrix)) {
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

        ++stats_.solves;
        return solve_gmres(op, rhs, solution, restart, max_iterations,
                           tolerance, &preconditioner_, controls, &workspace_);
    }

    void invalidate() {
        rows_.clear();
        columns_.clear();
        values_.clear();
    }

    const LinearSolverContextStats& stats() const noexcept { return stats_; }
    std::size_t workspace_reallocations() const noexcept {
        return workspace_.reallocations;
    }

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
    GmresWorkspace workspace_;
    LinearSolverContextStats stats_;
    std::vector<std::uint32_t> rows_;
    std::vector<std::uint32_t> columns_;
    std::vector<double> values_;
};

} // namespace cfdx::core
