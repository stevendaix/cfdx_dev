#pragma once

#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/null_space.h"
#include "cfdx/core/linalg/preconditioner.h"
#include "cfdx/core/linalg/sparse_matrix.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace cfdx::core {

// Right-preconditioned GMRES/FGMRES on the quotient space of a singular
// operator. The right and left null spaces must BOTH be described by `space`.
// In particular, this is not a general nonsymmetric/left-null-space solver.
// An incompatible RHS, incorrect null-space claim, or failed setup leaves
// the caller's iterate untouched. The original unprojected residual must
// also be certified by the dispatch entry point.
class ProjectedPreconditioner final : public Preconditioner {
public:
    ProjectedPreconditioner(Preconditioner& inner,
                            const NullSpaceProjector& space)
        : inner_(inner), space_(space) {}

    bool setup(const SparseMatrix& matrix) override {
        return inner_.setup(matrix);
    }

    bool apply(const Vector& residual, Vector& correction) const override {
        if (residual.size() != space_.dimension() || !residual.is_valid())
            return false;
        Vector projected = residual;
        space_.remove(projected);
        if (correction.size() != projected.size())
            correction.resize(projected.size());
        if (!inner_.apply(projected, correction) || !correction.is_valid())
            return false;
        space_.remove(correction);
        return correction.is_valid();
    }

    const char* name() const override { return "projected-preconditioner"; }

private:
    Preconditioner& inner_;
    const NullSpaceProjector& space_;
};

// Required for right preconditioning: a supplied right null mode must also
// annihilate A^T. The relative test is independent of the matrix unit scale.
inline bool is_bilateral_null_space(const SparseMatrix& matrix,
                                     const NullSpaceProjector& space,
                                     double tolerance = 1e-12) {
    if (matrix.n_rows() != matrix.n_cols() ||
        matrix.n_rows() != space.dimension() ||
        !std::isfinite(tolerance) || tolerance < 0.0 ||
        !space.is_null_space(matrix, tolerance))
        return false;

    const auto* row = matrix.row_offsets_data();
    const auto* col = matrix.columns_data();
    const auto* val = matrix.values_data();
    double frobenius = 0.0;
    for (std::size_t k = 0; k < matrix.nnz(); ++k) {
        if (!std::isfinite(val[k])) return false;
        frobenius = std::hypot(frobenius, val[k]);
    }
    for (const auto& mode : space.basis()) {
        Vector transpose_product(matrix.n_cols(), 0.0);
        for (std::size_t i = 0; i < matrix.n_rows(); ++i) {
            for (std::size_t k = row[i]; k < row[i + 1]; ++k)
                transpose_product(col[k]) += val[k] * mode(i);
        }
        if (transpose_product.norm2() >
            tolerance * std::max(frobenius, std::numeric_limits<double>::min()))
            return false;
    }
    return true;
}

inline SolverResult solve_projected_gmres(
    const SparseMatrix& matrix, const Vector& rhs, Vector& solution,
    const NullSpaceProjector& space, int restart, std::size_t max_iterations,
    double tolerance, Preconditioner& inner_preconditioner) {
    SolverResult failed;
    failed.status = SolverStatus::NOT_APPLICABLE;
    if (matrix.n_rows() == 0 || matrix.n_cols() != matrix.n_rows() ||
        rhs.size() != matrix.n_rows() || solution.size() != matrix.n_rows() ||
        !rhs.is_valid() || !solution.is_valid() ||
        restart <= 0 || max_iterations == 0 ||
        !(tolerance > 0.0) || !std::isfinite(tolerance) ||
        !space.is_compatible(rhs) || !is_bilateral_null_space(matrix, space))
        return failed;

    // Work on a copy so rejected inputs do not destroy an existing pressure
    // field. Projection removes the arbitrary gauge from the initial guess.
    Vector gauge_free_solution = solution;
    space.remove(gauge_free_solution);
    const auto ax0 = matrix.matvec(gauge_free_solution);
    double r0_norm = 0.0;
    for (std::size_t i = 0; i < rhs.size(); ++i)
        r0_norm = std::hypot(r0_norm, rhs(i) - ax0[i]);
    const double b_norm = rhs.norm2();
    if (r0_norm == 0.0 || (b_norm > 0.0 && r0_norm <= tolerance * b_norm)) {
        // Already satisfied in physical coordinates: do not build an AMG
        // hierarchy just to report zero Krylov iterations.
        SolverResult converged;
        converged.status = SolverStatus::CONVERGED;
        converged.residual = r0_norm;
        converged.residual_relative = b_norm > 0.0 ? r0_norm / b_norm : 0.0;
        solution = std::move(gauge_free_solution);
        return converged;
    }
    ProjectedPreconditioner projected_pc(inner_preconditioner, space);
    if (!projected_pc.setup(matrix)) return failed;

    LinearOperator projected_operator;
    projected_operator.size = matrix.n_rows();
    projected_operator.apply = [&matrix, &space](const Vector& input, Vector& output) {
        Vector gauge_free = input;
        space.remove(gauge_free);
        const auto product = matrix.matvec(gauge_free);
        if (output.size() != product.size()) output.resize(product.size());
        for (std::size_t i = 0; i < product.size(); ++i) output(i) = product[i];
        space.remove(output);
    };
    projected_operator.apply_fp64 = projected_operator.apply;

    const auto result = solve_fgmres(
        projected_operator, rhs, gauge_free_solution, restart,
        max_iterations, tolerance, &projected_pc);
    if (result.status != SolverStatus::NOT_APPLICABLE &&
        gauge_free_solution.is_valid()) {
        space.remove(gauge_free_solution);
        solution = std::move(gauge_free_solution);
    }
    return result;
}

} // namespace cfdx::core