#pragma once

#include "cfdx/core/linalg/advanced_preconditioners.h"
#include "cfdx/core/linalg/algebraic_line_preconditioner.h"
#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/linear_solver_models.h"
#include "cfdx/core/linalg/matrix_diagnostics.h"
#include "cfdx/core/linalg/null_space.h"
#include "cfdx/core/linalg/projected_gmres.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

namespace cfdx::core {

struct LinearSolveReport {
    SolverResult result;
    LinearSolverPlan plan;
    double physical_residual_relative = std::numeric_limits<double>::quiet_NaN();
    LinearScalingModel applied_scaling = LinearScalingModel::None;
};

inline double physical_residual_norm(
    const SparseMatrix& matrix,
    const Vector& solution,
    const Vector& rhs) {
    const auto ax = matrix.matvec(solution);
    double norm = 0.0;
    for (std::size_t i = 0; i < ax.size(); ++i)
        norm = std::hypot(norm, ax[i] - rhs(i));
    return norm;
}

inline double physical_relative_residual(
    const SparseMatrix& matrix,
    const Vector& solution,
    const Vector& rhs) {
    return physical_residual_norm(matrix, solution, rhs) /
        std::max(rhs.norm2(), 1e-300);
}

inline MatrixScaling to_matrix_scaling(LinearScalingModel model) {
    switch (model) {
        case LinearScalingModel::None: return MatrixScaling::None;
        case LinearScalingModel::Row: return MatrixScaling::Row;
        case LinearScalingModel::Column: return MatrixScaling::Column;
        case LinearScalingModel::RowColumn: return MatrixScaling::RowColumn;
        case LinearScalingModel::SymmetricDiagonal: return MatrixScaling::SymmetricDiagonal;
    }
    throw std::invalid_argument("unknown linear scaling model");
}

inline std::unique_ptr<Preconditioner> make_scalar_preconditioner(
    PreconditionerModel model,
    bool constant_null_space = false,
    bool require_smoother_certificate = true) {
    switch (model) {
        case PreconditionerModel::None:
            // CG's legacy null-preconditioner overload installs Jacobi.
            // An explicit identity object preserves the requested "none"
            // semantics consistently for every Krylov method.
            return std::make_unique<IdentityPreconditioner>();
        case PreconditionerModel::Jacobi:
            return std::make_unique<JacobiPreconditioner>();
        case PreconditionerModel::GaussSeidel:
            return std::make_unique<GaussSeidelPreconditioner>();
        case PreconditionerModel::ILU0:
            return std::make_unique<ILU0Preconditioner>();
        case PreconditionerModel::AnisotropicLine:
            return std::make_unique<AlgebraicLinePreconditioner>();
        case PreconditionerModel::NativeAMG:
            return std::make_unique<NativeBoomerAMGPreconditioner>(
                constant_null_space, require_smoother_certificate);
        case PreconditionerModel::AnisotropicBoomerAMG:
            return std::make_unique<NativeAnisotropicBoomerAMGPreconditioner>(
                constant_null_space);
        case PreconditionerModel::AnisotropicSmoothedAggregationAMG:
            return std::make_unique<NativeAnisotropicSmoothedAggregationAMGPreconditioner>(
                constant_null_space);
        case PreconditionerModel::SmoothedAggregationAMG:
            return std::make_unique<
                NativeSmoothedAggregationAMGPreconditioner>(
                    constant_null_space, require_smoother_certificate);
        case PreconditionerModel::NativeFieldSplit:
        case PreconditionerModel::CoupledBlockSchur:
            throw std::invalid_argument(
                "block preconditioning needs field metadata; use the coupled solver API");
        case PreconditionerModel::Auto:
        case PreconditionerModel::ILUT:
        case PreconditionerModel::FSAI:
        case PreconditionerModel::RAS:
        case PreconditionerModel::LSC:
        case PreconditionerModel::MGR:
        case PreconditionerModel::PCD:
            throw std::invalid_argument(
                "PCD requires the coupled pressure-velocity API with explicit pressure operators");
    }
    throw std::invalid_argument("unknown preconditioner model");
}

inline LinearSolveReport solve_linear_system(
    const SparseMatrix& matrix,
    const Vector& rhs,
    Vector& solution,
    LinearProblemKind problem,
    const LinearSolverRequest& request = {},
    std::size_t max_iterations = 1000,
    double tolerance = 1e-12) {
    LinearSolveReport report;
    // Route against actual CSR coefficients. An explicitly forced CG solve
    // with a non-SPD matrix fails closed without mutating the solution; keep
    // the longstanding NOT_APPLICABLE runtime contract for this case.
    const auto characteristics = measure_matrix_characteristics(
        matrix, problem == LinearProblemKind::CoupledPressureVelocity);
    // For b=0, the meaningful relative reference is the initial true
    // residual ||A*x_initial||. Using ||b||=0 or skipping certification
    // would either report infinity or accept an unverified homogeneous solve.
    const double physical_rhs_norm = rhs.norm2();
    const double homogeneous_initial_norm = physical_rhs_norm == 0.0
        ? physical_residual_norm(matrix, solution, rhs) : 0.0;
    if (request.krylov == KrylovModel::CG &&
        (!characteristics.numerically_symmetric ||
         !characteristics.positive_diagonal)) {
        report.plan = select_linear_solver(problem, matrix.n_rows(), request);
        report.plan.reason = "CG rejected: measured matrix is not symmetric with a positive diagonal";
        report.result.status = SolverStatus::NOT_APPLICABLE;
        report.physical_residual_relative =
            physical_relative_residual(matrix, solution, rhs);
        return report;
    }
    report.plan = select_linear_solver(problem, characteristics, request);

    if (request.scaling != LinearScalingModel::None &&
        request.null_space != NullSpaceModel::None)
        throw std::invalid_argument(
            "explicit algebraic scaling is not currently combined with a projected null-space solve");

    if (problem == LinearProblemKind::CoupledPressureVelocity &&
        request.scaling != LinearScalingModel::None)
        throw std::invalid_argument(
            "explicit scaling is currently exposed for scalar linear systems only");

    const SparseMatrix* solve_matrix = &matrix;
    MatrixScalingResult scaling_result;
    Vector scaled_rhs = rhs;
    Vector scaled_solution = solution;
    std::vector<double> column_scale(matrix.n_cols(), 1.0);
    if (request.scaling != LinearScalingModel::None) {
        scaling_result = scale_matrix(matrix, to_matrix_scaling(request.scaling));
        if (scaling_result.zero_rows != 0 || scaling_result.zero_columns != 0) {
            report.applied_scaling = request.scaling;
            report.result.status = SolverStatus::NOT_APPLICABLE;
            report.physical_residual_relative =
                std::numeric_limits<double>::infinity();
            return report;
        }
        solve_matrix = &scaling_result.matrix;
        scaled_rhs = Vector(rhs.size(), 0.0);
        scaled_solution = Vector(solution.size(), 0.0);
        column_scale = scaling_result.column_scale;
        for (std::size_t i = 0; i < rhs.size(); ++i) {
            scaled_rhs(i) = scaling_result.row_scale[i] * rhs(i);
            scaled_solution(i) = solution(i) / scaling_result.column_scale[i];
        }
        report.applied_scaling = request.scaling;
    }

    // Equilibration changes the residual norm. Converging only to the
    // requested tolerance in scaled coordinates may fail that same tolerance
    // in physical coordinates. Request one decade more inner accuracy and
    // still gate exclusively on the original physical true residual.
    const double inner_tolerance = request.scaling == LinearScalingModel::None
        ? tolerance : 0.1 * tolerance;

    std::unique_ptr<Preconditioner> preconditioner;
    if (problem == LinearProblemKind::CoupledPressureVelocity) {
        if (matrix.n_rows() % 4 != 0)
            throw std::invalid_argument("coupled pressure-velocity system must have 4N unknowns");
        if (report.plan.preconditioner == PreconditionerModel::CoupledBlockSchur)
            preconditioner = std::make_unique<CoupledBlockSchurAMGPreconditioner>(
                matrix.n_rows() / 4);
        else if (report.plan.preconditioner == PreconditionerModel::PCD)
            throw std::invalid_argument(
                "PCD requires explicit pressure operators; use the coupled physics solver API");
        else
            throw std::invalid_argument("coupled dispatch requires an explicit coupled block preconditioner");
    } else {
        // The SPD/CG certification path keeps the fail-closed smoother
        // certificate; a general nonsymmetric Krylov solve (GMRES, FGMRES,
        // BiCGStab) uses the relaxed preconditioner policy whose outer
        // true-residual checks remain authoritative.
        preconditioner = make_scalar_preconditioner(
            report.plan.preconditioner,
            report.plan.null_space == NullSpaceModel::Constant,
            report.plan.krylov == KrylovModel::CG);
    }
    Preconditioner* pc = preconditioner.get();
    const Vector& solve_rhs =
        request.scaling == LinearScalingModel::None ? rhs : scaled_rhs;
    Vector& solve_solution =
        request.scaling == LinearScalingModel::None ? solution : scaled_solution;
    std::unique_ptr<NullSpaceProjector> null_space;
    if (report.plan.null_space == NullSpaceModel::Constant)
        null_space = std::make_unique<NullSpaceProjector>(
            NullSpaceProjector::constant(matrix.n_rows()));

    switch (report.plan.krylov) {
        case KrylovModel::CG:
            if (null_space) {
                report.result = pc
                    ? solve_cg(*solve_matrix, solve_rhs, solve_solution, *pc, *null_space,
                               max_iterations, inner_tolerance)
                    : solve_cg(*solve_matrix, solve_rhs, solve_solution, *null_space,
                               max_iterations, inner_tolerance);
            } else {
                report.result = pc
                    ? solve_cg(*solve_matrix, solve_rhs, solve_solution, *pc, max_iterations, inner_tolerance)
                    : solve_cg(*solve_matrix, solve_rhs, solve_solution, max_iterations, inner_tolerance);
            }
            break;
        case KrylovModel::BiCGStab:
            report.result = solve_bicgstab(
                *solve_matrix, solve_rhs, solve_solution, max_iterations, inner_tolerance, pc);
            break;
        case KrylovModel::GMRES:
        case KrylovModel::FGMRES:
            if (null_space) {
                // The projected right-preconditioned Arnoldi recurrence also
                // covers changing multilevel corrections (FGMRES semantics).
                report.result = solve_projected_gmres(
                    *solve_matrix, solve_rhs, solve_solution, *null_space,
                    request.gmres_restart, max_iterations, inner_tolerance, *pc);
            } else if (report.plan.krylov == KrylovModel::GMRES) {
                report.result = solve_gmres(
                    *solve_matrix, solve_rhs, solve_solution, request.gmres_restart,
                    max_iterations, inner_tolerance, pc);
            } else {
                report.result = solve_fgmres(
                    *solve_matrix, solve_rhs, solve_solution, request.gmres_restart,
                    max_iterations, inner_tolerance, pc);
            }
            break;
        case KrylovModel::Auto:
        case KrylovModel::LGMRES:
        case KrylovModel::PipeCG:
        case KrylovModel::PipeFGMRES:
            report.result.status = SolverStatus::NOT_APPLICABLE;
            break;
    }
    if (request.scaling != LinearScalingModel::None) {
        for (std::size_t i = 0; i < solution.size(); ++i)
            solution(i) = column_scale[i] * scaled_solution(i);
    }
    // Every dispatch path uses the ORIGINAL unscaled physical system for
    // its final convergence certificate, including the homogeneous b=0 case.
    const double final_residual = physical_residual_norm(matrix, solution, rhs);
    const double reference = physical_rhs_norm > 0.0
        ? physical_rhs_norm : homogeneous_initial_norm;
    report.physical_residual_relative = reference > 0.0
        ? final_residual / reference
        : (final_residual == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
    if (report.result.status == SolverStatus::CONVERGED &&
        (!std::isfinite(report.physical_residual_relative) ||
         report.physical_residual_relative > tolerance))
        report.result.status = SolverStatus::MAX_ITER_REACHED;
    return report;
}

} // namespace cfdx::core
