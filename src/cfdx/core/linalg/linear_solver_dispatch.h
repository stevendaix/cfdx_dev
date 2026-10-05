#pragma once

#include "cfdx/core/linalg/advanced_preconditioners.h"
#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/linear_solver_models.h"
#include "cfdx/core/linalg/null_space.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

namespace cfdx::core {

struct LinearSolveReport {
    SolverResult result;
    LinearSolverPlan plan;
    MatrixScaling scaling = MatrixScaling::None;
    double physical_residual_relative = 0.0;
};

inline double relative_physical_residual(
    const SparseMatrix& matrix,
    const Vector& solution,
    const Vector& rhs) {
    const auto ax = matrix.matvec(solution);
    double residual_squared = 0.0;
    double rhs_squared = 0.0;
    for (std::size_t i = 0; i < ax.size(); ++i) {
        const double residual = ax[i] - rhs(i);
        residual_squared += residual * residual;
        rhs_squared += rhs(i) * rhs(i);
    }
    return std::sqrt(residual_squared) /
           std::max(std::sqrt(rhs_squared), 1.0e-30);
}

inline std::unique_ptr<Preconditioner> make_scalar_preconditioner(
    PreconditionerModel model,
    bool constant_null_space = false) {
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
        case PreconditionerModel::NativeAMG:
            return std::make_unique<NativeBoomerAMGPreconditioner>(
                constant_null_space);
        case PreconditionerModel::SmoothedAggregationAMG:
            return std::make_unique<
                NativeSmoothedAggregationAMGPreconditioner>(
                    constant_null_space);
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
    report.plan = select_linear_solver(problem, matrix.n_rows(), request);
    report.scaling = request.scaling;

    if (request.scaling != MatrixScaling::None &&
        report.plan.null_space == NullSpaceModel::Constant) {
        throw std::invalid_argument(
            "matrix scaling with a constant null-space projection is not yet supported");
    }

    const SparseMatrix* solve_matrix = &matrix;
    const Vector* solve_rhs = &rhs;
    Vector scaled_rhs;
    Vector physical_solution;
    MatrixScalingResult scaled;

    if (request.scaling != MatrixScaling::None) {
        scaled = scale_matrix(matrix, request.scaling);
        if (scaled.zero_rows != 0 || scaled.zero_columns != 0) {
            report.result.status = SolverStatus::NOT_APPLICABLE;
            report.physical_residual_relative =
                std::numeric_limits<double>::infinity();
            return report;
        }

        scaled_rhs = Vector(rhs.size());
        for (std::size_t i = 0; i < rhs.size(); ++i)
            scaled_rhs(i) = scaled.row_scale[i] * rhs(i);

        physical_solution = Vector(solution.size());
        solve_matrix = &scaled.matrix;
        solve_rhs = &scaled_rhs;
    }

    std::unique_ptr<Preconditioner> preconditioner;
    if (problem == LinearProblemKind::CoupledPressureVelocity) {
        if (solve_matrix->n_rows() % 4 != 0)
            throw std::invalid_argument("coupled pressure-velocity system must have 4N unknowns");
        if (report.plan.preconditioner == PreconditionerModel::CoupledBlockSchur)
            preconditioner = std::make_unique<CoupledBlockSchurAMGPreconditioner>(
                solve_matrix->n_rows() / 4);
        else if (report.plan.preconditioner == PreconditionerModel::PCD)
            throw std::invalid_argument(
                "PCD requires explicit pressure operators; use the coupled physics solver API");
        else
            throw std::invalid_argument("coupled dispatch requires an explicit coupled block preconditioner");
    } else {
        preconditioner = make_scalar_preconditioner(
            report.plan.preconditioner,
            report.plan.null_space == NullSpaceModel::Constant);
    }
    Preconditioner* pc = preconditioner.get();
    std::unique_ptr<NullSpaceProjector> null_space;
    if (report.plan.null_space == NullSpaceModel::Constant)
        null_space = std::make_unique<NullSpaceProjector>(
            NullSpaceProjector::constant(solve_matrix->n_rows()));

    switch (report.plan.krylov) {
        case KrylovModel::CG:
            if (null_space) {
                report.result = pc
                    ? solve_cg(*solve_matrix, *solve_rhs, solution, *pc, *null_space,
                               max_iterations, tolerance)
                    : solve_cg(*solve_matrix, *solve_rhs, solution, *null_space,
                               max_iterations, tolerance);
            } else {
                report.result = pc
                    ? solve_cg(matrix, rhs, solution, *pc, max_iterations, tolerance)
                    : solve_cg(matrix, rhs, solution, max_iterations, tolerance);
            }
            break;
        case KrylovModel::BiCGStab:
            report.result = solve_bicgstab(
                *solve_matrix, *solve_rhs, solution, max_iterations, tolerance, pc);
            break;
        case KrylovModel::GMRES:
            report.result = solve_gmres(
                *solve_matrix, *solve_rhs, solution, request.gmres_restart,
                max_iterations, tolerance, pc);
            break;
        case KrylovModel::FGMRES:
            report.result = solve_fgmres(
                *solve_matrix, *solve_rhs, solution, request.gmres_restart,
                max_iterations, tolerance, pc);
            break;
        case KrylovModel::Auto:
        case KrylovModel::LGMRES:
        case KrylovModel::PipeCG:
        case KrylovModel::PipeFGMRES:
            report.result.status = SolverStatus::NOT_APPLICABLE;
            break;
    }
    if (request.scaling != MatrixScaling::None) {
        for (std::size_t i = 0; i < solution.size(); ++i)
            physical_solution(i) = scaled.column_scale[i] * solution(i);
        solution = physical_solution;
    }

    report.physical_residual_relative =
        relative_physical_residual(matrix, solution, rhs);

    if (report.result.status == SolverStatus::CONVERGED &&
        report.physical_residual_relative > tolerance) {
        report.result.status = SolverStatus::MAX_ITER_REACHED;
    }
    return report;
}

} // namespace cfdx::core
