#pragma once

#include "cfdx/core/linalg/advanced_preconditioners.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/linear_solver_models.h"

#include <memory>
#include <stdexcept>

namespace cfdx::core {

struct LinearSolveReport {
    SolverResult result;
    LinearSolverPlan plan;
};

inline std::unique_ptr<Preconditioner> make_scalar_preconditioner(
    PreconditionerModel model) {
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
            return std::make_unique<NativeBoomerAMGPreconditioner>();
        case PreconditionerModel::NativeFieldSplit:
        case PreconditionerModel::CoupledBlockSchur:
            throw std::invalid_argument(
                "block preconditioning needs field metadata; use the coupled solver API");
        case PreconditionerModel::Auto:
        case PreconditionerModel::ILUT:
        case PreconditionerModel::SmoothedAggregationAMG:
        case PreconditionerModel::FSAI:
        case PreconditionerModel::RAS:
        case PreconditionerModel::LSC:
        case PreconditionerModel::MGR:
            throw std::invalid_argument("preconditioner is not available in scalar dispatch");
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
    auto preconditioner = make_scalar_preconditioner(report.plan.preconditioner);
    Preconditioner* pc = preconditioner.get();

    switch (report.plan.krylov) {
        case KrylovModel::CG:
            report.result = pc
                ? solve_cg(matrix, rhs, solution, *pc, max_iterations, tolerance)
                : solve_cg(matrix, rhs, solution, max_iterations, tolerance);
            break;
        case KrylovModel::BiCGStab:
            report.result = solve_bicgstab(
                matrix, rhs, solution, max_iterations, tolerance, pc);
            break;
        case KrylovModel::GMRES:
        case KrylovModel::FGMRES:
            report.result = solve_gmres(
                matrix, rhs, solution, request.gmres_restart,
                max_iterations, tolerance, pc);
            break;
        case KrylovModel::Auto:
        case KrylovModel::LGMRES:
        case KrylovModel::PipeCG:
        case KrylovModel::PipeFGMRES:
            report.result.status = SolverStatus::NOT_APPLICABLE;
            break;
    }
    return report;
}

} // namespace cfdx::core
