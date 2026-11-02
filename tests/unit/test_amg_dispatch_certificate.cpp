#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "common/test_harness.h"

#include <cmath>
#include <cstddef>

using namespace cfdx::core;
using namespace cfdx::testing;

static SparseMatrix make_matrix(bool negative_diagonal) {
    SparseMatrix A(18, 18);
    for (std::size_t i = 0; i < 18; ++i) {
        A.push_back(i, i, negative_diagonal && i == 9 ? -2.0 : 2.0);
        if (i > 0) A.push_back(i, i - 1, -1.0);
        if (i + 1 < 18) A.push_back(i, i + 1, -1.0);
    }
    A.finalize();
    return A;
}

int main() {
    run_case("cg_amg_dispatch_requires_spd_smoother_certificate", [] {
        auto A = make_matrix(true);
        Vector rhs(18, 1.0), solution(18, 0.0);
        LinearSolverRequest request;
        request.krylov = KrylovModel::CG;
        request.preconditioner = PreconditionerModel::NativeAMG;
        const auto report = solve_linear_system(
            A, rhs, solution, LinearProblemKind::PressurePoisson,
            request, 100, 1e-10);
        EXPECT_TRUE(report.result.status == SolverStatus::NOT_APPLICABLE);
        EXPECT_TRUE(report.result.iterations == 0);
    });
    run_case("gmres_amg_dispatch_uses_relaxed_smoother_policy", [] {
        auto A = make_matrix(true);
        Vector rhs(18, 1.0), solution(18, 0.0);
        LinearSolverRequest request;
        request.krylov = KrylovModel::GMRES;
        request.preconditioner = PreconditionerModel::NativeAMG;
        const auto report = solve_linear_system(
            A, rhs, solution, LinearProblemKind::PressurePoisson,
            request, 100, 1e-8);
        EXPECT_TRUE(report.result.status != SolverStatus::NOT_APPLICABLE);
        EXPECT_TRUE(std::isfinite(report.physical_residual_relative));
    });
    run_case("cg_amg_dispatch_still_solves_certified_spd", [] {
        auto A = make_matrix(false);
        Vector rhs(18, 1.0), solution(18, 0.0);
        LinearSolverRequest request;
        request.krylov = KrylovModel::CG;
        request.preconditioner = PreconditionerModel::NativeAMG;
        const auto report = solve_linear_system(
            A, rhs, solution, LinearProblemKind::PressurePoisson,
            request, 500, 1e-9);
        EXPECT_TRUE(report.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(report.physical_residual_relative <= 1e-9);
    });
    return run_all();
}
