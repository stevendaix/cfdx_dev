#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "common/test_harness.h"
#include <cmath>
using namespace cfdx::core;
using namespace cfdx::testing;
static SparseMatrix poisson(std::size_t n) {
    SparseMatrix a(n,n);
    for(std::size_t i=0;i<n;++i) {
        a.push_back(i,i,2.0);
        if (i) a.push_back(i,i-1,-1.0);
        if (i+1<n) a.push_back(i,i+1,-1.0);
    }
    a.finalize();return a;
}
int main() {
    run_case("v26_refinement_stays_off_by_default",[]{
        LinearSolverRequest req;
        EXPECT_TRUE(req.residual_refinement_steps==0);
        EXPECT_TRUE(req.scaling==LinearScalingModel::None);
    });
    run_case("v26_refinement_rejects_incompatible_requests",[]{
        auto a=poisson(20);Vector b(20,1.0),x(20,0.0);
        LinearSolverRequest req;req.krylov=KrylovModel::CG;
        req.residual_refinement_steps=1;
        EXPECT_THROW(solve_linear_system(a,b,x,LinearProblemKind::PressurePoisson,req,10,1e-9),std::invalid_argument);
        req.krylov=KrylovModel::FGMRES;
        req.scaling=LinearScalingModel::SymmetricDiagonal;
        EXPECT_THROW(solve_linear_system(a,b,x,LinearProblemKind::PressurePoisson,req,10,1e-9),std::invalid_argument);
        req.scaling=LinearScalingModel::None;req.residual_refinement_steps=9;
        EXPECT_THROW(solve_linear_system(a,b,x,LinearProblemKind::PressurePoisson,req,10,1e-9),std::invalid_argument);
    });
    run_case("v26_extended_residual_reconstructs_original_equation",[]{
        SparseMatrix a(1,1);a.push_back(0,0,1.0e12);a.push_back(0,0,-1.0e12);a.push_back(0,0,1.0);a.finalize();
        Vector b(1,1.),x(1,1.),r(1,0.0);
        EXPECT_TRUE(extended_precision_original_residual(a,x,b,r));
        EXPECT_NEAR(r(0),0.0,1.e-14);
        // Candidate correction direction must never be certified merely
        // because the componentwise/backward error is small.
        EXPECT_TRUE(physical_relative_residual(a,x,b)<=1.e-14);
    });
    run_case("v26_refinement_does_not_waive_original_residual_gate",[]{
        auto a=poisson(30);Vector b(30,1.0),x(30,0.0);
        LinearSolverRequest req;req.krylov=KrylovModel::FGMRES;
        req.preconditioner=PreconditionerModel::Jacobi;req.gmres_restart=8;
        req.residual_refinement_steps=2;
        auto result=solve_linear_system(a,b,x,LinearProblemKind::PressurePoisson,req,4,1e-14);
        auto r=physical_relative_residual(a,x,b);
        EXPECT_TRUE(std::isfinite(r));
        EXPECT_NEAR(result.physical_residual_relative,r,1e-13);
        EXPECT_TRUE(result.result.status!=SolverStatus::CONVERGED || r<=1e-14);
        EXPECT_TRUE(result.result.residual_relative<=1e-14 ||
                    result.result.status!=SolverStatus::CONVERGED);
    });
    run_case("v26_v18_nonsymmetric_general_operator_accepts_refinement",[]{
        SparseMatrix a(3,3);
        a.push_back(0,0,7.0);a.push_back(0,1,-2.0);
        a.push_back(1,0,-0.5);a.push_back(1,1,6.0);a.push_back(1,2,-1.0);
        a.push_back(2,1,-3.0);a.push_back(2,2,5.0);
        a.finalize();
        Vector b(3,1.0),x(3,0.0);
        LinearSolverRequest req;req.krylov=KrylovModel::FGMRES;
        req.preconditioner=PreconditionerModel::ILU0;
        req.residual_refinement_steps=2;
        // V18 sparse flux-pressure Jacobian is deliberately nonsymmetric;
        // the production path labels it General rather than PressurePoisson.
        const auto report=solve_linear_system(
            a,b,x,LinearProblemKind::General,req,3,1e-13);
        const double real=physical_relative_residual(a,x,b);
        EXPECT_TRUE(std::isfinite(real));
        EXPECT_NEAR(report.physical_residual_relative,real,1e-12);
        EXPECT_TRUE(report.result.status!=SolverStatus::CONVERGED || real<=1e-13);
    });
    run_case("v26_refinement_reduces_original_residual_from_capped_krylov",[]{
        auto a=poisson(30);Vector b(30,1.0),x0(30,0.0),x1(30,0.0);
        LinearSolverRequest req;req.krylov=KrylovModel::FGMRES;
        req.preconditioner=PreconditionerModel::Jacobi;req.gmres_restart=8;
        auto base=solve_linear_system(a,b,x0,LinearProblemKind::PressurePoisson,req,4,1e-9);
        req.residual_refinement_steps=4;
        auto improved=solve_linear_system(a,b,x1,LinearProblemKind::PressurePoisson,req,4,1e-9);
        const auto r0=physical_relative_residual(a,x0,b);
        const auto r1=physical_relative_residual(a,x1,b);
        EXPECT_TRUE(std::isfinite(r0)&&std::isfinite(r1));
        EXPECT_TRUE(r1<=r0);
        EXPECT_NEAR(improved.physical_residual_relative,r1,1e-12);
        EXPECT_TRUE(improved.result.status!=SolverStatus::CONVERGED || r1<=1e-9);
    });
    return run_all();
}
