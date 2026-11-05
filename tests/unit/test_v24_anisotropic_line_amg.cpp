#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "common/test_harness.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>

using namespace cfdx::core;
using namespace cfdx::testing;

static SparseMatrix anisotropic_spd_matrix(std::size_t nx, std::size_t ny)
{
    SparseMatrix A(nx*ny,nx*ny);
    constexpr double ax=1000.0, ay=1.0;
    for (std::size_t y=0;y<ny;++y) for (std::size_t x=0;x<nx;++x) {
        const auto row=y*nx+x;
        double diagonal=0.1;
        if (x) { A.push_back(row,row-1,-ax); diagonal+=ax; }
        if (x+1<nx) { A.push_back(row,row+1,-ax); diagonal+=ax; }
        if (y) { A.push_back(row,row-nx,-ay); diagonal+=ay; }
        if (y+1<ny) { A.push_back(row,row+nx,-ay); diagonal+=ay; }
        A.push_back(row,row,diagonal);
    }
    A.finalize();
    return A;
}

int main()
{
    run_case("v24_line_amg_auto_krylov_selects_fgmres", [] {
        LinearSolverRequest req;
        req.preconditioner=PreconditionerModel::AnisotropicLineAMG;
        const auto plan=select_linear_solver(LinearProblemKind::PressurePoisson,128,req);
        EXPECT_TRUE(plan.krylov==KrylovModel::FGMRES);
        EXPECT_TRUE(plan.preconditioner==PreconditionerModel::AnisotropicLineAMG);
    });
    run_case("v24_line_amg_forbids_explicit_cg", [] {
        LinearSolverRequest req;
        req.krylov=KrylovModel::CG;
        req.preconditioner=PreconditionerModel::AnisotropicLineAMG;
        EXPECT_THROW(select_linear_solver(LinearProblemKind::PressurePoisson,128,req),std::invalid_argument);
    });
    run_case("v24_line_amg_true_residual_anisotropic", [] {
        const auto A=anisotropic_spd_matrix(16,12);
        Vector rhs(A.n_rows(),0.0), solution(A.n_rows(),0.0);
        for (std::size_t i=0;i<rhs.size();++i)
            rhs(i)=1.0+0.3*std::sin(0.3*static_cast<double>(i));
        LinearSolverRequest req;
        req.preconditioner=PreconditionerModel::AnisotropicLineAMG;
        req.gmres_restart=80;
        const auto report=solve_linear_system(
            A,rhs,solution,LinearProblemKind::PressurePoisson,req,200,1e-9);
        EXPECT_TRUE(report.result.status==SolverStatus::CONVERGED);
        EXPECT_TRUE(std::isfinite(report.physical_residual_relative));
        EXPECT_TRUE(report.physical_residual_relative<=1e-9);
    });
    return run_all();
}
