#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "cfdx/core/linalg/vanka_patch_preconditioner.h"
#include "common/test_harness.h"

#include <cmath>
#include <stdexcept>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
SparseMatrix coupled(std::size_t nc, double anisotropy=1.0) {
    SparseMatrix A(4*nc,4*nc);
    const std::size_t nv=3*nc;
    for(std::size_t i=0;i<nc;++i) {
        const double scale = i%2==0 ? anisotropy : 1.0;
        A.push_back(i,i,5.0*scale);
        if(i) A.push_back(i,i-1,-0.25);
        if(i+1<nc) A.push_back(i,i+1,-0.25);
        A.push_back(i,nv+i,1.0);
        if(i+1<nc) A.push_back(i,nv+i+1,-1.0);
        A.push_back(nc+i,nc+i,4.0*scale);
        A.push_back(2*nc+i,2*nc+i,4.0*scale);
        if(i==0) {
            A.push_back(nv, nv, 1.0);
        } else {
            A.push_back(nv+i,i,-1.0);
            A.push_back(nv+i,i-1,1.0);
            A.push_back(nv+i,nv+i,0.1);
        }
    }
    A.finalize();
    return A;
}
Vector exact_solution(std::size_t n) {
    Vector x(n,0.0);
    for(std::size_t i=0;i<n;++i)
        x(i)=0.7+std::sin(0.091*static_cast<double>(i));
    return x;
}
Vector matmul(const SparseMatrix& A,const Vector& x) {
    Vector b(x.size(),0.0);
    A.matvec_into(x,b);
    return b;
}
}

int main() {
    run_case("pressure_patch_builds_overlapping_G_D_stencils", [] {
        const auto A=coupled(16,1000.0);
        VankaPressurePatchPreconditioner pc(16);
        EXPECT_TRUE(pc.setup(A));
        EXPECT_TRUE(pc.patch_count()==16);
        EXPECT_TRUE(pc.largest_patch()>=3);
        EXPECT_TRUE(pc.largest_patch()<=24);
        Vector r(A.n_rows(),1.0), x(A.n_rows(),0.0), alias=r;
        EXPECT_TRUE(pc.apply(r,x));
        EXPECT_TRUE(pc.apply(alias,alias));
        for(std::size_t i=0;i<x.size();++i)
            EXPECT_NEAR(x(i),alias(i),1e-12);
        EXPECT_TRUE(x.is_valid());
    });
    run_case("rejected_setup_and_numeric_update_preserve_active_patches", [] {
        const auto A=coupled(12,100.0);
        VankaPressurePatchPreconditioner pc(12);
        EXPECT_TRUE(pc.setup(A));
        Vector r(A.n_rows(),1.0), before(A.n_rows(),0.0), after(A.n_rows(),0.0);
        EXPECT_TRUE(pc.apply(r,before));
        EXPECT_FALSE(pc.setup(SparseMatrix{}));
        EXPECT_TRUE(pc.update_values(coupled(12,100.0))); // identical graph
        EXPECT_FALSE(pc.update_values(coupled(11,100.0)));
        EXPECT_TRUE(pc.apply(r,after));
        for(std::size_t i=0;i<before.size();++i)
            EXPECT_NEAR(before(i),after(i),1e-12);
    });
    run_case("vanka_coupled_FGMRES_true_residual_and_solution", [] {
        for(double stretch : {1.0,1e3}) {
            const auto A=coupled(20,stretch);
            const auto exact=exact_solution(A.n_rows());
            const auto b=matmul(A,exact);
            Vector x(A.n_rows(),0.0);
            LinearSolverRequest req;
            req.preconditioner=PreconditionerModel::VankaPressurePatch;
            req.gmres_restart=40;
            const auto report=solve_linear_system(A,b,x,
                LinearProblemKind::CoupledPressureVelocity,req,240,1e-12);
            EXPECT_TRUE(report.plan.krylov==KrylovModel::FGMRES);
            EXPECT_TRUE(report.result.status==SolverStatus::CONVERGED);
            EXPECT_TRUE(physical_relative_residual(A,x,b)<1e-12);
            for(std::size_t i=0;i<x.size();++i)
                EXPECT_NEAR(x(i),exact(i),2e-7);
        }
    });
    run_case("vanka_rejects_scalar_and_CG_dispatch", [] {
        const auto A=coupled(8);
        LinearSolverRequest req;
        req.preconditioner=PreconditionerModel::VankaPressurePatch;
        EXPECT_THROW(select_linear_solver(LinearProblemKind::PressurePoisson,
            measure_matrix_characteristics(A),req),std::invalid_argument);
        req.krylov=KrylovModel::CG;
        EXPECT_THROW(select_linear_solver(LinearProblemKind::CoupledPressureVelocity,
            measure_matrix_characteristics(A,true),req),std::invalid_argument);
    });
    run_case("singular_local_pivot_fails_closed", [] {
        // A completely unconstrained and disconnected pressure equation
        // must not yield a silently regularized inverse.
        SparseMatrix A(4,4);
        A.push_back(0,0,1); A.push_back(1,1,1); A.push_back(2,2,1);
        A.finalize();
        VankaPressurePatchPreconditioner pc(1);
        EXPECT_FALSE(pc.setup(A));
    });
    return run_all();
}