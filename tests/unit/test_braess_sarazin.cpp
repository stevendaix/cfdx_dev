// Algebraic Braess-Sarazin Schur identity + original-system FGMRES certificate.
#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "common/test_harness.h"
#include <algorithm>
#include <cmath>
using namespace cfdx::core;
using namespace cfdx::testing;
namespace {
SparseMatrix coupled(std::size_t nc,double stretch) {
    SparseMatrix A(4*nc,4*nc);
    const auto nv=3*nc;
    for(std::size_t i=0;i<nc;++i) {
        const double k=(i%2?1.0:stretch);
        A.push_back(i,i,6*k);
        A.push_back(i,nc+i,0.2);
        A.push_back(nc+i,i,0.15);
        A.push_back(nc+i,nc+i,5*k);
        A.push_back(2*nc+i,2*nc+i,4*k);
        if(i>0) A.push_back(i,i-1,-0.2);
        if(i+1<nc) A.push_back(i,i+1,-0.2);
        A.push_back(i,nv+i,1.0);
        if(i+1<nc) A.push_back(i,nv+i+1,-0.35);
        if(i==0) A.push_back(nv,nv,1.0);
        else {
            A.push_back(nv+i,i,-1.0);
            A.push_back(nv+i,i-1,1.0);
            A.push_back(nv+i,nv+i,0.5);
            A.push_back(nv+i,nv+i-1,-0.12);
        }
    }
    A.finalize(); return A;
}
double a(const SparseMatrix& A,std::size_t r,std::size_t c) {
    double total=0.0;
    for(std::size_t k=A.row_offsets_data()[r]; k<A.row_offsets_data()[r+1]; ++k)
        if(A.columns_data()[k]==c) total+=A.values_data()[k];
    return total;
}
}
int main() {
    run_case("braess_diagonal_inverse_matches_schur_definition", [] {
        constexpr std::size_t nc=14;
        const auto A=coupled(nc,1000.0);
        for (double alpha : {0.5, 1.0, 1.8}) {
            CoupledBlockSchurOptions options;
            options.velocity_approximation=CoupledSchurVelocityApproximation::Diagonal;
            options.velocity_inverse_scale=1.0/alpha;
            CoupledBlockSchurAMGPreconditioner pc(nc,options);
            EXPECT_TRUE(pc.setup(A));
            const auto& S=pc.schur_matrix();
            for(std::size_t i=0;i<nc;++i) for(std::size_t j=0;j<nc;++j) {
                double reference=a(A,3*nc+i,3*nc+j);
                for(std::size_t u=0;u<3*nc;++u) {
                    reference-=a(A,3*nc+i,u)*(1.0/alpha)/a(A,u,u)*a(A,u,3*nc+j);
                }
                EXPECT_NEAR(a(S,i,j),reference,1e-12*std::max(1.0,std::abs(reference)));
            }
        }
    });
    run_case("braess_fgmres_converges_true_physical_residual", [] {
        for(double alpha : {0.75,1.0,1.5}) for(double anisotropy : {1.0,1000.0}) {
            const auto A=coupled(25,anisotropy);
            Vector exact(A.n_rows(),0.0), rhs(A.n_rows(),0.0), x(A.n_rows(),0.0);
            for(std::size_t i=0;i<A.n_rows();++i)
                exact(i)=0.3+0.5*std::sin(0.1*static_cast<double>(i));
            A.matvec_into(exact,rhs);
            LinearSolverRequest req;
            req.krylov=KrylovModel::FGMRES;
            req.preconditioner=PreconditionerModel::BraessSarazin;
            req.braess_sarazin_alpha=alpha;
            req.gmres_restart=50;
            const auto result=solve_linear_system(A,rhs,x,
                LinearProblemKind::CoupledPressureVelocity,req,400,1e-11);
            EXPECT_TRUE(result.result.status==SolverStatus::CONVERGED);
            EXPECT_TRUE(physical_relative_residual(A,x,rhs)<=1e-11);
            double error2=0.0, norm2=0.0;
            for(std::size_t i=0;i<x.size();++i) {
                error2=std::hypot(error2,x(i)-exact(i));
                norm2=std::hypot(norm2,exact(i));
            }
            EXPECT_TRUE(error2/norm2 < 1e-7);
        }
    });
    run_case("braess_rejects_non_coupled_cg_and_unbounded_relaxation", [] {
        const auto A=coupled(10,1.0);
        LinearSolverRequest req;
        req.preconditioner=PreconditionerModel::BraessSarazin;
        req.braess_sarazin_alpha=0.0;
        EXPECT_THROW(select_linear_solver(LinearProblemKind::CoupledPressureVelocity,
            measure_matrix_characteristics(A,true),req),std::invalid_argument);
        req.braess_sarazin_alpha=1.0;
        EXPECT_THROW(select_linear_solver(LinearProblemKind::PressurePoisson,
            measure_matrix_characteristics(A),req),std::invalid_argument);
        req.krylov=KrylovModel::CG;
        EXPECT_THROW(select_linear_solver(LinearProblemKind::CoupledPressureVelocity,
            measure_matrix_characteristics(A,true),req),std::invalid_argument);
    });
    return run_all();
}
