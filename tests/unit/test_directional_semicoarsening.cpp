// Regression and physical-residual contracts for opt-in anisotropic multigrid.
// Directional semi-coarsening matches WEAK edges, preserving strong lines;
// line-GS performs true multiplicative forward/backward block solves.
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
SparseMatrix operator_matrix(std::size_t nx, std::size_t ny, double horizontal,
                             double vertical, double diagonal, double shift) {
    SparseMatrix A(nx * ny, nx * ny);
    for (std::size_t y = 0; y < ny; ++y)
        for (std::size_t x = 0; x < nx; ++x) {
            const auto i = y * nx + x;
            double diag = shift;
            const auto add = [&](std::size_t j, double a) {
                A.push_back(i, j, -a);
                diag += a;
            };
            if (x) add(i-1, horizontal);
            if (x+1<nx) add(i+1, horizontal);
            if (y) add(i-nx, vertical);
            if (y+1<ny) add(i+nx, vertical);
            if (x && y) add(i-nx-1, diagonal);
            if (x+1<nx && y+1<ny) add(i+nx+1, diagonal);
            if (x+1<nx && y) add(i-nx+1, diagonal);
            if (x && y+1<ny) add(i+nx-1, diagonal);
            A.push_back(i,i,diag);
        }
    A.finalize();
    return A;
}

Vector exact_solution(std::size_t n) {
    Vector x(n,0.0);
    for (std::size_t i=0;i<n;++i)
        x(i) = std::sin(0.29 * static_cast<double>(i)) +
               0.4 * std::cos(0.07 * static_cast<double>(i));
    return x;
}
Vector rhs_for(const SparseMatrix& A, const Vector& exact) {
    Vector rhs(A.n_rows(),0.0);
    A.matvec_into(exact,rhs);
    return rhs;
}
double true_relative(const SparseMatrix& A,const Vector& b,const Vector& x) {
    Vector Ax(b.size(),0.0);
    A.matvec_into(x,Ax);
    double r=0.0, bnorm=0.0;
    for(std::size_t i=0;i<b.size();++i) {
        r = std::hypot(r,b(i)-Ax(i));
        bnorm = std::hypot(bnorm,b(i));
    }
    return r/std::max(bnorm,1e-300);
}
double rel_solution_error(const Vector& x,const Vector& ref) {
    double err=0.0, norm=0.0;
    for(std::size_t i=0;i<x.size();++i) {
        err=std::hypot(err,x(i)-ref(i));
        norm=std::hypot(norm,ref(i));
    }
    return err/std::max(norm,1e-300);
}
}

int main() {
    run_case("multiplicative_line_GS_uses_updated_neighbour_lines", [] {
        const auto A=operator_matrix(16,12,1000.0,1.0,0.01,1.0);
        AlgebraicLinePreconditioner pc(32.0, 0.1, true);
        EXPECT_TRUE(pc.setup(A));
        EXPECT_TRUE(pc.longest_line()>=8);
        const auto exact=exact_solution(A.n_rows());
        const auto rhs=rhs_for(A,exact);
        Vector jacobi(rhs.size(),0.0), x(rhs.size(),0.0);
        EXPECT_TRUE(pc.apply(rhs,jacobi));
        std::vector<double> scratch;
        EXPECT_TRUE(pc.sweep(rhs,x,false,1.0,scratch));
        EXPECT_TRUE(pc.sweep(rhs,x,true,1.0,scratch));
        EXPECT_TRUE(true_relative(A,rhs,x) < true_relative(A,rhs,jacobi));
        EXPECT_TRUE(true_relative(A,rhs,x)<0.1);
    });
    run_case("directional_coarsening_preserves_strong_lines_and_constants", [] {
        const auto A=operator_matrix(24,16,10000.0,1.0,0.02,1.0);
        NativeDirectionalSemiCoarseningAMGPreconditioner pc;
        EXPECT_TRUE(pc.setup(A));
        EXPECT_TRUE(pc.finest_line_coverage() > A.n_rows()/2);
        EXPECT_TRUE(pc.finest_line_coupling_share() > 0.95);
        EXPECT_TRUE(pc.finest_directional_pairs() > A.n_rows()/5);
        EXPECT_NEAR(pc.prolongation_row_sum_min(),1.0,3e-13);
        EXPECT_NEAR(pc.prolongation_row_sum_max(),1.0,3e-13);
        EXPECT_TRUE(pc.operator_complexity()>1.0);
        EXPECT_TRUE(pc.grid_complexity()>1.0);
    });
    run_case("directional_amg_explicit_FGMRES_dispatch_and_CG_rejection", [] {
        const auto A=operator_matrix(16,12,1000.0,1.0,0.01,1.0);
        LinearSolverRequest req;
        req.preconditioner=PreconditionerModel::DirectionalSemiCoarseningAMG;
        const auto plan=select_linear_solver(
            LinearProblemKind::PressurePoisson,measure_matrix_characteristics(A),req);
        EXPECT_TRUE(plan.krylov==KrylovModel::FGMRES);
        req.krylov=KrylovModel::CG;
        EXPECT_THROW(select_linear_solver(
            LinearProblemKind::PressurePoisson,measure_matrix_characteristics(A),req),
            std::invalid_argument);
    });
    run_case("isotropic_fallback_keeps_hierarchy_and_true_residual", [] {
        const auto A=operator_matrix(12,12,1.0,1.0,0.01,1.0);
        NativeDirectionalSemiCoarseningAMGPreconditioner pc;
        EXPECT_TRUE(pc.setup(A));
        EXPECT_TRUE(pc.finest_directional_pairs()==0);
        const auto x=exact_solution(A.n_rows());
        const auto b=rhs_for(A,x);
        Vector z(A.n_rows(),0.0);
        EXPECT_TRUE(pc.apply(b,z));
        EXPECT_TRUE(std::isfinite(true_relative(A,b,z)));
    });
    run_case("transactional_rejected_hierarchy_keeps_previous_solution", [] {
        const auto A=operator_matrix(16,12,1000.0,1.0,0.01,1.0);
        NativeDirectionalSemiCoarseningAMGPreconditioner pc;
        EXPECT_TRUE(pc.setup(A));
        Vector b(A.n_rows(),1.0), old(A.n_rows(),0.0), fresh(A.n_rows(),0.0);
        EXPECT_TRUE(pc.apply(b,old));
        EXPECT_FALSE(pc.setup(SparseMatrix{}));
        EXPECT_TRUE(pc.apply(b,fresh));
        EXPECT_TRUE(rel_solution_error(fresh,old)<1e-13);
    });
    run_case("directional_numeric_refresh_rejects_dominant_rotation", [] {
        const auto A=operator_matrix(16,12,1000.0,1.0,0.01,1.0);
        const auto B=operator_matrix(16,12,1.0,1000.0,0.01,1.0);
        NativeDirectionalSemiCoarseningAMGPreconditioner pc;
        EXPECT_TRUE(pc.setup(A));
        EXPECT_FALSE(pc.update_values(B));
        Vector b(A.n_rows(),1.0), x(A.n_rows(),0.0);
        EXPECT_TRUE(pc.apply(b,x));
    });
    run_case("anisotropic_elliptic_true_residual_and_solution_x_y_diagonal", [] {
        for (const auto& [ax,ay,ad] :
             {std::tuple{1000.0,1.0,0.01},std::tuple{1.0,1000.0,0.01},
              std::tuple{0.5,0.5,1000.0},std::tuple{10000.0,1.0,0.01}}) {
            const auto A=operator_matrix(20,16,ax,ay,ad,0.5);
            const auto exact=exact_solution(A.n_rows());
            const auto b=rhs_for(A,exact);
            Vector x(A.n_rows(),0.0);
            LinearSolverRequest req;
            req.preconditioner=PreconditionerModel::DirectionalSemiCoarseningAMG;
            req.gmres_restart=60;
            const auto result=solve_linear_system(A,b,x,
                LinearProblemKind::PressurePoisson,req,350,1e-11);
            EXPECT_TRUE(result.result.status==SolverStatus::CONVERGED);
            EXPECT_TRUE(true_relative(A,b,x)<=1e-11);
            const auto err=rel_solution_error(x,exact);
            if (!(err<=1e-7))
                std::cerr << "physical solution err=" << err << " (ax=" << ax
                          << " ay=" << ay << " ad=" << ad
                          << ", true residual=" << true_relative(A,b,x) << ")\n";
            EXPECT_TRUE(err<=1e-7);
        }
    });
    run_case("directional_neumann_pressure_projected_true_residual", [] {
        const auto A=operator_matrix(20,16,1000.0,1.0,0.01,0.0);
        Vector exact=exact_solution(A.n_rows());
        double mean=0.0;
        for (std::size_t i=0;i<exact.size();++i) mean+=exact(i);
        mean/=static_cast<double>(exact.size());
        for (std::size_t i=0;i<exact.size();++i) exact(i)-=mean;
        const auto b=rhs_for(A,exact);
        Vector x(A.n_rows(),0.0);
        LinearSolverRequest request;
        request.preconditioner=PreconditionerModel::DirectionalSemiCoarseningAMG;
        request.null_space=NullSpaceModel::Constant;
        request.gmres_restart=60;
        const auto result=solve_linear_system(A,b,x,LinearProblemKind::PressurePoisson,
                                               request,300,1e-9);
        EXPECT_TRUE(result.result.status==SolverStatus::CONVERGED);
        EXPECT_TRUE(true_relative(A,b,x)<=1e-9);
        EXPECT_TRUE(rel_solution_error(x,exact)<=1e-5);
    });
    run_case("compatible_numeric_amg_update_is_not_silently_rebuilt", [] {
        const auto A=operator_matrix(20,16,1000.0,1.0,0.01,1.0);
        const auto B=operator_matrix(20,16,1100.0,1.0,0.01,1.0);
        NativeDirectionalSemiCoarseningAMGPreconditioner pc;
        EXPECT_TRUE(pc.setup(A));
        const auto pairs=pc.finest_directional_pairs();
        EXPECT_TRUE(pc.update_values(B));
        EXPECT_TRUE(pc.numeric_updates()==1);
        EXPECT_TRUE(pc.hierarchy_builds()==1);
        EXPECT_TRUE(pc.finest_directional_pairs()==pairs);
        const auto exact=exact_solution(B.n_rows());
        const auto rhs=rhs_for(B,exact);
        Vector correction(B.n_rows(),0.0);
        EXPECT_TRUE(pc.apply(rhs,correction));
        EXPECT_TRUE(std::isfinite(true_relative(B,rhs,correction)));
    });
    run_case("strong_line_selector_rejects_crossing_equal_strength_directions", [] {
        // Four diagonal couplings are equally strong: a path cannot capture
        // both crossing physical directions. Must use SA/Jacobi fallback.
        const auto A=operator_matrix(32,32,0.5,0.5,1000.0,1.0);
        NativeDirectionalSemiCoarseningAMGPreconditioner pc;
        EXPECT_TRUE(pc.setup(A));
        EXPECT_TRUE(pc.finest_line_coverage()==0);
        EXPECT_TRUE(pc.finest_line_coupling_share()<0.85);
        EXPECT_TRUE(pc.finest_directional_pairs()==0);
    });
    run_case("large_1024_directional_multilevel_true_residual_gate", [] {
        // Reproduces the first benchmark failure: a smoothed weak-edge P
        // caused enormous non-M-matrix coarse-grid amplification. Injection
        // and the strong-coupling-mass guard must keep this case convergent.
        const auto A=operator_matrix(32,32,1000.0,1.0,0.01,1.0);
        Vector exact(A.n_rows(),0.0);
        for (std::size_t i=0;i<exact.size();++i)
            exact(i)=1.0+std::sin(0.018*static_cast<double>(i)) +
                     0.2*std::cos(0.13*static_cast<double>(i));
        const auto b=rhs_for(A,exact);
        Vector x(A.n_rows(),0.0);
        LinearSolverRequest request;
        request.preconditioner=PreconditionerModel::DirectionalSemiCoarseningAMG;
        request.gmres_restart=50;
        const auto result=solve_linear_system(A,b,x,
            LinearProblemKind::PressurePoisson,request,300,1e-10);
        EXPECT_TRUE(result.result.status==SolverStatus::CONVERGED);
        EXPECT_TRUE(true_relative(A,b,x)<=1e-10);
        EXPECT_TRUE(rel_solution_error(x,exact)<=1e-7);
    });
    return run_all();
}