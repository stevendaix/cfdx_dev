// Phase A/B: anisotropy diagnostics and bounded-distance energy interpolation.
// Each test asserts in Release (the common harness uses live assertions).
#include "cfdx/core/linalg/anisotropy_analysis.h"
#include "cfdx/core/linalg/amg_preconditioner.h"
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "common/test_harness.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
SparseMatrix matrix(std::size_t nx, std::size_t ny,
                    double horizontal, double vertical, double oblique,
                    double shift = 0.5) {
    SparseMatrix A(nx*ny,nx*ny);
    for(std::size_t y=0;y<ny;++y) for(std::size_t x=0;x<nx;++x) {
        const std::size_t i = y*nx+x;
        double diag = shift;
        auto connect=[&](std::size_t j,double a) { A.push_back(i,j,-a); diag+=a; };
        if (x) connect(i-1,horizontal);
        if (x+1<nx) connect(i+1,horizontal);
        if (y) connect(i-nx,vertical);
        if (y+1<ny) connect(i+nx,vertical);
        if (x && y) connect(i-nx-1,oblique);
        if (x+1<nx && y+1<ny) connect(i+nx+1,oblique);
        if (x+1<nx && y) connect(i-nx+1,oblique);
        if (x && y+1<ny) connect(i+nx-1,oblique);
        A.push_back(i,i,diag);
    }
    A.finalize();
    return A;
}
Vector manufactured(std::size_t n) {
    Vector x(n,0.0);
    for(std::size_t i=0;i<n;++i)
        x(i)=0.6*std::sin(0.011*static_cast<double>(i))+
             0.25*std::cos(0.23*static_cast<double>(i))+1.0;
    return x;
}
Vector product(const SparseMatrix& A,const Vector& x) {
    Vector result(x.size(),0.0);
    A.matvec_into(x,result);
    return result;
}
double error(const Vector& x,const Vector& exact) {
    double r=0.0, n=0.0;
    for(std::size_t i=0;i<x.size();++i) {
        r=std::hypot(r,x(i)-exact(i));
        n=std::hypot(n,exact(i));
    }
    return r/n;
}
}

int main() {
    run_case("matrix_and_mesh_anisotropy_report", [] {
        constexpr std::size_t nx=20, ny=14;
        const auto A=matrix(nx,ny,1000.0,1.0,0.001);
        AnisotropyGeometry mesh;
        for(std::size_t y=0;y<ny;++y) for(std::size_t x=0;x<nx;++x) {
            mesh.cell_centers.push_back({static_cast<double>(x),
                static_cast<double>(y),0.0});
            mesh.cell_aspect_ratios.push_back(1000.0);
        }
        mesh.face_nonorthogonality_degrees={0.0,30.0};
        mesh.face_skewness={0.0,0.2};
        const auto report=analyze_anisotropy(A,&mesh);
        EXPECT_TRUE(report.rows==nx*ny);
        EXPECT_TRUE(report.nonzeros==A.nnz());
        EXPECT_TRUE(report.anisotropic_rows>nx*ny/2);
        EXPECT_TRUE(report.geometrically_directional_rows>nx*ny/2);
        EXPECT_TRUE(report.strong_directed_edges>nx*ny);
        EXPECT_TRUE(report.mean_neighbor_direction_alignment>0.9);
        EXPECT_TRUE(report.nonpositive_diagonals==0);
        EXPECT_TRUE(report.asymmetric_strong_edges==0);
        EXPECT_NEAR(report.max_cell_aspect_ratio,1000.0,1e-12);
        EXPECT_NEAR(report.max_face_nonorthogonality_degrees,30.0,1e-12);
        EXPECT_TRUE(report.json("pressure").find("\"equation\":\"pressure\"")!=std::string::npos);
    });
    run_case("isotropy_and_tiny_connections_are_not_false_anisotropy", [] {
        const auto A=matrix(12,12,1.0,1.0,0.0);
        const auto report=analyze_anisotropy(A);
        EXPECT_TRUE(report.anisotropic_rows==0);
        EXPECT_TRUE(report.worst_offdiagonal_contrast<=1.01);
        EXPECT_TRUE(!report.has_geometry);
        EXPECT_TRUE(report.geometric_dimension==0);
    });
    run_case("face_quality_only_is_recognized_as_geometry", [] {
        const auto A=matrix(4,4,1.0,1.0,0.0);
        AnisotropyGeometry g;
        g.face_nonorthogonality_degrees={5.0,65.0};
        g.face_skewness={0.1,0.9};
        const auto report=analyze_anisotropy(A,&g);
        EXPECT_TRUE(report.has_geometry);
        EXPECT_TRUE(report.geometric_dimension==0);
        EXPECT_NEAR(report.max_face_nonorthogonality_degrees,65.0,1e-14);
        EXPECT_NEAR(report.max_face_skewness,0.9,1e-14);
    });
    run_case("malformed_geometry_rejected_without_guessing", [] {
        const auto A=matrix(6,6,1.0,10.0,0.0);
        AnisotropyGeometry g;
        g.cell_centers.push_back({1.0,2.0,0.0});
        EXPECT_THROW(analyze_anisotropy(A,&g),std::invalid_argument);
    });
    run_case("energy_interpolation_preserves_constant_and_reduces_energy", [] {
        const auto A=matrix(16,16,100.0,1.0,20.0);
        FunctionalLinearOperator op(A.n_rows(), [&A](const Vector& x,Vector& y) {
            if (y.size()!=x.size()) y.resize(x.size());
            A.matvec_into(x,y);
        });
        MatrixFreeVcyclePreconditioner regular(op,0.7,4,4,0.25,25,
            AMGInterpolationPolicy::SmoothedAggregation,false,false);
        MatrixFreeVcyclePreconditioner optimized(op,0.7,4,4,0.25,25,
            AMGInterpolationPolicy::EnergyMinimizingAggregation,false,false);
        EXPECT_TRUE(regular.setup(A));
        EXPECT_TRUE(optimized.setup(A));
        EXPECT_NEAR(optimized.prolongation_row_sum_min(),1.0,1e-12);
        EXPECT_NEAR(optimized.prolongation_row_sum_max(),1.0,1e-12);
        EXPECT_TRUE(std::isfinite(optimized.first_prolongation_energy()));
        EXPECT_TRUE(optimized.first_prolongation_energy() <
                    regular.first_prolongation_energy());
        EXPECT_TRUE(optimized.first_prolongation_row_nnz_max()>1);
    });
    run_case("energy_minimized_amg_fgmres_true_residual_x_y_oblique", [] {
        for(const auto& configuration :
            {std::array<double,3>{1000.0,1.0,0.01},
             std::array<double,3>{1.0,1000.0,0.01},
             std::array<double,3>{1.0,1.0,1000.0}}) {
            const auto A=matrix(16,16,configuration[0],configuration[1],configuration[2]);
            const auto exact=manufactured(A.n_rows());
            const auto b=product(A,exact);
            Vector x(A.n_rows(),0.0);
            LinearSolverRequest req;
            req.preconditioner=PreconditionerModel::EnergyMinimizingAMG;
            req.gmres_restart=50;
            const auto plan=select_linear_solver(
                LinearProblemKind::PressurePoisson,measure_matrix_characteristics(A),req);
            EXPECT_TRUE(plan.krylov==KrylovModel::FGMRES);
            const auto result=solve_linear_system(
                A,b,x,LinearProblemKind::PressurePoisson,req,350,1e-10);
            EXPECT_TRUE(result.result.status==SolverStatus::CONVERGED);
            EXPECT_TRUE(physical_relative_residual(A,x,b)<=1e-10);
            EXPECT_TRUE(error(x,exact)<1e-7);
        }
    });
    run_case("energy_AMG_projected_Neumann_anisotropic_nullspace", [] {
        const auto A=matrix(16,16,1000.0,1.0,0.02,0.0);
        Vector exact=manufactured(A.n_rows());
        double mean=0.0;
        for(std::size_t i=0;i<exact.size();++i) mean+=exact(i);
        mean/=static_cast<double>(exact.size());
        for(std::size_t i=0;i<exact.size();++i) exact(i)-=mean;
        const Vector b=product(A,exact);
        Vector x(A.n_rows(),0.0);
        LinearSolverRequest req;
        req.preconditioner=PreconditionerModel::EnergyMinimizingAMG;
        req.null_space=NullSpaceModel::Constant;
        req.gmres_restart=45;
        const auto result=solve_linear_system(A,b,x,
            LinearProblemKind::PressurePoisson,req,350,1e-10);
        EXPECT_TRUE(result.result.status==SolverStatus::CONVERGED);
        EXPECT_TRUE(physical_relative_residual(A,x,b)<1e-10);
        EXPECT_TRUE(error(x,exact)<1e-7);
    });
    run_case("energy_amg_rejects_CG_and_nonsymmetric_operator", [] {
        const auto A=matrix(8,8,1000.0,1.0,0.0);
        LinearSolverRequest req;
        req.preconditioner=PreconditionerModel::EnergyMinimizingAMG;
        req.krylov=KrylovModel::CG;
        EXPECT_THROW(select_linear_solver(LinearProblemKind::PressurePoisson,
            measure_matrix_characteristics(A),req),std::invalid_argument);
        req.krylov=KrylovModel::Auto;
        EXPECT_THROW(select_linear_solver(LinearProblemKind::ScalarTransport,
            measure_matrix_characteristics(A),req),std::invalid_argument);
    });
    return run_all();
}