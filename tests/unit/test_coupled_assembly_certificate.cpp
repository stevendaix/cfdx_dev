#include "cfdx/physics/coupled_assembly_certificate.h"
#include "common/test_harness.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    // Two velocity and two pressure unknowns, pressure row 1 pinned.
    SparseMatrix A(4,4);
    A.push_back(0,0,4.0); A.push_back(0,1,0.2);
    A.push_back(0,2,1.0);
    A.push_back(1,1,3.0); A.push_back(1,2,-0.4);
    A.push_back(2,0,-1.0); A.push_back(2,1,0.5);
    A.push_back(2,2,0.1);
    A.push_back(3,3,1.0);
    A.finalize();
    Vector rhs(4,0.0), guess(4,0.0);
    rhs(0)=1.0; rhs(1)=-1.0; rhs(2)=0.2; rhs(3)=0.0;
    guess(0)=0.3; guess(1)=-0.4; guess(2)=0.5; guess(3)=0.0;
    const std::vector<double> adjusted{4.0,1.0};
    const std::vector<double> raw{4.0,2.0};
    CoupledAssemblyCertificateOptions options;
    options.pressure_reference_cell=1;
    options.pressure_schur_diagonal=&adjusted;
    options.pressure_schur_diagonal_raw=&raw;
    run_case("coupled_assembly_certificate_captures_nonzero_initial_residual", [&] {
        const auto c=certify_coupled_assembly_structure(A,rhs,guess,2,2,options);
        EXPECT_TRUE(c.passed);
        EXPECT_TRUE(c.nnz_momentum == 3);
        EXPECT_TRUE(c.nnz_gradient == 2);
        EXPECT_TRUE(c.nnz_divergence == 2);
        EXPECT_TRUE(c.nnz_pressure == 2);
        EXPECT_TRUE(c.momentum_linear_residual_linf > 0.5);
        EXPECT_TRUE(std::isfinite(c.continuity_linear_residual_linf));
        std::cout << "COUPLED_CSR initial_momentum_linear_residual="
                  << c.momentum_linear_residual_linf << '\n';
    });
    run_case("coupled_assembly_certificate_rejects_wrong_gauge_row", [&] {
        SparseMatrix wrong(4,4);
        for(std::size_t i=0;i<4;++i)
            for(std::size_t k=A.row_offsets_data()[i];k<A.row_offsets_data()[i+1];++k)
                wrong.push_back(i,A.columns_data()[k],A.values_data()[k]);
        wrong.push_back(3,0,0.01);
        wrong.finalize();
        EXPECT_TRUE(!certify_coupled_assembly_structure(
            wrong,rhs,guess,2,2,options).passed);
    });
    run_case("coupled_assembly_certificate_rejects_inconsistent_schur_diagonal", [&] {
        const std::vector<double> wrong_adjusted{4.1,1.0};
        auto bad=options; bad.pressure_schur_diagonal=&wrong_adjusted;
        EXPECT_TRUE(!certify_coupled_assembly_structure(A,rhs,guess,2,2,bad).passed);
    });
    run_case("coupled_assembly_certificate_rejects_zero_and_nonfinite_rows", [&] {
        SparseMatrix zero(4,4);
        zero.push_back(0,0,4.0); zero.push_back(1,1,3.0);
        zero.push_back(2,2,0.1); zero.finalize();
        EXPECT_TRUE(!certify_coupled_assembly_structure(
            zero,rhs,guess,2,2,options).passed);
        Vector invalid_rhs = rhs;
        invalid_rhs(0)=std::numeric_limits<double>::quiet_NaN();
        EXPECT_TRUE(!certify_coupled_assembly_structure(
            A,invalid_rhs,guess,2,2,options).passed);
    });
    return run_all();
}
