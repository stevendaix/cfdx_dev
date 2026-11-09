#include "cfdx/core/linalg/compensated_pressure_balance.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("two_word_pressure_difference_survives_fp64_rounding", [] {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 2.0); A.push_back(0, 1, -2.0);
        A.push_back(1, 0, -2.0); A.push_back(1, 1, 2.0);
        A.finalize();
        Vector rhs(2, 0.0), hi(2, 1e8), lo(2, 0.0);
        lo(1) = 2e-9;
        rhs(0) = -4e-9; rhs(1) = 4e-9;
        const auto accurate = audit_compensated_pressure_balance(A, rhs, hi, lo);
        EXPECT_TRUE(accurate.valid);
        EXPECT_TRUE(accurate.symmetric_m_matrix);
        EXPECT_TRUE(accurate.original_relative_l2 < 1e-13);
        EXPECT_TRUE(accurate.algebraic_face_relative_l2 < 1e-13);
        EXPECT_TRUE(accurate.max_row_face_identity < 1e-16);
        EXPECT_TRUE(!accurate.physical_face_flux_certified);
        const auto rounded = audit_compensated_pressure_balance(
            A, rhs, hi, Vector(2, 0.0));
        EXPECT_TRUE(rounded.original_relative_l2 > 0.99);
        EXPECT_TRUE(rounded.algebraic_face_relative_l2 > 0.99);
    });
    run_case("two_sum_accumulates_lost_pressure_increment", [] {
        double hi = 1e8, lo = 0.0;
        EXPECT_TRUE(compensated_add(hi, lo, 2e-9));
        EXPECT_TRUE(hi == 1e8);
        EXPECT_TRUE(lo == 2e-9);
        EXPECT_TRUE(compensated_add(hi, lo, -2e-9));
        EXPECT_TRUE(hi == 1e8);
        EXPECT_TRUE(lo == 0.0);
        EXPECT_TRUE(!compensated_add(hi, lo, INFINITY));
    });
    run_case("dirichlet_leakage_and_internal_pairs_are_balanced", [] {
        SparseMatrix A(3, 3);
        A.push_back(0,0,3); A.push_back(0,1,-2);
        A.push_back(1,0,-2); A.push_back(1,1,4); A.push_back(1,2,-2);
        A.push_back(2,1,-2); A.push_back(2,2,3);
        A.finalize();
        Vector rhs(3,0.0), hi(3,0.0), lo(3,0.0);
        hi(0)=2; hi(1)=3; hi(2)=1;
        rhs(0)=0; rhs(1)=6; rhs(2)=-3;
        const auto proof = audit_compensated_pressure_balance(A,rhs,hi,lo);
        EXPECT_TRUE(proof.valid);
        EXPECT_TRUE(proof.symmetric_m_matrix);
        EXPECT_TRUE(proof.original_relative_l2 < 1e-14);
        EXPECT_TRUE(proof.algebraic_face_relative_l2 < 1e-14);
    });
    run_case("asymmetric_operator_does_not_forge_face_certificate", [] {
        SparseMatrix A(2,2);
        A.push_back(0,0,2); A.push_back(0,1,-1);
        A.push_back(1,0,-2); A.push_back(1,1,3);
        A.finalize();
        Vector rhs(2,1), hi(2,1), lo(2,0);
        const auto proof = audit_compensated_pressure_balance(A,rhs,hi,lo);
        EXPECT_TRUE(proof.valid);
        EXPECT_TRUE(!proof.symmetric_m_matrix);
        EXPECT_TRUE(std::isinf(proof.algebraic_face_relative_l2));
    });
    run_case("zero_rhs_and_nonfinite_fields_are_refused", [] {
        SparseMatrix A(1,1); A.push_back(0,0,1); A.finalize();
        Vector rhs(1,0), hi(1,0), lo(1,0);
        EXPECT_TRUE(!audit_compensated_pressure_balance(A,rhs,hi,lo).valid);
        rhs(0)=1; hi(0)=INFINITY;
        EXPECT_TRUE(!audit_compensated_pressure_balance(A,rhs,hi,lo).valid);
    });
    return run_all();
}
