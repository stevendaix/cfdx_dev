#include "cfdx/core/linalg/matrix_diagnostics.h"
#include "cfdx/core/linalg/null_space.h"
#include "common/test_harness.h"

#include <cmath>
#include <limits>
#include <string>

using namespace cfdx::core;
using namespace cfdx::testing;

static SparseMatrix make_pathological_matrix() {
    SparseMatrix A(5, 5);
    A.push_back(0, 0, 1.0);
    A.push_back(0, 1, -0.5);
    A.push_back(1, 0, -0.5);
    A.push_back(1, 1, 2.0);
    A.push_back(2, 2, 1e-18);
    A.push_back(3, 3, 3.0);
    A.push_back(4, 4, 4.0);
    A.push_back(4, 4, 5.0);
    A.finalize();
    return A;
}

int main() {
    run_case("diagnostics_detect_structural_pathologies", [] {
        const auto d = diagnose_matrix(make_pathological_matrix(), 1e-12);
        EXPECT_TRUE(d.rows == 5);
        EXPECT_TRUE(d.columns == 5);
        EXPECT_TRUE(d.empty_rows == 0);
        EXPECT_TRUE(d.missing_diagonal == 0);
        EXPECT_TRUE(d.near_zero_diagonal >= 1);
        EXPECT_TRUE(d.isolated_dofs >= 3);
        EXPECT_TRUE(d.duplicate_entries == 1);
        EXPECT_TRUE(d.connected_components >= 3);
        EXPECT_TRUE(d.diagonal_dynamic_range > 1e17);

        const auto pathologies = classify_matrix_pathologies(d);
        EXPECT_TRUE(pathologies.size() >= 2);
        EXPECT_TRUE(std::string(to_string(pathologies.front())) != "none");
    });

    run_case("diagnostics_report_clean_connected_matrix", [] {
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 2.0);
        A.push_back(0, 1, -1.0);
        A.push_back(1, 0, -1.0);
        A.push_back(1, 1, 2.0);
        A.push_back(1, 2, -1.0);
        A.push_back(2, 1, -1.0);
        A.push_back(2, 2, 2.0);
        A.finalize();

        const auto d = diagnose_matrix(A);
        EXPECT_TRUE(d.finite);
        EXPECT_TRUE(d.empty_rows == 0);
        EXPECT_TRUE(d.empty_columns == 0);
        EXPECT_TRUE(d.missing_diagonal == 0);
        EXPECT_TRUE(d.isolated_dofs == 0);
        EXPECT_TRUE(d.connected_components == 1);
        EXPECT_TRUE(d.minimum_gershgorin_margin >= 0.0);
    });

    run_case("row_and_column_scaling_are_explicit", [] {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1e-6);
        A.push_back(0, 1, 2e-6);
        A.push_back(1, 0, 3e3);
        A.push_back(1, 1, 4e3);
        A.finalize();

        const auto scaled = scale_matrix(A, MatrixScaling::RowColumn);
        EXPECT_TRUE(scaled.applied);
        EXPECT_NEAR(scaled.matrix(0, 0), 1e-6 / 2e-6 / 3e3, 1e-20);
        EXPECT_NEAR(scaled.row_scale[0], 5e5, 1e-8);
        EXPECT_NEAR(scaled.row_scale[1], 2.5e-4, 1e-16);
        EXPECT_NEAR(scaled.column_scale[0], 1.0 / 3e3, 1e-16);
        EXPECT_NEAR(scaled.column_scale[1], 1.0 / 4e3, 1e-16);
    });

    run_case("symmetric_diagonal_scaling_preserves_symmetry", [] {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 4.0);
        A.push_back(0, 1, 2.0);
        A.push_back(1, 0, 2.0);
        A.push_back(1, 1, 9.0);
        A.finalize();

        const auto scaled = scale_matrix(A, MatrixScaling::SymmetricDiagonal);
        EXPECT_NEAR(scaled.matrix(0, 1), scaled.matrix(1, 0), 1e-14);
        EXPECT_NEAR(scaled.matrix(0, 0), 1.0, 1e-14);
        EXPECT_NEAR(scaled.matrix(1, 1), 1.0, 1e-14);
    });

    run_case("scaling_does_not_repair_zero_rows", [] {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 0.0);
        A.push_back(1, 1, 2.0);
        A.finalize();

        const auto scaled = scale_matrix(A, MatrixScaling::Row);
        EXPECT_TRUE(scaled.zero_rows == 1);
        EXPECT_NEAR(scaled.matrix(0, 0), 0.0, 0.0);
        EXPECT_NEAR(scaled.matrix(1, 1), 1.0, 1e-14);
    });

    run_case("nullspace_compatibility_remains_explicit", [] {
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 1.0);
        A.push_back(0, 1, -1.0);
        A.push_back(1, 0, -1.0);
        A.push_back(1, 1, 2.0);
        A.push_back(1, 2, -1.0);
        A.push_back(2, 1, -1.0);
        A.push_back(2, 2, 1.0);
        A.finalize();

        const auto null_space = NullSpaceProjector::constant(3);
        Vector compatible(3);
        compatible(0) = 1.0;
        compatible(1) = -2.0;
        compatible(2) = 1.0;
        Vector incompatible(3, 1.0);

        EXPECT_TRUE(null_space.is_compatible(compatible));
        EXPECT_TRUE(!null_space.is_compatible(incompatible));
        EXPECT_TRUE(rhs_compatibility_norm(null_space, incompatible) > 0.0);
        EXPECT_TRUE(null_space.is_null_space(A));
    });

    return run_all();
}
