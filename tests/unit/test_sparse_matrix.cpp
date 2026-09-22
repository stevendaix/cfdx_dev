// M0.8-T01 — Tests for SparseMatrix (CSR)\n// Includes assembly-order and storage-invariant regressions.

#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"
#include <limits>

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("sparse_matrix_default", []() {
        SparseMatrix A;
        EXPECT_TRUE(A.n_rows() == 0);
        EXPECT_TRUE(A.n_cols() == 0);
        EXPECT_TRUE(A.nnz() == 0);
    });

    run_case("sparse_matrix_explicit", []() {
        SparseMatrix A(3, 3);
        EXPECT_TRUE(A.n_rows() == 3);
        EXPECT_TRUE(A.n_cols() == 3);
    });

    run_case("sparse_matrix_push_back", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(0, 1, 2.0);
        A.push_back(1, 1, 3.0);
        A.finalize();
        EXPECT_TRUE(A.nnz() == 3);
        EXPECT_TRUE(A(0, 0) == 1.0);
        EXPECT_TRUE(A(0, 1) == 2.0);
        EXPECT_TRUE(A(1, 1) == 3.0);
        EXPECT_TRUE(A(1, 0) == 0.0);
    });

    run_case("sparse_matrix_consistent", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 2.0);
        A.finalize();
        EXPECT_TRUE(A.is_consistent());
    });

    run_case("sparse_matrix_row_out_of_range", []() {
        SparseMatrix A(2, 2);
        EXPECT_THROW(A.push_back(5, 0, 1.0), std::out_of_range);
    });

    run_case("sparse_matrix_col_out_of_range", []() {
        SparseMatrix A(2, 2);
        EXPECT_THROW(A.push_back(0, 5, 1.0), std::out_of_range);
    });

    run_case("sparse_matrix_matvec", []() {
        // A = [1 2 0]
        //     [0 0 3]
        SparseMatrix A(2, 3);
        A.push_back(0, 0, 1.0);
        A.push_back(0, 1, 2.0);
        A.push_back(1, 2, 3.0);
        A.finalize();

        Vector x(3);
        x(0) = 1.0;
        x(1) = 2.0;
        x(2) = 3.0;

        auto y = A.matvec(x);
        EXPECT_TRUE(y.size() == 2);
        EXPECT_TRUE(y[0] == 1.0 * 1.0 + 2.0 * 2.0);  // 5.0
        EXPECT_TRUE(y[1] == 3.0 * 3.0);               // 9.0
    });

    run_case("sparse_matrix_matvec_transpose", []() {
        SparseMatrix A(2, 3);
        A.push_back(0, 0, 1.0);
        A.push_back(0, 1, 2.0);
        A.push_back(1, 2, 3.0);
        A.finalize();

        Vector x(2);
        x(0) = 1.0;
        x(1) = 2.0;

        auto y = A.matvec_transpose(x);
        EXPECT_TRUE(y.size() == 3);
        EXPECT_TRUE(y[0] == 1.0 * 1.0);  // 1.0
        EXPECT_TRUE(y[1] == 2.0 * 1.0);  // 2.0
        EXPECT_TRUE(y[2] == 3.0 * 2.0);  // 6.0
    });

    run_case("sparse_matrix_clear", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 2.0);
        A.finalize();
        A.clear();
        EXPECT_TRUE(A.nnz() == 0);
        EXPECT_TRUE(A.is_consistent());
    });

    run_case("sparse_matrix_transpose_dimension_mismatch", []() {
        SparseMatrix A(2, 3);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 2, 2.0);
        A.finalize();

        Vector wrong(3);
        EXPECT_THROW(A.matvec_transpose(wrong), std::runtime_error);
    });

    run_case("sparse_matrix_csr_row_offsets_and_transpose_are_exact", []() {
        SparseMatrix A(3, 2);
        A.push_back(0, 0, 2.0);
        A.push_back(0, 1, -1.0);
        A.push_back(1, 0, 4.0);
        A.push_back(2, 1, 3.0);
        A.finalize();

        EXPECT_TRUE(A.is_consistent());
        EXPECT_TRUE(A.nnz() == 4);
        EXPECT_TRUE(A.row_offsets_data()[0] == 0);
        EXPECT_TRUE(A.row_offsets_data()[1] == 2);
        EXPECT_TRUE(A.row_offsets_data()[2] == 3);
        EXPECT_TRUE(A.row_offsets_data()[3] == 4);

        Vector x(2);
        x(0) = 5.0;
        x(1) = 7.0;
        const auto y = A.matvec(x);
        EXPECT_NEAR(y[0], 3.0, 1e-14);
        EXPECT_NEAR(y[1], 20.0, 1e-14);
        EXPECT_NEAR(y[2], 21.0, 1e-14);

        Vector xt(3);
        xt(0) = 11.0;
        xt(1) = -2.0;
        xt(2) = 3.0;
        const auto yt = A.matvec_transpose(xt);
        EXPECT_NEAR(yt[0], 14.0, 1e-14);
        EXPECT_NEAR(yt[1], -2.0, 1e-14);
    });

    run_case("sparse_matrix_finalize_reorders_out_of_order_rows_and_columns", []() {
        SparseMatrix A(3, 4);
        A.push_back(2, 3, 30.0);
        A.push_back(0, 2, 2.0);
        A.push_back(1, 3, 13.0);
        A.push_back(0, 0, 1.0);
        A.push_back(2, 1, 21.0);
        A.push_back(1, 0, 10.0);
        A.finalize();
        EXPECT_TRUE(A.is_consistent());
        EXPECT_TRUE(A.row_offsets_data()[1] == 2);
        EXPECT_TRUE(A.row_offsets_data()[2] == 4);
        EXPECT_TRUE(A.row_offsets_data()[3] == 6);
        EXPECT_TRUE(A(0, 0) == 1.0);
        EXPECT_TRUE(A(0, 2) == 2.0);
        EXPECT_TRUE(A(1, 0) == 10.0);
        EXPECT_TRUE(A(1, 3) == 13.0);
    });

    run_case("sparse_matrix_duplicate_entries_accumulate_consistently", []() {
        SparseMatrix A(1, 2);
        A.push_back(0, 0, 1.5);
        A.push_back(0, 0, 2.5);
        A.push_back(0, 1, -1.0);
        A.finalize();
        EXPECT_TRUE(A.is_consistent());
        EXPECT_NEAR(A(0, 0), 4.0, 1e-14);
        Vector x(2);
        x(0) = 2.0;
        x(1) = 3.0;
        EXPECT_NEAR(A.matvec(x)[0], 5.0, 1e-14);
    });

    run_case("sparse_matrix_finalize_is_idempotent", []() {
        SparseMatrix A(1, 2);
        A.push_back(0, 1, 2.0);
        A.push_back(0, 0, 1.0);
        A.finalize();
        A.finalize();
        EXPECT_TRUE(A.is_consistent());
        EXPECT_TRUE(A(0, 0) == 1.0);
    });

    run_case("sparse_matrix_rejects_nonfinite_assembly_value", []() {
        SparseMatrix A(1, 1);
        EXPECT_THROW(A.push_back(0, 0, std::nan("")), std::invalid_argument);
        EXPECT_THROW(A.push_back(0, 0, std::numeric_limits<double>::infinity()), std::invalid_argument);
    });

    run_case("sparse_matrix_rejects_push_after_finalize", []() {
        SparseMatrix A(1, 1);
        A.push_back(0, 0, 1.0);
        A.finalize();
        EXPECT_THROW(A.push_back(0, 0, 2.0), std::logic_error);
    });

    run_case("sparse_matrix_detects_mutated_nonfinite_storage", []() {
        SparseMatrix A(1, 1);
        A.push_back(0, 0, 1.0);
        A.finalize();
        A.values_data()[0] = std::nan("");
        EXPECT_FALSE(A.is_consistent());
    });

    run_case("sparse_matrix_clear_restores_empty_csr", []() {
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 2.0);
        A.push_back(2, 2, 3.0);
        A.finalize();
        A.clear();

        EXPECT_TRUE(A.n_rows() == 3);
        EXPECT_TRUE(A.n_cols() == 3);
        EXPECT_TRUE(A.nnz() == 0);
        EXPECT_TRUE(A.row_offsets_data()[0] == 0);
        EXPECT_TRUE(A.row_offsets_data()[1] == 0);
        EXPECT_TRUE(A.row_offsets_data()[2] == 0);
        EXPECT_TRUE(A.row_offsets_data()[3] == 0);
        EXPECT_TRUE(A.is_consistent());
    });

    return run_all();
}