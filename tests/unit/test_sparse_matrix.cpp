// M0.8-T01 — Tests for SparseMatrix (CSR)

#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "test_harness.h"

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

    return run_all();
}