// M0.8-T03 — Linear system tests
#include "cfdx/core/linalg/linear_system.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "common/test_harness.h"
#include <cstdio>

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("linear_system_construct", []() {
        // Test basic construction with empty matrix/vectors
        SparseMatrix A(3, 3);
        Vector b(3, 0.0);
        Vector x(3, 0.0);
        
        LinearSystem ls(A, b, x);
        
        EXPECT_TRUE(ls.n_rows() == 3);
        EXPECT_TRUE(ls.n_cols() == 3);
        EXPECT_TRUE(ls.A().n_rows() == 3);
        EXPECT_TRUE(ls.b().size() == 3);
        EXPECT_TRUE(ls.x().size() == 3);
        EXPECT_TRUE(ls.is_consistent());
    });

    run_case("linear_system_resize_vectors", []() {
        SparseMatrix A(5, 5);
        Vector b(5, 0.0);
        Vector x(5, 0.0);
        
        LinearSystem ls(A, b, x);
        
        b.resize(10);
        x.resize(10);
        
        EXPECT_TRUE(ls.n_rows() == 5);
        EXPECT_TRUE(ls.b().size() == 10);
        EXPECT_TRUE(ls.x().size() == 10);
    });

    run_case("linear_system_set_get", []() {
        SparseMatrix A(3, 3);
        Vector b(3, 0.0);
        Vector x(3, 0.0);
        
        A.push_back(0, 0, 2.0);
        A.push_back(0, 1, -1.0);
        A.push_back(1, 0, -1.0);
        A.push_back(1, 1, 2.0);
        A.push_back(1, 2, -1.0);
        A.push_back(2, 1, -1.0);
        A.push_back(2, 2, 2.0);
        A.finalize();
        
        b(0) = 1.0;
        b(1) = 0.0;
        b(2) = 1.0;
        
        LinearSystem ls(A, b, x);
        
        EXPECT_TRUE(A(0, 0) == 2.0);
        EXPECT_TRUE(A(0, 1) == -1.0);
        EXPECT_TRUE(A(1, 0) == -1.0);
        EXPECT_TRUE(A(1, 1) == 2.0);
        EXPECT_TRUE(A(1, 2) == -1.0);
        EXPECT_TRUE(A(2, 1) == -1.0);
        EXPECT_TRUE(A(2, 2) == 2.0);
        EXPECT_TRUE(b(0) == 1.0);
        EXPECT_TRUE(b(1) == 0.0);
        EXPECT_TRUE(b(2) == 1.0);
    });

    run_case("linear_system_matvec", []() {
        SparseMatrix A(3, 3);
        Vector x(3, 0.0);
        Vector b(3, 0.0);
        
        // A = [[2, -1, 0], [-1, 2, -1], [0, -1, 2]]
        A.push_back(0, 0, 2.0);
        A.push_back(0, 1, -1.0);
        A.push_back(1, 0, -1.0);
        A.push_back(1, 1, 2.0);
        A.push_back(1, 2, -1.0);
        A.push_back(2, 1, -1.0);
        A.push_back(2, 2, 2.0);
        A.finalize();
        
        x(0) = 1.0; x(1) = 2.0; x(2) = 3.0;
        
        LinearSystem ls(A, b, x);
        
        // b = A * x
        auto b_calc = A.matvec(x);
        
        EXPECT_TRUE(b_calc.size() == 3);
        EXPECT_TRUE(b_calc[0] == 0.0);   // 2*1 + (-1)*2 = 0
        EXPECT_TRUE(b_calc[1] == 0.0);   // -1*1 + 2*2 + (-1)*3 = 0
        EXPECT_TRUE(b_calc[2] == 4.0);   // (-1)*2 + 2*3 = 4... wait: -1*2 + 2*3 = 4
        
        // Verify linear system consistency
        EXPECT_TRUE(ls.is_consistent());
    });

    return run_all();
}
