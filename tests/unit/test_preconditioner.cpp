// M0.8-T07 — Tests for preconditioners

#include "cfdx/core/linalg/preconditioner.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

SparseMatrix make_tridiagonal() {
    SparseMatrix A(3, 3);
    A.push_back(0, 0, 4.0);
    A.push_back(0, 1, -1.0);
    A.push_back(1, 0, -2.0);
    A.push_back(1, 1, 4.0);
    A.push_back(1, 2, -1.0);
    A.push_back(2, 1, -1.0);
    A.push_back(2, 2, 3.0);
    A.finalize();
    return A;
}

int main() {
    run_case("ilu0_exact_tridiagonal", []() {
        const auto A = make_tridiagonal();
        Vector expected(3);
        expected(0) = 1.0;
        expected(1) = 2.0;
        expected(2) = 3.0;

        Vector r(3);
        const auto rhs = A.matvec(expected);
        for (std::size_t i = 0; i < 3; ++i) r(i) = rhs[i];

        ILU0Preconditioner ilu;
        EXPECT_TRUE(ilu.setup(A));

        Vector z;
        EXPECT_TRUE(ilu.apply(r, z));
        for (std::size_t i = 0; i < 3; ++i) {
            EXPECT_NEAR(z(i), expected(i), 1e-12);
        }
    });

    run_case("ilu0_requires_diagonal", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 1, 1.0);
        A.push_back(1, 0, 1.0);
        A.finalize();

        ILU0Preconditioner ilu;
        EXPECT_TRUE(!ilu.setup(A));
    });

    run_case("jacobi_and_gauss_seidel", []() {
        const auto A = make_tridiagonal();
        Vector r(3);
        r(0) = 4.0;
        r(1) = 2.0;
        r(2) = 3.0;

        JacobiPreconditioner jacobi;
        EXPECT_TRUE(jacobi.setup(A));
        Vector z;
        EXPECT_TRUE(jacobi.apply(r, z));
        EXPECT_NEAR(z(0), 1.0, 1e-12);
        EXPECT_NEAR(z(1), 0.5, 1e-12);
        EXPECT_NEAR(z(2), 1.0, 1e-12);

        GaussSeidelPreconditioner gs;
        EXPECT_TRUE(gs.setup(A));
        EXPECT_TRUE(gs.apply(r, z));
        EXPECT_NEAR(z(0), 1.0, 1e-12);
        EXPECT_NEAR(z(1), 1.0, 1e-12);
        EXPECT_NEAR(z(2), 1.0, 1e-12);
    });

    return run_all();
}
