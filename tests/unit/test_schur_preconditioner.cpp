#include "cfdx/core/linalg/block_preconditioner.h"
#include "cfdx/core/linalg/preconditioner.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::core;
using namespace cfdx::testing;

static SparseMatrix make_block_system() {
    SparseMatrix A(3, 3);
    A.push_back(0,0,4.0); A.push_back(0,1,1.0); A.push_back(0,2,1.0);
    A.push_back(1,0,1.0); A.push_back(1,1,3.0); A.push_back(1,2,2.0);
    A.push_back(2,0,2.0); A.push_back(2,1,1.0); A.push_back(2,2,5.0);
    A.finalize();
    return A;
}

int main() {
    run_case("schur_complement_setup_and_apply", [] {
        const auto A = make_block_system();
        SchurComplementPreconditioner p({0,1}, {2});
        EXPECT_TRUE(p.setup(A));

        Vector r(3);
        r(0) = 9.0; r(1) = 10.0; r(2) = 19.0;
        Vector z(3);
        EXPECT_TRUE(p.apply(r, z));
        EXPECT_TRUE(z.is_valid());

        EXPECT_NEAR(z(0), 17.0 / 11.0, 1e-12);
        EXPECT_NEAR(z(1), 31.0 / 11.0, 1e-12);
        EXPECT_NEAR(z(2), 72.0 / 23.0, 1e-12);
        EXPECT_TRUE(std::isfinite(z.norm_inf()));
    });

    run_case("schur_rejects_invalid_partition", [] {
        SparseMatrix A(3,3);
        A.push_back(0,0,1.0); A.push_back(1,1,1.0); A.push_back(2,2,1.0);
        A.finalize();

        SchurComplementPreconditioner overlap({0,1}, {1,2});
        EXPECT_TRUE(!overlap.setup(A));

        SchurComplementPreconditioner incomplete({0}, {2});
        EXPECT_TRUE(!incomplete.setup(A));
    });

    run_case("coupled_block_schur_setup_and_apply", [] {
        // One-cell saddle-point system:
        //
        // [ 4  0  0  1 ] [ux]   [4]
        // [ 0  5  0  2 ] [uy] = [5]
        // [ 0  0  6  3 ] [uz]   [6]
        // [-1 -2 -3  0 ] [ p]   [2]
        //
        // S = -D M^-1 G = 1/4 + 4/5 + 9/6 = 2.55.
        SparseMatrix A(4, 4);
        A.push_back(0, 0, 4.0); A.push_back(0, 3, 1.0);
        A.push_back(1, 1, 5.0); A.push_back(1, 3, 2.0);
        A.push_back(2, 2, 6.0); A.push_back(2, 3, 3.0);
        A.push_back(3, 0, -1.0); A.push_back(3, 1, -2.0);
        A.push_back(3, 2, -3.0);
        A.finalize();

        CoupledBlockSchurPreconditioner p(1);
        EXPECT_TRUE(p.setup(A));

        Vector r(4);
        r(0) = 4.0; r(1) = 5.0; r(2) = 6.0; r(3) = 2.0;
        Vector z(4);
        EXPECT_TRUE(p.apply(r, z));

        const double schur = 2.55;
        const double zp = 8.0 / schur;
        EXPECT_NEAR(z(3), zp, 1e-12);
        EXPECT_NEAR(z(0), 1.0 - zp / 4.0, 1e-12);
        EXPECT_NEAR(z(1), 1.0 - 2.0 * zp / 5.0, 1e-12);
        EXPECT_NEAR(z(2), 1.0 - 3.0 * zp / 6.0, 1e-12);
    });

    run_case("coupled_block_schur_uses_sparse_pressure_coupling", [] {
        // Two-cell pressure block with non-zero off-diagonal Schur terms.
        // M = 4 I, D = -I, G = [[1,0.5],[0.5,1]], C = 2 I,
        // so S = C - D M^-1 G = [[2.25,0.125],[0.125,2.25]].
        SparseMatrix A(8, 8);
        for (std::size_t c = 0; c < 2; ++c) {
            A.push_back(c, c, 4.0);
            A.push_back(2 + c, 2 + c, 4.0);
            A.push_back(4 + c, 4 + c, 4.0);
            A.push_back(6 + c, 6 + c, 2.0);
            A.push_back(6 + c, c, -1.0);
        }
        A.push_back(0, 6, 1.0);
        A.push_back(0, 7, 0.5);
        A.push_back(1, 6, 0.5);
        A.push_back(1, 7, 1.0);
        A.push_back(2, 6, 0.0);
        A.push_back(2, 7, 0.0);
        A.push_back(3, 6, 0.0);
        A.push_back(3, 7, 0.0);
        A.push_back(4, 6, 0.0);
        A.push_back(4, 7, 0.0);
        A.push_back(5, 6, 0.0);
        A.push_back(5, 7, 0.0);
        A.finalize();

        CoupledBlockSchurPreconditioner p(2);
        EXPECT_TRUE(p.setup(A));

        Vector r(8, 0.0);
        r(6) = 1.0;
        Vector z(8, 0.0);
        EXPECT_TRUE(p.apply(r, z));

        const double det = 2.25 * 2.25 - 0.125 * 0.125;
        const double p0 = 2.25 / det;
        const double p1 = -0.125 / det;
        EXPECT_NEAR(z(6), p0, 1e-12);
        EXPECT_NEAR(z(7), p1, 1e-12);
        EXPECT_NEAR(z(0), -(p0 + 0.5 * p1) / 4.0, 1e-12);
        EXPECT_NEAR(z(1), -(0.5 * p0 + p1) / 4.0, 1e-12);
        EXPECT_TRUE(std::isfinite(z.norm_inf()));
    });

    run_case("coupled_block_schur_accepts_pressure_gauge_row", [] {
        SparseMatrix A(4, 4);
        A.push_back(0, 0, 2.0);
        A.push_back(1, 1, 3.0);
        A.push_back(2, 2, 4.0);
        A.push_back(3, 3, 1.0);
        A.finalize();

        CoupledBlockSchurPreconditioner p(1);
        EXPECT_TRUE(p.setup(A));

        Vector r(4);
        r(0) = 2.0; r(1) = 3.0; r(2) = 4.0; r(3) = 7.0;
        Vector z(4);
        EXPECT_TRUE(p.apply(r, z));
        EXPECT_NEAR(z(0), 1.0, 1e-12);
        EXPECT_NEAR(z(1), 1.0, 1e-12);
        EXPECT_NEAR(z(2), 1.0, 1e-12);
        EXPECT_NEAR(z(3), 7.0, 1e-12);
    });

    return run_all();
}
