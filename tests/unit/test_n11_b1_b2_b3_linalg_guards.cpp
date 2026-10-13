#include "cfdx/core/linalg/advanced_preconditioners.h"
#include "cfdx/core/linalg/block_preconditioner.h"
#include "cfdx/core/linalg/preconditioner.h"
#include "common/test_harness.h"
#include <algorithm>
#include <cmath>
#include <limits>

using namespace cfdx::core;
using namespace cfdx::testing;

// N11 B1/B2/B3 linalg guard contracts (#480):
// - B1: BlockJacobi pivot rejection is relative to the block scale
//   (1e-14 * max |a_ij|), matching invert_dense, so uniformly rescaled
//   blocks stay invertible.
// - B2: the Schur preconditioner fails closed on a zero or non-finite
//   diagonal_schur_scale_ when the Diagonal factorization is requested.
// - B3: Gauss-Seidel diagonal and ILU(0) pivot/diagonal guards are relative
//   to the appropriate scale (GS: raw row scale; ILU(0): the factorized U
//   scale of the pivot's own row, never the dimensionless L multipliers)
//   instead of an absolute epsilon cut.

static SparseMatrix make_1d_laplacian(std::size_t n, double scale) {
    SparseMatrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        A.push_back(i, i, 2.0 * scale);
        if (i > 0) A.push_back(i, i - 1, -1.0 * scale);
        if (i + 1 < n) A.push_back(i, i + 1, -1.0 * scale);
    }
    A.finalize();
    return A;
}

static double residual_inf(const SparseMatrix& A, const Vector& z,
                           const Vector& r) {
    const auto Az = A.matvec(z);
    double max_err = 0.0;
    for (std::size_t i = 0; i < r.size(); ++i)
        max_err = std::max(max_err, std::abs(Az[i] - r(i)));
    return max_err;
}

int main() {
    run_case("b1_block_jacobi_accepts_uniformly_rescaled_block", [] {
        // A regular 2x2 block rescaled by 1e-31: the old absolute 1e-30
        // pivot cut rejected it; the scale-relative floor must invert it.
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 4.0e-31); A.push_back(0, 1, 1.0e-31);
        A.push_back(1, 0, 1.0e-31); A.push_back(1, 1, 3.0e-31);
        A.finalize();

        BlockDiagonalPreconditioner p({{0, 1}});
        EXPECT_TRUE(p.setup(A));

        Vector r(2), z(2, 0.0);
        r(0) = 4.0e-31;
        r(1) = 4.0e-31;
        EXPECT_TRUE(p.apply(r, z));
        EXPECT_NEAR(z(0), 8.0 / 11.0, 1e-12);
        EXPECT_NEAR(z(1), 12.0 / 11.0, 1e-12);
    });

    run_case("b1_block_jacobi_rejects_relative_pivot_floor", [] {
        // [[1,0],[1,1e-16]]: the 1e-16 pivot sits above the old absolute
        // 1e-30 cut but below 1e-14 * block_scale, so it must be rejected.
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 0, 1.0); A.push_back(1, 1, 1e-16);
        A.finalize();

        BlockDiagonalPreconditioner p({{0, 1}});
        EXPECT_TRUE(!p.setup(A));
    });

    run_case("b1_block_jacobi_exact_on_two_blocks", [] {
        SparseMatrix A(4, 4);
        A.push_back(0, 0, 4.0); A.push_back(0, 1, 1.0);
        A.push_back(1, 0, 1.0); A.push_back(1, 1, 3.0);
        A.push_back(2, 2, 2.0); A.push_back(2, 3, -1.0);
        A.push_back(3, 2, -1.0); A.push_back(3, 3, 2.0);
        A.finalize();

        BlockDiagonalPreconditioner p({{0, 1}, {2, 3}});
        EXPECT_TRUE(p.setup(A));

        Vector r(4), z(4, 0.0);
        r(0) = 5.0; r(1) = 11.0; r(2) = 3.0; r(3) = 6.0;
        EXPECT_TRUE(p.apply(r, z));
        EXPECT_NEAR(z(0), 4.0 / 11.0, 1e-12);
        EXPECT_NEAR(z(1), 39.0 / 11.0, 1e-12);
        EXPECT_NEAR(z(2), 4.0, 1e-12);
        EXPECT_NEAR(z(3), 5.0, 1e-12);
    });

    run_case("b1_block_jacobi_rejects_nonfinite_block_coefficient", [] {
        // SparseMatrix::push_back already rejects non-finite values at
        // assembly time, so the fail-closed B1 guards are exercised the
        // way an external numeric update would produce them: a finite
        // matrix is finalized, then one raw CSR coefficient is mutated
        // through values_data(). std::max then ignores the NaN during
        // the scale scan, so the per-coefficient finiteness check plus
        // the inverse finiteness gate are required to fail closed.
        SparseMatrix nan_block(2, 2);
        nan_block.push_back(0, 0, 1.0);
        nan_block.push_back(0, 1, 0.5);
        nan_block.push_back(1, 0, 0.0);
        nan_block.push_back(1, 1, 1.0);
        nan_block.finalize();
        nan_block.values_data()[1] =
            std::numeric_limits<double>::quiet_NaN();
        BlockDiagonalPreconditioner p({{0, 1}});
        EXPECT_TRUE(!p.setup(nan_block));

        SparseMatrix inf_block(2, 2);
        inf_block.push_back(0, 0, 1.0);
        inf_block.push_back(0, 1, 0.5);
        inf_block.push_back(1, 0, 0.0);
        inf_block.push_back(1, 1, 1.0);
        inf_block.finalize();
        inf_block.values_data()[1] =
            std::numeric_limits<double>::infinity();
        BlockDiagonalPreconditioner q({{0, 1}});
        EXPECT_TRUE(!q.setup(inf_block));

        // A failed setup must leave the preconditioner unusable rather
        // than applying a stale or partially built inverse.
        Vector r(2, 1.0), z(2, 0.0);
        EXPECT_TRUE(!p.apply(r, z));
    });

    run_case("b2_schur_diagonal_scale_fails_closed", [] {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 4.0); A.push_back(0, 1, 1.0);
        A.push_back(1, 0, 2.0); A.push_back(1, 1, 3.0);
        A.finalize();

        SchurComplementPreconditioner zero(
            {0}, {1}, SchurFactorization::Diagonal, 0.0);
        EXPECT_TRUE(!zero.setup(A));

        SchurComplementPreconditioner nonfinite(
            {0}, {1}, SchurFactorization::Diagonal,
            std::numeric_limits<double>::quiet_NaN());
        EXPECT_TRUE(!nonfinite.setup(A));

        // Non-Diagonal factorizations never read the scale, so scale 0.0
        // must stay accepted there (no behaviour change beyond B2).
        SchurComplementPreconditioner lower(
            {0}, {1}, SchurFactorization::Lower, 0.0);
        EXPECT_TRUE(lower.setup(A));

        // A finite non-zero Diagonal scale keeps its documented action:
        // y2 = scale * inv(S) r2 with S = A22 - A21 inv(A11) A12 = 2.5.
        SchurComplementPreconditioner minus(
            {0}, {1}, SchurFactorization::Diagonal, -1.0);
        EXPECT_TRUE(minus.setup(A));
        Vector r(2), z(2, 0.0);
        r(0) = 2.0; r(1) = 5.0;
        EXPECT_TRUE(minus.apply(r, z));
        EXPECT_NEAR(z(0), 0.5, 1e-12);
        EXPECT_NEAR(z(1), -2.0, 1e-12);
    });

    run_case("b3_gs_ilu0_accept_uniformly_rescaled_laplacian", [] {
        // The rescaled Laplacian is perfectly regular but every entry
        // sits far below the old absolute epsilon cut; both preconditioners
        // must accept it now that the guards are row-relative.
        const auto A = make_1d_laplacian(8, 1.0e-18);

        GaussSeidelPreconditioner gs;
        EXPECT_TRUE(gs.setup(A));
        ILU0Preconditioner ilu0;
        EXPECT_TRUE(ilu0.setup(A));

        Vector r(8, 1.0e-18), z(8, 0.0);
        EXPECT_TRUE(gs.apply(r, z));
        EXPECT_TRUE(z.is_valid());
        EXPECT_TRUE(ilu0.apply(r, z));
        EXPECT_TRUE(z.is_valid());
    });

    run_case("b3_ilu0_accepts_multiple_uniform_rescalings", [] {
        // The ILU(0) guards must be invariant to a uniform rescaling of a
        // regular matrix: the L multipliers stay O(1) while the U pivots
        // follow the matrix scale, so the reference scale must come from
        // the U part only.
        for (const double scale : {1.0, 1.0e-8, 1.0e-16, 1.0e-18}) {
            const auto A = make_1d_laplacian(8, scale);
            ILU0Preconditioner ilu0;
            EXPECT_TRUE(ilu0.setup(A));

            Vector r(8, 0.0), z(8, 0.0);
            for (std::size_t i = 0; i < r.size(); ++i)
                r(i) = scale * (1.0 + static_cast<double>(i));
            EXPECT_TRUE(ilu0.apply(r, z));
            EXPECT_TRUE(z.is_valid());
            // Zero fill on a tridiagonal pattern: ILU(0) is an exact LU
            // factorization, so the solve must satisfy A z = r at the
            // matrix scale.
            EXPECT_TRUE(residual_inf(A, z, r) <= 1.0e-12 * scale);
        }
    });

    run_case("b3_ilu0_accepts_contrasting_row_scales", [] {
        // Row 0 lives at 1e-18 while rows 1-2 live at O(1): the pivot of
        // row 0 must be judged against row 0's own factorized U scale, not
        // against the much larger scale of the row being eliminated.
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 2.0e-18); A.push_back(0, 1, -1.0e-18);
        A.push_back(1, 0, -1.0e-18); A.push_back(1, 1, 2.0); A.push_back(1, 2, -1.0);
        A.push_back(2, 1, -1.0); A.push_back(2, 2, 2.0);
        A.finalize();

        ILU0Preconditioner ilu0;
        EXPECT_TRUE(ilu0.setup(A));
        Vector r(3, 1.0), z(3, 0.0);
        EXPECT_TRUE(ilu0.apply(r, z));
        EXPECT_TRUE(z.is_valid());

        GaussSeidelPreconditioner gs;
        EXPECT_TRUE(gs.setup(A));
        EXPECT_TRUE(gs.apply(r, z));
        EXPECT_TRUE(z.is_valid());
    });

    run_case("b3_gs_forward_substitution_is_exact", [] {
        const auto A = make_1d_laplacian(8, 1.0);
        GaussSeidelPreconditioner gs;
        EXPECT_TRUE(gs.setup(A));

        Vector r(8), z(8, 0.0);
        for (std::size_t i = 0; i < r.size(); ++i)
            r(i) = 1.0 + static_cast<double>(i);
        EXPECT_TRUE(gs.apply(r, z));

        // GS applies (D + L)^-1 r: diag_i z_i + sum_{j<i} a_ij z_j = r_i.
        const auto* ro = A.row_offsets_data();
        const auto* ci = A.columns_data();
        const auto* av = A.values_data();
        for (std::size_t i = 0; i < 8; ++i) {
            double partial = 0.0;
            for (std::size_t k = ro[i]; k < ro[i + 1]; ++k)
                if (ci[k] <= i) partial += av[k] * z(ci[k]);
            EXPECT_NEAR(partial, r(i), 1e-12);
        }
    });

    run_case("b3_ilu0_is_exact_lu_on_zero_fill_pattern", [] {
        // The tridiagonal pattern needs no fill, so ILU(0) is an exact LU
        // factorization and the preconditioned solve must satisfy A z = r.
        const auto A = make_1d_laplacian(8, 1.0);
        ILU0Preconditioner ilu0;
        EXPECT_TRUE(ilu0.setup(A));

        Vector r(8), z(8, 0.0);
        for (std::size_t i = 0; i < r.size(); ++i)
            r(i) = 1.0 + static_cast<double>(i);
        EXPECT_TRUE(ilu0.apply(r, z));
        EXPECT_TRUE(z.is_valid());
        EXPECT_TRUE(residual_inf(A, z, r) < 1e-12);
    });

    run_case("b3_gs_rejects_row_negligible_diagonal", [] {
        // diag = 1e-10 used to pass the absolute epsilon cut while being
        // negligible against its own row scale 1e10; reject it now.
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0e-10); A.push_back(0, 1, 1.0e10);
        A.push_back(1, 1, 1.0);
        A.finalize();

        GaussSeidelPreconditioner gs;
        EXPECT_TRUE(!gs.setup(A));
    });

    run_case("b3_ilu0_rejects_singular_dense_block", [] {
        // [[1,1],[1,1]] factorizes to a zero trailing diagonal; the
        // row-relative floor must reject it.
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0); A.push_back(0, 1, 1.0);
        A.push_back(1, 0, 1.0); A.push_back(1, 1, 1.0);
        A.finalize();

        ILU0Preconditioner ilu0;
        EXPECT_TRUE(!ilu0.setup(A));
    });

    run_case("b3_gs_ilu0_reject_rows_without_diagonal", [] {
        SparseMatrix missing(2, 2);
        missing.push_back(0, 1, 1.0);
        missing.push_back(1, 0, 1.0);
        missing.finalize();

        GaussSeidelPreconditioner gs;
        EXPECT_TRUE(!gs.setup(missing));
        ILU0Preconditioner ilu0;
        EXPECT_TRUE(!ilu0.setup(missing));
    });

    run_case("b3_failed_setup_clears_stale_factorization", [] {
        // A successful setup followed by a failing one must not leave the
        // previous factorization usable: apply() has to fail on the
        // cleared state instead of silently reusing stale factors.
        const auto good = make_1d_laplacian(8, 1.0);
        ILU0Preconditioner ilu0;
        EXPECT_TRUE(ilu0.setup(good));
        Vector r(8, 1.0), z(8, 0.0);
        EXPECT_TRUE(ilu0.apply(r, z));

        SparseMatrix singular(2, 2);
        singular.push_back(0, 0, 1.0); singular.push_back(0, 1, 1.0);
        singular.push_back(1, 0, 1.0); singular.push_back(1, 1, 1.0);
        singular.finalize();
        EXPECT_TRUE(!ilu0.setup(singular));
        EXPECT_TRUE(!ilu0.apply(r, z));

        GaussSeidelPreconditioner gs;
        EXPECT_TRUE(gs.setup(good));
        SparseMatrix no_diagonal(2, 2);
        no_diagonal.push_back(0, 1, 1.0);
        no_diagonal.push_back(1, 0, 1.0);
        no_diagonal.finalize();
        EXPECT_TRUE(!gs.setup(no_diagonal));
        EXPECT_TRUE(!gs.apply(r, z));
    });

    return run_all();
}
