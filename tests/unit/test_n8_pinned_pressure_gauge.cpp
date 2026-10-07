#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/exact_schur.h"
#include "cfdx/core/linalg/lsc_bfbt_schur.h"
#include "cfdx/core/linalg/null_space.h"
#include "cfdx/core/linalg/pcd_schur.h"
#include "cfdx/core/linalg/simplerc_schur.h"
#include "cfdx/core/linalg/schur_approximation.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

#include <initializer_list>
#include <optional>
#include <tuple>

using namespace cfdx::core;
using namespace cfdx::testing;

static SparseMatrix make_sparse(
    std::size_t rows,
    std::size_t cols,
    const std::initializer_list<std::tuple<std::size_t, std::size_t, double>>& entries)
{
    SparseMatrix matrix(rows, cols);
    for (const auto& [row, col, value] : entries)
        matrix.push_back(row, col, value);
    matrix.finalize();
    return matrix;
}

static void expect_gauge_component(
    const Vector& rhs,
    const Vector& out,
    double expected_sign)
{
    EXPECT_TRUE(rhs.size() == 2 && out.size() == 2);
    EXPECT_NEAR(out(0), expected_sign * rhs(0), 1e-12);
    EXPECT_TRUE(std::isfinite(out(0)) && std::isfinite(out(1)));
}

int main()
{
    // Reference cell 0 is pinned. The gauge row/column is represented in the
    // pressure block as C(0,0)=1, with no velocity-pressure coupling in that
    // row/column. The physical pressure equation remains non-trivial on cell 1.
    //
    // This is deliberately a pinned system, not a singular mean-zero system:
    // applying NullSpaceProjector::constant() to a RHS with a non-zero reference
    // component would be mathematically incompatible and must not be used.
    const auto Auu = make_sparse(2, 2, {
        {0, 0, 1.0}, {1, 1, 1.0}
    });
    const auto G = make_sparse(2, 2, {
        {1, 1, 1.0}
    });
    const auto D = make_sparse(2, 2, {
        {1, 1, 1.0}
    });
    const auto C = make_sparse(2, 2, {
        {0, 0, 1.0}, {1, 1, 2.0}
    });
    const BlockOperator blocks(Auu, G, D, C);
    blocks.validate();

    Vector rhs(2, 0.0);
    rhs(0) = 0.75;
    rhs(1) = -0.25;

    run_case("pinned_gauge_is_not_a_mean_zero_null_space", [&] {
        const auto projector = NullSpaceProjector::constant(2);
        EXPECT_TRUE(!projector.is_compatible(rhs));

        Vector projected = rhs;
        projector.remove(projected);
        EXPECT_NEAR(projected(0), 0.5, 1e-12);

        // Projection changes the requested pinned gauge component, so it cannot
        // be used as the policy for this nonsingular reference-cell system.
        EXPECT_TRUE(std::abs(projected(0) - rhs(0)) > 1e-12);
    });

    run_case("exact_schur_preserves_pinned_gauge_component", [&] {
        ExactSchurApproximation exact(
            [](const Vector& input, Vector& output) {
                output = input;
                return true;
            });
        EXPECT_TRUE(exact.setup(blocks));

        Vector out(2, 0.0);
        EXPECT_TRUE(exact.apply(rhs, out));
        expect_gauge_component(rhs, out, 1.0);

        Vector schur(2, 0.0);
        EXPECT_TRUE(exact.apply_schur(rhs, schur));
        expect_gauge_component(rhs, schur, 1.0);
    });

    run_case("simple_and_simplec_preserve_pinned_gauge_component", [&] {
        for (const auto mode : {SimplerSchurMode::SIMPLE, SimplerSchurMode::SIMPLEC}) {
            SimplerSchurApproximation schur(mode);
            EXPECT_TRUE(schur.setup(blocks));

            Vector out(2, 0.0);
            EXPECT_TRUE(schur.apply(rhs, out));
            expect_gauge_component(rhs, out, 1.0);

            const auto assembled = schur.assembled_operator();
            EXPECT_TRUE(assembled.n_rows() == 2 && assembled.n_cols() == 2);
            EXPECT_NEAR(assembled.matvec(rhs)(0), rhs(0), 1e-12);
        }
    });

    run_case("lsc_and_bfbt_preserve_pinned_gauge_component", [&] {
        auto solve_pressure = [](const Vector& input, Vector& output) {
            output = input;
            return true;
        };

        for (const auto mode : {
                 LscBfbtSchurApproximation::Mode::LSC,
                 LscBfbtSchurApproximation::Mode::BFBT}) {
            LscBfbtSchurApproximation schur(mode, solve_pressure);
            EXPECT_TRUE(schur.setup(blocks));

            Vector out(2, 0.0);
            EXPECT_TRUE(schur.apply(rhs, out));
            expect_gauge_component(rhs, out, 1.0);
        }
    });

    run_case("pcd_preserves_pinned_gauge_component_with_reference_row", [&] {
        // PCD's pressure-side operators are already gauge-pinned by imposing a
        // unit row at the reference cell. Its inverse action carries the sign
        // convention of S^{-1} ~= -Fp^{-1}, so the gauge component is -rhs_ref.
        const auto mass = make_sparse(2, 2, {
            {0, 0, 1.0}, {1, 1, 1.0}
        });
        const auto laplacian = make_sparse(2, 2, {
            {0, 0, 1.0}, {1, 1, 1.0}
        });
        const auto convection_diffusion = make_sparse(2, 2, {
            {0, 0, 1.0}, {1, 1, 1.0}
        });

        auto solve_pressure = [](const Vector& input, Vector& output) {
            output = input;
            return true;
        };

        PcdSchurApproximation pcd(
            mass,
            laplacian,
            convection_diffusion,
            solve_pressure,
            solve_pressure,
            std::nullopt,
            0);
        EXPECT_TRUE(pcd.setup(blocks));

        Vector out(2, 0.0);
        EXPECT_TRUE(pcd.apply(rhs, out));
        expect_gauge_component(rhs, out, -1.0);
    });

    return run_all();
}
