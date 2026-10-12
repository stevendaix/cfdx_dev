#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/lsc_bfbt_schur.h"
#include "cfdx/core/linalg/null_space.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

#include <cmath>
#include <initializer_list>
#include <tuple>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

static SparseMatrix make_sparse(
    std::size_t rows, std::size_t cols,
    const std::initializer_list<std::tuple<std::size_t, std::size_t, double>>& entries) {
    SparseMatrix A(rows, cols);
    for (const auto& [i, j, value] : entries) A.push_back(i, j, value);
    A.finalize();
    return A;
}

static double norm2(const Vector& x) {
    return std::sqrt(x.dot(x));
}

int main() {
    // The pressure operator P=D Q^-1 G has the constant mode in its null
    // space. The pressure solve callback represents a gauge-fixed solve on
    // the compatible mean-zero subspace.
    const auto Auu = make_sparse(2, 2, {
        {0, 0, 4.0}, {1, 1, 3.0}
    });
    const auto G = make_sparse(2, 2, {
        {0, 0, 1.0}, {0, 1, -1.0},
        {1, 0, -1.0}, {1, 1, 1.0}
    });
    const auto D = G;
    const auto C = make_sparse(2, 2, {});
    const BlockOperator blocks(Auu, G, D, C);
    blocks.validate();

    const auto null_space = NullSpaceProjector::constant(2);

    // Build the exact pressure-side solve on the compatible subspace for
    // P = D Q^-1 G. The constant mode is singular and is removed explicitly.
    const auto make_pressure_solve = [](double q0, double q1) {
        const double lambda = q0 + q1;
        return [lambda](const Vector& rhs, Vector& z) {
            if (rhs.size() != 2) return false;
            if (z.size() != 2) z = Vector(2, 0.0);
            const double mode = 0.5 * (rhs(0) - rhs(1));
            z(0) = mode / lambda;
            z(1) = -mode / lambda;
            return true;
        };
    };
    const auto lsc_pressure_solve = make_pressure_solve(0.25, 1.0 / 3.0);
    const auto bfbt_pressure_solve = make_pressure_solve(1.0, 1.0);

    run_case("lsc_declares_and_applies_pressure_null_space_policy", [&] {
        LscBfbtSchurApproximation lsc(
            LscBfbtSchurApproximation::Mode::LSC,
            lsc_pressure_solve,
            {},
            null_space);
        EXPECT_TRUE(lsc.setup(blocks));
        EXPECT_TRUE(lsc.has_pressure_null_space_policy());

        Vector rhs(2, 0.0);
        rhs(0) = 1.0;
        rhs(1) = -1.0;
        Vector pressure(2, 7.0);
        EXPECT_TRUE(lsc.apply(rhs, pressure));
        EXPECT_NEAR(pressure(0) + pressure(1), 0.0, 1e-14);
        EXPECT_TRUE(std::isfinite(pressure(0)) && std::isfinite(pressure(1)));
    });

    run_case("bfbt_declares_and_applies_pressure_null_space_policy", [&] {
        LscBfbtSchurApproximation bfbt(
            LscBfbtSchurApproximation::Mode::BFBT,
            bfbt_pressure_solve,
            {1.0, 1.0},
            null_space);
        EXPECT_TRUE(bfbt.setup(blocks));
        EXPECT_TRUE(bfbt.has_pressure_null_space_policy());

        Vector rhs(2, 0.0);
        rhs(0) = 0.75;
        rhs(1) = -0.75; // compatible with the constant pressure null-space
        Vector pressure(2, 0.0);
        EXPECT_TRUE(bfbt.apply(rhs, pressure));
        EXPECT_NEAR(pressure(0) + pressure(1), 0.0, 1e-14);
        EXPECT_TRUE(norm2(pressure) > 0.0);
    });

    run_case("null_space_policy_rejects_incompatible_pressure_rhs", [&] {
        LscBfbtSchurApproximation lsc(
            LscBfbtSchurApproximation::Mode::LSC,
            lsc_pressure_solve,
            {},
            null_space);
        EXPECT_TRUE(lsc.setup(blocks));

        Vector rhs(2, 1.0);
        Vector pressure(2, 0.0);
        EXPECT_TRUE(!lsc.apply(rhs, pressure));
        EXPECT_NEAR(norm2(pressure), 0.0, 0.0);
    });

    run_case("null_space_policy_rejects_dimension_mismatch", [&] {
        const auto wrong_null_space = NullSpaceProjector::constant(3);
        LscBfbtSchurApproximation lsc(
            LscBfbtSchurApproximation::Mode::LSC,
            lsc_pressure_solve,
            {},
            wrong_null_space);
        EXPECT_TRUE(!lsc.setup(blocks));
    });

    return run_all();
}
