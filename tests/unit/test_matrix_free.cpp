// M0.8-T10 — Tests for matrix-free operators

#include "cfdx/core/numerics/matrix_free.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("matrix_free_diagonal", []() {
        MatrixFreeOperator op(
            [](const Vector& x, Vector& y) {
                y(0) = 2.0 * x(0);
                y(1) = 3.0 * x(1);
                y(2) = 4.0 * x(2);
            },
            3, 3);

        EXPECT_TRUE(op.valid());
        EXPECT_TRUE(op.n_rows() == 3);
        EXPECT_TRUE(op.n_cols() == 3);

        Vector x(3);
        x(0) = 1.0; x(1) = 2.0; x(2) = 3.0;
        Vector y;
        op.apply(x, y);

        EXPECT_NEAR(y(0), 2.0, 1e-12);
        EXPECT_NEAR(y(1), 6.0, 1e-12);
        EXPECT_NEAR(y(2), 12.0, 1e-12);
    });

    run_case("matrix_free_invalid_operator", []() {
        MatrixFreeOperator op;
        EXPECT_TRUE(!op.valid());

        Vector x(1);
        Vector y;
        bool threw = false;
        try {
            op.apply(x, y);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        EXPECT_TRUE(threw);
    });

    return run_all();
}
