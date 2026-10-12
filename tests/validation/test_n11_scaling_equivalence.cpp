#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/matrix_diagnostics.h"
#include "cfdx/core/linalg/preconditioner.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

using namespace cfdx::core;
using namespace cfdx::testing;

static SparseMatrix make_scaled_system() {
    SparseMatrix A(4, 4);
    A.push_back(0, 0, 2.0e-8);
    A.push_back(0, 1, -1.0e-8);
    A.push_back(1, 0, -3.0e3);
    A.push_back(1, 1, 4.0e3);
    A.push_back(1, 2, -1.0e3);
    A.push_back(2, 1, -2.0e-2);
    A.push_back(2, 2, 5.0e-2);
    A.push_back(2, 3, -1.0e-2);
    A.push_back(3, 2, -7.0e4);
    A.push_back(3, 3, 9.0e4);
    A.finalize();
    return A;
}

static Vector make_reference_solution() {
    Vector x(4);
    x(0) = 0.75;
    x(1) = -1.25;
    x(2) = 2.0;
    x(3) = -0.5;
    return x;
}

static Vector rhs_from(const SparseMatrix& A, const Vector& x) {
    const auto values = A.matvec(x);
    Vector b(x.size());
    for (std::size_t i = 0; i < values.size(); ++i) b(i) = values[i];
    return b;
}

static double relative_residual(const SparseMatrix& A,
                                const Vector& x,
                                const Vector& b) {
    const auto Ax = A.matvec(x);
    double r2 = 0.0;
    double b2 = 0.0;
    for (std::size_t i = 0; i < Ax.size(); ++i) {
        const double r = Ax[i] - b(i);
        r2 += r * r;
        b2 += b(i) * b(i);
    }
    return std::sqrt(r2) / std::max(std::sqrt(b2), 1.0e-300);
}

int main() {
    run_case("row_column_scaling_is_an_exact_algebraic_transform", [] {
        const auto A = make_scaled_system();
        const auto x = make_reference_solution();
        const auto b = rhs_from(A, x);
        const auto scaled = scale_matrix(A, MatrixScaling::RowColumn);

        Vector x_hat(4);
        Vector b_hat(4);
        for (std::size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(std::isfinite(scaled.row_scale[i]));
            EXPECT_TRUE(std::isfinite(scaled.column_scale[i]));
            EXPECT_TRUE(scaled.row_scale[i] > 0.0);
            EXPECT_TRUE(scaled.column_scale[i] > 0.0);
            b_hat(i) = scaled.row_scale[i] * b(i);
            x_hat(i) = x(i) / scaled.column_scale[i];
        }

        const auto transformed = scaled.matrix.matvec(x_hat);
        for (std::size_t i = 0; i < 4; ++i)
            EXPECT_NEAR(transformed[i], b_hat(i), 1.0e-12 * std::max(1.0, std::abs(b_hat(i))));
    });

    run_case("scaled_gmres_solution_maps_back_to_physical_solution", [] {
        const auto A = make_scaled_system();
        const auto x_exact = make_reference_solution();
        const auto b = rhs_from(A, x_exact);
        const auto scaled = scale_matrix(A, MatrixScaling::RowColumn);

        Vector b_hat(4);
        for (std::size_t i = 0; i < 4; ++i)
            b_hat(i) = scaled.row_scale[i] * b(i);

        Vector x_hat(4, 0.0);
        JacobiPreconditioner preconditioner;
        const auto result = solve_gmres(
            scaled.matrix, b_hat, x_hat, 40, 200, 1.0e-11, &preconditioner);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);

        Vector x_physical(4);
        for (std::size_t i = 0; i < 4; ++i)
            x_physical(i) = scaled.column_scale[i] * x_hat(i);

        EXPECT_TRUE(relative_residual(A, x_physical, b) < 1.0e-10);
        for (std::size_t i = 0; i < 4; ++i)
            EXPECT_NEAR(x_physical(i), x_exact(i), 1.0e-9);
    });

    run_case("scaling_reports_zero_rows_without_repair", [] {
        SparseMatrix A(2, 2);
        A.push_back(1, 1, 3.0);
        A.finalize();

        const auto scaled = scale_matrix(A, MatrixScaling::RowColumn);
        EXPECT_TRUE(scaled.zero_rows == 1);
        EXPECT_TRUE(scaled.zero_columns == 1);
        EXPECT_NEAR(scaled.matrix(0, 0), 0.0, 0.0);
    });

    return run_all();
}
