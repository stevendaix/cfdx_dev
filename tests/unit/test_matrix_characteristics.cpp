#include "cfdx/core/linalg/matrix_diagnostics.h"
#include "cfdx/core/linalg/linear_solver_models.h"
#include "common/test_harness.h"

#include <cmath>

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("measured_characteristics_feed_coupled_dispatch", [] {
        SparseMatrix matrix(3, 3);
        matrix.push_back(0, 0, 4.0);
        matrix.push_back(0, 1, -1.0);
        matrix.push_back(1, 0, -1.0);
        matrix.push_back(1, 1, 4.0);
        matrix.push_back(1, 2, -1.0);
        matrix.push_back(2, 1, -1.0);
        matrix.push_back(2, 2, 4.0);
        matrix.finalize();

        const auto characteristics = measure_matrix_characteristics(matrix, true);
        EXPECT_TRUE(characteristics.square);
        EXPECT_TRUE(characteristics.saddle_point);
        EXPECT_TRUE(characteristics.numerically_symmetric);
        EXPECT_TRUE(characteristics.diagonally_dominant);
        EXPECT_TRUE(std::abs(characteristics.average_nnz_per_row - 7.0 / 3.0) < 1e-14);

        LinearSolverRequest request;
        const auto plan = select_linear_solver(
            LinearProblemKind::CoupledPressureVelocity,
            characteristics,
            request);
        EXPECT_TRUE(plan.krylov == KrylovModel::FGMRES);
        EXPECT_TRUE(plan.preconditioner == PreconditionerModel::CoupledBlockSchur);
    });

    run_case("measured_characteristics_detect_nonsymmetry_and_scaling", [] {
        SparseMatrix matrix(2, 2);
        matrix.push_back(0, 0, 1.0e9);
        matrix.push_back(0, 1, 1.0);
        matrix.push_back(1, 0, 0.5);
        matrix.push_back(1, 1, 1.0e-2);
        matrix.finalize();

        const auto characteristics = measure_matrix_characteristics(matrix);
        EXPECT_TRUE(characteristics.square);
        EXPECT_TRUE(!characteristics.numerically_symmetric);
        EXPECT_TRUE(characteristics.strongly_scaled);
        EXPECT_TRUE(characteristics.coefficient_range > 1.0e10);
    });

    return 0;
}
