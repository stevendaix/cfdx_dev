#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/native_fieldsplit.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

SparseMatrix make_poisson(std::size_t n) {
    SparseMatrix matrix(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        matrix.push_back(i, i, 2.0);
        if (i > 0) matrix.push_back(i, i - 1, -1.0);
        if (i + 1 < n) matrix.push_back(i, i + 1, -1.0);
    }
    matrix.finalize();
    return matrix;
}

SparseMatrix make_two_field_system() {
    SparseMatrix matrix(4, 4);
    matrix.push_back(0, 0, 4.0); matrix.push_back(0, 1, 1.0);
    matrix.push_back(0, 2, 1.0);
    matrix.push_back(1, 0, 1.0); matrix.push_back(1, 1, 3.0);
    matrix.push_back(1, 3, 1.0);
    matrix.push_back(2, 0, 2.0); matrix.push_back(2, 2, 5.0);
    matrix.push_back(2, 3, 1.0);
    matrix.push_back(3, 1, 1.0); matrix.push_back(3, 2, 1.0);
    matrix.push_back(3, 3, 4.0);
    matrix.finalize();
    return matrix;
}

Vector multiply(const SparseMatrix& matrix, const Vector& x) {
    const auto values = matrix.matvec(x);
    Vector result(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) result(i) = values[i];
    return result;
}

double relative_true_residual(const SparseMatrix& matrix,
                              const Vector& x,
                              const Vector& rhs) {
    const auto ax = matrix.matvec(x);
    double residual2 = 0.0;
    double rhs2 = 0.0;
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        const double residual = rhs(i) - ax[i];
        residual2 += residual * residual;
        rhs2 += rhs(i) * rhs(i);
    }
    return std::sqrt(residual2 / std::max(rhs2, 1e-300));
}

} // namespace

int main() {
    run_case("boomer_style_amg_reduces_true_residual", [] {
        const auto matrix = make_poisson(64);
        Vector rhs(64);
        for (std::size_t i = 0; i < rhs.size(); ++i)
            rhs(i) = (i % 3 == 0) ? 1.0 : -0.25;

        NativeBoomerAMGPreconditioner amg;
        Vector correction;
        EXPECT_TRUE(!amg.apply(rhs, correction));
        EXPECT_TRUE(!amg.last_error().empty());
        EXPECT_TRUE(amg.setup(matrix));
        EXPECT_TRUE(amg.is_ready());
        EXPECT_TRUE(amg.coarse_size() > 0);
        EXPECT_TRUE(amg.coarse_size() < matrix.n_rows());
        EXPECT_TRUE(amg.apply(rhs, correction));
        EXPECT_TRUE(correction.is_valid());
        EXPECT_TRUE(relative_true_residual(matrix, correction, rhs) < 1.0);
    });

    run_case("cg_accepts_native_boomer_style_amg", [] {
        const auto matrix = make_poisson(96);
        Vector exact(96);
        for (std::size_t i = 0; i < exact.size(); ++i)
            exact(i) = std::sin(0.05 * static_cast<double>(i + 1));
        const Vector rhs = multiply(matrix, exact);

        NativeBoomerAMGPreconditioner amg;
        amg.configure(AMGMemoryPolicy::Balanced);
        Vector solution(96, 0.0);
        const auto result = solve_cg(matrix, rhs, solution, amg, 500, 1e-10);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(matrix, solution, rhs) < 1e-8);

        Vector native_solution(96, 0.0);
        const auto native = solve_cg(matrix, rhs, native_solution, 500, 1e-10);
        EXPECT_TRUE(native.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(matrix, native_solution, rhs) < 1e-8);
        EXPECT_TRUE((solution - native_solution).norm_inf() < 1e-7);
    });

    run_case("fieldsplit_schur_variants_are_well_defined", [] {
        const auto matrix = make_two_field_system();
        Vector rhs(4);
        rhs(0) = 1.0; rhs(1) = -2.0; rhs(2) = 3.0; rhs(3) = 0.5;

        const std::vector<FieldSplitStrategy> strategies{
            FieldSplitStrategy::Additive,
            FieldSplitStrategy::DiagonalSchur,
            FieldSplitStrategy::LowerSchur,
            FieldSplitStrategy::UpperSchur,
            FieldSplitStrategy::FullSchur};
        for (const auto strategy : strategies) {
            NativeFieldSplitPreconditioner split({0, 1}, {2, 3}, strategy);
            Vector correction(4, 0.0);
            EXPECT_TRUE(split.setup(matrix));
            EXPECT_TRUE(split.apply(rhs, correction));
            EXPECT_TRUE(correction.is_valid());
            EXPECT_TRUE(std::string(split.name()).find("NativeFieldSplit") == 0);
            if (strategy == FieldSplitStrategy::FullSchur)
                EXPECT_TRUE(relative_true_residual(matrix, correction, rhs) < 1e-12);
        }
    });

    run_case("gmres_uses_full_schur_fieldsplit", [] {
        const auto matrix = make_two_field_system();
        Vector exact(4);
        exact(0) = 1.0; exact(1) = -2.0; exact(2) = 0.5; exact(3) = 3.0;
        const Vector rhs = multiply(matrix, exact);

        NativeFieldSplitPreconditioner split(
            {0, 1}, {2, 3}, FieldSplitStrategy::FullSchur);
        Vector solution(4, 0.0);
        const auto result = solve_gmres(
            matrix, rhs, solution, 4, 20, 1e-12, &split);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(matrix, solution, rhs) < 1e-11);
        EXPECT_TRUE((solution - exact).norm_inf() < 1e-10);
    });

    run_case("fieldsplit_rejects_overlap_without_stale_success", [] {
        const auto matrix = make_two_field_system();
        NativeFieldSplitPreconditioner split(
            {0, 1}, {2, 3}, FieldSplitStrategy::LowerSchur);
        EXPECT_TRUE(split.setup(matrix));

        SparseMatrix wrong_size(3, 3);
        wrong_size.push_back(0, 0, 1.0);
        wrong_size.push_back(1, 1, 1.0);
        wrong_size.push_back(2, 2, 1.0);
        wrong_size.finalize();
        EXPECT_TRUE(!split.setup(wrong_size));
        EXPECT_TRUE(!split.last_error().empty());
        Vector rhs(4, 1.0), correction(4, 0.0);
        EXPECT_TRUE(!split.apply(rhs, correction));
    });

    return run_all();
}
