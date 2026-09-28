#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "cfdx/core/linalg/advanced_preconditioners.h"
#include "cfdx/core/linalg/amg_preconditioner.h"
#include "common/test_harness.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include "cfdx/core/linalg/linear_operator.h"

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

SparseMatrix make_poisson(std::size_t n, double scale = 1.0) {
    SparseMatrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        A.push_back(i, i, 2.0 * scale);
        if (i > 0) A.push_back(i, i - 1, -scale);
        if (i + 1 < n) A.push_back(i, i + 1, -scale);
    }
    A.finalize();
    return A;
}

SparseMatrix make_poisson_changed_pattern(std::size_t n, double scale = 1.0) {
    SparseMatrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        A.push_back(i, i, 2.0 * scale);
        if (i > 0) A.push_back(i, i - 1, -scale);
        if (i + 1 < n) A.push_back(i, i + 1, -scale);
    }
    A.push_back(0, n - 1, 0.125);
    A.finalize();
    return A;
}

SparseMatrix make_rhs_matrix(std::size_t n) {
    // A diagonally dominant tridiagonal matrix representative of a
    // convection/diffusion-like nonsymmetric transport operator.
    SparseMatrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        A.push_back(i, i, 2.4);
        if (i > 0) A.push_back(i, i - 1, -1.1);
        if (i + 1 < n) A.push_back(i, i + 1, -0.9);
    }
    A.finalize();
    return A;
}

Vector make_rhs(std::size_t n) {
    Vector b(n);
    for (std::size_t i = 0; i < n; ++i)
        b(i) = 1.0 + 0.25 * std::sin(0.13 * static_cast<double>(i));
    return b;
}

double relative_true_residual(const SparseMatrix& A,
                              const Vector& x,
                              const Vector& b) {
    const auto ax = A.matvec(x);
    double r2 = 0.0;
    double b2 = 0.0;
    for (std::size_t i = 0; i < b.size(); ++i) {
        const double r = b(i) - ax[i];
        r2 += r * r;
        b2 += b(i) * b(i);
    }
    return std::sqrt(r2 / std::max(b2, 1e-300));
}

struct BenchmarkRow {
    std::size_t n;
    const char* method;
    SolverStatus status;
    std::size_t iterations;
    double residual;
    double milliseconds;
};

template <typename Solver>
BenchmarkRow run_case(const char* name,
                      const SparseMatrix& A,
                      const Vector& b,
                      Solver&& solve) {
    Vector x(A.n_rows(), 0.0);
    const auto t0 = std::chrono::steady_clock::now();
    const SolverResult result = solve(x);
    const auto t1 = std::chrono::steady_clock::now();
    const double ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double true_r = relative_true_residual(A, x, b);
    std::cout << name << " n=" << A.n_rows()
              << " status=" << static_cast<int>(result.status)
              << " iterations=" << result.iterations
              << " reported_residual=" << result.residual
              << " true_residual=" << true_r
              << " ms=" << ms << '\n';
    return {A.n_rows(), name, result.status, result.iterations, true_r, ms};
}

class TestSparseOperator final : public LinearOperatorBase {
public:
    explicit TestSparseOperator(const SparseMatrix& A) : A_(A) {}

    std::size_t rows() const noexcept override { return A_.n_rows(); }
    std::size_t cols() const noexcept override { return A_.n_cols(); }

    void apply(const Vector& x, Vector& y) const override {
        const auto values = A_.matvec(x);
        if (y.size() != values.size()) y.resize(values.size());
        for (std::size_t i = 0; i < values.size(); ++i) y(i) = values[i];
    }

    void update(const SparseMatrix& A) { A_ = A; }

private:
    SparseMatrix A_;
};

} // namespace

int main() {
    run_case("amg_vs_jacobi_poisson", make_poisson(256), make_rhs(256),
             [](Vector& x) {
                 const auto A = make_poisson(256);
                 const auto b = make_rhs(256);
                 JacobiPreconditioner jacobi;
                 return solve_cg(A, b, x, jacobi, 5000, 1e-10);
             });

    run_case("amg_poisson", make_poisson(256), make_rhs(256),
             [](Vector& x) {
                 const auto A = make_poisson(256);
                 const auto b = make_rhs(256);
                 NativeBoomerAMGPreconditioner amg;
                 return solve_cg(A, b, x, amg, 5000, 1e-10);
             });

    run_case("ilu0_transport", make_rhs_matrix(256), make_rhs(256),
             [](Vector& x) {
                 const auto A = make_rhs_matrix(256);
                 const auto b = make_rhs(256);
                 ILU0Preconditioner ilu;
                 return solve_bicgstab(A, b, x, 5000, 1e-10, &ilu);
             });

    run_case("jacobi_transport", make_rhs_matrix(256), make_rhs(256),
             [](Vector& x) {
                 const auto A = make_rhs_matrix(256);
                 const auto b = make_rhs(256);
                 JacobiPreconditioner jacobi;
                 return solve_bicgstab(A, b, x, 5000, 1e-10, &jacobi);
             });

    run_case("dispatch_pressure_amg", make_poisson(64), make_rhs(64),
             [](Vector& x) {
                 const auto A = make_poisson(64);
                 const auto b = make_rhs(64);
                 const auto report = solve_linear_system(
                     A, b, x, LinearProblemKind::PressurePoisson, {}, 5000, 1e-10);
                 EXPECT_TRUE(report.plan.preconditioner == PreconditionerModel::NativeAMG);
                 return report.result;
             });

    run_case("dispatch_momentum_ilu0", make_rhs_matrix(64), make_rhs(64),
             [](Vector& x) {
                 const auto A = make_rhs_matrix(64);
                 const auto b = make_rhs(64);
                 const auto report = solve_linear_system(
                     A, b, x, LinearProblemKind::Momentum, {}, 5000, 1e-10);
                 EXPECT_TRUE(report.plan.preconditioner == PreconditionerModel::ILU0);
                 return report.result;
             });

    run_case("explicit_unavailable_preconditioner_is_rejected",
             make_poisson(64), make_rhs(64),
             [](Vector&) {
                 SolverResult result;
                 try {
                     LinearSolverRequest request;
                     request.preconditioner = PreconditionerModel::ILUT;
                     (void)select_linear_solver(
                         LinearProblemKind::PressurePoisson, 64, request);
                     EXPECT_TRUE(false);
                 } catch (const std::invalid_argument&) {
                     result.status = SolverStatus::NOT_APPLICABLE;
                 }
                 return result;
             });

    {
        const auto A = make_poisson(64);
        auto A2 = make_poisson(64, 1.25);
        NativeBoomerAMGPreconditioner uninitialized;
        EXPECT_TRUE(!uninitialized.update_values(A));
        EXPECT_TRUE(!uninitialized.is_ready());

        NativeBoomerAMGPreconditioner amg;
        EXPECT_TRUE(amg.setup(A));
        EXPECT_TRUE(amg.hierarchy_builds() == 1);
        EXPECT_TRUE(amg.numeric_updates() == 0);
        EXPECT_TRUE(amg.update_values(A2));
        EXPECT_TRUE(amg.hierarchy_builds() == 1);
        EXPECT_TRUE(amg.numeric_updates() == 1);

        const auto changed_pattern = make_poisson_changed_pattern(64, 1.25);
        EXPECT_TRUE(!amg.update_values(changed_pattern));
        EXPECT_TRUE(amg.last_error().find("unchanged CSR pattern") != std::string::npos);
        EXPECT_TRUE(amg.is_ready());
        EXPECT_TRUE(amg.hierarchy_builds() == 1);
        EXPECT_TRUE(amg.numeric_updates() == 1);

        Vector b = make_rhs(64);
        Vector x;
        EXPECT_TRUE(amg.apply(b, x));
        EXPECT_TRUE(relative_true_residual(A2, x, b) < 1.0);
    }

    {
        const auto A = make_poisson(64);
        const auto A2 = make_poisson(64, 1.25);
        const auto changed = make_poisson_changed_pattern(64, 1.25);
        TestSparseOperator op(A);
        MatrixFreeVcyclePreconditioner vcycle(op);
        EXPECT_TRUE(vcycle.setup(A));

        EXPECT_TRUE(vcycle.update_values(A2));
        Vector b = make_rhs(64);
        Vector z;
        EXPECT_TRUE(vcycle.apply(b, z));
        EXPECT_TRUE(relative_true_residual(A2, z, b) < 1.0);

        EXPECT_TRUE(!vcycle.update_values(changed));
        op.update(A2);
        EXPECT_TRUE(vcycle.apply(b, z));
        EXPECT_TRUE(relative_true_residual(A2, z, b) < 1.0);
    }

    {
        const auto A = make_poisson(256);
        const auto b = make_rhs(256);
        Vector x_j(256, 0.0);
        Vector x_a(256, 0.0);
        JacobiPreconditioner jacobi;
        NativeBoomerAMGPreconditioner amg;
        const auto rj = solve_cg(A, b, x_j, jacobi, 5000, 1e-10);
        const auto ra = solve_cg(A, b, x_a, amg, 5000, 1e-10);
        EXPECT_TRUE(rj.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(ra.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(A, x_j, b) < 1e-9);
        EXPECT_TRUE(relative_true_residual(A, x_a, b) < 1e-9);
        EXPECT_TRUE(ra.iterations < rj.iterations);
    }

    std::cout << "AMG/preconditioner qualification: PASS\n";
    return 0;
}
