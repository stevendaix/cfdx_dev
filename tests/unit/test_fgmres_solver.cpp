#include "cfdx/core/linalg/gmres_solver.h"
#include "common/test_harness.h"

#include <cmath>
#include <cstddef>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

SparseMatrix make_nonsymmetric_system(std::size_t n) {
    SparseMatrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        A.push_back(i, i, 4.0 + 0.05 * static_cast<double>(i));
        if (i > 0) A.push_back(i, i - 1, -1.2);
        if (i + 1 < n) A.push_back(i, i + 1, -0.4);
    }
    A.finalize();
    return A;
}

Vector matvec(const SparseMatrix& A, const Vector& x) {
    const auto values = A.matvec(x);
    Vector y(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) y(i) = values[i];
    return y;
}

double true_relative_residual(
    const SparseMatrix& A, const Vector& x, const Vector& b) {
    const auto ax = A.matvec(x);
    double rr2 = 0.0;
    double bb2 = 0.0;
    for (std::size_t i = 0; i < b.size(); ++i) {
        const double r = ax[i] - b(i);
        rr2 += r * r;
        bb2 += b(i) * b(i);
    }
    return std::sqrt(rr2 / std::max(bb2, 1e-300));
}

// Deliberately variable preconditioner: the map changes on every application.
// This is a contract test for FGMRES, not a performance benchmark.
class AlternatingPreconditioner final : public Preconditioner {
public:
    bool setup(const SparseMatrix& A) override {
        if (A.n_rows() != A.n_cols() || A.n_rows() == 0) return false;
        inv_diag_.assign(A.n_rows(), 0.0);
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            bool found = false;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                if (col[k] == i) {
                    if (val[k] == 0.0) return false;
                    inv_diag_[i] = 1.0 / val[k];
                    found = true;
                    break;
                }
            }
            if (!found) return false;
        }
        calls_ = 0;
        return true;
    }

    bool apply(const Vector& r, Vector& z) const override {
        if (r.size() != inv_diag_.size() || z.size() != r.size()) return false;
        ++calls_;
        if ((calls_ & 1U) != 0U) {
            for (std::size_t i = 0; i < r.size(); ++i) z(i) = r(i);
        } else {
            for (std::size_t i = 0; i < r.size(); ++i) z(i) = inv_diag_[i] * r(i);
        }
        return true;
    }

    const char* name() const override { return "Alternating"; }
    std::size_t calls() const { return calls_; }

private:
    std::vector<double> inv_diag_;
    mutable std::size_t calls_{0};
};

} // namespace

int main() {
    run_case("fgmres_accepts_variable_preconditioner_and_true_residual", [] {
        const auto A = make_nonsymmetric_system(32);
        Vector exact(32);
        for (std::size_t i = 0; i < exact.size(); ++i)
            exact(i) = std::sin(0.13 * static_cast<double>(i + 1));
        const auto b = matvec(A, exact);

        AlternatingPreconditioner pc;
        Vector x(32, 0.0);
        const auto result = solve_fgmres(A, b, x, 8, 200, 1e-11, &pc);

        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(pc.calls() > 1);
        EXPECT_TRUE(true_relative_residual(A, x, b) < 1e-9);
        EXPECT_TRUE(std::isfinite(result.residual));
        EXPECT_TRUE(std::abs(result.residual_relative -
                              true_relative_residual(A, x, b)) < 1e-10);
    });

    run_case("fgmres_matches_gmres_on_stationary_preconditioner", [] {
        const auto A = make_nonsymmetric_system(24);
        Vector exact(24);
        for (std::size_t i = 0; i < exact.size(); ++i)
            exact(i) = 1.0 - 0.02 * static_cast<double>(i);
        const auto b = matvec(A, exact);

        JacobiPreconditioner pc_gmres;
        JacobiPreconditioner pc_fgmres;
        Vector x_gmres(24, 0.0);
        Vector x_fgmres(24, 0.0);

        const auto gmres = solve_gmres(A, b, x_gmres, 8, 200, 1e-11, &pc_gmres);
        const auto fgmres = solve_fgmres(A, b, x_fgmres, 8, 200, 1e-11, &pc_fgmres);

        EXPECT_TRUE(gmres.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(fgmres.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(true_relative_residual(A, x_gmres, b) < 1e-9);
        EXPECT_TRUE(true_relative_residual(A, x_fgmres, b) < 1e-9);
        EXPECT_TRUE((x_gmres - x_fgmres).norm_inf() < 1e-10);
    });

    return run_all();
}
