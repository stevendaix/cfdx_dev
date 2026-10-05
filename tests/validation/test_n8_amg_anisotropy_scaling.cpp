#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/linear_operator.h"
#include "cfdx/core/linalg/cg_solver.h"
#include "common/test_harness.h"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

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
private:
    SparseMatrix A_;
};

SparseMatrix make_anisotropic_diffusion(std::size_t nx,
                                        std::size_t ny,
                                        double ax,
                                        double ay,
                                        double reaction) {
    const std::size_t n = nx * ny;
    SparseMatrix A(n, n);
    for (std::size_t y = 0; y < ny; ++y) {
        for (std::size_t x = 0; x < nx; ++x) {
            const std::size_t i = y * nx + x;
            double diagonal = reaction;
            if (x > 0) {
                A.push_back(i, i - 1, -ax);
                diagonal += ax;
            }
            if (x + 1 < nx) {
                A.push_back(i, i + 1, -ax);
                diagonal += ax;
            }
            if (y > 0) {
                A.push_back(i, i - nx, -ay);
                diagonal += ay;
            }
            if (y + 1 < ny) {
                A.push_back(i, i + nx, -ay);
                diagonal += ay;
            }
            A.push_back(i, i, diagonal);
        }
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

template <typename AMG>
void run_amg_case(const SparseMatrix& A,
                  const Vector& b,
                  const char* method,
                  double anisotropy,
                  double scale,
                  bool& structural_ok) {
    TestSparseOperator op(A);
    AMG amg;
    const bool setup_ok = amg.setup(A);
    const auto levels = amg.hierarchy_level_sizes();

    std::cout << "n8_amg_sweep"
              << " method=" << method
              << " anisotropy=" << anisotropy
              << " scale=" << scale
              << " setup=" << (setup_ok ? "PASS" : "FAIL")
              << " levels=" << levels.size()
              << " fine=" << A.n_rows();

    if (!setup_ok) {
        std::cout << " structural=FAIL" << '\n';
        structural_ok = false;
        return;
    }

    Vector z(A.n_rows(), 0.0);
    const bool apply_ok = amg.apply(b, z);
    const double apply_residual = relative_true_residual(A, z, b);

    Vector x(A.n_rows(), 0.0);
    const auto result = solve_cg(A, b, x, amg, 5000, 1e-9);
    const double true_residual = relative_true_residual(A, x, b);

    bool finite_evidence = std::isfinite(apply_residual) &&
                           std::isfinite(true_residual);
    const bool hierarchy_ok = levels.size() >= 2 &&
                              amg.first_prolongation_nnz() > 0;

    std::cout << " apply=" << (apply_ok ? "PASS" : "FAIL")
              << " hierarchy=" << (hierarchy_ok ? "PASS" : "FAIL")
              << " P_nnz=" << amg.first_prolongation_nnz()
              << " apply_true_residual=" << apply_residual
              << " cg_status=" << static_cast<int>(result.status)
              << " cg_iterations=" << result.iterations
              << " cg_true_residual=" << true_residual
              << " evidence=" << (finite_evidence ? "FINITE" : "NONFINITE")
              << '\n';

    // This campaign establishes the measured applicability envelope; it does
    // not invent a universal iteration-count or convergence threshold. Setup,
    // hierarchy construction and finite independently recomputed residuals are
    // hard structural/evidence gates. Convergence is retained as measured
    // evidence for later policy qualification.
    structural_ok = structural_ok && apply_ok && hierarchy_ok && finite_evidence;
}

} // namespace

int main() {
    std::cout << std::setprecision(17);

    constexpr std::size_t nx = 32;
    constexpr std::size_t ny = 32;
    constexpr double reaction = 1.0;

    // Deterministic anisotropy ratios and global coefficient scales. The
    // operator is a 5-point finite-volume-like diffusion stencil with a
    // positive reaction term, so every matrix is nonsingular and SPD.
    constexpr double anisotropies[] = {1.0, 10.0, 100.0, 1000.0};
    constexpr double scales[] = {1e-3, 1.0, 1e3};

    bool structural_ok = true;

    for (const double anisotropy : anisotropies) {
        for (const double scale : scales) {
            const double ax = scale * anisotropy;
            const double ay = scale;
            const auto A = make_anisotropic_diffusion(nx, ny, ax, ay, reaction * scale);
            const auto b = make_rhs(A.n_rows());

            std::cout << "n8_amg_operator"
                      << " anisotropy=" << anisotropy
                      << " scale=" << scale
                      << " ax=" << ax
                      << " ay=" << ay
                      << " rows=" << A.n_rows()
                      << " cols=" << A.n_cols() << '\n';

            run_amg_case<NativeBoomerAMGPreconditioner>(
                A, b, "NativeBoomerAMG", anisotropy, scale, structural_ok);
        }
    }

    EXPECT_TRUE(structural_ok);
    std::cout << "N8 AMG anisotropy/scaling qualification: "
              << (structural_ok ? "PASS" : "FAIL") << '\n';
    return 0;
}