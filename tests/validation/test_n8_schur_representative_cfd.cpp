#include "cfdx/core/linalg/exact_schur.h"
#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/lsc_bfbt_schur.h"
#include "cfdx/core/linalg/simplerc_schur.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <tuple>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

using Dense = std::vector<std::vector<double>>;

struct RepresentativeCase {
    std::size_t n = 0;
    double convection = 0.0;
};

struct Blocks {
    SparseMatrix Auu;
    SparseMatrix G;
    SparseMatrix D;
    SparseMatrix C;
};

SparseMatrix make_sparse(
    std::size_t rows, std::size_t cols,
    const std::vector<std::tuple<std::size_t, std::size_t, double>>& entries) {
    SparseMatrix A(rows, cols);
    for (const auto& [i, j, value] : entries) A.push_back(i, j, value);
    A.finalize();
    return A;
}

Blocks make_representative_case(const RepresentativeCase& cfg) {
    const std::size_t np = cfg.n * cfg.n;
    const std::size_t nu = 2 * np;
    std::vector<std::tuple<std::size_t, std::size_t, double>> auu_entries;
    std::vector<std::tuple<std::size_t, std::size_t, double>> g_entries;
    std::vector<std::tuple<std::size_t, std::size_t, double>> d_entries;
    std::vector<std::tuple<std::size_t, std::size_t, double>> c_entries;

    const auto cell = [n = cfg.n](std::size_t x, std::size_t y) {
        return y * n + x;
    };

    const double diffusion = 1.0;
    const double reaction = 0.5;
    const double beta = cfg.convection;

    // Two velocity components. The operator is a finite-volume-like
    // diffusion/reaction block with a controlled upwind convection bias.
    for (std::size_t y = 0; y < cfg.n; ++y) {
        for (std::size_t x = 0; x < cfg.n; ++x) {
            const std::size_t p = cell(x, y);
            for (std::size_t component = 0; component < 2; ++component) {
                const std::size_t row = component * np + p;
                double diagonal = reaction + 4.0 * diffusion + beta;
                auu_entries.emplace_back(row, row, diagonal);

                if (x > 0) auu_entries.emplace_back(
                    row, component * np + cell(x - 1, y), -diffusion - beta);
                if (x + 1 < cfg.n) auu_entries.emplace_back(
                    row, component * np + cell(x + 1, y), -diffusion);
                if (y > 0) auu_entries.emplace_back(
                    row, component * np + cell(x, y - 1), -diffusion);
                if (y + 1 < cfg.n) auu_entries.emplace_back(
                    row, component * np + cell(x, y + 1), -diffusion);
            }

            // Discrete pressure gradients: central in the interior and
            // one-sided at the boundary.
            const std::size_t ux = p;
            const std::size_t uy = np + p;
            if (x > 0)
                g_entries.emplace_back(ux, cell(x - 1, y), -0.5);
            else
                g_entries.emplace_back(ux, p, -1.0);
            if (x + 1 < cfg.n)
                g_entries.emplace_back(ux, cell(x + 1, y), 0.5);
            else
                g_entries.emplace_back(ux, p, 1.0);
            if (y > 0)
                g_entries.emplace_back(uy, cell(x, y - 1), -0.5);
            else
                g_entries.emplace_back(uy, p, -1.0);
            if (y + 1 < cfg.n)
                g_entries.emplace_back(uy, cell(x, y + 1), 0.5);
            else
                g_entries.emplace_back(uy, p, 1.0);

            // Deliberately non-transposed D: this models the algebraic
            // asymmetry that collocated finite-volume operators can exhibit.
            d_entries.emplace_back(p, ux, -0.95);
            d_entries.emplace_back(p, uy, -1.0);
            c_entries.emplace_back(p, p, 0.02);
        }
    }

    return {
        make_sparse(nu, nu, auu_entries),
        make_sparse(nu, np, g_entries),
        make_sparse(np, nu, d_entries),
        make_sparse(np, np, c_entries),
    };
}

Dense dense(const SparseMatrix& A) {
    Dense out(A.n_rows(), std::vector<double>(A.n_cols(), 0.0));
    for (std::size_t i = 0; i < A.n_rows(); ++i) {
        for (std::size_t k = A.row_offsets_data()[i];
             k < A.row_offsets_data()[i + 1]; ++k) {
            out[i][A.columns_data()[k]] += A.values_data()[k];
        }
    }
    return out;
}

Dense inverse(Dense A) {
    const std::size_t n = A.size();
    Dense I(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) I[i][i] = 1.0;

    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i)
            if (std::abs(A[i][k]) > std::abs(A[pivot][k])) pivot = i;
        EXPECT_TRUE(std::abs(A[pivot][k]) > 1e-12);
        std::swap(A[k], A[pivot]);
        std::swap(I[k], I[pivot]);

        const double d = A[k][k];
        for (std::size_t j = 0; j < n; ++j) {
            A[k][j] /= d;
            I[k][j] /= d;
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (i == k) continue;
            const double factor = A[i][k];
            for (std::size_t j = 0; j < n; ++j) {
                A[i][j] -= factor * A[k][j];
                I[i][j] -= factor * I[k][j];
            }
        }
    }
    return I;
}

Dense multiply(const Dense& A, const Dense& B) {
    Dense C(A.size(), std::vector<double>(B[0].size(), 0.0));
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t k = 0; k < B.size(); ++k)
            for (std::size_t j = 0; j < B[0].size(); ++j)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

Dense subtract(const Dense& A, const Dense& B) {
    Dense C = A;
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A[i].size(); ++j) C[i][j] -= B[i][j];
    return C;
}

Dense diagonal_inverse(const Dense& A) {
    Dense result(A.size(), std::vector<double>(A.size(), 0.0));
    for (std::size_t i = 0; i < A.size(); ++i) {
        EXPECT_TRUE(std::abs(A[i][i]) > 1e-12);
        result[i][i] = 1.0 / A[i][i];
    }
    return result;
}

double inf_norm(const Dense& A) {
    double result = 0.0;
    for (const auto& row : A) {
        double sum = 0.0;
        for (double value : row) sum += std::abs(value);
        result = std::max(result, sum);
    }
    return result;
}

Vector matvec(const Dense& A, const Vector& x) {
    Vector y(A.size(), 0.0);
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A[i].size(); ++j) y(i) += A[i][j] * x(j);
    return y;
}

double relative_error(const Vector& a, const Vector& b) {
    return (a - b).norm2() / std::max(1.0, b.norm2());
}

Vector make_rhs(std::size_t n, double phase) {
    Vector r(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        r(i) = std::sin(0.17 * static_cast<double>(i + 1) + phase)
             + 0.3 * std::cos(0.11 * static_cast<double>(i + 1) - phase);
    return r;
}

} // namespace

int main() {
    run_case("n8_schur_representative_cfd_matrix_qualification", [] {
        // Controlled finite-volume-like pressure/velocity matrices: 2-D
        // diffusion + convection velocity blocks, discrete pressure
        // gradients/divergence, and pressure stabilization. These are
        // representative algebraic CFD operators, not physical case oracles.
        const std::vector<RepresentativeCase> cases{
            {4, 0.0}, {4, 1.0}, {4, 4.0},
            {6, 0.0}, {6, 1.0}, {6, 4.0},
        };

        for (const auto& cfg : cases) {
            const auto data = make_representative_case(cfg);
            const BlockOperator blocks(data.Auu, data.G, data.D, data.C);
            blocks.validate();

            const Dense Auu = dense(blocks.Auu());
            const Dense G = dense(blocks.G());
            const Dense D = dense(blocks.D());
            const Dense C = dense(blocks.C());
            const Dense Auu_inv = inverse(Auu);
            const Dense Q_inv = diagonal_inverse(Auu);
            const Dense P = multiply(multiply(D, Q_inv), G);
            const Dense S = subtract(C, multiply(multiply(D, Auu_inv), G));
            const Dense P_inv = inverse(P);
            const Dense S_inv = inverse(S);

            const double cond_inf = inf_norm(Auu) * inf_norm(Auu_inv);
            EXPECT_TRUE(std::isfinite(cond_inf) && cond_inf > 1.0);

            auto pressure_solve = [P_inv](const Vector& rhs, Vector& solution) {
                solution = matvec(P_inv, rhs);
                return true;
            };
            ExactSchurApproximation exact(
                [Auu_inv](const Vector& rhs, Vector& solution) {
                    solution = matvec(Auu_inv, rhs);
                    return true;
                });
            EXPECT_TRUE(exact.setup(blocks));

            SimplerSchurApproximation simple(SimplerSchurMode::SIMPLE);
            SimplerSchurApproximation simplec(SimplerSchurMode::SIMPLEC);
            EXPECT_TRUE(simple.setup(blocks));
            EXPECT_TRUE(simplec.setup(blocks));

            double max_simple_error = 0.0;
            double max_simplec_error = 0.0;
            double max_lsc_error = 0.0;
            double max_bfbt_error = 0.0;

            for (const double phase : {0.0, 0.7}) {
                const Vector pressure = make_rhs(cfg.n * cfg.n, phase);
                const Vector rhs = make_rhs(cfg.n * cfg.n, phase + 0.31);

                Vector exact_operator(rhs.size(), 0.0);
                EXPECT_TRUE(exact.apply_schur(pressure, exact_operator));

                Vector simple_operator(rhs.size(), 0.0);
                Vector simplec_operator(rhs.size(), 0.0);
                EXPECT_TRUE(simple.apply(pressure, simple_operator));
                EXPECT_TRUE(simplec.apply(pressure, simplec_operator));

                max_simple_error = std::max(
                    max_simple_error, relative_error(simple_operator, exact_operator));
                max_simplec_error = std::max(
                    max_simplec_error, relative_error(simplec_operator, exact_operator));

                for (const auto mode : {LscBfbtSchurApproximation::Mode::LSC,
                                        LscBfbtSchurApproximation::Mode::BFBT}) {
                    const std::vector<double> qdiag =
                        mode == LscBfbtSchurApproximation::Mode::LSC
                            ? std::vector<double>{}
                            : [&] {
                                  std::vector<double> q(Auu.size(), 0.0);
                                  for (std::size_t i = 0; i < Auu.size(); ++i)
                                      q[i] = Auu[i][i];
                                  return q;
                              }();
                    LscBfbtSchurApproximation approx(mode, pressure_solve, qdiag);
                    EXPECT_TRUE(approx.setup(blocks));
                    Vector estimated(rhs.size(), 0.0);
                    EXPECT_TRUE(approx.apply(rhs, estimated));
                    const Vector exact_inverse = matvec(S_inv, rhs);
                    const double error = relative_error(estimated, exact_inverse);
                    if (mode == LscBfbtSchurApproximation::Mode::LSC)
                        max_lsc_error = std::max(max_lsc_error, error);
                    else
                        max_bfbt_error = std::max(max_bfbt_error, error);
                }
            }

            std::cout
                << "n8_schur_cfd_matrix"
                << " cells=" << cfg.n * cfg.n
                << " velocity_unknowns=" << blocks.velocity_size()
                << " pressure_unknowns=" << blocks.pressure_size()
                << " convection=" << cfg.convection
                << " Auu_nnz=" << blocks.Auu().nnz()
                << " cond_inf_Auu=" << cond_inf
                << " simple_error=" << max_simple_error
                << " simplec_error=" << max_simplec_error
                << " lsc_error=" << max_lsc_error
                << " bfbt_error=" << max_bfbt_error
                << '\n';

            // Measurement only: no approximation-quality threshold is invented
            // from this synthetic family. The next qualification step must
            // combine this dataset with physical CFD matrices and conditioning.
            EXPECT_TRUE(std::isfinite(max_simple_error));
            EXPECT_TRUE(std::isfinite(max_simplec_error));
            EXPECT_TRUE(std::isfinite(max_lsc_error));
            EXPECT_TRUE(std::isfinite(max_bfbt_error));
        }
    });

    return run_all();
}
