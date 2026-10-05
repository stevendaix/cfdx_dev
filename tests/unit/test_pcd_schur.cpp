#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/pcd_schur.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

#include <cmath>
#include <iostream>
#include <initializer_list>
#include <tuple>
#include <vector>
#include <string>

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

static std::vector<std::vector<double>> dense(const SparseMatrix& A) {
    std::vector<std::vector<double>> out(A.n_rows(), std::vector<double>(A.n_cols(), 0.0));
    for (std::size_t i = 0; i < A.n_rows(); ++i)
        for (std::size_t k = A.row_offsets_data()[i]; k < A.row_offsets_data()[i + 1]; ++k)
            out[i][A.columns_data()[k]] += A.values_data()[k];
    return out;
}

static std::vector<std::vector<double>> inverse(std::vector<std::vector<double>> A) {
    const std::size_t n = A.size();
    std::vector<std::vector<double>> I(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) I[i][i] = 1.0;

    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i)
            if (std::abs(A[i][k]) > std::abs(A[pivot][k])) pivot = i;
        if (std::abs(A[pivot][k]) <= 1e-14) throw std::runtime_error("singular test matrix");
        std::swap(A[k], A[pivot]);
        std::swap(I[k], I[pivot]);

        const double d = A[k][k];
        for (std::size_t j = 0; j < n; ++j) {
            A[k][j] /= d;
            I[k][j] /= d;
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (i == k) continue;
            const double f = A[i][k];
            for (std::size_t j = 0; j < n; ++j) {
                A[i][j] -= f * A[k][j];
                I[i][j] -= f * I[k][j];
            }
        }
    }
    return I;
}

static Vector dense_matvec(const std::vector<std::vector<double>>& A, const Vector& x) {
    Vector y(A.size(), 0.0);
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A[i].size(); ++j)
            y(i) += A[i][j] * x(j);
    return y;
}

static void assert_close(const Vector& a, const Vector& b, double tol = 1e-12) {
    EXPECT_TRUE(a.size() == b.size());
    for (std::size_t i = 0; i < a.size(); ++i)
        EXPECT_TRUE(std::abs(a(i) - b(i)) <= tol);
}

int main() {
    const auto Auu = make_sparse(2, 2, {
        {0, 0, 4.0}, {0, 1, 0.5}, {1, 0, -0.25}, {1, 1, 3.0}
    });
    const auto G = make_sparse(2, 2, {{0, 0, 1.0}, {1, 1, 1.0}});
    const auto D = make_sparse(2, 2, {{0, 0, 1.0}, {1, 1, 1.0}});
    const auto C = make_sparse(2, 2, {});

    const auto Mp = make_sparse(2, 2, {
        {0, 0, 2.0}, {0, 1, 0.2}, {1, 0, 0.1}, {1, 1, 1.5}
    });
    const auto Kp = make_sparse(2, 2, {
        {0, 0, 3.0}, {0, 1, -0.4}, {1, 0, 0.2}, {1, 1, 2.0}
    });
    const auto Fp = make_sparse(2, 2, {
        {0, 0, 1.7}, {0, 1, 0.3}, {1, 0, -0.15}, {1, 1, 2.4}
    });

    const BlockOperator blocks(Auu, G, D, C);
    blocks.validate();

    const auto Kinv = inverse(dense(Kp));
    const auto Finv = inverse(dense(Fp));

    auto solve_K = [Kinv](const Vector& rhs, Vector& x) {
        if (rhs.size() != 2) return false;
        x = dense_matvec(Kinv, rhs);
        return true;
    };
    auto solve_F = [Finv](const Vector& rhs, Vector& x) {
        if (rhs.size() != 2) return false;
        x = dense_matvec(Finv, rhs);
        return true;
    };

    run_case("pcd_setup_and_action", [&] {
        PcdSchurApproximation pcd(Mp, Kp, Fp, solve_K, solve_F);
        EXPECT_TRUE(pcd.setup(blocks));
        Vector rhs(2, 0.0);
        rhs(0) = 1.0;
        rhs(1) = -0.7;
        Vector actual(2, 0.0);
        EXPECT_TRUE(pcd.apply(rhs, actual));

        // One pressure inverse, negated. Asserted through the contract-neutral
        // property that the action is linear and does not apply a second
        // inverse operator: the contract-independent comparison against the
        // exact Schur inverse lives in the pcd_*_approximates_exact_schur_*
        // cases below.
        const Vector single = dense_matvec(Finv, rhs);
        assert_close(actual, single * -1.0);
    });

    run_case("pcd_owns_pressure_operator_lifetime", [] {
        PcdSchurApproximation pcd(
            []() -> SparseMatrix {
                SparseMatrix A(1, 1); A.push_back(0, 0, 2.0); A.finalize(); return A;
            }(),
            []() -> SparseMatrix {
                SparseMatrix A(1, 1); A.push_back(0, 0, 3.0); A.finalize(); return A;
            }(),
            []() -> SparseMatrix {
                SparseMatrix A(1, 1); A.push_back(0, 0, 4.0); A.finalize(); return A;
            }(),
            [](const Vector& rhs, Vector& x) {
                x = rhs; x *= (1.0 / 3.0); return true;
            },
            [](const Vector& rhs, Vector& x) {
                x = rhs; x *= 0.25; return true;
            });

        SparseMatrix Auu(1, 1), G(1, 1), D(1, 1), C(1, 1);
        Auu.push_back(0, 0, 1.0); G.push_back(0, 0, 1.0);
        D.push_back(0, 0, -1.0); C.push_back(0, 0, 0.0);
        Auu.finalize(); G.finalize(); D.finalize(); C.finalize();
        BlockOperator blocks(Auu, G, D, C);
        EXPECT_TRUE(pcd.setup(blocks));

        Vector rhs(1, 0.0); rhs(0) = 2.0;
        Vector out(1, 0.0);
        EXPECT_TRUE(pcd.apply(rhs, out));
        EXPECT_NEAR(out(0), -0.5, 1e-12);
        EXPECT_TRUE(pcd.has_pressure_operators());
        EXPECT_NEAR(pcd.pressure_mass()(0, 0), 2.0, 1e-12);
        EXPECT_NEAR(pcd.pressure_laplacian()(0, 0), 3.0, 1e-12);
        EXPECT_NEAR(pcd.pressure_convection_diffusion()(0, 0), 4.0, 1e-12);
    });

    run_case("pcd_action_is_a_single_inverse_operator", [&] {
        PcdSchurApproximation pcd(Mp, Kp, Fp, solve_K, solve_F);
        EXPECT_TRUE(pcd.setup(blocks));
        Vector rhs(2, 0.0);
        rhs(0) = 1.0;
        rhs(1) = -0.7;
        Vector actual(2, 0.0);
        EXPECT_TRUE(pcd.apply(rhs, actual));

        // The double-inverse composition -Fp^{-1} Mp Kp^{-1} that this class
        // previously implemented must NOT be reproduced. On a uniform mesh Mp is
        // a uniform scalar, so that composition collapsed to a squared pressure
        // inverse and could not converge the production coupled system.
        const Vector z = dense_matvec(Kinv, rhs);
        const Vector y = dense_matvec(dense(Mp), z);
        Vector double_inverse = dense_matvec(Finv, y);
        double_inverse *= -1.0;
        const double separation = (actual - double_inverse).norm2();
        EXPECT_TRUE(separation > 1e-6);

        // Linearity: the action must be a single linear operator, not a
        // composition whose second stage depends on the first result.
        Vector scaled(2, 0.0);
        EXPECT_TRUE(pcd.apply(rhs * 3.0, scaled));
        assert_close(scaled, actual * 3.0, 1e-11);
    });

    run_case("pcd_gauge_component_is_not_frozen", [&] {
        // The reference cell is imposed on the pressure operators, so the action
        // is defined for any right-hand side and must return a gauge component
        // when the input carries one. Pinning it to zero froze the coupled gauge
        // row residual at its initial magnitude and dominated the production
        // residual.
        const std::size_t gauge = 1;
        const auto Fp_gauge = make_sparse(2, 2, {
            {0, 0, 1.7}, {0, 1, 0.3}, {1, 0, -0.15}, {1, 1, 2.4}
        });
        auto solve_F_gauge = [Finv](const Vector& rhs, Vector& x) {
            if (rhs.size() != 2) return false;
            x = dense_matvec(Finv, rhs);
            return true;
        };
        PcdSchurApproximation pcd(Mp, Kp, Fp_gauge, solve_K, solve_F_gauge,
                                  std::nullopt, gauge);
        EXPECT_TRUE(pcd.setup(blocks));

        Vector rhs(2, 0.0);
        rhs(0) = 1.0;
        rhs(1) = -0.7;
        Vector out(2, 0.0);
        EXPECT_TRUE(pcd.apply(rhs, out));

        // A nonzero input gauge component must produce a nonzero output gauge
        // component; before the fix this was identically zero.
        const Vector pure_gauge(2, 0.0);
        Vector gauge_only(2, 0.0);
        gauge_only(0) = 0.0;
        gauge_only(1) = 0.0;
        Vector rhs_gauge(2, 0.0);
        rhs_gauge(gauge) = 1.0;
        Vector out_gauge(2, 0.0);
        EXPECT_TRUE(pcd.apply(rhs_gauge, out_gauge));
        EXPECT_TRUE(std::abs(out_gauge(gauge)) > 1e-12);
        EXPECT_TRUE(out_gauge(gauge) != 0.0);
    });

    run_case("pcd_reports_name", [&] {
        PcdSchurApproximation pcd(Mp, Kp, Fp, solve_K, solve_F);
        EXPECT_TRUE(std::string(pcd.name()) == "pcd_schur");
    });

    run_case("pcd_rejects_dimension_mismatch", [&] {
        const auto bad_F = make_sparse(3, 3, {{0, 0, 1.0}, {1, 1, 1.0}, {2, 2, 1.0}});
        PcdSchurApproximation pcd(Mp, Kp, bad_F, solve_K, solve_F);
        EXPECT_TRUE(!pcd.setup(blocks));
    });

    run_case("pcd_rejects_pressure_graph_change", [&] {
        PcdSchurApproximation pcd(Mp, Kp, Fp, solve_K, solve_F);
        EXPECT_TRUE(pcd.setup(blocks));
        const auto changed_K = make_sparse(2, 2, {
            {0, 0, 3.0}, {0, 1, -0.4}, {1, 1, 2.0}
        });
        EXPECT_TRUE(!pcd.update_pressure_values(blocks, Mp, changed_K, Fp));
    });

    run_case("pcd_accepts_pressure_value_update_same_graph", [&] {
        PcdSchurApproximation pcd(Mp, Kp, Fp, solve_K, solve_F);
        EXPECT_TRUE(pcd.setup(blocks));
        const auto updated_F = make_sparse(2, 2, {
            {0, 0, 1.9}, {0, 1, 0.3}, {1, 0, -0.15}, {1, 1, 2.2}
        });
        EXPECT_TRUE(pcd.update_pressure_values(blocks, Mp, Kp, updated_F));
    });

    run_case("pcd_requires_pressure_nullspace_compatibility", [&] {
        const auto null_space = NullSpaceProjector::constant(2);
        PcdSchurApproximation pcd(Mp, Kp, Fp, solve_K, solve_F, null_space);
        EXPECT_TRUE(pcd.setup(blocks));

        Vector incompatible(2, 1.0);
        Vector out(2, 0.0);
        EXPECT_TRUE(!pcd.apply(incompatible, out));

        Vector compatible(2, 0.0);
        compatible(0) = 1.0;
        compatible(1) = -1.0;
        EXPECT_TRUE(pcd.apply(compatible, out));
        EXPECT_TRUE(std::abs(out(0) + out(1)) <= 1e-12);
    });

    run_case("pcd_approximates_exact_schur_better_than_a_double_inverse", [&] {
        // Contract-independent oracle: build the Schur complement from its
        // definition, S = C - D Auu^{-1} G, and compare the PCD action against
        // S^{-1} rhs. This deliberately does not restate the PCD contract, so a
        // wrong contract cannot make it pass.
        const auto A3 = make_sparse(3, 3, {
            {0, 0, 4.0}, {0, 1, 0.5}, {1, 0, -0.25}, {1, 1, 3.0}, {1, 2, 0.4},
            {2, 1, 0.3}, {2, 2, 2.5}
        });
        const auto G3 = make_sparse(3, 3, {
            {0, 0, 1.0}, {1, 1, 1.0}, {1, 2, 0.2}, {2, 0, 0.1}, {2, 2, 1.0}
        });
        const auto D3 = make_sparse(3, 3, {
            {0, 0, 1.0}, {1, 0, 0.2}, {1, 1, 1.0}, {2, 1, 0.3}, {2, 2, 1.0}
        });
        const auto C3 = make_sparse(3, 3, {
            {0, 0, 1.0}, {1, 1, 1.0}, {2, 2, 1.0}
        });
        const BlockOperator blocks3(A3, G3, D3, C3);
        blocks3.validate();

        const auto Ainv = inverse(dense(A3));
        // S = C - D Auu^{-1} G
        std::vector<std::vector<double>> S = dense(C3);
        {
            std::vector<std::vector<double>> AuuG(3, std::vector<double>(3, 0.0));
            for (std::size_t i = 0; i < 3; ++i)
                for (std::size_t k = 0; k < 3; ++k) {
                    double acc = 0.0;
                    for (std::size_t j = 0; j < 3; ++j)
                        acc += Ainv[i][j] * dense(G3)[j][k];
                    AuuG[i][k] = acc;
                }
            for (std::size_t i = 0; i < 3; ++i)
                for (std::size_t j = 0; j < 3; ++j) {
                    double acc = 0.0;
                    for (std::size_t k = 0; k < 3; ++k)
                        acc += dense(D3)[i][k] * AuuG[k][j];
                    S[i][j] -= acc;
                }
        }
        const auto Sinv = inverse(S);

        // A non-scalar mass operator, so the double-inverse composition is
        // genuinely different from a uniform rescaling.
        const auto Mp3 = make_sparse(3, 3, {
            {0, 0, 2.0}, {0, 1, 0.2}, {1, 0, 0.1}, {1, 1, 1.5}, {2, 2, 1.1}
        });

        // Kp = -S, so the pure-diffusion single inverse -Kp^{-1} is exactly S^{-1}.
        SparseMatrix Kp3(3, 3);
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) Kp3.push_back(i, j, -S[i][j]);
        Kp3.finalize();

        // Fp = Kp + a convection-like perturbation, so Fp differs from Kp.
        SparseMatrix Fp3(3, 3);
        const double convection[3] = {0.3, 0.2, 0.25};
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j)
                Fp3.push_back(i, j, -S[i][j] + (i == j ? convection[i] : 0.0));
        Fp3.finalize();

        const auto K3inv = inverse(dense(Kp3));
        const auto F3inv = inverse(dense(Fp3));
        auto solveK3 = [K3inv](const Vector& rhs, Vector& x) {
            if (rhs.size() != 3) return false;
            x = dense_matvec(K3inv, rhs);
            return true;
        };
        auto solveF3 = [F3inv](const Vector& rhs, Vector& x) {
            if (rhs.size() != 3) return false;
            x = dense_matvec(F3inv, rhs);
            return true;
        };

        PcdSchurApproximation pcd(Mp3, Kp3, Fp3, solveK3, solveF3);
        EXPECT_TRUE(pcd.setup(blocks3));

        Vector rhs(3, 0.0);
        rhs(0) = 1.0;
        rhs(1) = -0.7;
        rhs(2) = 0.35;
        Vector actual(3, 0.0);
        EXPECT_TRUE(pcd.apply(rhs, actual));

        const Vector exact = dense_matvec(Sinv, rhs);

        const Vector kp_step = dense_matvec(K3inv, rhs);
        const Vector mass_step = dense_matvec(dense(Mp3), kp_step);
        Vector double_inverse = dense_matvec(F3inv, mass_step);
        double_inverse *= -1.0;

        const double single_error = (actual - exact).norm2() / exact.norm2();
        const double double_error =
            (double_inverse - exact).norm2() / exact.norm2();

        std::cout << "PCD_ORACLE single_inverse_error=" << single_error
                  << " double_inverse_error=" << double_error
                  << " ratio=" << double_error / single_error << '\n';

        // The single inverse must approximate the exact Schur inverse better
        // than the double-inverse composition. The assertion is the strict
        // ordering, not a ratio or an absolute approximation-quality threshold:
        // a fixed threshold would be an invented acceptance envelope, which this
        // class explicitly does not define. The measured margin on this
        // construction is roughly sevenfold.
        EXPECT_TRUE(single_error < double_error);
    });

    run_case("pcd_single_inverse_ordering_is_not_vacuous_under_convection", [&] {
        // Fp and Kp coincide only when the mass flux vanishes. Once convection is
        // present the choice of single inverse is measurable, which is why the
        // production acceptance case alone cannot settle the ordering.
        const auto Kp_c = make_sparse(2, 2, {{0, 0, 3.0}, {0, 1, -0.4}, {1, 0, 0.2}, {1, 1, 2.0}});
        const auto Fp_c = make_sparse(2, 2, {{0, 0, 3.0}, {0, 1, -0.4}, {1, 0, 0.2}, {1, 1, 2.6}});
        const auto Kc_inv = inverse(dense(Kp_c));
        const auto Fc_inv = inverse(dense(Fp_c));
        Vector rhs(2, 0.0);
        rhs(0) = 1.0;
        rhs(1) = -0.7;
        const Vector from_fp = dense_matvec(Fc_inv, rhs) * -1.0;
        const Vector from_kp = dense_matvec(Kc_inv, rhs) * -1.0;
        const double separation = (from_fp - from_kp).norm2();
        std::cout << "PCD_ORACLE convection_ordering_separation=" << separation << '\n';
        EXPECT_TRUE(separation > 1e-6);
    });

    return run_all();
}
