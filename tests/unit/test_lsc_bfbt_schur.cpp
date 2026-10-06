#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/exact_schur.h"
#include "cfdx/core/linalg/lsc_bfbt_schur.h"
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
        std::swap(A[k], A[pivot]);
        std::swap(I[k], I[pivot]);
        for (std::size_t j = 0; j < n; ++j) {
            A[k][j] /= A[k][k];
            I[k][j] /= A[k][k];
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
        for (std::size_t j = 0; j < A[i].size(); ++j) y(i) += A[i][j] * x(j);
    return y;
}

static double norm2(const Vector& x) {
    return std::sqrt(x.dot(x));
}

int main() {
    // Small Stokes-like system with C=0 and D=G^T.  The pressure Schur is
    // negative definite; LSC/BFBt therefore return the negative inverse
    // approximation required by the CFDX Schur convention.
    const auto Auu = make_sparse(2, 2, {
        {0, 0, 4.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 1, 3.0}
    });
    const auto G = make_sparse(2, 2, {
        {0, 0, 1.0}, {0, 1, 0.5}, {1, 0, 0.25}, {1, 1, 1.0}
    });
    const auto D = make_sparse(2, 2, {
        {0, 0, 1.0}, {0, 1, 0.25}, {1, 0, 0.5}, {1, 1, 1.0}
    });
    const auto C = make_sparse(2, 2, {});
    const BlockOperator blocks(Auu, G, D, C);
    blocks.validate();

    const auto P = dense(D);
    (void)P;
    const auto Ad = dense(Auu);
    const auto Gd = dense(G);
    const auto Dd = dense(D);
    const auto Ainv = inverse(Ad);

    const auto pressure_operator = [&Dd, &Gd](const Vector& x) {
        Vector y(2, 0.0);
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j) y(i) += Dd[i][j] * Gd[j][0] * x(0) + Dd[i][j] * Gd[j][1] * x(1);
        return y;
    };
    (void)pressure_operator;

    // Build P = D Q^-1 G and provide an exact small-system solve callback.
    const auto make_pressure_inverse = [&Dd, &Gd](const std::vector<double>& qinv) {
        std::vector<std::vector<double>> M(2, std::vector<double>(2, 0.0));
        for (std::size_t j = 0; j < 2; ++j)
            for (std::size_t i = 0; i < 2; ++i)
                for (std::size_t k = 0; k < 2; ++k)
                    M[i][j] += Dd[i][k] * qinv[k] * Gd[k][j];
        return inverse(M);
    };

    const auto run_mode = [&](LscBfbtSchurApproximation::Mode mode,
                              const std::vector<double>& qdiag) {
        std::vector<double> qinv(2, 1.0);
        if (!qdiag.empty()) {
            qinv[0] = 1.0 / qdiag[0];
            qinv[1] = 1.0 / qdiag[1];
        } else if (mode == LscBfbtSchurApproximation::Mode::LSC) {
            qinv[0] = 0.25;
            qinv[1] = 1.0 / 3.0;
        }
        const auto Pinv = make_pressure_inverse(qinv);
        auto solve_P = [Pinv](const Vector& rhs, Vector& z) {
            if (z.size() != rhs.size()) z = Vector(rhs.size(), 0.0);
            for (std::size_t i = 0; i < rhs.size(); ++i) {
                z(i) = 0.0;
                for (std::size_t j = 0; j < rhs.size(); ++j) z(i) += Pinv[i][j] * rhs(j);
            }
            return true;
        };
        LscBfbtSchurApproximation pre(mode, solve_P, qdiag);
        EXPECT_TRUE(pre.setup(blocks));
        Vector rhs(2, 0.0);
        rhs(0) = 1.0;
        rhs(1) = -0.4;
        Vector z(2, 0.0);
        EXPECT_TRUE(pre.apply(rhs, z));
        EXPECT_TRUE(std::isfinite(z(0)) && std::isfinite(z(1)));
        return z;
    };

    run_case("pressure_operator_matches_lsc_bfbt_definition", [&] {
        const auto lsc_p = LscBfbtSchurApproximation::assemble_pressure_operator(
            blocks, LscBfbtSchurApproximation::Mode::LSC);
        const auto bfbt_p = LscBfbtSchurApproximation::assemble_pressure_operator(
            blocks, LscBfbtSchurApproximation::Mode::BFBT);

        const auto expected = [](const SparseMatrix& Dm,
                                 const SparseMatrix& Gm,
                                 const std::vector<double>& q_inverse) {
            SparseMatrix Pm(Dm.n_rows(), Gm.n_cols());
            for (std::size_t row = 0; row < Dm.n_rows(); ++row) {
                for (std::size_t dk = Dm.row_offsets_data()[row];
                     dk < Dm.row_offsets_data()[row + 1]; ++dk) {
                    const std::size_t v = Dm.columns_data()[dk];
                    const double d = Dm.values_data()[dk] * q_inverse[v];
                    for (std::size_t gk = Gm.row_offsets_data()[v];
                         gk < Gm.row_offsets_data()[v + 1]; ++gk) {
                        Pm.push_back(row, Gm.columns_data()[gk],
                                     d * Gm.values_data()[gk]);
                    }
                }
            }
            Pm.finalize();
            return Pm;
        };

        const auto expected_lsc = expected(D, G, {0.25, 1.0 / 3.0});
        const auto expected_bfbt = expected(D, G, {1.0, 1.0});
        const auto max_abs_diff = [](const SparseMatrix& a, const SparseMatrix& b) {
            double result = 0.0;
            EXPECT_TRUE(a.n_rows() == b.n_rows() && a.n_cols() == b.n_cols());
            for (std::size_t row = 0; row < a.n_rows(); ++row) {
                for (std::size_t col = 0; col < a.n_cols(); ++col) {
                    double av = 0.0;
                    double bv = 0.0;
                    for (std::size_t k = a.row_offsets_data()[row];
                         k < a.row_offsets_data()[row + 1]; ++k)
                        if (a.columns_data()[k] == col) av += a.values_data()[k];
                    for (std::size_t k = b.row_offsets_data()[row];
                         k < b.row_offsets_data()[row + 1]; ++k)
                        if (b.columns_data()[k] == col) bv += b.values_data()[k];
                    result = std::max(result, std::abs(av - bv));
                }
            }
            return result;
        };

        EXPECT_NEAR(max_abs_diff(lsc_p, expected_lsc), 0.0, 1e-15);
        EXPECT_NEAR(max_abs_diff(bfbt_p, expected_bfbt), 0.0, 1e-15);
        EXPECT_TRUE(lsc_p.n_rows() == blocks.pressure_size());
        EXPECT_TRUE(bfbt_p.n_rows() == blocks.pressure_size());
    });

    run_case("lsc_and_bfbt_setup_and_apply", [&] {
        (void)run_mode(LscBfbtSchurApproximation::Mode::LSC, {});
        (void)run_mode(LscBfbtSchurApproximation::Mode::BFBT, {});
    });

    run_case("lsc_matches_documented_algebra", [&] {
        const std::vector<double> qinv{0.25, 1.0 / 3.0};
        const auto Pinv = make_pressure_inverse(qinv);
        auto solve_P = [Pinv](const Vector& rhs, Vector& z) {
            if (z.size() != rhs.size()) z = Vector(rhs.size(), 0.0);
            for (std::size_t i = 0; i < rhs.size(); ++i)
                for (std::size_t j = 0; j < rhs.size(); ++j) z(i) += Pinv[i][j] * rhs(j);
            return true;
        };
        LscBfbtSchurApproximation lsc(
            LscBfbtSchurApproximation::Mode::LSC, solve_P, {});
        EXPECT_TRUE(lsc.setup(blocks));

        Vector rhs(2, 0.0);
        rhs(0) = 0.7;
        rhs(1) = -1.2;
        Vector z(2, 0.0);
        EXPECT_TRUE(lsc.apply(rhs, z));

        // For this small example Q=diag(Auu), and the implementation must
        // equal -P^-1 E P^-1 exactly to floating-point roundoff.
        const auto Esolve = [&Dd, &Ad, &Gd](const Vector& r) {
            const auto Qinv = std::vector<double>{0.25, 1.0 / 3.0};
            Vector g(2, 0.0), qg(2, 0.0), aqg(2, 0.0), qaqg(2, 0.0), e(2, 0.0);
            for (std::size_t i = 0; i < 2; ++i)
                for (std::size_t j = 0; j < 2; ++j) g(i) += Gd[i][j] * r(j);
            for (std::size_t i = 0; i < 2; ++i) qg(i) = Qinv[i] * g(i);
            for (std::size_t i = 0; i < 2; ++i)
                for (std::size_t j = 0; j < 2; ++j) aqg(i) += Ad[i][j] * qg(j);
            for (std::size_t i = 0; i < 2; ++i) qaqg(i) = Qinv[i] * aqg(i);
            for (std::size_t i = 0; i < 2; ++i)
                for (std::size_t j = 0; j < 2; ++j) e(i) += Dd[i][j] * qaqg(j);
            return e;
        };
        Vector y(2, 0.0);
        EXPECT_TRUE(solve_P(rhs, y));
        Vector e = Esolve(y);
        Vector ref(2, 0.0);
        EXPECT_TRUE(solve_P(e, ref));
        ref *= -1.0;
        EXPECT_NEAR(norm2(z - ref), 0.0, 1e-12);
    });

    run_case("numeric_refresh_rejects_graph_change", [&] {
        const std::vector<double> qdiag{2.0, 5.0};
        auto solve_P = [](const Vector& rhs, Vector& z) {
            if (z.size() != rhs.size()) z = Vector(rhs.size(), 0.0);
            for (std::size_t i = 0; i < rhs.size(); ++i) z(i) = rhs(i);
            return true;
        };
        LscBfbtSchurApproximation pre(
            LscBfbtSchurApproximation::Mode::BFBT, solve_P, qdiag);
        EXPECT_TRUE(pre.setup(blocks));

        const auto G_changed = make_sparse(2, 2, {
            {0, 0, 1.0}, {0, 1, 0.75}, {1, 0, 0.25}
        });
        const BlockOperator changed(Auu, G_changed, D, C);
        EXPECT_TRUE(!pre.update_values(changed));
    });

    run_case("numeric_refresh_updates_values_without_rebuild", [&] {
        const std::vector<double> qdiag{2.0, 5.0};
        auto solve_P = [](const Vector& rhs, Vector& z) {
            if (z.size() != rhs.size()) z = Vector(rhs.size(), 0.0);
            for (std::size_t i = 0; i < rhs.size(); ++i) z(i) = rhs(i);
            return true;
        };

        LscBfbtSchurApproximation pre(
            LscBfbtSchurApproximation::Mode::BFBT, solve_P, qdiag);
        EXPECT_TRUE(pre.setup(blocks));

        Vector rhs(2, 0.0);
        rhs(0) = 1.0;
        rhs(1) = -0.4;
        Vector before(2, 0.0);
        EXPECT_TRUE(pre.apply(rhs, before));

        const auto D_changed = make_sparse(2, 2, {
            {0, 0, 1.2}, {0, 1, 0.25}, {1, 0, 0.5}, {1, 1, 1.0}
        });
        const BlockOperator changed(Auu, G, D_changed, C);
        EXPECT_TRUE(pre.update_values(changed));

        Vector after(2, 0.0);
        EXPECT_TRUE(pre.apply(rhs, after));
        EXPECT_TRUE(norm2(after - before) > 1e-12);
    });

    run_case("exact_schur_numeric_refresh_contract", [&] {
        ExactSchurApproximation exact(
            [](const Vector& rhs, Vector& y) {
                if (y.size() != rhs.size()) y = Vector(rhs.size(), 0.0);
                for (std::size_t i = 0; i < rhs.size(); ++i) y(i) = rhs(i) / 2.0;
                return true;
            });

        EXPECT_TRUE(exact.setup(blocks));
        Vector p(2, 0.0);
        p(0) = 1.0;
        p(1) = -0.5;
        Vector before(2, 0.0);
        EXPECT_TRUE(exact.apply_schur(p, before));

        const auto D_changed = make_sparse(2, 2, {
            {0, 0, 1.2}, {0, 1, 0.25}, {1, 0, 0.5}, {1, 1, 1.0}
        });
        const BlockOperator changed(Auu, G, D_changed, C);
        EXPECT_TRUE(exact.update_values(changed));

        Vector after(2, 0.0);
        EXPECT_TRUE(exact.apply_schur(p, after));
        EXPECT_TRUE(norm2(after - before) > 1e-12);

        const auto G_changed = make_sparse(2, 2, {
            {0, 0, 1.0}, {0, 1, 0.5}, {1, 0, 0.25}
        });
        const BlockOperator graph_changed(Auu, G_changed, D_changed, C);
        EXPECT_TRUE(!exact.update_values(graph_changed));

        Vector preserved(2, 0.0);
        EXPECT_TRUE(exact.apply_schur(p, preserved));
        EXPECT_NEAR(norm2(preserved - after), 0.0, 1e-14);
    });

    run_case("exact_schur_explicit_pressure_null_space", [&] {
        // Auu = I, G = [1 -1], D = -G^T gives S = G^T G.
        // The constant pressure mode is therefore an exact null mode and the
        // compatible pressure operator is SPD on the projected subspace.
        const auto Auu_ns = make_sparse(2, 2, {
            {0, 0, 1.0}, {1, 1, 1.0}
        });
        const auto G_ns = make_sparse(2, 2, {
            {0, 0, 1.0}, {0, 1, -1.0}
        });
        const auto D_ns = make_sparse(2, 2, {
            {0, 0, -1.0}, {1, 0, 1.0}
        });
        const BlockOperator ns_blocks(Auu_ns, G_ns, D_ns, C);
        NullSpaceProjector pressure_null_space = NullSpaceProjector::constant(2);

        ExactSchurApproximation exact(
            [](const Vector& rhs, Vector& y) {
                if (y.size() != rhs.size()) y = Vector(rhs.size(), 0.0);
                for (std::size_t i = 0; i < rhs.size(); ++i) y(i) = rhs(i);
                return true;
            },
            ExactSchurApproximation::Controls{100, 1e-12},
            pressure_null_space);

        EXPECT_TRUE(exact.setup(ns_blocks));

        Vector compatible_rhs(2, 0.0);
        compatible_rhs(0) = 1.0;
        compatible_rhs(1) = -1.0;
        Vector pressure(2, 0.0);
        EXPECT_TRUE(exact.apply(compatible_rhs, pressure));
        EXPECT_NEAR(pressure(0) + pressure(1), 0.0, 1e-12);

        Vector incompatible_rhs(2, 1.0);
        Vector rejected_pressure(2, 0.0);
        EXPECT_TRUE(!exact.apply(incompatible_rhs, rejected_pressure));

        Vector null_mode(2, 1.0);
        Vector projected(2, 0.0);
        EXPECT_TRUE(exact.apply_schur(null_mode, projected));
        EXPECT_NEAR(norm2(projected), 0.0, 1e-14);
    });

    run_case("bfbt_uses_explicit_velocity_scaling", [&] {
        const std::vector<double> qdiag{2.0, 5.0};
        const auto z = run_mode(LscBfbtSchurApproximation::Mode::BFBT, qdiag);
        EXPECT_TRUE(std::isfinite(z(0)) && std::isfinite(z(1)));
    });

    run_case("lsc_bfbt_declares_the_inverse_action", [&] {
        // apply() returns -P^{-1} E P^{-1} r, which approximates S^{-1} r, so both
        // modes are inverse actions. An adapter building a preconditioner needs
        // this action; one assembling a residual would be silently wrong.
        const auto Pinv = make_pressure_inverse({1.0, 1.0});
        auto solve_P = [Pinv](const Vector& rhs, Vector& z) {
            if (z.size() != rhs.size()) z = Vector(rhs.size(), 0.0);
            for (std::size_t i = 0; i < rhs.size(); ++i) {
                z(i) = 0.0;
                for (std::size_t j = 0; j < rhs.size(); ++j) z(i) += Pinv[i][j] * rhs(j);
            }
            return true;
        };

        LscBfbtSchurApproximation lsc(LscBfbtSchurApproximation::Mode::LSC, solve_P, {});
        LscBfbtSchurApproximation bfbt(LscBfbtSchurApproximation::Mode::BFBT, solve_P, {});
        EXPECT_TRUE(lsc.action() == SchurAction::InverseOperator);
        EXPECT_TRUE(bfbt.action() == SchurAction::InverseOperator);
        EXPECT_TRUE(provides_action(lsc, SchurAction::InverseOperator));
        EXPECT_TRUE(!provides_action(lsc, SchurAction::Operator));

        // The declared action must hold through setup, not only on a bare instance.
        EXPECT_TRUE(lsc.setup(blocks));
        EXPECT_TRUE(provides_action(lsc, SchurAction::InverseOperator));
    });


    run_case("lsc_bfbt_pressure_nullspace_projection_is_explicit", [&] {
        // A singular pressure-side operator is the natural situation before
        // pinning a pressure reference. The Schur layer must either receive an
        // explicit null-space policy or reject an incompatible RHS; it must
        // never silently pretend that the singular inverse is ordinary.
        const auto auu = make_sparse(2, 2, {
            {0, 0, 2.0}, {1, 1, 3.0}
        });
        const auto g = make_sparse(2, 2, {
            {0, 0, 1.0}, {1, 1, 1.0}
        });
        const auto d = make_sparse(2, 2, {
            {0, 0, 1.0}, {1, 1, 1.0}
        });
        const auto c = make_sparse(2, 2, {});
        const BlockOperator singular_blocks(auu, g, d, c);

        const auto solve_singular = [](const Vector& rhs, Vector& z) {
            // This callback deliberately refuses an incompatible RHS. The
            // Schur implementation must propagate that failure.
            double mean = 0.0;
            for (std::size_t i = 0; i < rhs.size(); ++i) mean += rhs(i);
            if (std::abs(mean) > 1e-12) return false;
            z = rhs;
            return true;
        };

        LscBfbtSchurApproximation lsc(
            LscBfbtSchurApproximation::Mode::LSC, solve_singular);
        EXPECT_TRUE(lsc.setup(singular_blocks));

        Vector incompatible(2, 0.0);
        incompatible(0) = 1.0;
        Vector out(2, 0.0);
        EXPECT_TRUE(!lsc.apply(incompatible, out));
    });

    run_case("lsc_bfbt_apply_propagates_pressure_solve_failure", [&] {
        const auto auu = make_sparse(2, 2, {
            {0, 0, 2.0}, {1, 1, 3.0}
        });
        const auto g = make_sparse(2, 2, {
            {0, 0, 1.0}, {1, 1, 1.0}
        });
        const auto d = make_sparse(2, 2, {
            {0, 0, 1.0}, {1, 1, 1.0}
        });
        const auto c = make_sparse(2, 2, {});
        const BlockOperator blocks(auu, g, d, c);

        std::size_t calls = 0;
        const auto failing_solver = [&calls](const Vector&, Vector&) {
            ++calls;
            return false;
        };
        LscBfbtSchurApproximation lsc(
            LscBfbtSchurApproximation::Mode::LSC, failing_solver);
        EXPECT_TRUE(lsc.setup(blocks));

        Vector rhs(2, 1.0);
        Vector out(2, 0.0);
        EXPECT_TRUE(!lsc.apply(rhs, out));
        EXPECT_TRUE(calls == 1);
    });

    return run_all();
}
