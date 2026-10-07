#include "cfdx/core/linalg/exact_schur.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/lsc_bfbt_schur.h"
#include "cfdx/core/linalg/pcd_schur.h"
#include "cfdx/core/linalg/simplerc_schur.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <string>
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

    for (std::size_t y = 0; y < cfg.n; ++y) {
        for (std::size_t x = 0; x < cfg.n; ++x) {
            const std::size_t p = cell(x, y);
            for (std::size_t component = 0; component < 2; ++component) {
                const std::size_t row = component * np + p;
                const double diagonal = reaction + 4.0 * diffusion + beta;
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

            const std::size_t ux = p;
            const std::size_t uy = np + p;
            if (x > 0) g_entries.emplace_back(ux, cell(x - 1, y), -0.5);
            else g_entries.emplace_back(ux, p, -1.0);
            if (x + 1 < cfg.n) g_entries.emplace_back(ux, cell(x + 1, y), 0.5);
            else g_entries.emplace_back(ux, p, 1.0);
            if (y > 0) g_entries.emplace_back(uy, cell(x, y - 1), -0.5);
            else g_entries.emplace_back(uy, p, -1.0);
            if (y + 1 < cfg.n) g_entries.emplace_back(uy, cell(x, y + 1), 0.5);
            else g_entries.emplace_back(uy, p, 1.0);
            if (p == 0) g_entries.emplace_back(ux, p, 0.25);

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
    for (std::size_t i = 0; i < A.n_rows(); ++i)
        for (std::size_t k = A.row_offsets_data()[i];
             k < A.row_offsets_data()[i + 1]; ++k)
            out[i][A.columns_data()[k]] += A.values_data()[k];
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

SparseMatrix make_pcd_pressure_operator(std::size_t n, double convection) {
    const std::size_t np = n * n;
    std::vector<std::tuple<std::size_t, std::size_t, double>> entries;
    const auto cell = [n](std::size_t x, std::size_t y) { return y * n + x; };
    for (std::size_t y = 0; y < n; ++y) {
        for (std::size_t x = 0; x < n; ++x) {
            const std::size_t p = cell(x, y);
            entries.emplace_back(p, p, 5.0 + convection);
            if (x > 0) entries.emplace_back(p, cell(x - 1, y), -1.0 - convection);
            if (x + 1 < n) entries.emplace_back(p, cell(x + 1, y), -1.0);
            if (y > 0) entries.emplace_back(p, cell(x, y - 1), -1.0);
            if (y + 1 < n) entries.emplace_back(p, cell(x, y + 1), -1.0);
        }
    }
    return make_sparse(np, np, entries);
}

Vector matvec(const Dense& A, const Vector& x) {
    Vector y(A.size(), 0.0);
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A[i].size(); ++j) y(i) += A[i][j] * x(j);
    return y;
}

double relative_residual(const Dense& A, const Vector& x, const Vector& b) {
    const Vector r = b - matvec(A, x);
    return r.norm2() / std::max(b.norm2(), 1e-300);
}

Vector make_rhs(std::size_t n, double phase) {
    Vector r(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        r(i) = std::sin(0.17 * static_cast<double>(i + 1) + phase)
             + 0.3 * std::cos(0.11 * static_cast<double>(i + 1) - phase);
    return r;
}

class CountingCallbackPreconditioner final : public Preconditioner {
public:
    using Apply = std::function<bool(const Vector&, Vector&)>;

    CountingCallbackPreconditioner(std::string name, Apply apply)
        : name_(std::move(name)), apply_(std::move(apply)) {}

    bool setup(const SparseMatrix& A) override {
        size_ = A.n_rows();
        return A.n_rows() == A.n_cols() && size_ > 0;
    }

    bool apply(const Vector& r, Vector& z) const override {
        ++applications_;
        return apply_(r, z);
    }

    const char* name() const override { return name_.c_str(); }
    std::size_t applications() const noexcept { return applications_; }

private:
    std::string name_;
    Apply apply_;
    std::size_t size_ = 0;
    mutable std::size_t applications_ = 0;
};

struct MethodResult {
    const char* name;
    SolverResult solver;
    double true_residual = 0.0;
    std::size_t preconditioner_applications = 0;
    long long setup_us = 0;
    long long solve_us = 0;
};

} // namespace

int main() {
    run_case("n8_schur_complete_krylov_qualification", [] {
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
            const Dense S =
                subtract(C, multiply(multiply(D, Auu_inv), G));
            const Dense S_inv = inverse(S);
            const Dense Q_inv = [&] {
                Dense q(Auu.size(), std::vector<double>(Auu.size(), 0.0));
                for (std::size_t i = 0; i < Auu.size(); ++i)
                    q[i][i] = 1.0 / Auu[i][i];
                return q;
            }();
            const Dense P = multiply(multiply(D, Q_inv), G);
            const Dense P_inv = inverse(P);

            const auto Mp = make_sparse(
                cfg.n * cfg.n, cfg.n * cfg.n,
                [&] {
                    std::vector<std::tuple<std::size_t, std::size_t, double>> e;
                    for (std::size_t i = 0; i < cfg.n * cfg.n; ++i)
                        e.emplace_back(i, i, 1.0);
                    return e;
                }());
            const auto Kp = make_pcd_pressure_operator(cfg.n, 0.0);
            const auto Fp = make_pcd_pressure_operator(cfg.n, cfg.convection);
            const Dense Kp_inv = inverse(dense(Kp));
            const Dense Fp_inv = inverse(dense(Fp));

            auto exact_action = [S](const Vector& rhs, Vector& out) {
                out = matvec(S, rhs);
            };
            LinearOperator schur_operator{cfg.n * cfg.n, exact_action};

            std::size_t method_records = 0;
            std::size_t converged_approximation_records = 0;
            double max_true_residual = 0.0;
            double max_reported_true_residual_gap = 0.0;

            auto run_method = [&](const char* name,
                                  long long setup_us,
                                  CountingCallbackPreconditioner& pc,
                                  const Vector& rhs) {
                EXPECT_TRUE(pc.setup(Kp));
                Vector x(rhs.size(), 0.0);
                const auto begin = std::chrono::steady_clock::now();
                const int restart = static_cast<int>(std::min<std::size_t>(cfg.n * cfg.n, 64));
                const SolverResult result = solve_fgmres(
                    schur_operator, rhs, x, restart, 200, 1e-10, &pc);
                const auto end = std::chrono::steady_clock::now();
                MethodResult out{
                    name,
                    result,
                    relative_residual(S, x, rhs),
                    pc.applications(),
                    setup_us,
                    std::chrono::duration_cast<std::chrono::microseconds>(
                        end - begin).count()};
                ++method_records;
                max_true_residual = std::max(max_true_residual, out.true_residual);
                max_reported_true_residual_gap = std::max(
                    max_reported_true_residual_gap,
                    std::abs(out.solver.residual_relative - out.true_residual));
                if (out.solver.status == SolverStatus::CONVERGED &&
                    std::string(name) != "exact")
                    ++converged_approximation_records;
                EXPECT_TRUE(std::isfinite(out.true_residual));
                EXPECT_TRUE(out.solver.status != SolverStatus::NOT_APPLICABLE);
                std::cout
                    << "n8_schur_complete"
                    << " cells=" << cfg.n * cfg.n
                    << " convection=" << cfg.convection
                    << " method=" << out.name
                    << " krylov=FGMRES"
                    << " requested_method=" << out.name
                    << " resolved_method=" << pc.name()
                    << " status=" << static_cast<int>(out.solver.status)
                    << " converged="
                    << (out.solver.status == SolverStatus::CONVERGED ? 1 : 0)
                    << " iterations=" << out.solver.iterations
                    << " reported_residual=" << out.solver.residual
                    << " reported_relative=" << out.solver.residual_relative
                    << " true_residual=" << out.true_residual
                    << " preconditioner_applications="
                    << out.preconditioner_applications
                    << " setup_us=" << out.setup_us
                    << " solve_us=" << out.solve_us
                    << '\n';
                return out;
            };

            const Vector rhs = make_rhs(cfg.n * cfg.n, 0.23);

            const auto exact_setup_begin = std::chrono::steady_clock::now();
            CountingCallbackPreconditioner exact_pc(
                "exact_schur_inverse",
                [S_inv](const Vector& r, Vector& z) {
                    z = matvec(S_inv, r);
                    return true;
                });
            const auto exact_setup_end = std::chrono::steady_clock::now();
            const auto exact_setup_us =
                std::chrono::duration_cast<std::chrono::microseconds>(
                    exact_setup_end - exact_setup_begin).count();
            const auto exact = run_method(
                "exact", exact_setup_us, exact_pc, rhs);
            EXPECT_TRUE(exact.solver.status == SolverStatus::CONVERGED);

            SimplerSchurApproximation simple(SimplerSchurMode::SIMPLE);
            auto begin = std::chrono::steady_clock::now();
            EXPECT_TRUE(simple.setup(blocks));
            const auto simple_op = simple.assembled_operator();
            const Dense simple_inv = inverse(dense(simple_op));
            auto end = std::chrono::steady_clock::now();
            const auto simple_setup_us =
                std::chrono::duration_cast<std::chrono::microseconds>(
                    end - begin).count();
            CountingCallbackPreconditioner simple_pc(
                "simple_schur_inverse",
                [simple_inv](const Vector& r, Vector& z) {
                    z = matvec(simple_inv, r);
                    return true;
                });
            const auto simple_result =
                run_method("simple", simple_setup_us, simple_pc, rhs);

            SimplerSchurApproximation simplec(SimplerSchurMode::SIMPLEC);
            begin = std::chrono::steady_clock::now();
            EXPECT_TRUE(simplec.setup(blocks));
            const auto simplec_op = simplec.assembled_operator();
            const Dense simplec_inv = inverse(dense(simplec_op));
            end = std::chrono::steady_clock::now();
            const auto simplec_setup_us =
                std::chrono::duration_cast<std::chrono::microseconds>(
                    end - begin).count();
            CountingCallbackPreconditioner simplec_pc(
                "simplec_schur_inverse",
                [simplec_inv](const Vector& r, Vector& z) {
                    z = matvec(simplec_inv, r);
                    return true;
                });
            const auto simplec_result =
                run_method("simplec", simplec_setup_us, simplec_pc, rhs);

            const auto pressure_solve =
                [P_inv](const Vector& r, Vector& z) {
                    z = matvec(P_inv, r);
                    return true;
                };
            const std::vector<double> qdiag = [&] {
                std::vector<double> q(Auu.size(), 0.0);
                for (std::size_t i = 0; i < Auu.size(); ++i) q[i] = Auu[i][i];
                return q;
            }();

            for (const auto mode :
                 {LscBfbtSchurApproximation::Mode::LSC,
                  LscBfbtSchurApproximation::Mode::BFBT}) {
                const char* name =
                    mode == LscBfbtSchurApproximation::Mode::LSC ? "lsc" : "bfbt";
                begin = std::chrono::steady_clock::now();
                LscBfbtSchurApproximation approx(mode, pressure_solve, qdiag);
                EXPECT_TRUE(approx.setup(blocks));
                end = std::chrono::steady_clock::now();
                const auto setup_us =
                    std::chrono::duration_cast<std::chrono::microseconds>(
                        end - begin).count();
                CountingCallbackPreconditioner pc(
                    name,
                    [&approx](const Vector& r, Vector& z) {
                        return approx.apply(r, z);
                    });
                run_method(name, setup_us, pc, rhs);
            }

            begin = std::chrono::steady_clock::now();
            auto solve_kp = [Kp_inv](const Vector& r, Vector& z) {
                z = matvec(Kp_inv, r);
                return true;
            };
            auto solve_fp = [Fp_inv](const Vector& r, Vector& z) {
                z = matvec(Fp_inv, r);
                return true;
            };
            PcdSchurApproximation pcd(Mp, Kp, Fp, solve_kp, solve_fp);
            EXPECT_TRUE(pcd.setup(blocks));
            end = std::chrono::steady_clock::now();
            const auto pcd_setup_us =
                std::chrono::duration_cast<std::chrono::microseconds>(
                    end - begin).count();
            CountingCallbackPreconditioner pcd_pc(
                "pcd_schur",
                [&pcd](const Vector& r, Vector& z) {
                    return pcd.apply(r, z);
                });
            const auto pcd_result = run_method("pcd", pcd_setup_us, pcd_pc, rhs);

            // A non-zero convection point must actually distinguish the two
            // pressure operators. This prevents the PCD campaign from passing
            // only because Fp and Kp collapse to the same operator.
            if (cfg.convection > 0.0) {
                const Vector fp_action = matvec(Fp_inv, rhs);
                const Vector kp_action = matvec(Kp_inv, rhs);
                EXPECT_TRUE((fp_action - kp_action).norm2() > 1e-8);
                EXPECT_TRUE(pcd_result.solver.status != SolverStatus::NOT_APPLICABLE);
            }

            // Non-convergence is deliberately reported in the machine-readable
            // record rather than hidden behind a tolerance relaxation or a
            // method substitution. The exact Schur oracle is the only mandatory
            // convergence gate at this stage; the approximation data define the
            // empirical envelope used by the subsequent N8 acceptance decision.
            EXPECT_TRUE(std::isfinite(simple_result.true_residual));
            EXPECT_TRUE(std::isfinite(simplec_result.true_residual));

            // The empirical envelope is only meaningful if every declared
            // method is exercised on every representative matrix. This is a
            // coverage gate, not an approximation-quality threshold.
            EXPECT_TRUE(method_records == 6);
            EXPECT_TRUE(std::isfinite(max_true_residual));
            EXPECT_TRUE(std::isfinite(max_reported_true_residual_gap));

            std::cout
                << "n8_schur_envelope"
                << " cells=" << cfg.n * cfg.n
                << " convection=" << cfg.convection
                << " method_records=" << method_records
                << " converged_approximation_records="
                << converged_approximation_records
                << " max_true_residual=" << max_true_residual
                << " max_reported_true_residual_gap="
                << max_reported_true_residual_gap
                << '\n';
        }
    });

    return run_all();
}