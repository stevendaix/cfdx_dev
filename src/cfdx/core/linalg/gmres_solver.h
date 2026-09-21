#pragma once

#include "sparse_matrix.h"
#include "cg_solver.h"
#include "vector.h"
#include "preconditioner.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

namespace cfdx::core {

struct LinearOperator {
    std::size_t size = 0;
    std::function<void(const Vector&, Vector&)> apply;
};

inline SolverResult solve_gmres_operator(
    const LinearOperator& op,
    const Vector& b,
    Vector& x,
    int restart = 30,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12,
    const Preconditioner* preconditioner = nullptr)
{
    SolverResult result;
    if (op.size == 0 || !op.apply || b.size() != op.size || x.size() != op.size ||
        restart <= 0 || max_iter == 0 || tolerance <= 0.0) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }

    const std::size_t n = op.size;
    const double b_norm = b.norm2();
    const double tol = tolerance * std::max(b_norm, 1.0);
    const int m = std::min<int>(restart, static_cast<int>(n));

        auto apply_operator = [&](const std::vector<double>& in, std::vector<double>& out) {
        Vector vin(n), vout(n);
        for(std::size_t i=0;i<n;++i) vin(i)=in[i];
        op.apply(vin,vout);
        out.resize(n);
        for(std::size_t i=0;i<n;++i) out[i]=vout(i);
    };

    auto residual = [&](std::vector<double>& r) {
        std::vector<double> ax;
        std::vector<double> xv(n);
        for(std::size_t i=0;i<n;++i) xv[i]=x(i);
        apply_operator(xv,ax);
        r.resize(n);
        for(std::size_t i=0;i<n;++i) r[i]=b(i)-ax[i];
    };

    std::vector<double> r, w(n), z(n);
    residual(r);
    auto norm = [](const std::vector<double>& v) {
        double s = 0.0;
        for (double a : v) s += a * a;
        return std::sqrt(s);
    };

    double beta = norm(r);
    if (beta <= tol) {
        result.status = SolverStatus::CONVERGED;
        result.iterations = 0;
        result.residual = beta;
        result.residual_relative = beta / std::max(b_norm, 1.0);
        return result;
    }

    std::vector<std::vector<double>> V(static_cast<std::size_t>(m + 1), std::vector<double>(n));
    std::vector<std::vector<double>> Z(static_cast<std::size_t>(m), std::vector<double>(n));
    std::vector<std::vector<double>> H(static_cast<std::size_t>(m + 1), std::vector<double>(m, 0.0));
    std::vector<double> cs(m, 0.0), sn(m, 0.0), g(m + 1, 0.0), y(m, 0.0);

    std::size_t iterations = 0;
    while (iterations < max_iter) {
        residual(r);
        beta = norm(r);
        if (beta <= tol) break;

        for (std::size_t i = 0; i < n; ++i) V[0][i] = r[i] / beta;
        for (auto& row : H) std::fill(row.begin(), row.end(), 0.0);
        std::fill(cs.begin(), cs.end(), 0.0);
        std::fill(sn.begin(), sn.end(), 0.0);
        std::fill(g.begin(), g.end(), 0.0);
        g[0] = beta;

        int used = 0;
        double estimated_residual = beta;

        for (int j = 0; j < m && iterations < max_iter; ++j) {
            if (preconditioner) {
                for (std::size_t i = 0; i < n; ++i) z[i] = V[j][i];
                Vector rv(n), zv(n);
                for (std::size_t i = 0; i < n; ++i) rv(i) = z[i];
                if (!preconditioner->apply(rv, zv)) {
                    result.status = SolverStatus::NOT_APPLICABLE;
                    result.iterations = iterations;
                    return result;
                }
                for (std::size_t i = 0; i < n; ++i) Z[j][i] = zv(i);
            } else {
                Z[j] = V[j];
            }

            apply_operator(Z[j], w);

            for (int i = 0; i <= j; ++i) {
                double h = 0.0;
                for (std::size_t k = 0; k < n; ++k) h += V[i][k] * w[k];
                H[i][j] = h;
                for (std::size_t k = 0; k < n; ++k) w[k] -= h * V[i][k];
            }

            H[j + 1][j] = norm(w);
            if (H[j + 1][j] > 0.0) {
                for (std::size_t k = 0; k < n; ++k) V[j + 1][k] = w[k] / H[j + 1][j];
            }

            for (int i = 0; i < j; ++i) {
                const double h0 = H[i][j];
                const double h1 = H[i + 1][j];
                H[i][j] = cs[i] * h0 + sn[i] * h1;
                H[i + 1][j] = -sn[i] * h0 + cs[i] * h1;
            }

            const double a = H[j][j];
            const double b2 = H[j + 1][j];
            const double rho = std::hypot(a, b2);
            if (rho > 0.0) {
                cs[j] = a / rho;
                sn[j] = b2 / rho;
                H[j][j] = rho;
                H[j + 1][j] = 0.0;
                const double gj = g[j];
                const double gj1 = g[j + 1];
                g[j] = cs[j] * gj + sn[j] * gj1;
                g[j + 1] = -sn[j] * gj + cs[j] * gj1;
            }

            ++iterations;
            estimated_residual = std::abs(g[j + 1]);
            used = j + 1;
            if (estimated_residual <= tol || H[j + 1][j] == 0.0) break;
        }

        if (used == 0) break;

        std::fill(y.begin(), y.end(), 0.0);
        for (int i = used - 1; i >= 0; --i) {
            double sum = g[i];
            for (int j = i + 1; j < used; ++j) sum -= H[i][j] * y[j];
            if (std::abs(H[i][i]) <= 1e-30) {
                result.status = SolverStatus::DIVERGED;
                result.iterations = iterations;
                result.residual = estimated_residual;
                result.residual_relative = estimated_residual / std::max(b_norm, 1.0);
                return result;
            }
            y[i] = sum / H[i][i];
        }

        for (int j = 0; j < used; ++j)
            for (std::size_t i = 0; i < n; ++i)
                x(i) += Z[j][i] * y[j];

        residual(r);
        beta = norm(r);
        if (beta <= tol) {
            result.status = SolverStatus::CONVERGED;
            result.iterations = iterations;
            result.residual = beta;
            result.residual_relative = beta / std::max(b_norm, 1.0);
            return result;
        }
    }

    result.status = SolverStatus::MAX_ITER_REACHED;
    result.iterations = iterations;
    result.residual = beta;
    result.residual_relative = beta / std::max(b_norm, 1.0);
    return result;
}

inline SolverResult solve_gmres(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    int restart = 30,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12,
    const Preconditioner* preconditioner = nullptr)
{
    LinearOperator op;
    op.size=A.n_rows();
    op.apply=[&A](const Vector& in, Vector& out){ A.matvec(in,out); };
    return solve_gmres_operator(op,b,x,restart,max_iter,tolerance,preconditioner);
}

} // namespace cfdx::core
