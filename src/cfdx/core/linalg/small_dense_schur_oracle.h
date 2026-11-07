#pragma once

// Diagnostic reference for small, gauge-fixed, possibly NON-symmetric FV
// saddle-point systems. Dense LU is a verification oracle, not a scalable
// production replacement for AMG/PCD or a guarantee of CFD convergence.

#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/schur_approximation.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace cfdx::core {

class SmallDenseCertifiedLU {
public:
    explicit SmallDenseCertifiedLU(std::size_t size_limit = 64)
        : size_limit_(size_limit) {}

    bool setup(const SparseMatrix& matrix) {
        const std::size_t n = matrix.n_rows();
        if (n == 0 || n > size_limit_ || matrix.n_cols() != n ||
            n > std::numeric_limits<std::size_t>::max() / n)
            return false;
        std::vector<double> scaled(n*n, 0.0), row_scale(n, 0.0);
        std::vector<std::size_t> pivots(n);
        const auto* rows = matrix.row_offsets_data();
        const auto* cols = matrix.columns_data();
        const auto* values = matrix.values_data();
        if (rows[0] != 0 || rows[n] != matrix.nnz()) return false;
        for (std::size_t i = 0; i < n; ++i) {
            if (rows[i] > rows[i+1] || rows[i+1] > matrix.nnz()) return false;
            for (std::size_t k = rows[i]; k < rows[i+1]; ++k) {
                if (cols[k] >= n || !std::isfinite(values[k])) return false;
                scaled[i*n+cols[k]] += values[k];
            }
            double row_norm = 0.0;
            for (std::size_t j = 0; j < n; ++j)
                row_norm += std::abs(scaled[i*n+j]);
            if (!std::isfinite(row_norm) || !(row_norm > 0.0)) return false;
            row_scale[i] = row_norm;
            for (std::size_t j = 0; j < n; ++j)
                scaled[i*n+j] /= row_norm;
        }
        // Partial pivoting on the row-scaled matrix. The pivot history is
        // replayed on a separately row-scaled right hand side by solve().
        for (std::size_t k=0; k<n; ++k) {
            std::size_t pivot = k;
            double best = std::abs(scaled[k*n+k]);
            for (std::size_t i=k+1; i<n; ++i) {
                const double value = std::abs(scaled[i*n+k]);
                if (value > best) { best = value; pivot = i; }
            }
            if (!std::isfinite(best) || best <= 1e-13) return false;
            pivots[k] = pivot;
            if (pivot != k)
                for (std::size_t j=0; j<n; ++j)
                    std::swap(scaled[k*n+j], scaled[pivot*n+j]);
            for (std::size_t i=k+1; i<n; ++i) {
                scaled[i*n+k] /= scaled[k*n+k];
                if (!std::isfinite(scaled[i*n+k])) return false;
                for (std::size_t j=k+1; j<n; ++j)
                    scaled[i*n+j] -= scaled[i*n+k] * scaled[k*n+j];
            }
        }
        for (double value:scaled)
            if (!std::isfinite(value)) return false;
        original_ = matrix;
        factors_ = std::move(scaled);
        row_scale_ = std::move(row_scale);
        pivots_ = std::move(pivots);
        return true;
    }

    bool solve(const Vector& rhs, Vector& x) const {
        const std::size_t n = pivots_.size();
        if (n == 0 || rhs.size() != n || !rhs.is_valid()) return false;
        std::vector<double> y(n, 0.0);
        for (std::size_t i=0; i<n; ++i) y[i] = rhs(i) / row_scale_[i];
        for (std::size_t k=0; k<n; ++k) std::swap(y[k],y[pivots_[k]]);
        for (std::size_t i=0; i<n; ++i)
            for (std::size_t j=0; j<i; ++j)
                y[i] -= factors_[i*n+j]*y[j];
        for (std::size_t i=n; i-- > 0;) {
            for (std::size_t j=i+1; j<n; ++j)
                y[i] -= factors_[i*n+j]*y[j];
            y[i] /= factors_[i*n+i];
            if (!std::isfinite(y[i])) return false;
        }
        Vector candidate(n, 0.0), actual(n, 0.0);
        for (std::size_t i=0; i<n; ++i) candidate(i) = y[i];
        original_.matvec_into(candidate,actual);
        if (!actual.is_valid()) return false;
        double residual_norm=0.0, rhs_norm=0.0;
        for (std::size_t i=0; i<n; ++i) {
            residual_norm=std::hypot(residual_norm,actual(i)-rhs(i));
            rhs_norm=std::hypot(rhs_norm,rhs(i));
        }
        // This separate true-residual guard avoids promoting a numerically
        // singular reference solve or unstable local LU as ground truth.
        if (!std::isfinite(residual_norm) || !std::isfinite(rhs_norm) ||
            residual_norm > 1e-12 + 1e-10*rhs_norm) return false;
        x=std::move(candidate);
        return true;
    }

private:
    std::size_t size_limit_;
    SparseMatrix original_;
    std::vector<double> factors_, row_scale_;
    std::vector<std::size_t> pivots_;
};

// Accurate small-matrix Schur oracle: S = C - D * Auu^{-1} * G, with no
// assumption that G = -D^T, that S is SPD, or that C is zero. Replaces the
// incorrect practice of requiring CG/positive p^TSp on an arbitrary
// nonsymmetric saddle-point Schur. Both LU factorizations are bounded.
class SmallDenseExactSchurInverse final : public SchurApproximation {
public:
    explicit SmallDenseExactSchurInverse(std::size_t maximum_unknowns = 64)
        : size_limit_(maximum_unknowns),
          velocity_lu_(maximum_unknowns), pressure_lu_(maximum_unknowns) {}

    const char* name() const noexcept override { return "small_dense_exact_schur_oracle"; }
    SchurAction action() const noexcept override { return SchurAction::InverseOperator; }

    bool setup(const BlockOperator& blocks) override {
        if (!blocks.is_valid() || blocks.velocity_size() == 0 ||
            blocks.pressure_size() == 0 || blocks.velocity_size() > size_limit_ ||
            blocks.pressure_size() > size_limit_) return false;
        SmallDenseCertifiedLU velocity(size_limit_), pressure(size_limit_);
        if (!velocity.setup(blocks.Auu())) return false;
        const std::size_t np=blocks.pressure_size();
        const std::size_t nu=blocks.velocity_size();
        SparseMatrix S(np,np);
        for (std::size_t col=0; col<np; ++col) {
            Vector basis(np,0.0), gp(nu,0.0), y(nu,0.0);
            Vector cp(np,0.0), dy(np,0.0);
            basis(col)=1.0;
            blocks.G().matvec_into(basis,gp);
            if (!gp.is_valid() || !velocity.solve(gp,y)) return false;
            blocks.C().matvec_into(basis,cp);
            blocks.D().matvec_into(y,dy);
            if (!cp.is_valid() || !dy.is_valid()) return false;
            for (std::size_t row=0; row<np; ++row) {
                const double value = cp(row)-dy(row);
                if (!std::isfinite(value)) return false;
                if (value != 0.0) S.push_back(row,col,value);
            }
        }
        S.finalize();
        if (!pressure.setup(S)) return false;
        // Full success before mutating any observable state.
        velocity_lu_=std::move(velocity);
        pressure_lu_=std::move(pressure);
        Auu_=blocks.Auu(); G_=blocks.G(); D_=blocks.D(); C_=blocks.C();
        schur_matrix_=std::move(S);
        ready_=true;
        return true;
    }

    bool update_values(const BlockOperator& blocks) override {
        if (!ready_ || !blocks.is_valid() ||
            !same_pattern(Auu_,blocks.Auu()) ||
            !same_pattern(G_,blocks.G()) ||
            !same_pattern(D_,blocks.D()) ||
            !same_pattern(C_,blocks.C())) return false;
        return setup(blocks); // transactional on singular/invalid updates
    }

    bool apply(const Vector& rhs, Vector& output) const override {
        return ready_ && pressure_lu_.solve(rhs,output);
    }

    bool solve_velocity(const Vector& rhs, Vector& output) const {
        return ready_ && velocity_lu_.solve(rhs,output);
    }

    const SparseMatrix& schur_matrix() const noexcept { return schur_matrix_; }

private:
    static bool same_pattern(const SparseMatrix& a, const SparseMatrix& b) {
        if (a.n_rows()!=b.n_rows() || a.n_cols()!=b.n_cols() || a.nnz()!=b.nnz())
            return false;
        for (std::size_t i=0; i<=a.n_rows(); ++i)
            if (a.row_offsets_data()[i]!=b.row_offsets_data()[i]) return false;
        for (std::size_t k=0; k<a.nnz(); ++k)
            if (a.columns_data()[k]!=b.columns_data()[k]) return false;
        return true;
    }

    std::size_t size_limit_;
    SmallDenseCertifiedLU velocity_lu_, pressure_lu_;
    SparseMatrix Auu_, G_, D_, C_, schur_matrix_;
    bool ready_{false};
};

} // namespace cfdx::core
