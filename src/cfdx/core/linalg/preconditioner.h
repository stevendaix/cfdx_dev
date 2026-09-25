#pragma once

#include "sparse_matrix.h"
#include "vector.h"
#include "mixed_precision.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace cfdx::core {

class Preconditioner {
public:
    virtual ~Preconditioner() = default;
    virtual bool setup(const SparseMatrix& A) = 0;
    virtual bool apply(const Vector& r, Vector& z) const = 0;
    virtual const char* name() const = 0;
};

class JacobiPreconditioner final : public Preconditioner {
public:
    bool setup(const SparseMatrix& A) override;
    bool apply(const Vector& r, Vector& z) const override;
    const char* name() const override { return "Jacobi"; }
private:
    std::vector<double> inv_diag_;
};

inline bool JacobiPreconditioner::setup(const SparseMatrix& A) {
    if (A.n_rows() != A.n_cols()) return false;
    const auto* row = A.row_offsets_data();
    const auto* col = A.columns_data();
    const auto* val = A.values_data();
    inv_diag_.assign(A.n_rows(), 0.0);
    for (std::size_t i = 0; i < A.n_rows(); ++i) {
        bool found = false;
        for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
            if (col[k] == i) {
                if (val[k] == 0.0) { inv_diag_.clear(); return false; }
                inv_diag_[i] = 1.0 / val[k];
                found = true;
                break;
            }
        }
        if (!found) { inv_diag_.clear(); return false; }
    }
    return true;
}

inline bool JacobiPreconditioner::apply(const Vector& r, Vector& z) const {
    if (r.size() != inv_diag_.size() || z.size() != r.size()) return false;
    for (std::size_t i = 0; i < r.size(); ++i) z(i) = inv_diag_[i] * r(i);
    return true;
}

// Cell-wise block Jacobi for the 4-variable pressure-based coupled system.
//
// The coupled unknown ordering is [Ux(0..N-1), Uy, Uz, p].  A scalar Jacobi
// preconditioner ignores the strong intra-cell velocity/pressure coupling.
// This preconditioner extracts the 4x4 diagonal block belonging to each cell
// and applies its exact dense inverse. It is deliberately local: it does not
// hide any global solve or approximate the Schur complement.
class CellBlockJacobiPreconditioner final : public Preconditioner {
public:
    explicit CellBlockJacobiPreconditioner(std::size_t n_cells) : n_cells_(n_cells) {}

    bool setup(const SparseMatrix& A) override {
        if (n_cells_ == 0 || A.n_rows() != A.n_cols() ||
            A.n_rows() != 4 * n_cells_)
            return false;

        inv_blocks_.assign(n_cells_, {});
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();

        for (std::size_t c = 0; c < n_cells_; ++c) {
            std::array<double, 16> block{};
            for (std::size_t r = 0; r < 4; ++r) {
                const std::size_t global_row = r * n_cells_ + c;
                for (std::size_t k = row[global_row]; k < row[global_row + 1]; ++k) {
                    const std::size_t global_col = col[k];
                    for (std::size_t q = 0; q < 4; ++q) {
                        if (global_col == q * n_cells_ + c) {
                            block[4 * r + q] += val[k];
                            break;
                        }
                    }
                }
            }

            // Gauss-Jordan inversion with partial pivoting. The block is only
            // 4x4, so an explicit local factorization is cheaper and clearer
            // than routing it through the global sparse solver.
            std::array<double, 32> aug{};
            for (std::size_t r = 0; r < 4; ++r) {
                for (std::size_t q = 0; q < 4; ++q)
                    aug[8 * r + q] = block[4 * r + q];
                aug[8 * r + 4 + r] = 1.0;
            }

            // Momentum and pressure rows can differ by several orders of
            // magnitude. Row-scale the local block before pivoting so the
            // factorization is insensitive to those physical units.
            std::array<double, 4> row_scale{};
            double block_norm_inf = 0.0;
            for (std::size_t r = 0; r < 4; ++r) {
                double row_sum = 0.0;
                for (std::size_t q = 0; q < 4; ++q)
                    row_sum += std::abs(block[4 * r + q]);
                row_scale[r] = row_sum;
                block_norm_inf = std::max(block_norm_inf, row_sum);
            }
            if (!std::isfinite(block_norm_inf) || block_norm_inf == 0.0) {
                inv_blocks_.clear();
                return false;
            }
            for (std::size_t r = 0; r < 4; ++r) {
                if (!(row_scale[r] > 0.0) || !std::isfinite(row_scale[r])) {
                    inv_blocks_.clear();
                    return false;
                }
                const double scale = 1.0 / row_scale[r];
                for (std::size_t q = 0; q < 4; ++q)
                    aug[8 * r + q] = block[4 * r + q] * scale;
            }

            for (std::size_t k = 0; k < 4; ++k) {
                std::size_t pivot = k;
                double pivot_abs = std::abs(aug[8 * k + k]);
                for (std::size_t r = k + 1; r < 4; ++r) {
                    const double candidate = std::abs(aug[8 * r + k]);
                    if (candidate > pivot_abs) {
                        pivot = r;
                        pivot_abs = candidate;
                    }
                }
                // The matrix being pivoted has already been row-scaled.
                // pivot_abs is therefore dimensionless and O(1). Comparing it
                // with the *unscaled* block norm is incorrect and can reject a
                // perfectly invertible block whenever pressure and momentum
                // coefficients have different units/scales.
                //
                // Use the norm of the scaled matrix instead. Each scaled row
                // has infinity norm <= 1 by construction.
                const double scaled_norm_inf = 1.0;
                const double pivot_floor =
                    64.0 * std::numeric_limits<double>::epsilon() * scaled_norm_inf;
                if (!std::isfinite(pivot_abs) || pivot_abs <= pivot_floor) {
                    inv_blocks_.clear();
                    return false;
                }
                if (pivot != k) {
                    for (std::size_t q = 0; q < 8; ++q)
                        std::swap(aug[8 * k + q], aug[8 * pivot + q]);
                }
                const double diagonal = aug[8 * k + k];
                for (std::size_t q = 0; q < 8; ++q)
                    aug[8 * k + q] /= diagonal;
                for (std::size_t r = 0; r < 4; ++r) {
                    if (r == k) continue;
                    const double factor = aug[8 * r + k];
                    if (factor == 0.0) continue;
                    for (std::size_t q = 0; q < 8; ++q)
                        aug[8 * r + q] -= factor * aug[8 * k + q];
                }
            }

            // The inverse above is B^{-1} for B = S A. Hence
            // A^{-1} = B^{-1} S.
            for (std::size_t r = 0; r < 4; ++r)
                for (std::size_t q = 0; q < 4; ++q)
                    inv_blocks_[c][4 * r + q] =
                        aug[8 * r + 4 + q] * row_scale[q];

            double inverse_residual = 0.0;
            for (std::size_t r = 0; r < 4; ++r) {
                double row_error = 0.0;
                double row_norm = 0.0;
                for (std::size_t q = 0; q < 4; ++q) {
                    double value = 0.0;
                    for (std::size_t k = 0; k < 4; ++k)
                        value += block[4 * r + k] * inv_blocks_[c][4 * k + q];
                    const double target = r == q ? 1.0 : 0.0;
                    row_error = std::max(row_error, std::abs(value - target));
                    row_norm = std::max(row_norm, std::abs(block[4 * r + q]));
                }
                inverse_residual = std::max(
                    inverse_residual, row_error / std::max(row_norm, 1e-300));
            }
            if (!std::isfinite(inverse_residual) || inverse_residual > 1e-10) {
                inv_blocks_.clear();
                return false;
            }
        }
        return true;
    }

    bool apply(const Vector& r, Vector& z) const override {
        if (r.size() != 4 * n_cells_ || z.size() != r.size() ||
            inv_blocks_.size() != n_cells_)
            return false;
        for (std::size_t c = 0; c < n_cells_; ++c) {
            const double in[4] = {
                r(c), r(n_cells_ + c), r(2 * n_cells_ + c), r(3 * n_cells_ + c)};
            double out[4]{};
            for (std::size_t row = 0; row < 4; ++row)
                for (std::size_t col_idx = 0; col_idx < 4; ++col_idx)
                    out[row] += inv_blocks_[c][4 * row + col_idx] * in[col_idx];
            z(c) = out[0];
            z(n_cells_ + c) = out[1];
            z(2 * n_cells_ + c) = out[2];
            z(3 * n_cells_ + c) = out[3];
        }
        return true;
    }

    const char* name() const override { return "CellBlockJacobi"; }

private:
    std::size_t n_cells_ = 0;
    std::vector<std::array<double, 16>> inv_blocks_;
};


/**
 * SIMPLE-type block Schur preconditioner for the monolithic
 * [ M  G ] [u] = [fu]
 * [ D  C ] [p]   [fp]
 * pressure-velocity system.
 *
 * The velocity block is approximated by its diagonal and the pressure
 * Schur complement by the diagonal of D diag(M)^-1 G - C.  This is the
 * algebraic analogue of the SIMPLE pressure correction and, unlike scalar
 * Jacobi, explicitly represents the velocity-pressure coupling.
 *
 * Application:
 *   zu = Mdiag^-1 * ru
 *   zp = (D zu - rp) / (D Mdiag^-1 G - C)_diag
 *   zu = Mdiag^-1 * (ru - G zp)
 *
 * The reference-pressure row is treated as an identity equation.
 */
class CoupledBlockSchurPreconditioner final : public Preconditioner {
public:
    explicit CoupledBlockSchurPreconditioner(std::size_t n_cells)
        : n_cells_(n_cells), nv_(3 * n_cells) {}

    bool setup(const SparseMatrix& A) override {
        const std::size_t n = 4 * n_cells_;
        if (n_cells_ == 0 || A.n_rows() != n || A.n_cols() != n)
            return false;

        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();

        inv_velocity_diag_.assign(nv_, 0.0);
        for (std::size_t i = 0; i < nv_; ++i) {
            double d = 0.0;
            bool found = false;
            for (std::uint32_t k = row[i]; k < row[i + 1]; ++k) {
                if (col[k] == i) {
                    d += val[k];
                    found = true;
                }
            }
            if (!found || !std::isfinite(d) || d == 0.0)
                return false;
            inv_velocity_diag_[i] = 1.0 / d;
        }

        pressure_velocity_rows_.assign(n_cells_, {});
        velocity_pressure_rows_.assign(nv_, {});
        for (std::size_t p = 0; p < n_cells_; ++p) {
            const std::size_t prow = nv_ + p;
            for (std::uint32_t k = row[prow]; k < row[prow + 1]; ++k) {
                if (col[k] < nv_) pressure_velocity_rows_[p].push_back({col[k], val[k]});
            }
        }
        for (std::size_t i = 0; i < nv_; ++i) {
            for (std::uint32_t k = row[i]; k < row[i + 1]; ++k) {
                if (col[k] >= nv_) velocity_pressure_rows_[i].push_back({col[k] - nv_, val[k]});
            }
        }

        schur_diag_.assign(n_cells_, 0.0);
        for (std::size_t p = 0; p < n_cells_; ++p) {
            const std::size_t prow = nv_ + p;

            // The pressure reference row is an identity row after gauge fixing.
            double cdiag = 0.0;
            bool pressure_identity = false;
            for (std::uint32_t k = row[prow]; k < row[prow + 1]; ++k) {
                if (col[k] == prow) {
                    cdiag += val[k];
                    if (std::abs(cdiag - 1.0) <= 64.0 * std::numeric_limits<double>::epsilon())
                        pressure_identity = true;
                }
            }
            if (pressure_identity) {
                schur_diag_[p] = 1.0;
                continue;
            }

            // Form the diagonal of D diag(M)^-1 G - C.
            // The pressure row provides D entries. For each coupled velocity
            // DOF, find the matching pressure coefficient in its momentum row.
            double s = -cdiag;
            for (std::uint32_t dk = row[prow]; dk < row[prow + 1]; ++dk) {
                const std::size_t u = col[dk];
                if (u >= nv_) continue;
                const double d = val[dk];

                double g = 0.0;
                for (std::uint32_t gk = row[u]; gk < row[u + 1]; ++gk) {
                    if (col[gk] == prow) {
                        g += val[gk];
                    }
                }
                s += d * inv_velocity_diag_[u] * g;
            }

            if (!std::isfinite(s) || std::abs(s) <=
                    64.0 * std::numeric_limits<double>::epsilon()) {
                return false;
            }
            schur_diag_[p] = s;
        }
        return true;
    }

    bool apply(const Vector& r, Vector& z) const override {
        if (r.size() != 4 * n_cells_ || z.size() != r.size() ||
            inv_velocity_diag_.size() != nv_ ||
            schur_diag_.size() != n_cells_)
            return false;

        // First velocity predictor: diag(M)^-1 r_u.
        for (std::size_t i = 0; i < nv_; ++i)
            z(i) = inv_velocity_diag_[i] * r(i);

        // Pressure correction from the approximate Schur complement.
        // rhs_p = D M^-1 r_u - r_p.
        for (std::size_t p = 0; p < n_cells_; ++p) {
            double rhs = -r(nv_ + p);
            // D is stored in the pressure row of the original matrix. The
            // Schur application is reconstructed from the matrix-independent
            // sparsity captured during setup below.
            for (const auto& e : pressure_velocity_rows_[p])
                rhs += e.second * z(e.first);
            z(nv_ + p) = rhs / schur_diag_[p];
        }

        // Velocity correction: M^-1 (r_u - G z_p).
        for (std::size_t i = 0; i < nv_; ++i) {
            double rhs = r(i);
            for (const auto& e : velocity_pressure_rows_[i])
                rhs -= e.second * z(nv_ + e.first);
            z(i) = inv_velocity_diag_[i] * rhs;
        }
        return true;
    }

    const char* name() const override { return "CoupledBlockSchur"; }

private:
    std::size_t n_cells_ = 0;
    std::size_t nv_ = 0;
    std::vector<double> inv_velocity_diag_;
    std::vector<double> schur_diag_;
    // Sparse coupling maps cached at setup so apply() is O(nnz) and does not
    // rescan the matrix.
    std::vector<std::vector<std::pair<std::size_t, double>>> pressure_velocity_rows_;
    std::vector<std::vector<std::pair<std::size_t, double>>> velocity_pressure_rows_;
};

class MixedPrecisionJacobiPreconditioner final : public Preconditioner {
public:
    explicit MixedPrecisionJacobiPreconditioner(SolverPrecision precision=SolverPrecision::FP32)
        : precision_(precision) {}
    bool setup(const SparseMatrix& A) override {
        if (A.n_rows()!=A.n_cols()) return false;
        const auto* row=A.row_offsets_data(); const auto* col=A.columns_data(); const auto* val=A.values_data();
        inv_diag_.assign(A.n_rows(),0.0);
        for(std::size_t i=0;i<A.n_rows();++i){
            bool found=false;
            for(std::size_t k=row[i];k<row[i+1];++k) if(col[k]==i){
                if(!std::isfinite(val[k]) || val[k]==0.0) { inv_diag_.clear(); return false; }
                inv_diag_[i]=1.0/val[k]; found=true; break;
            }
            if(!found) { inv_diag_.clear(); return false; }
        }
        return true;
    }
    bool apply(const Vector& r, Vector& z) const override {
        if(r.size()!=inv_diag_.size() || z.size()!=r.size()) return false;
        for(std::size_t i=0;i<r.size();++i){
            if(precision_==SolverPrecision::FP32)
                z(i)=static_cast<double>(static_cast<float>(r(i))*static_cast<float>(inv_diag_[i]));
            else
                z(i)=static_cast<double>(r(i))*static_cast<double>(inv_diag_[i]);
        }
        return true;
    }
    const char* name() const override { return "MixedPrecisionJacobi"; }
private:
    SolverPrecision precision_;
    std::vector<double> inv_diag_;
};

} // namespace cfdx::core
