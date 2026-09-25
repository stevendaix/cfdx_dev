#pragma once

#include "sparse_matrix.h"
#include "vector.h"
#include "mixed_precision.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <iostream>
#include <cstdlib>
#include <utility>
#include <vector>

namespace cfdx::core {

class Preconditioner {
public:
    virtual ~Preconditioner() = default;
    virtual bool setup(const SparseMatrix& A) = 0;
    virtual bool apply(const Vector& r, Vector& z) const = 0;
    virtual const char* name() const = 0;
};

class IdentityPreconditioner final : public Preconditioner {
public:
    bool setup(const SparseMatrix& A) override {
        if (A.n_rows() != A.n_cols() || A.n_rows() == 0)
            return false;
        n_ = A.n_rows();
        return true;
    }

    bool apply(const Vector& r, Vector& z) const override {
        if (r.size() != n_ || z.size() != n_)
            return false;
        for (std::size_t i = 0; i < n_; ++i)
            z(i) = r(i);
        return true;
    }

    const char* name() const override { return "Identity"; }

private:
    std::size_t n_{0};
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
                // Scale the singularity test with the block itself. A fixed
                // absolute threshold is incorrect because momentum and pressure
                // coefficients have different physical units and magnitudes.
                const double pivot_floor =
                    64.0 * std::numeric_limits<double>::epsilon() * block_norm_inf;
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
 * Approximate block-LU preconditioner for the pressure-based coupled system
 *
 *     [ M  G ] [u] = [ru]
 *     [ D  0 ] [p]   [rp]
 *
 * with CFDX ordering [Ux, Uy, Uz, p]. The velocity block is inverted with
 * cell-local scalar diagonal factors and the pressure Schur complement is
 * approximated by the diagonal of -D M^{-1} G.
 */
class CoupledBlockSchurPreconditioner final : public Preconditioner {
public:
    explicit CoupledBlockSchurPreconditioner(std::size_t n_cells)
        : n_cells_(n_cells), nv_(3 * n_cells) {}

    bool setup(const SparseMatrix& A) override {
        if (n_cells_ == 0 || A.n_rows() != A.n_cols() ||
            A.n_rows() != 4 * n_cells_)
            return false;

        inv_velocity_diag_.assign(nv_, 0.0);
        inv_schur_diag_.assign(n_cells_, 0.0);

        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        row_offsets_.assign(row, row + A.n_rows() + 1);
        columns_.assign(col, col + A.nnz());
        values_.assign(val, val + A.nnz());

        auto diagonal = [&](std::size_t r) {
            double d = 0.0;
            bool found = false;
            for (std::size_t k = row[r]; k < row[r + 1]; ++k) {
                if (col[k] == r) {
                    d += val[k];
                    found = true;
                }
            }
            return std::pair<bool, double>{found, d};
        };

        for (std::size_t j = 0; j < nv_; ++j) {
            const auto [found, d] = diagonal(j);
            if (found && std::isfinite(d) &&
                std::abs(d) > 64.0 * std::numeric_limits<double>::epsilon()) {
                inv_velocity_diag_[j] = 1.0 / d;
                continue;
            }

            // Some CFDX transport rows are intentionally represented without
            // a stored diagonal (for example an inactive/degenerate component).
            // Do not make that representation detail disable the whole
            // saddle-point preconditioner. Use a conservative row scaling:
            // the row infinity norm, or unity for an entirely empty row.
            double row_scale = 0.0;
            for (std::size_t k = row[j]; k < row[j + 1]; ++k)
                row_scale = std::max(row_scale, std::abs(val[k]));
            if (!std::isfinite(row_scale))
                return false;
            inv_velocity_diag_[j] = 1.0 / std::max(row_scale, 1.0);
        }

        // The coupled FV continuity assembly already contains the
        // pressure Schur operator used by Rhie-Chow: its A_pp diagonal is
        // the sum of the positive face pressure-response coefficients D.
        // Use that discrete pressure operator directly as the Schur
        // approximation. This is consistent with the segregated pressure
        // equation and, unlike a purely local G_jc product, does not suffer
        // cancellation of opposing face-area contributions on Cartesian cells.
        for (std::size_t c = 0; c < n_cells_; ++c) {
            const std::size_t pressure_row = nv_ + c;
            const auto [has_pressure_diagonal, pressure_diagonal] =
                diagonal(pressure_row);

            double schur = 0.0;
            if (has_pressure_diagonal && std::isfinite(pressure_diagonal) &&
                std::abs(pressure_diagonal) >
                    64.0 * std::numeric_limits<double>::epsilon()) {
                schur = pressure_diagonal;
            } else {
                // Conservative algebraic fallback for matrices that expose
                // no explicit pressure diagonal: approximate -D M^-1 G.
                const std::size_t pressure_col = nv_ + c;
                for (std::size_t kd = row[pressure_row];
                     kd < row[pressure_row + 1]; ++kd) {
                    const std::size_t velocity_col = col[kd];
                    if (velocity_col >= nv_)
                        continue;
                    double g = 0.0;
                    for (std::size_t kg = row[velocity_col];
                         kg < row[velocity_col + 1]; ++kg) {
                        if (col[kg] == pressure_col)
                            g += val[kg];
                    }
                    schur -= val[kd] * inv_velocity_diag_[velocity_col] * g;
                }
            }

            if (!std::isfinite(schur) ||
                std::abs(schur) <=
                    64.0 * std::numeric_limits<double>::epsilon()) {
                if (std::getenv("CFDX_DEBUG_COUPLED"))
                    std::cerr << "COUPLED_SCHUR_SETUP pressure_failure cell=" << c
                              << " has_diag=" << has_pressure_diagonal
                              << " pressure_diag=" << pressure_diagonal
                              << " schur=" << schur << "\\n";
                clear();
                return false;
            }
            inv_schur_diag_[c] = 1.0 / schur;
        }
        return true;
    }

    bool apply(const Vector& r, Vector& z) const override {
        if (r.size() != 4 * n_cells_ || z.size() != r.size() ||
            inv_velocity_diag_.size() != nv_ ||
            inv_schur_diag_.size() != n_cells_)
            return false;

        Vector velocity_predictor(nv_);
        for (std::size_t i = 0; i < nv_; ++i)
            velocity_predictor(i) = inv_velocity_diag_[i] * r(i);

        // Pressure Schur RHS: rp - D M^{-1} ru.
        for (std::size_t c = 0; c < n_cells_; ++c) {
            double rhs_p = r(nv_ + c);
            const std::size_t rb = row_offsets_[nv_ + c];
            const std::size_t re = row_offsets_[nv_ + c + 1];
            for (std::size_t k = rb; k < re; ++k) {
                const std::size_t j = columns_[k];
                if (j < nv_)
                    rhs_p -= values_[k] * velocity_predictor(j);
            }
            z(nv_ + c) = inv_schur_diag_[c] * rhs_p;
        }

        // Velocity back-substitution: M^{-1}(ru - G zp).
        for (std::size_t j = 0; j < nv_; ++j) {
            double value = velocity_predictor(j);
            const std::size_t rb = row_offsets_[j];
            const std::size_t re = row_offsets_[j + 1];
            for (std::size_t k = rb; k < re; ++k) {
                const std::size_t pcol = columns_[k];
                if (pcol >= nv_)
                    value -= inv_velocity_diag_[j] * values_[k] * z(pcol);
            }
            z(j) = value;
        }
        return true;
    }

    const char* name() const override { return "CoupledBlockSchur"; }

private:
    void clear() {
        inv_velocity_diag_.clear();
        inv_schur_diag_.clear();
        row_offsets_.clear();
        columns_.clear();
        values_.clear();
    }

    std::size_t n_cells_ = 0;
    std::size_t nv_ = 0;
    std::vector<double> inv_velocity_diag_;
    std::vector<double> inv_schur_diag_;
    std::vector<std::uint32_t> row_offsets_;
    std::vector<std::uint32_t> columns_;
    std::vector<double> values_;
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
