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

            double block_norm_inf = 0.0;
            for (std::size_t r = 0; r < 4; ++r) {
                double row_sum = 0.0;
                for (std::size_t q = 0; q < 4; ++q)
                    row_sum += std::abs(block[4 * r + q]);
                block_norm_inf = std::max(block_norm_inf, row_sum);
            }
            if (!std::isfinite(block_norm_inf) || block_norm_inf == 0.0) {
                inv_blocks_.clear();
                return false;
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

            for (std::size_t r = 0; r < 4; ++r)
                for (std::size_t q = 0; q < 4; ++q)
                    inv_blocks_[c][4 * r + q] = aug[8 * r + 4 + q];

            // Verify the computed inverse instead of relying only on pivot
            // magnitudes. This catches badly scaled blocks whose pivots survive
            // the threshold but whose inverse is numerically unreliable.
            double inverse_residual = 0.0;
            for (std::size_t r = 0; r < 4; ++r) {
                for (std::size_t q = 0; q < 4; ++q) {
                    double value = 0.0;
                    for (std::size_t k = 0; k < 4; ++k)
                        value += block[4 * r + k] * inv_blocks_[c][4 * k + q];
                    const double target = r == q ? 1.0 : 0.0;
                    inverse_residual = std::max(
                        inverse_residual, std::abs(value - target));
                }
            }
            if (!std::isfinite(inverse_residual) ||
                inverse_residual > 1024.0 * std::numeric_limits<double>::epsilon()) {
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
