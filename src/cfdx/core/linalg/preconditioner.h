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

    void set_pressure_schur_diagonal(const std::vector<double>& diagonal) {
        pressure_schur_diagonal_ = diagonal;
    }

    bool setup(const SparseMatrix& A) override {
        if (n_cells_ == 0 || A.n_rows() != A.n_cols() ||
            A.n_rows() != 4 * n_cells_)
            return false;

        inv_velocity_diag_.assign(nv_, 0.0);
        inv_schur_diag_.assign(n_cells_, 0.0);
        schur_row_offsets_.clear();
        schur_columns_.clear();
        schur_values_.clear();
        schur_ilu_valid_ = false;

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

            double row_scale = 0.0;
            for (std::size_t k = row[j]; k < row[j + 1]; ++k)
                row_scale = std::max(row_scale, std::abs(val[k]));
            if (!std::isfinite(row_scale))
                return false;
            inv_velocity_diag_[j] = 1.0 / std::max(row_scale, 1.0);
        }

        // Form the pressure Schur approximation explicitly:
        //
        //   S_tilde = C - D diag(M)^-1 G .
        //
        // Unlike the former diagonal-only implementation, retain the complete
        // pressure-pressure sparsity pattern. This lets the pressure block
        // carry neighbouring pressure coupling instead of collapsing it to
        // one scalar per cell.
        schur_row_offsets_.resize(n_cells_ + 1, 0);
        std::vector<std::vector<std::pair<std::size_t, double>>> schur_rows(n_cells_);

        for (std::size_t c = 0; c < n_cells_; ++c) {
            const std::size_t pr = nv_ + c;
            auto& entries = schur_rows[c];

            auto add_entry = [&](std::size_t q, double contribution) {
                const std::size_t pc = q - nv_;
                for (auto& item : entries) {
                    if (item.first == pc) {
                        item.second += contribution;
                        return;
                    }
                }
                entries.emplace_back(pc, contribution);
            };

            for (std::size_t k = row[pr]; k < row[pr + 1]; ++k) {
                const std::size_t q = col[k];
                if (q >= nv_)
                    add_entry(q, val[k]);
            }

            // For each pressure-row D coefficient, multiply by the matching
            // momentum-row G coefficients. Only diagonal M^{-1} is used here;
            // a future block-M option can replace this without changing the
            // Schur interface.
            for (std::size_t k = row[pr]; k < row[pr + 1]; ++k) {
                const std::size_t u = col[k];
                if (u >= nv_)
                    continue;
                const double d = val[k];
                for (std::size_t gk = row[u]; gk < row[u + 1]; ++gk) {
                    const std::size_t q = col[gk];
                    if (q < nv_)
                        continue;
                    add_entry(q, -d * inv_velocity_diag_[u] * val[gk]);
                }
            }

            std::sort(entries.begin(), entries.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });

            double diag = 0.0;
            for (const auto& [q, value] : entries)
                if (q == c) diag = value;

            const double row_l1 = [&]() {
                double s = 0.0;
                for (const auto& [q, value] : entries) s += std::abs(value);
                return s;
            }();

            const double scale_floor =
                128.0 * std::numeric_limits<double>::epsilon() *
                std::max(row_l1, 1e-300);

            if (!std::isfinite(row_l1) || !std::isfinite(diag))
                return false;

            // A gauge identity row has no physical Schur operator.
            if (std::abs(diag) <= scale_floor && row_l1 <= scale_floor) {
                if (pressure_schur_diagonal_.size() == n_cells_ &&
                    std::isfinite(pressure_schur_diagonal_[c]) &&
                    pressure_schur_diagonal_[c] > 0.0) {
                    entries.emplace_back(c, pressure_schur_diagonal_[c]);
                    std::sort(entries.begin(), entries.end(),
                              [](const auto& a, const auto& b) { return a.first < b.first; });
                } else {
                    entries.emplace_back(c, 1.0);
                    std::sort(entries.begin(), entries.end(),
                              [](const auto& a, const auto& b) { return a.first < b.first; });
                }
            }

            schur_row_offsets_[c + 1] =
                schur_row_offsets_[c] + entries.size();
        }

        schur_columns_.reserve(schur_row_offsets_.back());
        schur_values_.reserve(schur_row_offsets_.back());
        for (const auto& entries : schur_rows) {
            for (const auto& [q, value] : entries) {
                if (!std::isfinite(value))
                    return false;
                schur_columns_.push_back(static_cast<std::uint32_t>(q));
                schur_values_.push_back(value);
            }
        }

        // Factor the explicit sparse Schur matrix in its own sparsity pattern.
        // ILU(0) is used deliberately: no hidden dense pressure solve and no
        // extra fill are introduced into the coupled preconditioner.
        schur_ilu_values_ = schur_values_;
        schur_ilu_valid_ = factor_schur_ilu0();
        if (!schur_ilu_valid_) {
            // The diagonal is still a valid algebraic fallback. Keep it
            // available so an ILU pivot failure does not make the entire
            // coupled system "not applicable".
            for (std::size_t c = 0; c < n_cells_; ++c) {
                const std::size_t rb = schur_row_offsets_[c];
                const std::size_t re = schur_row_offsets_[c + 1];
                double d = 0.0;
                for (std::size_t k = rb; k < re; ++k)
                    if (schur_columns_[k] == c) d = schur_values_[k];
                const double floor =
                    128.0 * std::numeric_limits<double>::epsilon() *
                    std::max(std::abs(d), 1e-300);
                if (!std::isfinite(d) || std::abs(d) <= floor)
                    return false;
                inv_schur_diag_[c] = 1.0 / d;
            }
        }
        return true;
    }

    bool apply(const Vector& r, Vector& z) const override {
        if (r.size() != 4 * n_cells_ || z.size() != r.size() ||
            inv_velocity_diag_.size() != nv_ ||
            schur_row_offsets_.size() != n_cells_ + 1)
            return false;

        Vector velocity_predictor(nv_);
        for (std::size_t i = 0; i < nv_; ++i)
            velocity_predictor(i) = inv_velocity_diag_[i] * r(i);

        Vector pressure_rhs(n_cells_);
        for (std::size_t c = 0; c < n_cells_; ++c) {
            double rhs_p = r(nv_ + c);
            const std::size_t rb = row_offsets_[nv_ + c];
            const std::size_t re = row_offsets_[nv_ + c + 1];
            for (std::size_t k = rb; k < re; ++k) {
                const std::size_t j = columns_[k];
                if (j < nv_)
                    rhs_p -= values_[k] * velocity_predictor(j);
            }
            pressure_rhs(c) = rhs_p;
        }

        if (schur_ilu_valid_) {
            if (!solve_schur_ilu0(pressure_rhs, z)) return false;
        } else {
            if (inv_schur_diag_.size() != n_cells_) return false;
            for (std::size_t c = 0; c < n_cells_; ++c)
                z(nv_ + c) = inv_schur_diag_[c] * pressure_rhs(c);
        }

        // Velocity back-substitution: M^-1(ru - G zp).
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

    const char* name() const override { return "CoupledBlockSchur-ILU0"; }

private:
    bool factor_schur_ilu0() {
        const std::size_t n = n_cells_;
        if (schur_row_offsets_.size() != n + 1 ||
            schur_ilu_values_.size() != schur_columns_.size())
            return false;

        auto find_position = [&](std::size_t row, std::size_t col) -> std::size_t {
            const auto first = schur_columns_.begin() + schur_row_offsets_[row];
            const auto last = schur_columns_.begin() + schur_row_offsets_[row + 1];
            const auto it = std::lower_bound(first, last, static_cast<std::uint32_t>(col));
            if (it == last || *it != static_cast<std::uint32_t>(col))
                return schur_columns_.size();
            return static_cast<std::size_t>(it - schur_columns_.begin());
        };

        double matrix_scale = 0.0;
        for (const double value : schur_ilu_values_)
            matrix_scale = std::max(matrix_scale, std::abs(value));
        if (!std::isfinite(matrix_scale) || matrix_scale == 0.0)
            return false;

        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t rb = schur_row_offsets_[i];
            const std::size_t re = schur_row_offsets_[i + 1];

            for (std::size_t pos = rb; pos < re; ++pos) {
                const std::size_t j = schur_columns_[pos];
                if (j >= i) break;

                const std::size_t pivot_pos = find_position(j, j);
                if (pivot_pos == schur_columns_.size()) return false;
                const double pivot = schur_ilu_values_[pivot_pos];
                const double pivot_floor =
                    128.0 * std::numeric_limits<double>::epsilon() *
                    std::max(matrix_scale, std::abs(pivot));
                if (!std::isfinite(pivot) || std::abs(pivot) <= pivot_floor)
                    return false;

                schur_ilu_values_[pos] /= pivot;
                const double lij = schur_ilu_values_[pos];

                for (std::size_t qpos = schur_row_offsets_[j];
                     qpos < schur_row_offsets_[j + 1]; ++qpos) {
                    const std::size_t q = schur_columns_[qpos];
                    if (q <= j) continue;
                    const std::size_t target = find_position(i, q);
                    if (target == schur_columns_.size()) continue;
                    schur_ilu_values_[target] -=
                        lij * schur_ilu_values_[qpos];
                }
            }

            const std::size_t diag_pos = find_position(i, i);
            if (diag_pos == schur_columns_.size()) return false;
            const double pivot = schur_ilu_values_[diag_pos];
            const double pivot_floor =
                128.0 * std::numeric_limits<double>::epsilon() *
                std::max(matrix_scale, std::abs(pivot));
            if (!std::isfinite(pivot) || std::abs(pivot) <= pivot_floor)
                return false;
        }
        return true;
    }

    bool solve_schur_ilu0(const Vector& rhs, Vector& x) const {
        if (rhs.size() != n_cells_ || x.size() != 4 * n_cells_)
            return false;
        Vector y(n_cells_, 0.0);

        // L has unit diagonal.
        for (std::size_t i = 0; i < n_cells_; ++i) {
            double value = rhs(i);
            for (std::size_t k = schur_row_offsets_[i];
                 k < schur_row_offsets_[i + 1]; ++k) {
                const std::size_t j = schur_columns_[k];
                if (j >= i) break;
                value -= schur_ilu_values_[k] * y(j);
            }
            y(i) = value;
        }

        // U back solve.
        for (std::size_t ii = 0; ii < n_cells_; ++ii) {
            const std::size_t i = n_cells_ - 1 - ii;
            double value = y(i);
            double diag = 0.0;
            for (std::size_t k = schur_row_offsets_[i];
                 k < schur_row_offsets_[i + 1]; ++k) {
                const std::size_t j = schur_columns_[k];
                if (j == i) diag = schur_ilu_values_[k];
                else if (j > i) value -= schur_ilu_values_[k] * x(nv_ + j);
            }
            if (!std::isfinite(diag) || diag == 0.0) return false;
            x(nv_ + i) = value / diag;
        }
        return true;
    }

    void clear() {
        inv_velocity_diag_.clear();
        inv_schur_diag_.clear();
        row_offsets_.clear();
        columns_.clear();
        values_.clear();
        schur_row_offsets_.clear();
        schur_columns_.clear();
        schur_values_.clear();
        schur_ilu_values_.clear();
        schur_ilu_valid_ = false;
    }

    std::size_t n_cells_ = 0;
    std::size_t nv_ = 0;
    std::vector<double> inv_velocity_diag_;
    std::vector<double> inv_schur_diag_;
    std::vector<double> pressure_schur_diagonal_;
    std::vector<std::uint32_t> row_offsets_;
    std::vector<std::uint32_t> columns_;
    std::vector<double> values_;
    std::vector<std::uint32_t> schur_row_offsets_;
    std::vector<std::uint32_t> schur_columns_;
    std::vector<double> schur_values_;
    std::vector<double> schur_ilu_values_;
    bool schur_ilu_valid_ = false;
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
