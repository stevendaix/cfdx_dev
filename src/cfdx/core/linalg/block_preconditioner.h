#pragma once

#include "preconditioner.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace cfdx::core {

// Dense block-Jacobi preconditioner. Each declared diagonal block is inverted
// during setup and applied independently; off-block couplings are ignored.
class BlockDiagonalPreconditioner final : public Preconditioner {
public:
    explicit BlockDiagonalPreconditioner(std::vector<std::vector<std::size_t>> blocks)
        : blocks_(std::move(blocks)) {}

    bool setup(const SparseMatrix& A) override {
        if (A.n_rows() != A.n_cols()) return false;
        n_ = A.n_rows();
        if (!validate()) return false;

        inverse_blocks_.clear();
        inverse_blocks_.reserve(blocks_.size());
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();

        for (const auto& block : blocks_) {
            const std::size_t m = block.size();
            std::vector<double> aug(m * 2 * m, 0.0);
            auto at = [m, &aug](std::size_t i, std::size_t j) -> double& {
                return aug[i * (2 * m) + j];
            };

            for (std::size_t i = 0; i < m; ++i) {
                const std::size_t global_row = block[i];
                for (std::size_t k = row[global_row]; k < row[global_row + 1]; ++k) {
                    for (std::size_t j = 0; j < m; ++j) {
                        if (col[k] == block[j]) {
                            at(i, j) = val[k];
                            break;
                        }
                    }
                }
                at(i, m + i) = 1.0;
            }

            for (std::size_t pivot = 0; pivot < m; ++pivot) {
                std::size_t best = pivot;
                double best_abs = std::abs(at(pivot, pivot));
                for (std::size_t i = pivot + 1; i < m; ++i) {
                    const double a = std::abs(at(i, pivot));
                    if (a > best_abs) { best_abs = a; best = i; }
                }
                if (best_abs <= 1e-30) return false;

                if (best != pivot) {
                    for (std::size_t j = 0; j < 2 * m; ++j)
                        std::swap(at(pivot, j), at(best, j));
                }

                const double inv_pivot = 1.0 / at(pivot, pivot);
                for (std::size_t j = 0; j < 2 * m; ++j) at(pivot, j) *= inv_pivot;

                for (std::size_t i = 0; i < m; ++i) {
                    if (i == pivot) continue;
                    const double factor = at(i, pivot);
                    if (factor == 0.0) continue;
                    for (std::size_t j = 0; j < 2 * m; ++j)
                        at(i, j) -= factor * at(pivot, j);
                }
            }

            std::vector<double> inv(m * m);
            for (std::size_t i = 0; i < m; ++i)
                for (std::size_t j = 0; j < m; ++j)
                    inv[i * m + j] = at(i, m + j);
            inverse_blocks_.push_back(std::move(inv));
        }
        return true;
    }

    bool apply(const Vector& r, Vector& z) const override {
        if (r.size() != n_ || z.size() != n_ || inverse_blocks_.size() != blocks_.size()) return false;
        for (std::size_t b = 0; b < blocks_.size(); ++b) {
            const auto& block = blocks_[b];
            const auto& inv = inverse_blocks_[b];
            for (std::size_t i = 0; i < block.size(); ++i) {
                double value = 0.0;
                for (std::size_t j = 0; j < block.size(); ++j)
                    value += inv[i * block.size() + j] * r(block[j]);
                z(block[i]) = value;
            }
        }
        return true;
    }

    const char* name() const override { return "BlockJacobi"; }

private:
    bool validate() const {
        if (n_ == 0 || blocks_.empty()) return false;
        std::vector<bool> covered(n_, false);
        for (const auto& block : blocks_) {
            if (block.empty()) return false;
            for (std::size_t i : block) {
                if (i >= n_ || covered[i]) return false;
                covered[i] = true;
            }
        }
        return std::all_of(covered.begin(), covered.end(), [](bool v) { return v; });
    }

    std::vector<std::vector<std::size_t>> blocks_;
    std::vector<std::vector<double>> inverse_blocks_;
    std::size_t n_{0};
};

} // namespace cfdx::core
