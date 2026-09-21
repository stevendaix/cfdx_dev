#pragma once

#include "preconditioner.h"
#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace cfdx::core {

// Block-Jacobi diagonal scaling. The declared blocks are validated for
// complete, non-overlapping coverage; dense block inversion is not claimed.
class BlockDiagonalPreconditioner final : public Preconditioner {
public:
    explicit BlockDiagonalPreconditioner(std::vector<std::vector<std::size_t>> blocks)
        : blocks_(std::move(blocks)) {}

    bool setup(const SparseMatrix& A) override {
        if (A.n_rows() != A.n_cols()) return false;
        n_ = A.n_rows();
        inv_diag_.assign(n_, 0.0);
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        for (std::size_t i = 0; i < n_; ++i) {
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
        return validate();
    }

    bool apply(const Vector& r, Vector& z) const override {
        if (r.size() != n_ || z.size() != n_ || inv_diag_.size() != n_ || !validate()) return false;
        for (std::size_t i = 0; i < n_; ++i) z(i) = inv_diag_[i] * r(i);
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
    std::vector<double> inv_diag_;
    std::size_t n_{0};
};

} // namespace cfdx::core
