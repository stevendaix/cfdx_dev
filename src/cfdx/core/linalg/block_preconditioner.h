#pragma once

#include "sparse_matrix.h"
#include "vector.h"
#include <cstddef>
#include <stdexcept>
#include <algorithm>
#include <vector>

namespace cfdx {
namespace core {

/**
 * Block-Jacobi diagonal scaling preconditioner.
 *
 * This class deliberately does NOT invert the dense/coupled sub-block.
 * The name is retained for API compatibility; use a Schur-complement
 * implementation for coupled pressure/velocity systems.
 * Blocks are contiguous ranges; this implementation applies diagonal scaling
 * inside each block and deliberately does not invert the coupled block.
 *
 * This is intentionally independent of CFD physics. For pressure-velocity
 * coupling (e.g. SIMPLE), a Schur-complement or block-factorization
 * preconditioner is required for stronger coupling treatment.
 */
class BlockDiagonalPreconditioner {
public:
    struct Block { std::size_t begin, size; };

    void add_block(std::size_t begin, std::size_t size) {
        blocks_.push_back({begin, size});
    }

    bool setup(const SparseMatrix& A) {
        if (A.n_rows() != A.n_cols()) return false;
        A_ = &A;
        return validate_blocks(A.n_rows());
    }

    bool apply(const Vector& r, Vector& z) const {
        if (!A_ || r.size() != A_->n_rows()) return false;
        if (!validate_blocks(r.size())) return false;
        z.resize(r.size());
        z.fill(0.0);
        for (const auto& b : blocks_) {
            if (b.begin + b.size > r.size()) return false;
            for (std::size_t i = b.begin; i < b.begin + b.size; ++i) {
                const double d = (*A_)(i, i);
                if (d == 0.0) return false;
                z(i) = r(i) / d;
            }
        }
        return true;
    }
private:
    bool validate_blocks(std::size_t n) const {
        if (blocks_.empty()) return false;
        std::vector<bool> covered(n, false);
        for (const auto& b : blocks_) {
            if (b.size == 0 || b.begin > n || b.size > n - b.begin) return false;
            for (std::size_t i = b.begin; i < b.begin + b.size; ++i) {
                if (covered[i]) return false;
                covered[i] = true;
            }
        }
        return std::all_of(covered.begin(), covered.end(), [](bool v) { return v; });
    }

    const SparseMatrix* A_ = nullptr;
    std::vector<Block> blocks_;
};

} // namespace core
} // namespace cfdx
