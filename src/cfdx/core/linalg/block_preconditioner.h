#pragma once

#include "sparse_matrix.h"
#include "vector.h"
#include <cstddef>
#include <stdexcept>

namespace cfdx {
namespace core {

/**
 * Generic block-Jacobi baseline preconditioner.
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
        return true;
    }

    bool apply(const Vector& r, Vector& z) const {
        if (!A_ || r.size() != A_->n_rows()) return false;
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
    const SparseMatrix* A_ = nullptr;
    std::vector<Block> blocks_;
};

} // namespace core
} // namespace cfdx
