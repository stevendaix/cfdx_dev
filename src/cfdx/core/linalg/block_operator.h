#pragma once

#include "cfdx/core/linalg/sparse_matrix.h"

#include <cstddef>
#include <stdexcept>

namespace cfdx::core {

// Algebraic block view of a coupled pressure-velocity operator:
//
//     A = [ Auu  G ]
//         [ D    C ]
//
// The blocks are deliberately exposed without assuming G = -D^T.  This is
// required for collocated finite-volume discretisations, where the discrete
// gradient/divergence pair need not be exact transposes.
class BlockOperator {
public:
    BlockOperator(const SparseMatrix& Auu,
                  const SparseMatrix& G,
                  const SparseMatrix& D,
                  const SparseMatrix& C)
        : Auu_(&Auu), G_(&G), D_(&D), C_(&C) {}

    const SparseMatrix& Auu() const noexcept { return *Auu_; }
    const SparseMatrix& G() const noexcept { return *G_; }
    const SparseMatrix& D() const noexcept { return *D_; }
    const SparseMatrix& C() const noexcept { return *C_; }

    std::size_t velocity_size() const noexcept { return Auu_->n_rows(); }
    std::size_t pressure_size() const noexcept { return C_->n_rows(); }

    bool is_valid() const noexcept {
        return Auu_ != nullptr && G_ != nullptr && D_ != nullptr && C_ != nullptr &&
               Auu_->n_rows() == Auu_->n_cols() &&
               G_->n_rows() == Auu_->n_rows() &&
               D_->n_cols() == Auu_->n_cols() &&
               G_->n_cols() == C_->n_rows() &&
               D_->n_rows() == C_->n_rows() &&
               C_->n_rows() == C_->n_cols();
    }

    void validate() const {
        if (!is_valid())
            throw std::invalid_argument("BlockOperator: inconsistent saddle-point block dimensions");
    }

private:
    const SparseMatrix* Auu_;
    const SparseMatrix* G_;
    const SparseMatrix* D_;
    const SparseMatrix* C_;
};

} // namespace cfdx::core
