// M0.8-T03 — Linear system
//
// Spécification CFDX v0.7 §33 :
//   LinearSystem(A, b, x)
//
//   L'API ne doit pas exposer la représentation CSR directement.

#pragma once

#include "sparse_matrix.h"
#include "vector.h"
#include <cstddef>
#include <stdexcept>

namespace cfdx {
namespace core {

class LinearSystem {
public:
    LinearSystem(const SparseMatrix& A, const Vector& b, Vector& x)
        : A_(A), b_(b), x_(x) {}

    const SparseMatrix& A() const noexcept { return A_; }
    const Vector& b() const noexcept { return b_; }
    Vector& x() noexcept { return x_; }
    const Vector& x() const noexcept { return x_; }

    std::size_t n_rows() const noexcept { return A_.n_rows(); }
    std::size_t n_cols() const noexcept { return A_.n_cols(); }

    // Vérifie la cohérence des dimensions.
    bool is_consistent() const noexcept {
        return b_.size() == A_.n_rows() && x_.size() == A_.n_cols();
    }

private:
    const SparseMatrix& A_;
    const Vector& b_;
    Vector& x_;
};

}  // namespace core
}  // namespace cfdx