#pragma once

#include "cfdx/core/linalg/sparse_matrix.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::core {

// Sparse CSR matrix product C = A B.
//
// The implementation is intentionally independent of the storage internals:
// it uses the public CSR accessors and assembles each result row through a
// sorted map.  This gives deterministic column ordering and combines duplicate
// contributions without a dense intermediate matrix.
//
// This is infrastructure for algebraic Schur operators such as DG, DAuuG
// and pressure-side operators used by LSC/BFBt/PCD.
inline SparseMatrix sparse_matmul(const SparseMatrix& A,
                                  const SparseMatrix& B,
                                  double drop_tolerance = 0.0) {
    if (A.n_cols() != B.n_rows())
        throw std::invalid_argument("sparse_matmul: dimension mismatch");
    if (drop_tolerance < 0.0)
        throw std::invalid_argument("sparse_matmul: negative drop tolerance");

    SparseMatrix C(A.n_rows(), B.n_cols());

    for (std::size_t i = 0; i < A.n_rows(); ++i) {
        std::map<std::uint32_t, double> row;
        const auto a_begin = A.row_offsets_data()[i];
        const auto a_end = A.row_offsets_data()[i + 1];

        for (auto ak = a_begin; ak < a_end; ++ak) {
            const std::size_t k = A.columns_data()[ak];
            const double aik = A.values_data()[ak];
            const auto b_begin = B.row_offsets_data()[k];
            const auto b_end = B.row_offsets_data()[k + 1];

            for (auto bk = b_begin; bk < b_end; ++bk)
                row[B.columns_data()[bk]] += aik * B.values_data()[bk];
        }

        for (const auto& [j, value] : row) {
            if (drop_tolerance == 0.0 || std::abs(value) > drop_tolerance)
                C.push_back(i, j, value);
        }
    }

    C.finalize();
    return C;
}

// Compare only the symbolic CSR graph (dimensions, row offsets and column
// indices), deliberately ignoring numerical values.  Schur/AMG update paths
// can use this predicate to distinguish a value refresh from a symbolic
// rebuild.
inline bool same_sparse_pattern(const SparseMatrix& A,
                                const SparseMatrix& B) noexcept {
    if (A.n_rows() != B.n_rows() ||
        A.n_cols() != B.n_cols() ||
        A.nnz() != B.nnz())
        return false;

    const auto* a_rows = A.row_offsets_data();
    const auto* b_rows = B.row_offsets_data();
    for (std::size_t i = 0; i <= A.n_rows(); ++i) {
        if (a_rows[i] != b_rows[i]) return false;
    }

    const auto* a_cols = A.columns_data();
    const auto* b_cols = B.columns_data();
    for (std::size_t k = 0; k < A.nnz(); ++k) {
        if (a_cols[k] != b_cols[k]) return false;
    }
    return true;
}

} // namespace cfdx::core
