#pragma once

#include "cfdx/core/linalg/sparse_matrix.h"

#include <cstddef>
#include <cstdint>

namespace cfdx::core {

// Extract one block from the fixed production coupled ordering
// [Ux, Uy, Uz, p] without changing the assembled CSR operator.
//
// row_block/col_block are 0 for the 3N velocity field and 1 for pressure.
// The returned matrix owns its CSR storage; callers may therefore retain the
// extracted blocks while the source matrix is updated or destroyed.
inline SparseMatrix extract_coupled_block(const SparseMatrix& A,
                                          std::size_t row_block,
                                          std::size_t col_block,
                                          std::size_t row_size,
                                          std::size_t col_size) {
    if (A.n_rows() != A.n_cols() || A.n_rows() == 0 ||
        A.n_rows() % 4 != 0 || row_block > 1 || col_block > 1)
        return SparseMatrix();

    const std::size_t velocity_size = A.n_rows() * 3 / 4;
    const std::size_t roffset = row_block == 0 ? 0 : velocity_size;
    const std::size_t coffset = col_block == 0 ? 0 : velocity_size;

    if (roffset + row_size > A.n_rows() ||
        coffset + col_size > A.n_cols())
        return SparseMatrix();

    SparseMatrix block(row_size, col_size);
    const auto* rows = A.row_offsets_data();
    const auto* cols = A.columns_data();
    const auto* values = A.values_data();

    for (std::size_t r = 0; r < row_size; ++r) {
        const std::size_t gr = roffset + r;
        for (std::uint32_t k = rows[gr]; k < rows[gr + 1]; ++k) {
            const std::size_t gc = cols[k];
            if (gc < coffset || gc >= coffset + col_size)
                continue;
            block.push_back(r, gc - coffset, values[k]);
        }
    }
    block.finalize();
    return block;
}

} // namespace cfdx::core
