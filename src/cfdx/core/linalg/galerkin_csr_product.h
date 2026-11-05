#pragma once

#include "cfdx/core/linalg/sparse_matrix.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::core {

// Dependency-free, two-phase Galerkin product C = P^T A P.
//
// symbolic(): construct the AP and P^T(AP) graphs exactly once.
// numeric():  reuse those graphs with marker-indexed sparse accumulators,
//             without per-entry allocation, map lookups or COO sorting.
//
// The input A pattern and the ordering of transfer entries are part of the
// symbolic contract. Near-zero coarse entries are pruned as in the legacy
// Galerkin implementation. If pruning changes a subsequent level's graph,
// numeric() must reject that level and the AMG owner must rebuild completely.
// This is deliberately NOT a drop-tolerance/non-Galerkin approximation.
class GalerkinCSRProduct final {
public:
    using TransferRow = std::vector<std::pair<std::size_t, double>>;
    using TransferRows = std::vector<TransferRow>;
    using Index = SparseMatrix::Index;

    bool symbolic(const SparseMatrix& A, const TransferRows& P,
                  std::size_t coarse_size) {
        if (!A.is_consistent() || A.n_rows() == 0 ||
            A.n_rows() != A.n_cols() || P.size() != A.n_rows() ||
            coarse_size == 0 ||
            coarse_size > std::numeric_limits<Index>::max())
            return false;

        GalerkinCSRProduct candidate;
        candidate.n_fine_ = A.n_rows();
        candidate.n_coarse_ = coarse_size;
        candidate.input_offsets_.assign(A.row_offsets_data(),
                                        A.row_offsets_data() + A.n_rows() + 1);
        candidate.input_columns_.assign(A.columns_data(),
                                        A.columns_data() + A.nnz());
        // Build P^T as a contiguous CSR adjacency rather than allocating one
        // std::vector per coarse row.  This reduces allocation pressure and
        // improves locality on the large NACA AMG levels.  Preserve the
        // original stable (fine row, entry) order for numeric equivalence.
        std::vector<Index> transpose_counts(coarse_size, 0);
        for (std::size_t i = 0; i < P.size(); ++i) {
            for (std::size_t j = 0; j < P[i].size(); ++j) {
                const auto [ci, weight] = P[i][j];
                if (ci >= coarse_size || !std::isfinite(weight) ||
                    j > std::numeric_limits<Index>::max() ||
                    transpose_counts[ci] == std::numeric_limits<Index>::max())
                    return false;
                ++transpose_counts[ci];
            }
        }
        candidate.transpose_offsets_.resize(coarse_size + 1, 0);
        for (std::size_t ci = 0; ci < coarse_size; ++ci) {
            const std::uint64_t next =
                static_cast<std::uint64_t>(candidate.transpose_offsets_[ci]) +
                transpose_counts[ci];
            if (next > std::numeric_limits<Index>::max()) return false;
            candidate.transpose_offsets_[ci + 1] = static_cast<Index>(next);
        }
        const std::size_t trans_nnz = candidate.transpose_offsets_.back();
        candidate.transpose_fine_rows_.resize(trans_nnz);
        candidate.transpose_p_positions_.resize(trans_nnz);
        auto cursor = candidate.transpose_offsets_;
        for (std::size_t i = 0; i < P.size(); ++i) {
            for (std::size_t j = 0; j < P[i].size(); ++j) {
                const std::size_t ci = P[i][j].first;
                const std::size_t position = cursor[ci]++;
                candidate.transpose_fine_rows_[position] = static_cast<Index>(i);
                candidate.transpose_p_positions_[position] = static_cast<Index>(j);
            }
        }

        // Symbolic sparse A*P. A row only needs the coarse columns in the
        // prolongation rows of its nonzero fine neighbours.
        candidate.ap_offsets_.reserve(A.n_rows() + 1);
        candidate.ap_offsets_.push_back(0);
        std::vector<Index> scratch;
        // Sparse row marker prevents the O(stencil fan-out) duplicates that
        // otherwise dominate symbolic setup for smoothed aggregation.
        std::vector<std::size_t> seen(coarse_size,
                                      std::numeric_limits<std::size_t>::max());
        const auto* a_row = A.row_offsets_data();
        const auto* a_col = A.columns_data();
        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            scratch.clear();
            for (std::size_t k = a_row[i]; k < a_row[i + 1]; ++k)
                for (const auto& [cj, weight] : P[a_col[k]]) {
                    (void)weight;
                    if (seen[cj] != i) {
                        seen[cj] = i;
                        scratch.push_back(static_cast<Index>(cj));
                    }
                }
            std::sort(scratch.begin(), scratch.end());
            if (scratch.size() > std::numeric_limits<Index>::max() -
                                 candidate.ap_columns_.size()) return false;
            candidate.ap_columns_.insert(candidate.ap_columns_.end(),
                                         scratch.begin(), scratch.end());
            candidate.ap_offsets_.push_back(
                static_cast<Index>(candidate.ap_columns_.size()));
        }

        std::fill(seen.begin(), seen.end(),
                  std::numeric_limits<std::size_t>::max());
        // Symbolic P^T(AP). Each row ci collects the AP columns of the fine
        // points to which interpolation column ci connects.
        candidate.coarse_offsets_.reserve(coarse_size + 1);
        candidate.coarse_offsets_.push_back(0);
        for (std::size_t ci = 0; ci < coarse_size; ++ci) {
            scratch.clear();
            for (std::size_t entry = candidate.transpose_offsets_[ci];
                 entry < candidate.transpose_offsets_[ci + 1]; ++entry) {
                const std::size_t i = candidate.transpose_fine_rows_[entry];
                for (std::size_t p = candidate.ap_offsets_[i];
                     p < candidate.ap_offsets_[i + 1]; ++p) {
                    const auto cj = candidate.ap_columns_[p];
                    if (seen[cj] != ci) {
                        seen[cj] = ci;
                        scratch.push_back(cj);
                    }
                }
            }
            std::sort(scratch.begin(), scratch.end());
            if (scratch.size() > std::numeric_limits<Index>::max() -
                                 candidate.coarse_columns_.size()) return false;
            candidate.coarse_columns_.insert(candidate.coarse_columns_.end(),
                                             scratch.begin(), scratch.end());
            candidate.coarse_offsets_.push_back(
                static_cast<Index>(candidate.coarse_columns_.size()));
        }
        *this = std::move(candidate);
        return true;
    }

    // Assemble a new CSR coarse matrix during full AMG setup.  Historic
    // pruning of tiny entries is preserved (no non-Galerkin dropping).
    bool numeric(const SparseMatrix& A, const TransferRows& P,
                 SparseMatrix& out, std::size_t* dropped_entries = nullptr) const {
        std::vector<double> c_values;
        if (!accumulate(A, P, c_values)) return false;
        std::vector<Index> offsets;
        std::vector<Index> columns;
        std::vector<double> values;
        offsets.reserve(n_coarse_ + 1);
        columns.reserve(coarse_columns_.size());
        values.reserve(coarse_columns_.size());
        offsets.push_back(0);
        std::size_t dropped = 0;
        for (std::size_t ci = 0; ci < n_coarse_; ++ci) {
            for (std::size_t p = coarse_offsets_[ci]; p < coarse_offsets_[ci + 1]; ++p) {
                if (std::abs(c_values[p]) > 1e-30) {
                    columns.push_back(coarse_columns_[p]);
                    values.push_back(c_values[p]);
                } else {
                    ++dropped;
                }
            }
            offsets.push_back(static_cast<Index>(columns.size()));
        }
        out = SparseMatrix::from_sorted_csr(
            n_coarse_, n_coarse_, std::move(offsets),
            std::move(columns), std::move(values));
        if (dropped_entries) *dropped_entries += dropped;
        return true;
    }

    // Faster numeric-only AMG refresh. `out` is a staged copy of the old
    // coarse matrix. Require exactly the same pruned CSR graph, even when
    // coefficients cross the 1e-30 pruning threshold. Validate all values
    // and all rows BEFORE writing any coefficient: on failure, `out` remains
    // byte-for-byte unchanged and the owner can safely trigger full setup.
    bool numeric_into_existing(const SparseMatrix& A, const TransferRows& P,
                               SparseMatrix& out,
                               std::size_t* dropped_entries = nullptr) const {
        if (&A == &out || out.n_rows() != n_coarse_ ||
            out.n_cols() != n_coarse_ || !out.is_consistent()) return false;
        std::vector<double> c_values;
        if (!accumulate(A, P, c_values)) return false;
        const auto* old_offsets = out.row_offsets_data();
        const auto* old_columns = out.columns_data();
        std::size_t dropped = 0;
        for (std::size_t ci = 0; ci < n_coarse_; ++ci) {
            std::size_t cursor = old_offsets[ci];
            for (std::size_t p = coarse_offsets_[ci]; p < coarse_offsets_[ci + 1]; ++p) {
                if (std::abs(c_values[p]) > 1e-30) {
                    if (cursor == old_offsets[ci + 1] ||
                        old_columns[cursor] != coarse_columns_[p]) return false;
                    ++cursor;
                } else {
                    ++dropped;
                }
            }
            if (cursor != old_offsets[ci + 1]) return false;
        }
        std::size_t out_pos = 0;
        auto* out_values = out.values_data();
        for (const double value : c_values)
            if (std::abs(value) > 1e-30) out_values[out_pos++] = value;
        if (dropped_entries) *dropped_entries += dropped;
        return true;
    }

    std::size_t symbolic_ap_nnz() const noexcept { return ap_columns_.size(); }
    std::size_t symbolic_coarse_nnz() const noexcept { return coarse_columns_.size(); }

private:
    bool accumulate(const SparseMatrix& A, const TransferRows& P,
                    std::vector<double>& c_values) const {
        if (A.n_rows() != n_fine_ || A.n_cols() != n_fine_ ||
            P.size() != n_fine_ || !A.is_consistent() ||
            A.nnz() != input_columns_.size() || n_coarse_ == 0 ||
            !std::equal(input_offsets_.begin(), input_offsets_.end(),
                        A.row_offsets_data()) ||
            !std::equal(input_columns_.begin(), input_columns_.end(),
                        A.columns_data())) return false;

        constexpr std::size_t missing = std::numeric_limits<std::size_t>::max();
        std::vector<std::size_t> lookup(n_coarse_, missing);
        std::vector<double> ap_values(ap_columns_.size(), 0.0);
        const auto* a_row = A.row_offsets_data();
        const auto* a_col = A.columns_data();
        const auto* a_val = A.values_data();
        for (std::size_t i = 0; i < n_fine_; ++i) {
            for (std::size_t p = ap_offsets_[i]; p < ap_offsets_[i + 1]; ++p)
                lookup[ap_columns_[p]] = p;
            for (std::size_t k = a_row[i]; k < a_row[i + 1]; ++k) {
                for (const auto& [cj, weight] : P[a_col[k]]) {
                    if (cj >= n_coarse_ || !std::isfinite(weight) ||
                        lookup[cj] == missing) return false;
                    const std::size_t pos = lookup[cj];
                    ap_values[pos] += a_val[k] * weight;
                    if (!std::isfinite(ap_values[pos])) return false;
                }
            }
            for (std::size_t p = ap_offsets_[i]; p < ap_offsets_[i + 1]; ++p)
                lookup[ap_columns_[p]] = missing;
        }

        c_values.assign(coarse_columns_.size(), 0.0);
        for (std::size_t ci = 0; ci < n_coarse_; ++ci) {
            for (std::size_t p = coarse_offsets_[ci]; p < coarse_offsets_[ci + 1]; ++p)
                lookup[coarse_columns_[p]] = p;
            for (std::size_t t = transpose_offsets_[ci];
                 t < transpose_offsets_[ci + 1]; ++t) {
                const std::size_t i = transpose_fine_rows_[t];
                const std::size_t entry = transpose_p_positions_[t];
                if (entry >= P[i].size() || P[i][entry].first != ci ||
                    !std::isfinite(P[i][entry].second)) return false;
                const double weight = P[i][entry].second;
                for (std::size_t p = ap_offsets_[i]; p < ap_offsets_[i + 1]; ++p) {
                    const std::size_t pos = lookup[ap_columns_[p]];
                    if (pos == missing) return false;
                    c_values[pos] += weight * ap_values[p];
                    if (!std::isfinite(c_values[pos])) return false;
                }
            }
            for (std::size_t p = coarse_offsets_[ci]; p < coarse_offsets_[ci + 1]; ++p)
                lookup[coarse_columns_[p]] = missing;
        }
        return true;
    }

    std::size_t n_fine_ = 0;
    std::size_t n_coarse_ = 0;
    std::vector<Index> input_offsets_;
    std::vector<Index> input_columns_;
    std::vector<Index> ap_offsets_;
    std::vector<Index> ap_columns_;
    std::vector<Index> coarse_offsets_;
    std::vector<Index> coarse_columns_;
    // CSR P^T adjacency: for each coarse column, fine rows and transfer
    // positions in stable input order.  No per-column heap allocations.
    std::vector<Index> transpose_offsets_;
    std::vector<Index> transpose_fine_rows_;
    std::vector<Index> transpose_p_positions_;
};

} // namespace cfdx::core
