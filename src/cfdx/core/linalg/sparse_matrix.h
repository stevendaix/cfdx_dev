// M0.8-T01 — Sparse matrix (CSR)
//
// Spécification CFDX v0.7 §34 :
//   La représentation interne initiale est CSR.
//   L'API ne doit pas exposer directement cette représentation.
//
// CSR (Compressed Sparse Row) :
//   values[]      : valeurs non nulles, une ligne après l'autre
//   columns[]     : indice de colonne de chaque valeur
//   row_offsets[i] : début de la ligne i dans values[]
//   row_offsets[n_rows+1] = nnz

#pragma once

#include "cfdx/core/linalg/vector.h"
#include <vector>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <limits>

namespace cfdx {
namespace core {

class SparseMatrix {
public:
    using Value = double;
    using Index = std::uint32_t;  // Local index for CSR columns

    SparseMatrix() = default;

    explicit SparseMatrix(std::size_t n_rows, std::size_t n_cols)
        : n_rows_(n_rows), n_cols_(n_cols),
          row_offsets_(n_rows + 1, 0),
          finalized_(false) {
        if (n_rows > std::numeric_limits<Index>::max() ||
            n_cols > std::numeric_limits<Index>::max()) {
            throw std::overflow_error("SparseMatrix: dimensions exceed CSR index range");
        }
    }

    // --- Dimensions ---

    std::size_t n_rows() const noexcept { return n_rows_; }
    std::size_t n_cols() const noexcept { return n_cols_; }
    std::size_t nnz() const noexcept { return finalized_ ? values_.size() : pending_values_.size(); }

    // --- Construction ---

    // Add an entry (row, col, value). Assembly order is unrestricted;
    // finalize() canonicalizes the staged entries into CSR.
    void push_back(std::size_t row, std::size_t col, Value value) {
        if (finalized_) {
            throw std::logic_error("SparseMatrix: push_back after finalize");
        }
        if (row >= n_rows_) {
            throw std::out_of_range("SparseMatrix: row out of range");
        }
        if (col >= n_cols_) {
            throw std::out_of_range("SparseMatrix: column out of range");
        }
        if (!std::isfinite(value)) {
            throw std::invalid_argument("SparseMatrix: non-finite value");
        }
        pending_rows_.push_back(row);
        pending_columns_.push_back(static_cast<Index>(col));
        pending_values_.push_back(value);
    }

    // Finalize unordered COO contributions into UNIQUE, sorted CSR entries.
    // Assembly may add the same (row,column) more than once (e.g. implicit
    // diffusion + source + boundary flux). Before this consolidation, Jacobi
    // and other smoothers that read the first diagonal saw a different
    // operator from matvec(), which sums ALL duplicate entries. The resulting
    // preconditioner was inconsistent with A and could stall Krylov iteration.
    // Build the new CSR transactionally: an overflowing coefficient must not
    // publish a partially finalized matrix or silently become NaN/Inf.
    void finalize() {
        if (finalized_) return;
        const std::size_t count = pending_values_.size();
        if (count > std::numeric_limits<Index>::max()) {
            throw std::overflow_error("SparseMatrix: nnz exceeds CSR index range");
        }
        std::vector<std::size_t> order(count);
        for (std::size_t k = 0; k < count; ++k) order[k] = k;
        std::stable_sort(order.begin(), order.end(), [this](std::size_t a, std::size_t b) {
            if (pending_rows_[a] != pending_rows_[b]) return pending_rows_[a] < pending_rows_[b];
            return pending_columns_[a] < pending_columns_[b];
        });
        std::vector<Value> new_values;
        std::vector<Index> new_columns;
        std::vector<Index> new_offsets(n_rows_ + 1, 0);
        new_values.reserve(count);
        new_columns.reserve(count);
        std::size_t prior_row = n_rows_;
        Index prior_column = 0;
        long double accumulated = 0.0L;
        const auto flush = [&] {
            if (prior_row == n_rows_) return;
            if (!std::isfinite(accumulated) ||
                std::abs(accumulated) > std::numeric_limits<Value>::max()) {
                throw std::overflow_error("SparseMatrix: duplicate entry sum is not finite");
            }
            new_values.push_back(static_cast<Value>(accumulated));
            new_columns.push_back(prior_column);
            ++new_offsets[prior_row + 1];
        };
        for (const auto k : order) {
            const auto row = pending_rows_[k];
            const auto column = pending_columns_[k];
            if (row != prior_row || column != prior_column) {
                flush();
                prior_row = row;
                prior_column = column;
                accumulated = 0.0L;
            }
            accumulated += static_cast<long double>(pending_values_[k]);
        }
        flush();
        for (std::size_t i = 0; i < n_rows_; ++i)
            new_offsets[i + 1] += new_offsets[i];
        values_.swap(new_values);
        columns_.swap(new_columns);
        row_offsets_.swap(new_offsets);
        pending_rows_.clear();
        pending_columns_.clear();
        pending_values_.clear();
        finalized_ = true;
    }

    // Assemble already sorted and unique CSR without staging COO entries or
    // sorting them again. Used by the symbolic/numeric Galerkin product.
    // The factory rejects malformed graphs instead of silently merging them.
    static SparseMatrix from_sorted_csr(
        std::size_t n_rows, std::size_t n_cols,
        std::vector<Index> offsets, std::vector<Index> columns,
        std::vector<Value> values) {
        SparseMatrix result(n_rows, n_cols);
        if (offsets.size() != n_rows + 1 || offsets.empty() ||
            offsets.front() != 0 || offsets.back() != columns.size() ||
            columns.size() != values.size() ||
            columns.size() > std::numeric_limits<Index>::max()) {
            throw std::invalid_argument("SparseMatrix::from_sorted_csr: invalid dimensions");
        }
        for (std::size_t row = 0; row < n_rows; ++row) {
            if (offsets[row] > offsets[row + 1] ||
                offsets[row + 1] > columns.size())
                throw std::invalid_argument("SparseMatrix::from_sorted_csr: invalid offsets");
            for (std::size_t k = offsets[row]; k < offsets[row + 1]; ++k) {
                if (columns[k] >= n_cols || !std::isfinite(values[k]) ||
                    (k > offsets[row] && columns[k - 1] >= columns[k])) {
                    throw std::invalid_argument("SparseMatrix::from_sorted_csr: invalid row");
                }
            }
        }
        result.row_offsets_ = std::move(offsets);
        result.columns_ = std::move(columns);
        result.values_ = std::move(values);
        result.finalized_ = true;
        return result;
    }

    // --- Accès ---

    // Accès à l'entrée (row, col). Retourne 0.0 si absente.
    Value operator()(std::size_t row, std::size_t col) const {
        check_row(row);
        const Index start = row_offsets_[row];
        const Index end = row_offsets_[row + 1];
        Value sum = 0.0;
        for (Index k = start; k < end; ++k) {
            if (columns_[k] == static_cast<Index>(col)) sum += values_[k];
        }
        return sum;
    }

    // --- Opérations ---

    // Produit matrice-vecteur : y = A * x (API std::vector).
    std::vector<Value> matvec(const std::vector<Value>& x) const {
        if (x.size() != n_cols_) {
            throw std::runtime_error("SparseMatrix::matvec: dimension mismatch");
        }
        std::vector<Value> y(n_rows_, 0.0);
        for (std::size_t i = 0; i < n_rows_; ++i) {
            const Index start = row_offsets_[i];
            const Index end = row_offsets_[i + 1];
            Value sum = 0.0;
            for (Index k = start; k < end; ++k) {
                sum += values_[k] * x[columns_[k]];
            }
            y[i] = sum;
        }
        return y;
    }

    // Produit matrice-vecteur : y = A^T * x (API std::vector).
    std::vector<Value> matvec_transpose(const std::vector<Value>& x) const {
        if (x.size() != n_rows_) {
            throw std::runtime_error("SparseMatrix::matvec_transpose: dimension mismatch");
        }
        std::vector<Value> y(n_cols_, 0.0);
        for (std::size_t i = 0; i < n_rows_; ++i) {
            const Index start = row_offsets_[i];
            const Index end = row_offsets_[i + 1];
            const Value xi = x[i];
            for (Index k = start; k < end; ++k) {
                y[columns_[k]] += values_[k] * xi;
            }
        }
        return y;
    }

    // Surcharge Vector (API principale).
    std::vector<Value> matvec(const Vector& x) const {
        return matvec(std::vector<Value>(x.data(), x.data() + x.size()));
    }
    std::vector<Value> matvec_transpose(const Vector& x) const {
        return matvec_transpose(std::vector<Value>(x.data(), x.data() + x.size()));
    }

    // Allocation-free matrix-vector product for solver/operator callbacks.
    // The output must be pre-sized; aliasing is supported for square matrices.
    void matvec_into(const Vector& x, Vector& y) const {
        if (x.size() != n_cols_ || y.size() != n_rows_) {
            throw std::runtime_error("SparseMatrix::matvec_into: dimension mismatch");
        }
        if (&x == &y) {
            const Vector input = x;
            matvec_into(input, y);
            return;
        }
        for (std::size_t i = 0; i < n_rows_; ++i) {
            Value sum = 0.0;
            for (Index k = row_offsets_[i]; k < row_offsets_[i + 1]; ++k) {
                sum += values_[k] * x(columns_[k]);
            }
            y(i) = sum;
        }
    }

    // --- Accès bulk ---

    const Value* values_data() const noexcept { return values_.data(); }
    const std::uint32_t* columns_data() const noexcept { return columns_.data(); }
    const std::uint32_t* row_offsets_data() const noexcept { return row_offsets_.data(); }

    Value* values_data() noexcept { return values_.data(); }
    std::uint32_t* columns_data() noexcept { return columns_.data(); }
    std::uint32_t* row_offsets_data() noexcept { return row_offsets_.data(); }

    // --- Validation ---

    bool is_consistent() const noexcept {
        if (!finalized_) {
            return pending_rows_.empty() && pending_columns_.empty() && pending_values_.empty();
        }
        if (row_offsets_.size() != n_rows_ + 1) return false;
        if (row_offsets_[0] != 0) return false;
        if (row_offsets_[n_rows_] != static_cast<Index>(values_.size())) return false;
        for (std::size_t i = 1; i <= n_rows_; ++i) {
            if (row_offsets_[i] < row_offsets_[i - 1]) return false;
        }
        for (std::size_t i = 0; i < n_rows_; ++i) {
            const Index begin = row_offsets_[i];
            const Index end = row_offsets_[i + 1];
            for (Index k = begin; k < end; ++k) {
                if (columns_[k] >= static_cast<Index>(n_cols_)) return false;
                if (!std::isfinite(values_[k])) return false;
                if (k > begin && columns_[k] <= columns_[k - 1]) return false;
            }
        }
        return true;
    }

    void clear() {
        values_.clear();
        columns_.clear();
        pending_rows_.clear();
        pending_columns_.clear();
        pending_values_.clear();
        row_offsets_.assign(n_rows_ + 1, 0);
        finalized_ = false;
    }

private:
    void check_row(std::size_t i) const {
        if (i >= n_rows_) {
            throw std::out_of_range("SparseMatrix: row out of range");
        }
    }

    std::size_t n_rows_ = 0;
    std::size_t n_cols_ = 0;
    std::vector<Value> values_;
    std::vector<std::uint32_t> columns_;
    std::vector<std::uint32_t> row_offsets_;
    std::vector<std::size_t> pending_rows_;
    std::vector<std::uint32_t> pending_columns_;
    std::vector<Value> pending_values_;
    bool finalized_ = true;
};

}  // namespace core
}  // namespace cfdx
