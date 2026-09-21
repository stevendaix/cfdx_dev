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

namespace cfdx {
namespace core {

class SparseMatrix {
public:
    using Value = double;
    using Index = std::uint32_t;  // Local index for CSR columns

    SparseMatrix() = default;

    explicit SparseMatrix(std::size_t n_rows, std::size_t n_cols)
        : n_rows_(n_rows), n_cols_(n_cols),
          row_offsets_(n_rows + 1, 0) {}

    // --- Dimensions ---

    std::size_t n_rows() const noexcept { return n_rows_; }
    std::size_t n_cols() const noexcept { return n_cols_; }
    std::size_t nnz() const noexcept { return values_.size(); }

    // --- Construction ---

    // Ajoute une entrée (row, col, value). Doit être appelé par ligne,
    // dans l'ordre des colonnes croissant.
    void push_back(std::size_t row, std::size_t col, Value value) {
        if (row >= n_rows_) {
            throw std::out_of_range("SparseMatrix: row out of range");
        }
        if (col >= n_cols_) {
            throw std::out_of_range("SparseMatrix: column out of range");
        }
        values_.push_back(value);
        columns_.push_back(static_cast<Index>(col));
        // row_offsets[row+1] sera mis à jour à la fin.
        // On utilise un compteur séparé pour les entrées par ligne.
        row_offsets_[row + 1]++;
    }

    // Finalise la structure CSR après l'ajout des entrées.
    // row_offsets devient un préfixe sum.
    void finalize() {
        for (std::size_t i = 0; i < n_rows_; ++i) {
            row_offsets_[i + 1] += row_offsets_[i];
        }
        finalized_ = true;
    }

    // --- Accès ---

    // Accès à l'entrée (row, col). Retourne 0.0 si absente.
    Value operator()(std::size_t row, std::size_t col) const {
        check_row(row);
        const Index start = row_offsets_[row];
        const Index end = row_offsets_[row + 1];
        // Les colonnes sont supposées triées.
        for (Index k = start; k < end; ++k) {
            if (columns_[k] == static_cast<Index>(col)) {
                return values_[k];
            }
        }
        return 0.0;
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

    // Return the diagonal in O(nnz) time. Missing diagonal entries are zero.
    // This is intentionally extracted on demand so matrix construction remains
    // lightweight while iterative solvers avoid repeated per-row searches.
    std::vector<Value> diagonal() const {
        std::vector<Value> diag(n_rows_, 0.0);
        for (std::size_t i = 0; i < n_rows_; ++i) {
            const Index start = row_offsets_[i];
            const Index end = row_offsets_[i + 1];
            for (Index k = start; k < end; ++k) {
                if (columns_[k] == static_cast<Index>(i)) {
                    diag[i] = values_[k];
                    break;
                }
            }
        }
        return diag;
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
        if (row_offsets_.size() != n_rows_ + 1) return false;
        if (row_offsets_[0] != 0) return false;
        if (row_offsets_[n_rows_] != static_cast<Index>(values_.size())) return false;
        for (std::size_t i = 1; i <= n_rows_; ++i) {
            if (row_offsets_[i] < row_offsets_[i - 1]) return false;
        }
        for (const auto c : columns_) {
            if (c >= static_cast<Index>(n_cols_)) return false;
        }
        return true;
    }

    void clear() {
        values_.clear();
        columns_.clear();
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
    bool finalized_ = true;
};

}  // namespace core
}  // namespace cfdx
