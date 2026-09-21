#pragma once

#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"

#include <cstddef>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::core {

class MatrixFreeOperator {
public:
    using ApplyFunction = std::function<void(const Vector&, Vector&)>;

    MatrixFreeOperator() = default;

    MatrixFreeOperator(std::size_t rows, std::size_t cols, ApplyFunction apply)
        : rows_(rows), cols_(cols), apply_(std::move(apply))
    {
        if (!apply_) throw std::invalid_argument("MatrixFreeOperator: apply function is empty");
    }

    std::size_t n_rows() const noexcept { return rows_; }
    std::size_t n_cols() const noexcept { return cols_; }

    void apply(const Vector& x, Vector& y) const
    {
        if (x.size() != cols_ || y.size() != rows_)
            throw std::invalid_argument("MatrixFreeOperator: vector dimension mismatch");
        apply_(x, y);
    }

    Vector operator()(const Vector& x) const
    {
        if (x.size() != cols_)
            throw std::invalid_argument("MatrixFreeOperator: vector dimension mismatch");
        Vector y(rows_, 0.0);
        apply(x, y);
        return y;
    }

private:
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    ApplyFunction apply_;
};

inline MatrixFreeOperator matrix_free_from_sparse(const SparseMatrix& A)
{
    if (!A.is_consistent())
        throw std::invalid_argument("matrix_free_from_sparse: inconsistent matrix");

    return MatrixFreeOperator(
        A.n_rows(), A.n_cols(),
        [&A](const Vector& x, Vector& y) {
            const auto result = A.matvec(x);
            for (std::size_t i = 0; i < result.size(); ++i) y(i) = result[i];
        });
}

} // namespace cfdx::core
