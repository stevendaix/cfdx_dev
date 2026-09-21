#pragma once

#include "cfdx/core/linalg/vector.h"
#include <cstddef>
#include <functional>
#include <stdexcept>

namespace cfdx {
namespace core {

/** Matrix-free linear operator y = A(x), with no stored matrix. */
class MatrixFreeOperator {
public:
    using ApplyFunction = std::function<void(const Vector&, Vector&)>;

    MatrixFreeOperator() = default;
    explicit MatrixFreeOperator(ApplyFunction fn, std::size_t rows, std::size_t cols)
        : apply_(std::move(fn)), rows_(rows), cols_(cols) {}

    std::size_t n_rows() const noexcept { return rows_; }
    std::size_t n_cols() const noexcept { return cols_; }

    void apply(const Vector& x, Vector& y) const {
        if (!apply_) throw std::runtime_error("MatrixFreeOperator: no operator");
        if (x.size() != cols_) throw std::runtime_error("MatrixFreeOperator: input dimension mismatch");
        y.resize(rows_);
        apply_(x, y);
    }

    bool valid() const noexcept { return static_cast<bool>(apply_); }

private:
    ApplyFunction apply_;
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
};

} // namespace core
} // namespace cfdx
