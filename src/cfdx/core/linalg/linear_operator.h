#pragma once
#include "cfdx/core/linalg/vector.h"
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <utility>

namespace cfdx::core {

class LinearOperatorBase {
public:
    virtual ~LinearOperatorBase() = default;
    virtual std::size_t rows() const noexcept = 0;
    virtual std::size_t cols() const noexcept = 0;
    virtual void apply(const Vector& x, Vector& y) const = 0;

    // Optional diagonal access used by matrix-free smoothers/preconditioners.
    virtual bool has_diagonal() const noexcept { return false; }
    virtual void diagonal(Vector&) const {
        throw std::logic_error("linear operator does not expose a diagonal");
    }
};

class FunctionalLinearOperator final : public LinearOperatorBase {
public:
    using Apply = std::function<void(const Vector&, Vector&)>;
    FunctionalLinearOperator(std::size_t n, Apply apply) : n_(n), apply_(std::move(apply)) {
        if (n_ == 0 || !apply_) throw std::invalid_argument("invalid linear operator");
    }
    std::size_t rows() const noexcept override { return n_; }
    std::size_t cols() const noexcept override { return n_; }
    void apply(const Vector& x, Vector& y) const override {
        if (x.size() != n_) throw std::invalid_argument("LinearOperator: input size mismatch");
        if (y.size() != n_) y.resize(n_, 0.0);
        apply_(x, y);
    }
private:
    std::size_t n_;
    Apply apply_;
};

template<class Operator>
class MatrixFreeOperator final : public LinearOperatorBase {
public:
    MatrixFreeOperator(std::size_t rows, std::size_t cols, const Operator& op)
        : rows_(rows), cols_(cols), op_(op) {}
    std::size_t rows() const noexcept override { return rows_; }
    std::size_t cols() const noexcept override { return cols_; }
    void apply(const Vector& x, Vector& y) const override { op_.apply(x, y); }
private:
    std::size_t rows_, cols_;
    const Operator& op_;
};

} // namespace cfdx::core
