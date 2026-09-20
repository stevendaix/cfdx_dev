// M0.8-T02 — Vector
//
// Spécification CFDX v0.7 §33 :
//   L'algèbre linéaire est indépendante de la physique.
//   Objets : Vector, SparseMatrix, LinearSystem, Preconditioner, LinearSolver.

#pragma once

#include <vector>
#include <cstddef>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace cfdx {
namespace core {

class Vector {
public:
    using Value = double;
    using Index = std::size_t;  // Local index for vector elements

    Vector() = default;

    explicit Vector(std::size_t n, Value init = 0.0)
        : data_(n, init) {}

    // --- Dimensions ---

    std::size_t size() const noexcept { return data_.size(); }
    bool empty() const noexcept { return data_.empty(); }

    // --- Accès ---

    Value operator()(std::size_t i) const {
        check_index(i);
        return data_[i];
    }
    Value& operator()(std::size_t i) {
        check_index(i);
        return data_[i];
    }

    // --- Opérations ---

    void fill(Value v) { std::fill(data_.begin(), data_.end(), v); }
    void resize(std::size_t n) { data_.resize(n, 0.0); }
    void clear() { data_.clear(); }

    // --- Accès bulk ---

    const Value* data() const noexcept { return data_.data(); }
    Value* data() noexcept { return data_.data(); }

    // --- Normes ---

    Value norm1() const {
        Value s = 0.0;
        for (const auto v : data_) s += std::abs(v);
        return s;
    }
    Value norm2() const {
        Value s = 0.0;
        for (const auto v : data_) s += v * v;
        return std::sqrt(s);
    }
    Value norm_inf() const {
        Value s = 0.0;
        for (const auto v : data_) s = std::max(s, std::abs(v));
        return s;
    }

    // --- Opérateurs ---

    Vector& operator*=(Value s) {
        for (auto& v : data_) v *= s;
        return *this;
    }
    Vector operator*(Value s) const {
        Vector r = *this;
        r *= s;
        return r;
    }

    Vector& operator+=(const Vector& o) {
        if (o.size() != size()) {
            throw std::runtime_error("Vector::operator+=: dimension mismatch");
        }
        for (std::size_t i = 0; i < size(); ++i) data_[i] += o.data_[i];
        return *this;
    }
    Vector operator+(const Vector& o) const {
        Vector r = *this;
        r += o;
        return r;
    }

    Vector& operator-=(const Vector& o) {
        if (o.size() != size()) {
            throw std::runtime_error("Vector::operator-=: dimension mismatch");
        }
        for (std::size_t i = 0; i < size(); ++i) data_[i] -= o.data_[i];
        return *this;
    }
    Vector operator-(const Vector& o) const {
        Vector r = *this;
        r -= o;
        return r;
    }

    Value dot(const Vector& o) const {
        if (o.size() != size()) {
            throw std::runtime_error("Vector::dot: dimension mismatch");
        }
        Value s = 0.0;
        for (std::size_t i = 0; i < size(); ++i) s += data_[i] * o.data_[i];
        return s;
    }

    // --- Validation ---

    bool is_valid() const {
        for (const auto v : data_) {
            if (!std::isfinite(v)) return false;
        }
        return true;
    }

private:
    void check_index(std::size_t i) const {
        if (i >= size()) {
            throw std::out_of_range("Vector: index out of range");
        }
    }

    std::vector<Value> data_;
};

}  // namespace core
}  // namespace cfdx
