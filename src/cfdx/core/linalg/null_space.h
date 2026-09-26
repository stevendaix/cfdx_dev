#pragma once

#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::core {

/// Orthonormal basis and projector for known null modes of a linear operator.
/// This mirrors the useful core of PETSc MatNullSpace without introducing a
/// PETSc dependency. Incompatible right-hand sides are detected, never fixed
/// silently.
class NullSpaceProjector {
public:
    explicit NullSpaceProjector(std::vector<Vector> basis)
        : dimension_(basis.empty() ? 0 : basis.front().size()) {
        if (basis.empty() || dimension_ == 0)
            throw std::invalid_argument("null space basis must not be empty");

        basis_.reserve(basis.size());
        for (auto vector : basis) {
            if (vector.size() != dimension_)
                throw std::invalid_argument("null space basis dimension mismatch");
            if (!vector.is_valid())
                throw std::invalid_argument("null space basis contains non-finite values");
            const double input_norm = stable_norm(vector);
            if (!(input_norm > 0.0))
                throw std::invalid_argument("null space basis contains a zero vector");

            // Two modified Gram-Schmidt passes limit loss of orthogonality for
            // user-provided near-null-space vectors.
            for (int pass = 0; pass < 2; ++pass) {
                for (const auto& accepted : basis_) {
                    const double coefficient = vector.dot(accepted);
                    for (std::size_t i = 0; i < dimension_; ++i)
                        vector(i) -= coefficient * accepted(i);
                }
            }
            const double norm = stable_norm(vector);
            const double threshold = 64.0 * std::numeric_limits<double>::epsilon() *
                                     std::sqrt(static_cast<double>(dimension_)) *
                                     input_norm;
            if (!(norm > threshold) || !std::isfinite(norm))
                throw std::invalid_argument("null space basis is linearly dependent");
            vector *= 1.0 / norm;
            basis_.push_back(std::move(vector));
        }
    }

    static NullSpaceProjector constant(std::size_t dimension) {
        if (dimension == 0)
            throw std::invalid_argument("constant null space requires a positive dimension");
        return NullSpaceProjector({Vector(dimension, 1.0)});
    }

    std::size_t dimension() const noexcept { return dimension_; }
    std::size_t basis_size() const noexcept { return basis_.size(); }
    const std::vector<Vector>& basis() const noexcept { return basis_; }

    void remove(Vector& vector) const {
        check_vector(vector);
        for (const auto& mode : basis_) {
            const double coefficient = vector.dot(mode);
            for (std::size_t i = 0; i < dimension_; ++i)
                vector(i) -= coefficient * mode(i);
        }
    }

    void remove(std::vector<double>& vector) const {
        if (vector.size() != dimension_)
            throw std::invalid_argument("null space projection dimension mismatch");
        for (const auto value : vector)
            if (!std::isfinite(value))
                throw std::invalid_argument("cannot project a non-finite vector");
        for (const auto& mode : basis_) {
            double coefficient = 0.0;
            for (std::size_t i = 0; i < dimension_; ++i)
                coefficient += vector[i] * mode(i);
            for (std::size_t i = 0; i < dimension_; ++i)
                vector[i] -= coefficient * mode(i);
        }
    }

    double component_norm(const Vector& vector) const {
        check_vector(vector);
        double squared = 0.0;
        for (const auto& mode : basis_) {
            const double coefficient = vector.dot(mode);
            squared += coefficient * coefficient;
        }
        return std::sqrt(std::max(0.0, squared));
    }

    bool is_compatible(const Vector& rhs,
                       double relative_tolerance = 1e-12,
                       double absolute_tolerance = 1e-14) const {
        if (!(relative_tolerance >= 0.0) || !std::isfinite(relative_tolerance) ||
            !(absolute_tolerance >= 0.0) || !std::isfinite(absolute_tolerance))
            throw std::invalid_argument("null space compatibility tolerances must be finite and non-negative");
        const double threshold = absolute_tolerance + relative_tolerance * rhs.norm2();
        return component_norm(rhs) <= threshold;
    }

    double operator_residual(const SparseMatrix& matrix) const {
        if (matrix.n_rows() != dimension_ || matrix.n_cols() != dimension_)
            throw std::invalid_argument("null space operator dimension mismatch");

        double matrix_norm = 0.0;
        for (std::size_t k = 0; k < matrix.nnz(); ++k) {
            const double value = matrix.values_data()[k];
            if (!std::isfinite(value))
                throw std::invalid_argument("null space operator contains non-finite values");
            matrix_norm = std::hypot(matrix_norm, value);
        }
        double worst = 0.0;
        for (const auto& mode : basis_) {
            const auto product = matrix.matvec(mode);
            double product_norm = 0.0;
            for (const double value : product)
                product_norm = std::hypot(product_norm, value);
            worst = std::max(
                worst,
                product_norm / std::max(
                    matrix_norm, std::numeric_limits<double>::min()));
        }
        return worst;
    }

    bool is_null_space(const SparseMatrix& matrix,
                       double relative_tolerance = 1e-12) const {
        if (!(relative_tolerance >= 0.0) || !std::isfinite(relative_tolerance))
            throw std::invalid_argument("null space operator tolerance must be finite and non-negative");
        return operator_residual(matrix) <= relative_tolerance;
    }

private:
    static double stable_norm(const Vector& vector) {
        double norm = 0.0;
        for (std::size_t i = 0; i < vector.size(); ++i)
            norm = std::hypot(norm, vector(i));
        return norm;
    }

    void check_vector(const Vector& vector) const {
        if (vector.size() != dimension_)
            throw std::invalid_argument("null space projection dimension mismatch");
        if (!vector.is_valid())
            throw std::invalid_argument("cannot project a non-finite vector");
    }

    std::size_t dimension_ = 0;
    std::vector<Vector> basis_;
};

} // namespace cfdx::core
