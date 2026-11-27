#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

template <std::size_t Rows, std::size_t Cols>
using Matrix = std::array<std::array<double, Cols>, Rows>;

template <std::size_t N>
using Square = Matrix<N, N>;

template <std::size_t N>
Square<N> identity() {
    Square<N> result{};
    for (std::size_t i = 0; i < N; ++i) result[i][i] = 1.0;
    return result;
}

template <std::size_t Rows, std::size_t Cols>
Matrix<Cols, Rows> transpose(const Matrix<Rows, Cols>& input) {
    Matrix<Cols, Rows> result{};
    for (std::size_t i = 0; i < Rows; ++i)
        for (std::size_t j = 0; j < Cols; ++j)
            result[j][i] = input[i][j];
    return result;
}

template <std::size_t Rows, std::size_t Inner, std::size_t Cols>
Matrix<Rows, Cols> multiply(const Matrix<Rows, Inner>& lhs,
                            const Matrix<Inner, Cols>& rhs) {
    Matrix<Rows, Cols> result{};
    for (std::size_t i = 0; i < Rows; ++i)
        for (std::size_t k = 0; k < Inner; ++k)
            for (std::size_t j = 0; j < Cols; ++j)
                result[i][j] += lhs[i][k] * rhs[k][j];
    return result;
}

template <std::size_t N>
Square<N> add(const Square<N>& lhs, const Square<N>& rhs) {
    Square<N> result{};
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j)
            result[i][j] = lhs[i][j] + rhs[i][j];
    return result;
}

template <std::size_t N>
Square<N> scaled(double factor, const Square<N>& input) {
    Square<N> result{};
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j)
            result[i][j] = factor * input[i][j];
    return result;
}

template <std::size_t N>
Square<N> inverse(const Square<N>& input) {
    std::array<std::array<double, 2 * N>, N> augmented{};
    double scale = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            if (!std::isfinite(input[i][j]))
                throw std::runtime_error("inverse: non-finite matrix");
            augmented[i][j] = input[i][j];
            scale = std::max(scale, std::abs(input[i][j]));
        }
        augmented[i][N + i] = 1.0;
    }
    if (!(scale > 0.0)) throw std::runtime_error("inverse: zero matrix");
    const double pivot_tolerance =
        64.0 * std::numeric_limits<double>::epsilon() * scale;

    for (std::size_t col = 0; col < N; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < N; ++row)
            if (std::abs(augmented[row][col]) >
                std::abs(augmented[pivot][col]))
                pivot = row;
        if (!(std::abs(augmented[pivot][col]) > pivot_tolerance))
            throw std::runtime_error("inverse: singular or numerically rank-deficient matrix");
        if (pivot != col) std::swap(augmented[pivot], augmented[col]);

        const double diagonal = augmented[col][col];
        for (double& value : augmented[col]) value /= diagonal;
        for (std::size_t row = 0; row < N; ++row) {
            if (row == col) continue;
            const double factor = augmented[row][col];
            for (std::size_t j = 0; j < 2 * N; ++j)
                augmented[row][j] -= factor * augmented[col][j];
        }
    }

    Square<N> result{};
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j)
            result[i][j] = augmented[i][N + j];
    return result;
}

template <std::size_t N>
void require_spd(const Square<N>& matrix, const std::string& name) {
    Square<N> lower{};
    double scale = 0.0;
    for (const auto& row : matrix)
        for (double value : row) scale = std::max(scale, std::abs(value));
    const double tolerance =
        128.0 * std::numeric_limits<double>::epsilon() * scale;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            if (std::abs(matrix[i][j] - matrix[j][i]) > tolerance)
                throw std::runtime_error(name + ": matrix is not symmetric");
            double value = matrix[i][j];
            for (std::size_t k = 0; k < j; ++k)
                value -= lower[i][k] * lower[j][k];
            if (i == j) {
                if (!(value > tolerance))
                    throw std::runtime_error(name + ": matrix is not positive definite");
                lower[i][j] = std::sqrt(value);
            } else {
                lower[i][j] = value / lower[j][j];
            }
        }
    }
}

template <std::size_t N>
double max_difference(const Square<N>& lhs, const Square<N>& rhs) {
    double error = 0.0;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j)
            error = std::max(error, std::abs(lhs[i][j] - rhs[i][j]));
    return error;
}

void check_schur_identity(const Matrix<3, 3>& A,
                          const Matrix<2, 3>& B,
                          const Matrix<2, 2>& W,
                          const std::string& name) {
    require_spd(A, name + " A");
    require_spd(W, name + " W");
    const auto A_inverse = inverse(A);
    const auto W_inverse = inverse(W);
    const auto B_transpose = transpose(B);
    const auto S = multiply(multiply(B, A_inverse), B_transpose);
    require_spd(S, name + " Schur");

    for (double gamma : {1.0e-3, 0.05, 0.7, 5.0, 40.0}) {
        const auto grad_div = multiply(
            multiply(B_transpose, W_inverse), B);
        const auto A_gamma = add(A, scaled(gamma, grad_div));
        require_spd(A_gamma, name + " augmented velocity block");

        const auto S_gamma = multiply(
            multiply(B, inverse(A_gamma)), B_transpose);
        require_spd(S_gamma, name + " augmented Schur");

        const auto direct_inverse = inverse(S_gamma);
        const auto woodbury_inverse =
            add(inverse(S), scaled(gamma, W_inverse));
        const double error = max_difference(direct_inverse, woodbury_inverse);
        const double scale = std::max(1.0, [&] {
            double maximum = 0.0;
            for (const auto& row : woodbury_inverse)
                for (double value : row)
                    maximum = std::max(maximum, std::abs(value));
            return maximum;
        }());
        if (error > 2.0e-10 * scale)
            throw std::runtime_error(name + ": Woodbury Schur identity mismatch");
        std::cout << name << " gamma=" << gamma
                  << " max_abs_error=" << error << '\n';
    }
}

void check_singular_pressure_rejected() {
    const Matrix<2, 3> rank_deficient_B{{
        {{1.0, -1.0, 0.0}},
        {{2.0, -2.0, 0.0}}
    }};
    const Matrix<3, 3> A{{
        {{4.0, 1.0, 0.0}},
        {{1.0, 3.0, 0.2}},
        {{0.0, 0.2, 2.0}}
    }};
    const auto S = multiply(
        multiply(rank_deficient_B, inverse(A)), transpose(rank_deficient_B));
    bool rejected = false;
    try {
        (void)inverse(S);
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    if (!rejected)
        throw std::runtime_error("rank-deficient pressure Schur was not rejected");
}

}  // namespace

int main() {
    try {
        const Matrix<3, 3> A{{
            {{4.0, 1.0, 0.0}},
            {{1.0, 3.0, 0.2}},
            {{0.0, 0.2, 2.0}}
        }};
        const Matrix<2, 3> B_pinned{{
            {{1.0, -1.0, 0.0}},
            {{0.0, 1.0, -1.0}}
        }};
        const Matrix<2, 2> W{{
            {{2.0, 0.2}},
            {{0.2, 1.5}}
        }};
        check_schur_identity(A, B_pinned, W, "pinned_pressure");

        // Periodic three-cell incidence has B^T * 1 = 0. Restrict pressure
        // to an orthonormal mean-zero basis before applying the inverse identity.
        const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
        const double inv_sqrt6 = 1.0 / std::sqrt(6.0);
        const Matrix<3, 2> Q{{
            {{inv_sqrt2, inv_sqrt6}},
            {{-inv_sqrt2, inv_sqrt6}},
            {{0.0, -2.0 * inv_sqrt6}}
        }};
        const Matrix<3, 3> B_periodic{{
            {{-1.0, 1.0, 0.0}},
            {{0.0, -1.0, 1.0}},
            {{1.0, 0.0, -1.0}}
        }};
        const auto B_mean_zero = multiply(transpose(Q), B_periodic);
        check_schur_identity(A, B_mean_zero, identity<2>(), "mean_zero_pressure");

        check_singular_pressure_rejected();
        std::cout << "Augmented-Lagrangian Schur oracle: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Augmented-Lagrangian Schur oracle: " << error.what() << '\n';
        return 1;
    }
}
