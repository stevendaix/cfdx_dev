#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
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


template <std::size_t Rows, std::size_t Cols>
std::array<double, Rows> multiply_vector(
    const Matrix<Rows, Cols>& matrix,
    const std::array<double, Cols>& vector) {
    std::array<double, Rows> result{};
    for (std::size_t i = 0; i < Rows; ++i)
        for (std::size_t j = 0; j < Cols; ++j)
            result[i] += matrix[i][j] * vector[j];
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


void check_general_cross_augmentation() {
    const Matrix<3, 3> A{{
        {{4.0, 1.0, 0.0}},
        {{1.0, 3.0, 0.2}},
        {{0.0, 0.2, 2.0}}
    }};
    const Matrix<2, 3> D{{
        {{1.0, -1.0, 0.0}},
        {{0.0, 1.0, -1.0}}
    }};
    // Deliberately non-adjoint: G != -D^T.
    const Matrix<3, 2> G{{
        {{-0.98, -0.01}},
        {{1.03, -0.985}},
        {{-0.02, 1.025}}
    }};
    const Matrix<2, 2> W{{
        {{2.0, 0.2}},
        {{0.2, 1.5}}
    }};
    require_spd(A, "general cross-augmentation A");
    require_spd(W, "general cross-augmentation W");

    const auto A_inverse = inverse(A);
    const auto W_inverse = inverse(W);
    const auto S = scaled(
        -1.0, multiply(multiply(D, A_inverse), G));
    const auto standard_grad_div = multiply(
        multiply(transpose(D), W_inverse), D);

    for (double gamma : {1.0e-3, 0.05, 0.7, 5.0, 40.0}) {
        // General cross augmentation. It preserves the exact Schur inverse
        // shift for C=0 without assuming G=-D^T, but A_gamma is generally
        // nonsymmetric and is deliberately not asserted SPD.
        const auto cross_update = multiply(multiply(G, W_inverse), D);
        const auto A_gamma = add(A, scaled(-gamma, cross_update));
        const auto A_gamma_inverse = inverse(A_gamma);
        const auto S_gamma = scaled(
            -1.0, multiply(multiply(D, A_gamma_inverse), G));
        const auto direct_inverse = inverse(S_gamma);
        const auto expected_inverse =
            add(inverse(S), scaled(gamma, W_inverse));
        const double error = max_difference(direct_inverse, expected_inverse);
        if (error > 2.0e-10)
            throw std::runtime_error(
                "general cross-augmentation Schur identity mismatch");
        std::cout << "general_cross gamma=" << gamma
                  << " max_abs_error=" << error << '\n';

        if (gamma == 0.7 &&
            max_difference(A_gamma, transpose(A_gamma)) < 1.0e-6)
            throw std::runtime_error(
                "non-adjoint cross-augmentation unexpectedly remained symmetric");

        // Negative control 1: the conventional D^T W^{-1} D update does not
        // inherit the exact identity for this non-adjoint G.
        const auto A_gamma_symmetric = add(
            A, scaled(gamma, standard_grad_div));
        const auto S_gamma_symmetric = scaled(
            -1.0, multiply(
                multiply(D, inverse(A_gamma_symmetric)), G));
        const double conventional_error = max_difference(
            inverse(S_gamma_symmetric), expected_inverse);
        if (gamma == 0.7 && conventional_error < 1.0e-3)
            throw std::runtime_error(
                "non-adjoint conventional grad-div negative control did not fail");
    }

    // Negative control 2: adding a nonzero C pressure block invalidates the
    // simple inverse-shift identity, even for the cross augmentation.
    const double gamma = 0.7;
    const auto A_gamma = add(
        A, scaled(-gamma, multiply(multiply(G, W_inverse), D)));
    const auto S_gamma = scaled(
        -1.0, multiply(multiply(D, inverse(A_gamma)), G));
    const auto C = scaled(0.1, identity<2>());
    const auto S_with_C = add(S, C);
    const auto S_gamma_with_C = add(S_gamma, C);
    const auto expected_with_C =
        add(inverse(S_with_C), scaled(gamma, W_inverse));
    const double nonzero_C_error = max_difference(
        inverse(S_gamma_with_C), expected_with_C);
    if (nonzero_C_error < 1.0e-3)
        throw std::runtime_error(
            "nonzero-C negative control did not expose the missing C correction");
}


void check_C_aware_row_elimination() {
    const Matrix<3, 3> A{{
        {{4.0, 1.0, 0.0}},
        {{1.0, 3.0, 0.2}},
        {{0.0, 0.2, 2.0}}
    }};
    const Matrix<2, 3> D{{
        {{1.0, -1.0, 0.0}},
        {{0.0, 1.0, -1.0}}
    }};
    const Matrix<3, 2> G{{
        {{-0.98, -0.01}},
        {{1.03, -0.985}},
        {{-0.02, 1.025}}
    }};
    const Matrix<2, 2> C = scaled(0.1, identity<2>());
    const Matrix<2, 2> W{{
        {{2.0, 0.2}},
        {{0.2, 1.5}}
    }};
    require_spd(A, "C-aware row elimination A");
    require_spd(W, "C-aware row elimination W");

    const auto GtG = multiply(transpose(G), G);
    require_spd(GtG, "C-aware row elimination GtG");
    const auto L = multiply(
        multiply(C, inverse(GtG)), transpose(G));
    const auto LG = multiply(L, G);
    if (max_difference(LG, C) > 2.0e-12)
        throw std::runtime_error("C-aware row elimination failed to factor C=L*G");

    const auto LA = multiply(L, A);
    Matrix<2, 3> D_tilde{};
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            D_tilde[i][j] = D[i][j] - LA[i][j];

    const auto A_inverse = inverse(A);
    const auto W_inverse = inverse(W);
    const auto S = add(
        C, scaled(-1.0, multiply(multiply(D, A_inverse), G)));
    const auto S_transformed = scaled(
        -1.0, multiply(multiply(D_tilde, A_inverse), G));
    if (max_difference(S, S_transformed) > 2.0e-12)
        throw std::runtime_error(
            "C-aware row elimination changed the original Schur operator");

    for (double gamma : {1.0e-3, 0.05, 0.7, 5.0, 40.0}) {
        const auto cross_update = multiply(
            multiply(G, W_inverse), D_tilde);
        const auto A_gamma = add(A, scaled(-gamma, cross_update));
        const auto S_gamma = scaled(
            -1.0, multiply(multiply(D_tilde, inverse(A_gamma)), G));
        const auto actual_inverse = inverse(S_gamma);
        const auto expected_inverse =
            add(inverse(S), scaled(gamma, W_inverse));
        if (max_difference(actual_inverse, expected_inverse) > 2.0e-10)
            throw std::runtime_error(
                "C-aware row-elimination Schur identity mismatch");
    }

    // Verify both block row operations, including their right-hand sides,
    // preserve a known solution of the original nonzero-C system.
    const std::array<double, 3> velocity{{0.5, -1.0, 2.0}};
    const std::array<double, 2> pressure{{0.25, -0.75}};
    const auto Au = multiply_vector(A, velocity);
    const auto Gp = multiply_vector(G, pressure);
    const auto Du = multiply_vector(D, velocity);
    const auto Cp = multiply_vector(C, pressure);
    std::array<double, 3> momentum_rhs{};
    std::array<double, 2> continuity_rhs{};
    for (std::size_t i = 0; i < 3; ++i)
        momentum_rhs[i] = Au[i] + Gp[i];
    for (std::size_t i = 0; i < 2; ++i)
        continuity_rhs[i] = Du[i] + Cp[i];

    const auto Lf = multiply_vector(L, momentum_rhs);
    std::array<double, 2> transformed_continuity_rhs{};
    for (std::size_t i = 0; i < 2; ++i)
        transformed_continuity_rhs[i] = continuity_rhs[i] - Lf[i];

    const double gamma = 0.7;
    const auto A_gamma = add(
        A, scaled(-gamma, multiply(
            multiply(G, W_inverse), D_tilde)));
    const auto Winv_g = multiply_vector(W_inverse, transformed_continuity_rhs);
    const auto G_Winv_g = multiply_vector(G, Winv_g);
    const auto A_gamma_u = multiply_vector(A_gamma, velocity);
    std::array<double, 3> transformed_momentum_rhs{};
    for (std::size_t i = 0; i < 3; ++i)
        transformed_momentum_rhs[i] =
            momentum_rhs[i] - gamma * G_Winv_g[i];
    const auto Gp_again = multiply_vector(G, pressure);
    const auto Dtilde_u = multiply_vector(D_tilde, velocity);
    for (std::size_t i = 0; i < 3; ++i)
        if (std::abs(A_gamma_u[i] + Gp_again[i] -
                     transformed_momentum_rhs[i]) > 2.0e-12)
            throw std::runtime_error(
                "C-aware top-row RHS transformation changed the solution");
    for (std::size_t i = 0; i < 2; ++i)
        if (std::abs(Dtilde_u[i] - transformed_continuity_rhs[i]) > 2.0e-12)
            throw std::runtime_error(
                "C-aware bottom-row RHS transformation changed the solution");
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
        const auto Q_gram = multiply(transpose(Q), Q);
        if (max_difference(Q_gram, identity<2>()) > 2.0e-15)
            throw std::runtime_error("mean-zero pressure basis is not orthonormal");
        const auto B_transpose_periodic = transpose(B_periodic);
        for (std::size_t pressure_dof = 0; pressure_dof < 3; ++pressure_dof) {
            double constant_mode_action = 0.0;
            for (std::size_t cell = 0; cell < 3; ++cell)
                constant_mode_action += B_transpose_periodic[pressure_dof][cell];
            if (std::abs(constant_mode_action) > 1.0e-15)
                throw std::runtime_error("periodic incidence does not annihilate constant pressure");
        }
        const auto B_mean_zero = multiply(transpose(Q), B_periodic);
        check_schur_identity(A, B_mean_zero, identity<2>(), "mean_zero_pressure");

        check_singular_pressure_rejected();
        check_general_cross_augmentation();
        check_C_aware_row_elimination();
        std::cout << "Augmented-Lagrangian Schur oracle: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Augmented-Lagrangian Schur oracle: " << error.what() << '\n';
        return 1;
    }
}
