#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace {

using CellVector = std::array<double, 3>;
using GradientVector = std::array<double, 2>;
using Matrix = std::array<std::array<double, 3>, 3>;
using Reconstruction = std::array<std::array<double, 3>, 2>;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

GradientVector apply_reconstruction(
    const Reconstruction& g, const CellVector& x) {
    GradientVector result{};
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t c = 0; c < 3; ++c)
            result[i] += g[i][c] * x[c];
    return result;
}

CellVector apply_transpose(
    const Reconstruction& g, const GradientVector& y) {
    CellVector result{};
    for (std::size_t c = 0; c < 3; ++c)
        for (std::size_t i = 0; i < 2; ++i)
            result[c] += g[i][c] * y[i];
    return result;
}

Matrix gram_matrix(const Reconstruction& g) {
    Matrix a{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 2; ++k)
                a[i][j] += g[k][i] * g[k][j];
    return a;
}

CellVector apply_explicit(const Matrix& a, const CellVector& x) {
    CellVector result{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            result[i] += a[i][j] * x[j];
    return result;
}

double dot(const CellVector& x, const CellVector& y) {
    double value = 0.0;
    for (std::size_t i = 0; i < 3; ++i) value += x[i] * y[i];
    return value;
}

void check_close(double a, double b, double tolerance, const char* message) {
    require(std::abs(a - b) <= tolerance, message);
}

}  // namespace

int main() {
    // Exact counterexample also formalized in
    // CFDX.FiniteVolume.EnergyDiffusion (Lean).
    const Reconstruction g{{
        {{1.0, 0.0, -1.0}},
        {{1.0 / 5.0, 1.0, -6.0 / 5.0}}
    }};
    const Matrix expected{{
        {{26.0 / 25.0, 1.0 / 5.0, -31.0 / 25.0}},
        {{1.0 / 5.0, 1.0, -6.0 / 5.0}},
        {{-31.0 / 25.0, -6.0 / 5.0, 61.0 / 25.0}}
    }};
    const Matrix a = gram_matrix(g);
    constexpr double tolerance = 2.0e-15;

    // Exact coefficient oracle for G^T G.
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            check_close(a[i][j], expected[i][j], tolerance,
                        "G^T G coefficient mismatch");
            check_close(a[i][j], a[j][i], tolerance,
                        "Gram matrix is not symmetric");
        }
    }

    // Conservation: G annihilates constants, hence A does too.
    const CellVector ones{{1.0, 1.0, 1.0}};
    const auto constant_gradient = apply_reconstruction(g, ones);
    check_close(constant_gradient[0], 0.0, tolerance,
                "first reconstruction row does not annihilate constants");
    check_close(constant_gradient[1], 0.0, tolerance,
                "second reconstruction row does not annihilate constants");
    const auto constant_action = apply_explicit(a, ones);
    for (double value : constant_action)
        check_close(value, 0.0, tolerance, "Gram operator does not conserve constants");

    // Factorized and explicit actions must agree for every basis vector.
    for (std::size_t basis = 0; basis < 3; ++basis) {
        CellVector e{};
        e[basis] = 1.0;
        const auto factorized = apply_transpose(g, apply_reconstruction(g, e));
        const auto explicit_action = apply_explicit(a, e);
        for (std::size_t i = 0; i < 3; ++i)
            check_close(factorized[i], explicit_action[i], tolerance,
                        "factorized action differs from explicit CSR oracle");
    }

    // Positive semidefiniteness follows from x^T G^T G x = ||Gx||^2.
    const std::array<CellVector, 5> samples{{
        {{0.0, 0.0, 0.0}},
        {{1.0, -2.0, 3.0}},
        {{-4.0, 0.5, 2.0}},
        {{1.0, 1.0, 1.0}},
        {{1.0e3, -2.0e3, 3.0e3}}
    }};
    for (const auto& x : samples) {
        const auto gx = apply_reconstruction(g, x);
        const double energy = dot(x, apply_explicit(a, x));
        const double factorized_energy =
            gx[0] * gx[0] + gx[1] * gx[1];
        check_close(energy, factorized_energy,
                    2.0e-9 * (1.0 + factorized_energy),
                    "energy identity failed");
        require(energy >= -2.0e-9 * (1.0 + factorized_energy),
                "Gram operator has negative energy");
    }

    // Crucial limitation: PSD and conservation do not imply a Z/M-matrix.
    require(a[0][1] > 0.0,
            "counterexample lost its positive off-diagonal coefficient");
    return 0;
}
