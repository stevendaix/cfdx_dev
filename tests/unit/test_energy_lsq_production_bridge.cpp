#include "cfdx/core/fvm/least_squares_gradient.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace {

constexpr std::size_t cell_count = 7;
constexpr double tolerance = 2.0e-13;

using CellVector = std::array<double, cell_count>;
using Matrix = std::array<std::array<double, cell_count>, cell_count>;
using Gradient = std::array<double, 3>;
using Reconstruction = std::array<Gradient, cell_count>;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void close(double actual, double expected, const char* message) {
    if (std::abs(actual - expected) > tolerance * (1.0 + std::abs(expected)))
        throw std::runtime_error(message);
}

cfdx::core::Vec3 evaluate_production_lsq(const CellVector& values) {
    using cfdx::core::Vec3;
    const Vec3 centre{0.0, 0.0, 0.0};
    const std::vector<Vec3> neighbours{
        {1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {0.0, -1.0, 0.0},
        {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}
    };
    const std::vector<double> neighbour_values(
        values.begin() + 1, values.end());
    cfdx::core::StencilQuality quality;
    const Vec3 gradient = cfdx::core::least_squares_gradient(
        centre, values[0], neighbours, neighbour_values, &quality);
    require(quality.rank == 3 && quality.full_rank,
            "production LSQ stencil must remain full rank");
    require(quality.finite && std::isfinite(quality.condition_estimate),
            "production LSQ stencil quality must be finite");
    return gradient;
}

Gradient apply_reconstruction(
    const Reconstruction& reconstruction, const CellVector& values) {
    Gradient result{};
    for (std::size_t c = 0; c < cell_count; ++c)
        for (std::size_t g = 0; g < 3; ++g)
            result[g] += reconstruction[c][g] * values[c];
    return result;
}

CellVector apply_transpose_metric(
    const Reconstruction& reconstruction, const Gradient& gradient) {
    // M = diag(2, 3, 5), an explicitly SPD gradient-space metric.
    constexpr std::array<double, 3> metric{{2.0, 3.0, 5.0}};
    CellVector result{};
    for (std::size_t c = 0; c < cell_count; ++c)
        for (std::size_t g = 0; g < 3; ++g)
            result[c] += reconstruction[c][g] * metric[g] * gradient[g];
    return result;
}

Matrix assemble_gram(const Reconstruction& reconstruction) {
    constexpr std::array<double, 3> metric{{2.0, 3.0, 5.0}};
    Matrix result{};
    for (std::size_t i = 0; i < cell_count; ++i)
        for (std::size_t j = 0; j < cell_count; ++j)
            for (std::size_t g = 0; g < 3; ++g)
                result[i][j] += reconstruction[i][g] *
                                metric[g] * reconstruction[j][g];
    return result;
}

CellVector apply_matrix(const Matrix& matrix, const CellVector& values) {
    CellVector result{};
    for (std::size_t i = 0; i < cell_count; ++i)
        for (std::size_t j = 0; j < cell_count; ++j)
            result[i] += matrix[i][j] * values[j];
    return result;
}

}  // namespace

int main() {
    // Freeze geometry, neighbour ordering, inverse-distance-squared weights,
    // automatic rank selection and the production one-sided-Jacobi solver.
    // Recover the actual linear reconstruction by applying the production
    // kernel to each cell basis vector; do not hard-code its coefficients.
    Reconstruction reconstruction{};
    for (std::size_t basis = 0; basis < cell_count; ++basis) {
        CellVector e{};
        e[basis] = 1.0;
        const Gradient column = evaluate_production_lsq(e);
        for (std::size_t g = 0; g < 3; ++g)
            reconstruction[basis][g] = column[g];
    }

    // The production reconstruction must annihilate constants.
    CellVector ones{};
    ones.fill(1.0);
    const Gradient constant_gradient = evaluate_production_lsq(ones);
    for (double value : constant_gradient)
        close(value, 0.0, "production LSQ must annihilate constants");

    // Check linearity on a nontrivial field independently of the basis probes.
    const CellVector field{{1.25, 2.0, -0.5, 3.0, 0.25, -2.0, 4.5}};
    const Gradient direct_gradient = evaluate_production_lsq(field);
    const Gradient recovered_gradient =
        apply_reconstruction(reconstruction, field);
    for (std::size_t g = 0; g < 3; ++g)
        close(direct_gradient[g], recovered_gradient[g],
              "basis-recovered map must reproduce production LSQ");

    // Build the energy operator from the production-derived G and SPD M.
    const Matrix gram = assemble_gram(reconstruction);
    for (std::size_t i = 0; i < cell_count; ++i) {
        for (std::size_t j = 0; j < cell_count; ++j) {
            close(gram[i][j], gram[j][i],
                  "G^T M G must be symmetric");
        }
    }

    // Compare factorized and assembled actions for a mixed field and every
    // basis vector. This is the exact operator bridge; it is not yet a
    // production diffusion discretization or a mesh-flux implementation.
    const auto factorized = apply_transpose_metric(
        reconstruction, recovered_gradient);
    const auto assembled = apply_matrix(gram, field);
    for (std::size_t i = 0; i < cell_count; ++i)
        close(factorized[i], assembled[i],
              "factorized action must equal assembled G^T M G");

    for (std::size_t basis = 0; basis < cell_count; ++basis) {
        CellVector e{};
        e[basis] = 1.0;
        const Gradient ge = evaluate_production_lsq(e);
        const auto factorized_basis = apply_transpose_metric(reconstruction, ge);
        const auto assembled_basis = apply_matrix(gram, e);
        for (std::size_t i = 0; i < cell_count; ++i)
            close(factorized_basis[i], assembled_basis[i],
                  "basis action mismatch for G^T M G");
    }

    // Energy identity and PSD follow from the production-derived factor G.
    const double energy = [&]() {
        const auto action = apply_matrix(gram, field);
        double value = 0.0;
        for (std::size_t i = 0; i < cell_count; ++i)
            value += field[i] * action[i];
        return value;
    }();
    const double factorized_energy =
        2.0 * recovered_gradient[0] * recovered_gradient[0] +
        3.0 * recovered_gradient[1] * recovered_gradient[1] +
        5.0 * recovered_gradient[2] * recovered_gradient[2];
    close(energy, factorized_energy,
          "x^T G^T M G x must equal the weighted gradient energy");
    require(energy >= -tolerance,
            "SPD metric must produce nonnegative energy");

    const auto constant_action = apply_matrix(gram, ones);
    for (double value : constant_action)
        close(value, 0.0, "G^T M G must annihilate constants");
    return 0;
}
