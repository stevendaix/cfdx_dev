#include "cfdx/core/fvm/least_squares_gradient.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/core/mesh/mesh.h"

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
using TwoCellField = std::array<double, 2>;
using TwoCellGradients = std::array<Gradient, 2>;
using TwoCellGradientColumns = std::array<TwoCellGradients, 2>;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void close(double actual, double expected, const char* message) {
    if (std::abs(actual - expected) > tolerance * (1.0 + std::abs(expected)))
        throw std::runtime_error(message);
}

Gradient evaluate_production_lsq(const CellVector& values) {
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
    return {gradient.x, gradient.y, gradient.z};
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


cfdx::core::Mesh make_unit_cube_mesh() {
    using namespace cfdx::core;
    Mesh mesh;
    mesh.points().resize(8);
    mesh.points().set(0, 0.0, 0.0, 0.0);
    mesh.points().set(1, 1.0, 0.0, 0.0);
    mesh.points().set(2, 1.0, 1.0, 0.0);
    mesh.points().set(3, 0.0, 1.0, 0.0);
    mesh.points().set(4, 0.0, 0.0, 1.0);
    mesh.points().set(5, 1.0, 0.0, 1.0);
    mesh.points().set(6, 1.0, 1.0, 1.0);
    mesh.points().set(7, 0.0, 1.0, 1.0);
    mesh.faces().push_face({0, 3, 2, 1});
    mesh.faces().push_face({4, 5, 6, 7});
    mesh.faces().push_face({0, 1, 5, 4});
    mesh.faces().push_face({3, 7, 6, 2});
    mesh.faces().push_face({0, 4, 7, 3});
    mesh.faces().push_face({1, 2, 6, 5});
    mesh.ownership().resize(6);
    for (std::size_t face = 0; face < 6; ++face) {
        mesh.ownership().set_owner(face, 0);
        mesh.ownership().set_neighbour(face, FaceOwnership::BOUNDARY);
    }
    mesh.cells().push_cell({0, 1, 2, 3, 4, 5});
    return mesh;
}

Gradient evaluate_mesh_lsq(
    const cfdx::core::Mesh& mesh,
    double cell_value,
    const std::vector<cfdx::core::BoundaryGradientCondition>& conditions) {
    using namespace cfdx::core;
    Field<double, Location::CELL> field(1, "phi", "1", 1);
    field(0) = cell_value;
    const auto result = compute_gradient_least_squares(
        field, mesh, BoundaryGradientPolicy::ZERO_GRADIENT_GHOST, &conditions);
    return {result(0, 0), result(0, 1), result(0, 2)};
}


cfdx::core::Mesh make_two_cell_hex_mesh() {
    using namespace cfdx::core;
    Mesh mesh;
    mesh.points().resize(12);
    mesh.points().set(0, 0.0, 0.0, 0.0);
    mesh.points().set(1, 1.0, 0.0, 0.0);
    mesh.points().set(2, 1.0, 1.0, 0.0);
    mesh.points().set(3, 0.0, 1.0, 0.0);
    mesh.points().set(4, 0.0, 0.0, 1.0);
    mesh.points().set(5, 1.0, 0.0, 1.0);
    mesh.points().set(6, 1.0, 1.0, 1.0);
    mesh.points().set(7, 0.0, 1.0, 1.0);
    mesh.points().set(8, 2.0, 0.0, 0.0);
    mesh.points().set(9, 2.0, 1.0, 0.0);
    mesh.points().set(10, 2.0, 0.0, 1.0);
    mesh.points().set(11, 2.0, 1.0, 1.0);

    mesh.faces().push_face({0, 4, 7, 3});   // x = 0
    mesh.faces().push_face({8, 10, 11, 9}); // x = 2
    mesh.faces().push_face({1, 2, 6, 5});  // shared x = 1
    mesh.faces().push_face({0, 3, 2, 1});  // z = 0, cell 0
    mesh.faces().push_face({1, 2, 9, 8});  // z = 0, cell 1
    mesh.faces().push_face({4, 5, 6, 7});  // z = 1, cell 0
    mesh.faces().push_face({5, 10, 11, 6}); // z = 1, cell 1
    mesh.faces().push_face({0, 1, 5, 4});  // y = 0, cell 0
    mesh.faces().push_face({1, 8, 10, 5}); // y = 0, cell 1
    mesh.faces().push_face({3, 7, 6, 2});  // y = 1, cell 0
    mesh.faces().push_face({2, 6, 11, 9}); // y = 1, cell 1

    mesh.ownership().resize(11);
    for (std::size_t face = 0; face < 11; ++face) {
        const std::size_t owner = face == 1 || face == 4 || face == 6 ||
                                  face == 8 || face == 10 ? 1 : 0;
        mesh.ownership().set_owner(face, owner);
        mesh.ownership().set_neighbour(
            face, face == 2 ? 1 : FaceOwnership::BOUNDARY);
    }
    mesh.cells().push_cell({0, 2, 3, 5, 7, 9});
    mesh.cells().push_cell({2, 1, 4, 6, 8, 10});
    return mesh;
}

TwoCellGradients evaluate_two_cell_mesh_lsq(
    const cfdx::core::Mesh& mesh,
    const TwoCellField& cell_values,
    const std::vector<cfdx::core::BoundaryGradientCondition>& conditions) {
    using namespace cfdx::core;
    Field<double, Location::CELL> field(2, "phi", "1", 1);
    field(0) = cell_values[0];
    field(1) = cell_values[1];
    const auto result = compute_gradient_least_squares(
        field, mesh, BoundaryGradientPolicy::ZERO_GRADIENT_GHOST, &conditions);
    return TwoCellGradients{{
        Gradient{result(0, 0), result(0, 1), result(0, 2)},
        Gradient{result(1, 0), result(1, 1), result(1, 2)}
    }};
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

    // Mesh-level affine boundary contract. The deliberately mixed boundary
    // types avoid the degenerate all-Dirichlet symmetric one-cell case:
    // face 4 is Dirichlet, face 5 Neumann, and the remaining faces are
    // zero-gradient. Geometry, stencil membership and boundary types stay fixed.
    {
        using namespace cfdx::core;
        const Mesh mesh = make_unit_cube_mesh();
        std::vector<BoundaryGradientCondition> homogeneous(6);
        homogeneous[4] = {BoundaryGradientConditionType::DIRICHLET, 0.0};
        homogeneous[5] = {BoundaryGradientConditionType::NEUMANN, 0.0};
        const Gradient linear_column = evaluate_mesh_lsq(mesh, 1.0, homogeneous);

        std::vector<BoundaryGradientCondition> prescribed = homogeneous;
        prescribed[4].value = 2.25;
        prescribed[5].value = -0.75;
        prescribed[0] = {BoundaryGradientConditionType::NEUMANN, 0.4};
        const Gradient offset = evaluate_mesh_lsq(mesh, 0.0, prescribed);

        for (double cell_value : {-2.0, 0.0, 0.5, 3.0}) {
            const Gradient actual = evaluate_mesh_lsq(mesh, cell_value, prescribed);
            for (std::size_t component = 0; component < 3; ++component) {
                close(actual[component],
                      offset[component] + cell_value * linear_column[component],
                      "mesh LSQ with fixed nonhomogeneous BCs must be affine in cell values");
            }
        }

        // The offset must vanish when all prescribed boundary data vanish.
        const Gradient zero_offset = evaluate_mesh_lsq(mesh, 0.0, homogeneous);
        for (double value : zero_offset)
            close(value, 0.0, "homogeneous boundary data must have zero affine offset");
    }


    // Two-cell mesh-level affine decomposition. Unlike the preceding
    // one-cell fixture, this recovers a genuine 6x2 global cell-to-gradient
    // map from the production mesh wrapper. Boundary types and geometry are
    // frozen while nonhomogeneous data are represented by a separate offset.
    {
        using namespace cfdx::core;
        const Mesh mesh = make_two_cell_hex_mesh();
        std::vector<BoundaryGradientCondition> prescribed(11);
        for (auto& bc : prescribed)
            bc = {BoundaryGradientConditionType::ZERO_GRADIENT, 0.0};
        prescribed[0] = {BoundaryGradientConditionType::DIRICHLET, 1.25};
        prescribed[1] = {BoundaryGradientConditionType::NEUMANN, -0.5};
        prescribed[3] = {BoundaryGradientConditionType::NEUMANN, 0.4};

        auto homogeneous = prescribed;
        for (auto& bc : homogeneous) bc.value = 0.0;

        TwoCellGradientColumns columns{};
        for (std::size_t basis = 0; basis < 2; ++basis) {
            TwoCellField e{};
            e[basis] = 1.0;
            const auto gradient = evaluate_two_cell_mesh_lsq(mesh, e, homogeneous);
            for (std::size_t cell = 0; cell < 2; ++cell)
                columns[basis][cell] = gradient[cell];
        }

        const TwoCellField zero_field{};
        const auto boundary_offset =
            evaluate_two_cell_mesh_lsq(mesh, zero_field, prescribed);
        bool has_nonzero_offset = false;
        for (const auto& cell_gradient : boundary_offset)
            for (double value : cell_gradient)
                has_nonzero_offset |= std::abs(value) > tolerance;
        require(has_nonzero_offset,
                "nonhomogeneous boundary data must produce a nonzero affine offset");

        for (const TwoCellField values : {
                 TwoCellField{{-2.0, 0.5}},
                 TwoCellField{{0.0, 0.0}},
                 TwoCellField{{1.25, -3.0}},
                 TwoCellField{{4.0, 2.0}}}) {
            const auto actual =
                evaluate_two_cell_mesh_lsq(mesh, values, prescribed);
            for (std::size_t cell = 0; cell < 2; ++cell) {
                for (std::size_t component = 0; component < 3; ++component) {
                    const double homogeneous_action =
                        columns[0][cell][component] * values[0] +
                        columns[1][cell][component] * values[1];
                    close(actual[cell][component],
                          boundary_offset[cell][component] + homogeneous_action,
                          "two-cell production LSQ must equal G*phi + H(b)");
                }
            }
        }

        // Verify additivity across independent physical boundary sources,
        // not merely homogeneity along one combined data vector.
        auto dirichlet_only = homogeneous;
        dirichlet_only[0].value = prescribed[0].value;
        auto neumann_only = homogeneous;
        neumann_only[1].value = prescribed[1].value;
        neumann_only[3].value = prescribed[3].value;
        const auto dirichlet_offset =
            evaluate_two_cell_mesh_lsq(mesh, zero_field, dirichlet_only);
        const auto neumann_offset =
            evaluate_two_cell_mesh_lsq(mesh, zero_field, neumann_only);
        for (std::size_t cell = 0; cell < 2; ++cell)
            for (std::size_t component = 0; component < 3; ++component)
                close(boundary_offset[cell][component],
                      dirichlet_offset[cell][component] +
                          neumann_offset[cell][component],
                      "Dirichlet and Neumann affine offsets must superpose");

        // Also verify homogeneity by doubling the complete boundary vector.
        auto doubled_boundary = prescribed;
        for (auto& bc : doubled_boundary) bc.value *= 2.0;
        const auto doubled_offset =
            evaluate_two_cell_mesh_lsq(mesh, zero_field, doubled_boundary);
        for (std::size_t cell = 0; cell < 2; ++cell)
            for (std::size_t component = 0; component < 3; ++component)
                close(doubled_offset[cell][component],
                      2.0 * boundary_offset[cell][component],
                      "boundary offset must be homogeneous in prescribed BC data");
    }

    return 0;
}
