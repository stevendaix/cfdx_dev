#include "cfdx/core/numerics/laplacian_assembly.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

Mesh make_unit_cube() {
    Mesh m;
    m.points().resize(8);
    m.points().set(0, 0.0, 0.0, 0.0);
    m.points().set(1, 1.0, 0.0, 0.0);
    m.points().set(2, 1.0, 1.0, 0.0);
    m.points().set(3, 0.0, 1.0, 0.0);
    m.points().set(4, 0.0, 0.0, 1.0);
    m.points().set(5, 1.0, 0.0, 1.0);
    m.points().set(6, 1.0, 1.0, 1.0);
    m.points().set(7, 0.0, 1.0, 1.0);

    m.faces().push_face({0, 3, 2, 1});
    m.faces().push_face({4, 5, 6, 7});
    m.faces().push_face({0, 1, 5, 4});
    m.faces().push_face({3, 7, 6, 2});
    m.faces().push_face({0, 4, 7, 3});
    m.faces().push_face({1, 2, 6, 5});

    m.ownership().resize(6);
    for (std::size_t i = 0; i < 6; ++i) {
        m.ownership().set_owner(i, 0);
        m.ownership().set_neighbour(i, FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0, 1, 2, 3, 4, 5});
    return m;
}

} // namespace

int main() {
    run_case("laplacian_nonzero_dirichlet", []() {
        const Mesh mesh = make_unit_cube();
        SparseMatrix A;
        Vector b;
        const std::vector<bool> constrained{true};
        const std::vector<double> values{3.5};

        assemble_laplacian_csr(mesh, A, b, constrained, values);

        EXPECT_NEAR(A(0, 0), 1.0, 1e-14);
        EXPECT_NEAR(b(0), 3.5, 1e-14);
    });

    run_case("laplacian_homogeneous_dirichlet_compatibility", []() {
        const Mesh mesh = make_unit_cube();
        SparseMatrix A;
        Vector b;
        const std::vector<bool> constrained{true};

        assemble_laplacian_csr(mesh, A, b, constrained);

        EXPECT_NEAR(A(0, 0), 1.0, 1e-14);
        EXPECT_NEAR(b(0), 0.0, 1e-14);
    });

    run_case("laplacian_rejects_bad_dirichlet_values", []() {
        const Mesh mesh = make_unit_cube();
        SparseMatrix A;
        Vector b;
        const std::vector<bool> constrained{true};

        EXPECT_THROW(
            assemble_laplacian_csr(
                mesh, A, b, constrained, std::vector<double>{1.0, 2.0}),
            std::runtime_error);

        EXPECT_THROW(
            assemble_laplacian_csr(
                mesh, A, b, constrained,
                std::vector<double>{std::numeric_limits<double>::quiet_NaN()}),
            std::runtime_error);
    });

    return run_all();
}
