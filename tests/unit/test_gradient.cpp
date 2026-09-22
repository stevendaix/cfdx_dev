// M0.7-T01 — Tests for Gauss gradient

#include "cfdx/core/numerics/gradient.h"
#include "cfdx/core/fvm/least_squares_gradient.h"
#include "cfdx/core/mesh/mesh.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

// Construit un cube unité [0,1]^3 — 1 cellule, 6 faces, 8 sommets.
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

    // 6 faces (CCW vues de l'extérieur)
    m.faces().push_face({0, 3, 2, 1});  // 0 : bas
    m.faces().push_face({4, 5, 6, 7});  // 1 : haut
    m.faces().push_face({0, 1, 5, 4});  // 2 : avant
    m.faces().push_face({3, 7, 6, 2});  // 3 : arrière
    m.faces().push_face({0, 4, 7, 3});  // 4 : gauche
    m.faces().push_face({1, 2, 6, 5});  // 5 : droite

    m.ownership().resize(6);
    for (std::size_t i = 0; i < 6; ++i) {
        m.ownership().set_owner(i, 0);
        m.ownership().set_neighbour(i, FaceOwnership::BOUNDARY);
    }

    m.cells().push_cell({0, 1, 2, 3, 4, 5});

    BoundaryPatches bp;
    const char* names[6] = {"bottom", "top", "front", "back", "left", "right"};
    const PatchType types[6] = {PatchType::WALL, PatchType::WALL, PatchType::INLET,
                                 PatchType::OUTLET, PatchType::SYMMETRY, PatchType::WALL};
    for (int i = 0; i < 6; ++i) {
        Patch p;
        p.name = names[i];
        p.type = types[i];
        p.face_ids = {static_cast<std::uint32_t>(i)};
        bp.add_patch(p);
    }
    m.set_boundary(bp);

    return m;
}

int main() {
    run_case("gradient_dimension3", []() {
        Mesh m = make_unit_cube();
        ScalarCellField f(1, "p", "Pa", 1);
        f(0) = 0.5;
        auto grad = compute_gradient_gauss(f, m);
        EXPECT_TRUE(grad.size() == 1);
        EXPECT_TRUE(grad.dimension() == 3);
    });

    run_case("gradient_constant_field", []() {
        // Un champ constant a un gradient nul (par symétrie géométrique).
        Mesh m = make_unit_cube();
        ScalarCellField f(1, "p", "Pa", 1);
        f(0) = 42.0;

        auto grad = compute_gradient_gauss(f, m);
        double gx, gy, gz;
        grad.get(0, gx, gy, gz);
        EXPECT_NEAR(gx, 0.0, 1e-12);
        EXPECT_NEAR(gy, 0.0, 1e-12);
        EXPECT_NEAR(gz, 0.0, 1e-12);
    });

    run_case("gradient_reuses_geometry_cache", []() {
        Mesh m = make_unit_cube();
        ScalarCellField f(1, "p", "Pa", 1);
        f(0) = 42.0;
        const GeometryCache geometry = make_geometry_cache(m);
        auto grad = compute_gradient_gauss(f, m, geometry);
        EXPECT_NEAR(grad(0, 0), 0.0, 1e-12);
        EXPECT_NEAR(grad(0, 1), 0.0, 1e-12);
        EXPECT_NEAR(grad(0, 2), 0.0, 1e-12);
    });

    run_case("gradient_gauss_linear_field_exact_on_interior_cell", []() {
        // Three Cartesian cells in x.  The middle cell has neighbours on
        // both x-faces, so the zero-gradient boundary contract cannot
        // contaminate the linear-field oracle.
        Mesh m;
        m.points().resize(16);
        const double p[16][3] = {
            {0,0,0},{1,0,0},{2,0,0},{3,0,0},
            {0,1,0},{1,1,0},{2,1,0},{3,1,0},
            {0,0,1},{1,0,1},{2,0,1},{3,0,1},
            {0,1,1},{1,1,1},{2,1,1},{3,1,1}
        };
        for (std::size_t i = 0; i < 16; ++i)
            m.points().set(i, p[i][0], p[i][1], p[i][2]);

        const std::vector<std::vector<Index>> faces = {
            {0,8,12,4}, {0,1,9,8}, {4,12,13,5}, {0,4,5,1},
            {8,9,13,12}, {9,10,14,13}, {1,5,6,2}, {1,2,10,9},
            {5,13,14,6}, {1,5,6,2}, {10,11,15,14}, {2,3,11,10},
            {6,14,15,7}, {2,6,7,3}, {11,15,14,10}
        };
        // Replace the duplicated face above with the six faces per cell:
        // use an explicit, compact 3-cell topology below.
        m = Mesh{};
        m.points().resize(16);
        for (std::size_t i = 0; i < 16; ++i)
            m.points().set(i, p[i][0], p[i][1], p[i][2]);

        const std::vector<std::vector<Index>> fs = {
            {0,8,12,4}, {0,1,5,4}, {8,9,13,12}, {1,2,6,5},
            {9,10,14,13}, {2,3,7,6}, {10,11,15,14},
            {0,4,12,8}, {1,9,13,5}, {2,10,14,6},
            {3,11,15,7}, {4,5,6,7}, {8,12,13,9},
            {9,13,14,10}, {10,14,15,11}
        };
        for (const auto& face : fs) m.faces().push_face(face);

        m.ownership().resize(fs.size());
        for (std::size_t f = 0; f < fs.size(); ++f) {
            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        }
        // Face ordering: left boundary, x-internal 0-1, x-internal 1-2,
        // right boundary, then y/z boundary faces.
        m.ownership().set_owner(0, 0);
        m.ownership().set_owner(1, 0);
        m.ownership().set_neighbour(1, 1);
        m.ownership().set_owner(2, 1);
        m.ownership().set_neighbour(2, 2);
        m.ownership().set_owner(3, 2);
        for (std::size_t f = 4; f < fs.size(); ++f) {
            // y/z faces are assigned to the corresponding cell below.
            const std::size_t owner = (f <= 5) ? 0 : (f <= 8 ? 1 : 2);
            m.ownership().set_owner(f, owner);
        }
        m.cells().push_cell({0,1,4,5,7,8});
        m.cells().push_cell({1,2,6,9,12,13});
        m.cells().push_cell({2,3,10,11,14,15});

        const auto geometry = make_geometry_cache(m);
        ScalarCellField f(3, "phi", "1", 1);
        for (std::size_t c = 0; c < 3; ++c)
            f(c) = geometry.cell_centres[c].x;

        const auto grad = compute_gradient_gauss(f, m, geometry);
        EXPECT_NEAR(grad(1, 0), 1.0, 1e-12);
        EXPECT_NEAR(grad(1, 1), 0.0, 1e-12);
        EXPECT_NEAR(grad(1, 2), 0.0, 1e-12);
    });

    run_case("gradient_gauss_quadratic_field_exact_on_interior_cell", []() {
        Mesh m;
        m.points().resize(16);
        const double p[16][3] = {
            {0,0,0},{1,0,0},{2,0,0},{3,0,0},
            {0,1,0},{1,1,0},{2,1,0},{3,1,0},
            {0,0,1},{1,0,1},{2,0,1},{3,0,1},
            {0,1,1},{1,1,1},{2,1,1},{3,1,1}
        };
        for (std::size_t i = 0; i < 16; ++i)
            m.points().set(i, p[i][0], p[i][1], p[i][2]);

        const std::vector<std::vector<Index>> fs = {
            {0,8,12,4}, {0,1,5,4}, {8,9,13,12}, {1,2,6,5},
            {9,10,14,13}, {2,3,7,6}, {10,11,15,14},
            {0,4,12,8}, {1,9,13,5}, {2,10,14,6},
            {3,11,15,7}, {4,5,6,7}, {8,12,13,9},
            {9,13,14,10}, {10,14,15,11}
        };
        for (const auto& face : fs) m.faces().push_face(face);
        m.ownership().resize(fs.size());
        for (std::size_t face = 0; face < fs.size(); ++face)
            m.ownership().set_neighbour(face, FaceOwnership::BOUNDARY);
        m.ownership().set_owner(0, 0);
        m.ownership().set_owner(1, 0); m.ownership().set_neighbour(1, 1);
        m.ownership().set_owner(2, 1); m.ownership().set_neighbour(2, 2);
        m.ownership().set_owner(3, 2);
        for (std::size_t face = 4; face < fs.size(); ++face)
            m.ownership().set_owner(face, face <= 5 ? 0 : (face <= 8 ? 1 : 2));
        m.cells().push_cell({0,1,4,5,7,8});
        m.cells().push_cell({1,2,6,9,12,13});
        m.cells().push_cell({2,3,10,11,14,15});

        const auto geometry = make_geometry_cache(m);
        ScalarCellField f(3, "phi", "1", 1);
        for (std::size_t c = 0; c < 3; ++c) {
            const double x = geometry.cell_centres[c].x;
            f(c) = x * x;
        }

        const auto grad = compute_gradient_gauss(f, m, geometry);
        // Middle-cell centre is x=1.5, hence d(x^2)/dx = 3.
        EXPECT_NEAR(grad(1, 0), 3.0, 1e-12);
        EXPECT_NEAR(grad(1, 1), 0.0, 1e-12);
        EXPECT_NEAR(grad(1, 2), 0.0, 1e-12);
    });

    run_case("gradient_size_mismatch", []() {
        Mesh m = make_unit_cube();
        ScalarCellField f(2, "p", "Pa", 1);
        EXPECT_THROW(compute_gradient_gauss(f, m), std::runtime_error);
    });

    run_case("gradient_vector_field_rejected", []() {
        Mesh m = make_unit_cube();
        Field<double, Location::CELL> f(1, "U", "m/s", 3);
        f.set(0, 1.0, 2.0, 3.0);
        EXPECT_THROW(compute_gradient_gauss(f, m), std::runtime_error);
    });

    run_case("gradient_metadata", []() {
        Mesh m = make_unit_cube();
        ScalarCellField f(1, "p", "Pa", 1);
        f(0) = 0.5;
        auto grad = compute_gradient_gauss(f, m);
        EXPECT_TRUE(grad.name() == "p_grad");
        EXPECT_TRUE(grad.loc() == Location::CELL);
        EXPECT_TRUE(grad.metadata().unit == "Pa/m");
    });

    run_case("least_squares_exact_linear_2d", []() {
        const Vec3 centre{0.2, -0.3, 0.0};
        const Vec3 expected{2.0, -3.0, 0.0};
        const std::vector<Vec3> neighbours = {
            {1.1, -0.3, 0.0},
            {-0.4, 0.2, 0.0},
            {0.7, 1.4, 0.0},
            {-1.0, -1.1, 0.0}
        };
        std::vector<double> values;
        for (const auto& p : neighbours)
            values.push_back(1.7 + expected.x * p.x + expected.y * p.y);
        const double centre_value =
            1.7 + expected.x * centre.x + expected.y * centre.y;
        const Vec3 grad = least_squares_gradient(
            centre, centre_value, neighbours, values);
        EXPECT_NEAR(grad.x, expected.x, 1e-12);
        EXPECT_NEAR(grad.y, expected.y, 1e-12);
        EXPECT_NEAR(grad.z, 0.0, 1e-12);
    });

    run_case("least_squares_exact_linear_3d", []() {
        const Vec3 centre{0.1, -0.2, 0.3};
        const Vec3 expected{1.5, -2.0, 0.75};
        const std::vector<Vec3> neighbours = {
            {0.9, -0.1, 0.5},
            {-0.6, 0.7, 0.2},
            {0.4, -1.1, 1.2},
            {-0.8, -0.9, -0.4},
            {1.2, 0.4, -0.7}
        };
        std::vector<double> values;
        for (const auto& p : neighbours)
            values.push_back(-0.4 + expected.x * p.x +
                             expected.y * p.y + expected.z * p.z);
        const double centre_value =
            -0.4 + expected.x * centre.x +
            expected.y * centre.y + expected.z * centre.z;
        const Vec3 grad = least_squares_gradient(
            centre, centre_value, neighbours, values);
        EXPECT_NEAR(grad.x, expected.x, 1e-12);
        EXPECT_NEAR(grad.y, expected.y, 1e-12);
        EXPECT_NEAR(grad.z, expected.z, 1e-12);
    });

    run_case("least_squares_constant_field", []() {
        const Vec3 centre{0.0, 0.0, 0.0};
        const std::vector<Vec3> neighbours = {
            {1.0, 0.2, 0.0}, {-0.4, 0.9, 0.0},
            {0.3, -0.7, 0.0}, {-0.8, -0.2, 0.0}
        };
        const std::vector<double> values(neighbours.size(), 42.0);
        const Vec3 grad = least_squares_gradient(
            centre, 42.0, neighbours, values);
        EXPECT_NEAR(grad.x, 0.0, 1e-12);
        EXPECT_NEAR(grad.y, 0.0, 1e-12);
        EXPECT_NEAR(grad.z, 0.0, 1e-12);
    });

    run_case("least_squares_rank_deficient_2d", []() {
        const Vec3 centre{0.0, 0.0, 0.0};
        const std::vector<Vec3> neighbours = {
            {1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}
        };
        const std::vector<double> values = {3.0, -1.0};
        const Vec3 grad = least_squares_gradient(
            centre, 1.0, neighbours, values);
        EXPECT_NEAR(grad.x, 2.0, 1e-12);
        EXPECT_NEAR(grad.y, 0.0, 1e-12);
        EXPECT_NEAR(grad.z, 0.0, 1e-12);
    });

    run_case("gradient_empty_mesh", []() {
        Mesh m;
        ScalarCellField f(0, "p", "Pa", 1);
        auto grad = compute_gradient_gauss(f, m);
        EXPECT_TRUE(grad.size() == 0);
    });

    return run_all();
}