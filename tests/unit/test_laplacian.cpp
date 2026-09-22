// M0.7-T03 — Tests for Laplacian

#include "cfdx/core/numerics/laplacian.h"
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

Mesh make_two_cell_unit_cubes() {
    Mesh m;
    m.points().resize(12);
    const double p[12][3] = {
        {0,0,0},{1,0,0},{2,0,0},{0,1,0},{1,1,0},{2,1,0},
        {0,0,1},{1,0,1},{2,0,1},{0,1,1},{1,1,1},{2,1,1}};
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i, p[i][0], p[i][1], p[i][2]);

    m.faces().push_face({0,6,9,3});
    m.faces().push_face({0,1,7,6});
    m.faces().push_face({3,9,10,4});
    m.faces().push_face({0,3,4,1});
    m.faces().push_face({6,7,10,9});
    m.faces().push_face({7,10,4,1});
    m.faces().push_face({2,5,11,8});
    m.faces().push_face({1,2,8,7});
    m.faces().push_face({4,10,11,5});
    m.faces().push_face({1,4,5,2});
    m.faces().push_face({7,8,11,10});

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 5; ++f) {
        m.ownership().set_owner(f, 0);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(5, 0);
    m.ownership().set_neighbour(5, 1);
    for (std::size_t f = 6; f < 11; ++f) {
        m.ownership().set_owner(f, 1);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    m.cells().push_cell({5,6,7,8,9,10});
    return m;
}

int main() {
    run_case("laplacian_dimension1", []() {
        Mesh m = make_unit_cube();
        ScalarCellField f(1, "p", "Pa", 1);
        f(0) = 0.5;
        auto lap = compute_laplacian(f, m);
        EXPECT_TRUE(lap.size() == 1);
        EXPECT_TRUE(lap.dimension() == 1);
    });

    run_case("laplacian_constant_field", []() {
        // Un champ constant a un gradient nul, donc un laplacien nul.
        Mesh m = make_unit_cube();
        ScalarCellField f(1, "p", "Pa", 1);
        f(0) = 42.0;

        auto lap = compute_laplacian(f, m);
        EXPECT_NEAR(lap(0), 0.0, 1e-12);
    });

    run_case("laplacian_size_mismatch", []() {
        Mesh m = make_unit_cube();
        ScalarCellField f(2, "p", "Pa", 1);
        EXPECT_THROW(compute_laplacian(f, m), std::runtime_error);
    });

    run_case("laplacian_vector_field_rejected", []() {
        Mesh m = make_unit_cube();
        Field<double, Location::CELL> f(1, "U", "m/s", 3);
        f.set(0, 1.0, 2.0, 3.0);
        EXPECT_THROW(compute_laplacian(f, m), std::runtime_error);
    });

    run_case("laplacian_metadata", []() {
        Mesh m = make_unit_cube();
        ScalarCellField f(1, "p", "Pa", 1);
        f(0) = 0.5;
        auto lap = compute_laplacian(f, m);
        EXPECT_TRUE(lap.name() == "p_lap");
        EXPECT_TRUE(lap.loc() == Location::CELL);
        EXPECT_TRUE(lap.metadata().unit == "Pa/m^2");
    });

    run_case("laplacian_empty_mesh", []() {
        Mesh m;
        ScalarCellField f(0, "p", "Pa", 1);
        auto lap = compute_laplacian(f, m);
        EXPECT_TRUE(lap.size() == 0);
    });

    run_case("laplacian_corrected_matches_orthogonal_on_orthogonal_mesh", []() {
        Mesh m = make_two_cell_unit_cubes();
        ScalarCellField f(2, "p", "Pa", 1);
        f(0) = 0.5;
        f(1) = 1.5;
        auto orth = compute_laplacian(f, m, LaplacianScheme::ORTHOGONAL);
        auto corrected = compute_laplacian(f, m, LaplacianScheme::CORRECTED);
        EXPECT_NEAR(corrected(0), orth(0), 1e-12);
        EXPECT_NEAR(corrected(1), orth(1), 1e-12);
    });

    run_case("laplacian_corrected_preserves_internal_flux_on_skewed_face", []() {
        Mesh m = make_two_cell_unit_cubes();
        ScalarCellField f(2, "p", "Pa", 1);
        f(0) = 0.5;
        f(1) = 1.5;

        GeometryCache geometry = make_geometry_cache(m);
        geometry.face_Sf[5].y = 0.25;

        auto corrected = compute_laplacian(f, m, geometry, LaplacianScheme::CORRECTED);
        EXPECT_NEAR(corrected(0) + corrected(1), 0.0, 1e-12);
        EXPECT_TRUE(std::isfinite(corrected(0)));
        EXPECT_TRUE(std::isfinite(corrected(1)));
    });

    run_case("laplacian_limited_scheme_is_explicitly_unsupported", []() {
        Mesh m = make_unit_cube();
        ScalarCellField f(1, "p", "Pa", 1);
        f(0) = 42.0;
        EXPECT_THROW(compute_laplacian(f, m, LaplacianScheme::LIMITED), std::runtime_error);
    });

    run_case("laplacian_reuses_geometry_cache", []() {
        Mesh m = make_unit_cube();
        ScalarCellField f(1, "p", "Pa", 1);
        f(0) = 42.0;
        const GeometryCache geometry = make_geometry_cache(m);
        auto lap = compute_laplacian(f, m, geometry);
        EXPECT_NEAR(lap(0), 0.0, 1e-12);
    });

    run_case("laplacian_scheme_from_string", []() {
        EXPECT_TRUE(laplacian_scheme_from_string("orthogonal") == LaplacianScheme::ORTHOGONAL);
        EXPECT_TRUE(laplacian_scheme_from_string("corrected") == LaplacianScheme::CORRECTED);
        EXPECT_TRUE(laplacian_scheme_from_string("limited") == LaplacianScheme::LIMITED);
        EXPECT_TRUE(laplacian_scheme_from_string("uncorrected") == LaplacianScheme::UNCORRECTED);
        EXPECT_THROW(laplacian_scheme_from_string("bogus"), std::runtime_error);
    });

    run_case("laplacian_linear_field_two_cells", []() {
        Mesh m = make_two_cell_unit_cubes();
        ScalarCellField f(2, "p", "Pa", 1);
        // phi = x at the two cell centres: 0.5 and 1.5.
        f(0) = 0.5;
        f(1) = 1.5;
        auto lap = compute_laplacian(f, m, LaplacianScheme::ORTHOGONAL);
        // The shared face conductance is A/d = 1. With zero-gradient
        // boundaries, the two rows receive equal and opposite fluxes.
        EXPECT_NEAR(lap(0), 1.0, 1e-12);
        EXPECT_NEAR(lap(1), -1.0, 1e-12);
    });

    return run_all();
}