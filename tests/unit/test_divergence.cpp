// M0.7-T02 — Tests for divergence

#include "cfdx/core/numerics/divergence.h"
#include "cfdx/core/numerics/flux.h"
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
    run_case("divergence_dimension1", []() {
        Mesh m = make_unit_cube();
        ScalarFaceField phi(6, "phi", "m^3/s", 1);
        phi.fill(0.0);
        auto div = compute_divergence(phi, m);
        EXPECT_TRUE(div.size() == 1);
        EXPECT_TRUE(div.dimension() == 1);
    });

    run_case("divergence_zero_flux", []() {
        // Un flux nul donne une divergence nulle.
        Mesh m = make_unit_cube();
        ScalarFaceField phi(6, "phi", "m^3/s", 1);
        phi.fill(0.0);
        auto div = compute_divergence(phi, m);
        EXPECT_TRUE(div(0) == 0.0);
    });

    run_case("divergence_constant_flux", []() {
        // Un flux constant non nul donne une divergence non nulle
        // (la somme des contributions n'est pas nulle pour un champ constant).
        Mesh m = make_unit_cube();
        ScalarFaceField phi(6, "phi", "m^3/s", 1);
        phi.fill(1.0);
        auto div = compute_divergence(phi, m);
        // Toutes les faces appartiennent à la cellule 0 (owner).
        // div = Σ phi_f = 6 * 1.0 = 6.0
        EXPECT_TRUE(div(0) == 6.0);
    });

    run_case("divergence_reuses_geometry_cache", []() {
        Mesh m = make_unit_cube();
        ScalarFaceField phi(6, "phi", "m^3/s", 1);
        phi.fill(0.0);
        const GeometryCache geometry = make_geometry_cache(m);
        auto div = compute_divergence(phi, m, geometry);
        EXPECT_NEAR(div(0), 0.0, 1e-12);
    });

    run_case("divergence_uniform_vector_field_is_zero", []() {
        Mesh m = make_unit_cube();
        Field<double, Location::FACE> U(m.n_faces(), "U", "m/s", 3);
        U.fill(0.0);
        U.component_data(0)[0] = 1.0;
        U.component_data(0)[1] = 1.0;
        U.component_data(0)[2] = 1.0;
        U.component_data(0)[3] = 1.0;
        U.component_data(0)[4] = 1.0;
        U.component_data(0)[5] = 1.0;

        auto phi = compute_flux(U, m);
        auto div = compute_divergence(phi, m);
        // A constant vector field has zero analytical divergence.
        EXPECT_NEAR(div(0), 0.0, 1e-12);
    });

    run_case("divergence_internal_flux_is_conservative", []() {
        // The same internal-face flux must enter one cell and leave the other.
        // This is a discrete conservation invariant independent of the
        // particular cell values.
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

        ScalarFaceField phi(m.n_faces(), "phi", "m^3/s", 1);
        phi.fill(0.0);
        phi(5) = 2.0;
        auto div = compute_divergence(phi, m);
        EXPECT_NEAR(div(0), 2.0, 1e-12);
        EXPECT_NEAR(div(1), -2.0, 1e-12);
    });

    run_case("divergence_size_mismatch", []() {
        Mesh m = make_unit_cube();
        ScalarFaceField phi(5, "phi", "m^3/s", 1);
        EXPECT_THROW(compute_divergence(phi, m), std::runtime_error);
    });

    run_case("divergence_vector_field_rejected", []() {
        Mesh m = make_unit_cube();
        Field<double, Location::FACE> phi(6, "U", "m/s", 3);
        EXPECT_THROW(compute_divergence(phi, m), std::runtime_error);
    });

    run_case("divergence_metadata", []() {
        Mesh m = make_unit_cube();
        ScalarFaceField phi(6, "phi", "m^3/s", 1);
        phi.fill(0.0);
        auto div = compute_divergence(phi, m);
        EXPECT_TRUE(div.name() == "phi_div");
        EXPECT_TRUE(div.loc() == Location::CELL);
        EXPECT_TRUE(div.metadata().unit == "1/s");
    });

    run_case("divergence_empty_mesh", []() {
        Mesh m;
        ScalarFaceField phi(0, "phi", "m^3/s", 1);
        auto div = compute_divergence(phi, m);
        EXPECT_TRUE(div.size() == 0);
    });

    return run_all();
}