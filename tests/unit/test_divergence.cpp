// M0.7-T02 — Tests for divergence

#include "cfdx/core/numerics/divergence.h"
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