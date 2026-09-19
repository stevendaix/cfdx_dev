// M0.7-T05 — Tests for surface / volume integration

#include "cfdx/core/numerics/integrate.h"
#include "cfdx/core/mesh/mesh.h"
#include "test_harness.h"

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
    run_case("surface_integrate_dimension3", []() {
        Mesh m = make_unit_cube();
        ScalarFaceField phi(6, "phi", "m^3/s", 1);
        phi.fill(0.0);
        auto result = surface_integrate(phi, m);
        EXPECT_TRUE(result.size() == 1);
        EXPECT_TRUE(result.dimension() == 3);
    });

    run_case("surface_integrate_zero", []() {
        Mesh m = make_unit_cube();
        ScalarFaceField phi(6, "phi", "m^3/s", 1);
        phi.fill(0.0);
        auto result = surface_integrate(phi, m);
        double x, y, z;
        result.get(0, x, y, z);
        EXPECT_TRUE(x == 0.0);
        EXPECT_TRUE(y == 0.0);
        EXPECT_TRUE(z == 0.0);
    });

    run_case("surface_integrate_constant", []() {
        // φ_f = 1 pour toutes les faces → Σ Sf_f = 0 (fermeture).
        Mesh m = make_unit_cube();
        ScalarFaceField phi(6, "phi", "m^3/s", 1);
        phi.fill(1.0);
        auto result = surface_integrate(phi, m);
        double x, y, z;
        result.get(0, x, y, z);
        EXPECT_NEAR(x, 0.0, 1e-12);
        EXPECT_NEAR(y, 0.0, 1e-12);
        EXPECT_NEAR(z, 0.0, 1e-12);
    });

    run_case("surface_integrate_size_mismatch", []() {
        Mesh m = make_unit_cube();
        ScalarFaceField phi(5, "phi", "m^3/s", 1);
        EXPECT_THROW(surface_integrate(phi, m), std::runtime_error);
    });

    run_case("surface_integrate_metadata", []() {
        Mesh m = make_unit_cube();
        ScalarFaceField phi(6, "phi", "m^3/s", 1);
        phi.fill(0.0);
        auto result = surface_integrate(phi, m);
        EXPECT_TRUE(result.name() == "phi_surf_int");
        EXPECT_TRUE(result.loc() == Location::CELL);
        EXPECT_TRUE(result.metadata().unit == "m^3/s");
    });

    run_case("volume_integrate_unit_cube", []() {
        Mesh m = make_unit_cube();
        ScalarCellField f(1, "p", "Pa", 1);
        f(0) = 1.0;
        double vol = volume_integrate(f, m);
        // Volume du cube unité = 1.0.
        EXPECT_NEAR(vol, 1.0, 1e-12);
    });

    run_case("volume_integrate_constant", []() {
        Mesh m = make_unit_cube();
        ScalarCellField f(1, "p", "Pa", 1);
        f(0) = 3.0;
        double vol = volume_integrate(f, m);
        EXPECT_NEAR(vol, 3.0, 1e-12);
    });

    run_case("volume_integrate_size_mismatch", []() {
        Mesh m = make_unit_cube();
        ScalarCellField f(2, "p", "Pa", 1);
        EXPECT_THROW(volume_integrate(f, m), std::runtime_error);
    });

    run_case("volume_integrate_vector_rejected", []() {
        Mesh m = make_unit_cube();
        Field<double, Location::CELL> f(1, "U", "m/s", 3);
        f.set(0, 1.0, 2.0, 3.0);
        EXPECT_THROW(volume_integrate(f, m), std::runtime_error);
    });

    run_case("volume_integrate_empty_mesh", []() {
        Mesh m;
        ScalarCellField f(0, "p", "Pa", 1);
        double vol = volume_integrate(f, m);
        EXPECT_TRUE(vol == 0.0);
    });

    return run_all();
}