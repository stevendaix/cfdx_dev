// M0.7-T04 — Tests for flux

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
    run_case("flux_dimension1", []() {
        Mesh m = make_unit_cube();
        Field<double, Location::FACE> U(6, "U", "m/s", 3);
        U.fill(0.0);
        auto phi = compute_flux(U, m);
        EXPECT_TRUE(phi.size() == 6);
        EXPECT_TRUE(phi.dimension() == 1);
    });

    run_case("flux_zero_velocity", []() {
        Mesh m = make_unit_cube();
        Field<double, Location::FACE> U(6, "U", "m/s", 3);
        U.fill(0.0);
        auto phi = compute_flux(U, m);
        for (std::size_t f = 0; f < 6; ++f) {
            EXPECT_TRUE(phi(f) == 0.0);
        }
    });

    run_case("flux_uniform_x", []() {
        // U = (1, 0, 0) → flux = U · Sf = Sf_x.
        // Pour un cube unité, les faces x=0 et x=1 ont Sf_x = ±1.
        Mesh m = make_unit_cube();
        Field<double, Location::FACE> U(6, "U", "m/s", 3);
        for (std::size_t f = 0; f < 6; ++f) {
            U.set(f, 1.0, 0.0, 0.0);
        }
        auto phi = compute_flux(U, m);
        // Face 4 (gauche, x=0) : Sf_x = -1 → flux = -1.
        // Face 5 (droite, x=1) : Sf_x = +1 → flux = +1.
        // Les autres faces ont Sf_x = 0 → flux = 0.
        EXPECT_TRUE(phi(4) == -1.0);
        EXPECT_TRUE(phi(5) == 1.0);
        EXPECT_TRUE(phi(0) == 0.0);
        EXPECT_TRUE(phi(1) == 0.0);
    });

    run_case("flux_size_mismatch", []() {
        Mesh m = make_unit_cube();
        Field<double, Location::FACE> U(5, "U", "m/s", 3);
        EXPECT_THROW(compute_flux(U, m), std::runtime_error);
    });

    run_case("flux_wrong_dimension", []() {
        Mesh m = make_unit_cube();
        Field<double, Location::FACE> U(6, "U", "m/s", 2);  // dim=2
        EXPECT_THROW(compute_flux(U, m), std::runtime_error);
    });

    run_case("flux_metadata", []() {
        Mesh m = make_unit_cube();
        Field<double, Location::FACE> U(6, "U", "m/s", 3);
        U.fill(0.0);
        auto phi = compute_flux(U, m);
        EXPECT_TRUE(phi.name() == "U_flux");
        EXPECT_TRUE(phi.loc() == Location::FACE);
        EXPECT_TRUE(phi.metadata().unit == "m^3/s");
    });

    run_case("flux_empty_mesh", []() {
        Mesh m;
        Field<double, Location::FACE> U(0, "U", "m/s", 3);
        auto phi = compute_flux(U, m);
        EXPECT_TRUE(phi.size() == 0);
    });

    return run_all();
}