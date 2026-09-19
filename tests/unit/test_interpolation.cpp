// M0.6-T01 — Tests for cell→face interpolation

#include "cfdx/core/numerics/interpolation.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/ownership.h"
#include "test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

// Construit un maillage 1D : 2 cellules, 3 faces (2 internes + 2 frontières = 4 faces au total)
//   face 0 : frontière gauche, owner=0
//   face 1 : interne, owner=0, neighbour=1
//   face 2 : interne, owner=1, neighbour=0
//   face 3 : frontière droite, owner=1
Mesh make_two_cell_mesh() {
    Mesh m;

    m.points().resize(4);
    m.points().set(0, 0.0, 0.0, 0.0);
    m.points().set(1, 1.0, 0.0, 0.0);
    m.points().set(2, 2.0, 0.0, 0.0);
    m.points().set(3, 3.0, 0.0, 0.0);

    // 4 faces (segment 2 pts)
    m.faces().push_face({0, 1});
    m.faces().push_face({1, 2});
    m.faces().push_face({2, 1});
    m.faces().push_face({2, 3});

    m.ownership().resize(4);
    m.ownership().set_owner(0, 0);
    m.ownership().set_neighbour(0, FaceOwnership::BOUNDARY);
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);
    m.ownership().set_owner(2, 1);
    m.ownership().set_neighbour(2, 0);
    m.ownership().set_owner(3, 1);
    m.ownership().set_neighbour(3, FaceOwnership::BOUNDARY);

    m.cells().push_cell({0});
    m.cells().push_cell({1, 2, 3});

    BoundaryPatches bp;
    Patch left, right;
    left.name = "left";  left.type = PatchType::INLET;  left.face_ids = {0};
    right.name = "right"; right.type = PatchType::OUTLET; right.face_ids = {3};
    bp.add_patch(left);
    bp.add_patch(right);
    m.set_boundary(bp);

    return m;
}

int main() {
    run_case("interpolate_upwind_scalar", []() {
        Mesh m = make_two_cell_mesh();
        ScalarCellField f(2, "p", "Pa", 1);
        f(0) = 10.0;
        f(1) = 20.0;

        auto face_field = interpolate_cell_to_face(f, m, InterpScheme::UPWIND);
        EXPECT_TRUE(face_field.size() == 4);
        EXPECT_TRUE(face_field.dimension() == 1);
        // Upwind = owner partout
        EXPECT_TRUE(face_field(0) == 10.0);
        EXPECT_TRUE(face_field(1) == 10.0);
        EXPECT_TRUE(face_field(2) == 20.0);
        EXPECT_TRUE(face_field(3) == 20.0);
    });

    run_case("interpolate_linear_scalar", []() {
        Mesh m = make_two_cell_mesh();
        ScalarCellField f(2, "p", "Pa", 1);
        f(0) = 10.0;
        f(1) = 30.0;

        auto face_field = interpolate_cell_to_face(f, m, InterpScheme::LINEAR);
        // face 0 : owner=0, boundary → 10.0
        EXPECT_TRUE(face_field(0) == 10.0);
        // face 1 : 0.5*(10+30) = 20.0
        EXPECT_TRUE(face_field(1) == 20.0);
        // face 2 : 0.5*(30+10) = 20.0
        EXPECT_TRUE(face_field(2) == 20.0);
        // face 3 : owner=1, boundary → 30.0
        EXPECT_TRUE(face_field(3) == 30.0);
    });

    run_case("interpolate_limited_scalar", []() {
        Mesh m = make_two_cell_mesh();
        ScalarCellField f(2, "p", "Pa", 1);
        f(0) = 10.0;
        f(1) = 30.0;

        auto face_field = interpolate_cell_to_face(f, m, InterpScheme::LIMITED);
        // Limited = clamp(linear, min, max) = linear ici
        EXPECT_TRUE(face_field(0) == 10.0);
        EXPECT_TRUE(face_field(1) == 20.0);
        EXPECT_TRUE(face_field(2) == 20.0);
        EXPECT_TRUE(face_field(3) == 30.0);
    });

    run_case("interpolate_limited_clamps", []() {
        // Limited doit borner la valeur linéaire dans [min, max].
        // On utilise un champ scalaire avec des valeurs extrêmes.
        Mesh m = make_two_cell_mesh();
        ScalarCellField f(2, "p", "Pa", 1);
        f(0) = 0.0;
        f(1) = 100.0;

        auto face_field = interpolate_cell_to_face(f, m, InterpScheme::LIMITED);
        // linear = 50.0, borné par [0, 100] → 50.0
        EXPECT_TRUE(face_field(1) == 50.0);
        EXPECT_TRUE(face_field(2) == 50.0);
    });

    run_case("interpolate_vector_dim3", []() {
        Mesh m = make_two_cell_mesh();
        Field<double, Location::CELL> f(2, "U", "m/s", 3);
        f.set(0, 1.0, 2.0, 3.0);
        f.set(1, 7.0, 8.0, 9.0);

        auto face_field = interpolate_cell_to_face(f, m, InterpScheme::LINEAR);
        EXPECT_TRUE(face_field.dimension() == 3);
        double x, y, z;
        face_field.get(1, x, y, z);
        EXPECT_TRUE(x == 4.0);
        EXPECT_TRUE(y == 5.0);
        EXPECT_TRUE(z == 6.0);
    });

    run_case("interpolate_metadata", []() {
        Mesh m = make_two_cell_mesh();
        ScalarCellField f(2, "p", "Pa", 1);
        f(0) = 1.0; f(1) = 2.0;

        auto face_field = interpolate_cell_to_face(f, m, InterpScheme::LINEAR);
        EXPECT_TRUE(face_field.name() == "p");
        EXPECT_TRUE(face_field.loc() == Location::FACE);
        EXPECT_TRUE(face_field.metadata().unit == "Pa");
    });

    run_case("interp_scheme_from_string", []() {
        EXPECT_TRUE(interp_scheme_from_string("linear") == InterpScheme::LINEAR);
        EXPECT_TRUE(interp_scheme_from_string("upwind") == InterpScheme::UPWIND);
        EXPECT_TRUE(interp_scheme_from_string("limited") == InterpScheme::LIMITED);
        EXPECT_THROW(interp_scheme_from_string("bogus"), std::runtime_error);
    });

    run_case("interpolate_empty_mesh", []() {
        Mesh m;
        ScalarCellField f(0, "p", "Pa", 1);
        auto face_field = interpolate_cell_to_face(f, m, InterpScheme::UPWIND);
        EXPECT_TRUE(face_field.size() == 0);
    });

    run_case("interpolate_owner_out_of_range", []() {
        Mesh m = make_two_cell_mesh();
        ScalarCellField f(2, "p", "Pa", 1);
        m.ownership().set_owner(0, 99);  // owner hors de la portée
        EXPECT_THROW(interpolate_cell_to_face(f, m, InterpScheme::UPWIND),
                     std::runtime_error);
    });

    return run_all();
}