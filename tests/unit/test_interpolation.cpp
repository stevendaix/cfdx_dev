// M0.6-T01 — Tests for cell→face interpolation

#include "cfdx/core/numerics/interpolation.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/ownership.h"
#include "cfdx/core/numerics/flux.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

// Construit un maillage 1D : 2 cellules, 4 faces (2 internes + 2 frontières)
//   face 0 : frontière gauche, owner=0, vertices {0, 1}
//   face 1 : interne, owner=0, neighbour=1, vertices {1, 2}
//   face 2 : interne, owner=1, neighbour=0, vertices {2, 1}
//   face 3 : frontière droite, owner=1, vertices {2, 3}
// But we need at least 3 vertices per face for 2D/3D geometry.
// So let's use triangular faces with z=0 (flat in XY plane).
Mesh make_two_cell_mesh() {
    Mesh m;

    // 4 points forming two triangles
    m.points().resize(6);
    m.points().set(0, 0.0, 0.0, 0.0);
    m.points().set(1, 1.0, 0.0, 0.0);
    m.points().set(2, 2.0, 0.0, 0.0);
    m.points().set(3, 0.0, 1.0, 0.0);
    m.points().set(4, 1.0, 1.0, 0.0);
    m.points().set(5, 2.0, 1.0, 0.0);

    // 4 faces (triangles - 3 vertices each)
    // Face 0: left boundary triangle
    m.faces().push_face({0, 3, 1});
    // Face 1: internal face between cell 0 and 1
    m.faces().push_face({1, 4, 2});
    // Face 2: internal face between cell 1 and 0 (reversed)
    m.faces().push_face({2, 4, 1});
    // Face 3: right boundary triangle
    m.faces().push_face({2, 5, 4});

    m.ownership().resize(4);
    m.ownership().set_owner(0, 0);
    m.ownership().set_neighbour(0, FaceOwnership::BOUNDARY);
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);
    m.ownership().set_owner(2, 1);
    m.ownership().set_neighbour(2, 0);
    m.ownership().set_owner(3, 1);
    m.ownership().set_neighbour(3, FaceOwnership::BOUNDARY);

    // Cell 0: faces 0, 1
    m.cells().push_cell({0, 1});
    // Cell 1: faces 2, 3
    m.cells().push_cell({2, 3});

    BoundaryPatches bp;
    Patch left, right;
    left.name = "left";  left.type = PatchType::INLET;  left.face_ids = {0};
    right.name = "right"; right.type = PatchType::OUTLET; right.face_ids = {3};
    bp.add_patch(left);
    bp.add_patch(right);
    m.set_boundary(bp);

    return m;
}

// Helper to create a simple face velocity field (constant 1.0 in x direction)
Field<double, Location::FACE> make_constant_velocity(const Mesh& m) {
    Field<double, Location::FACE> U(m.n_faces(), "U", "m/s", 3);
    double* ux = U.component_data(0);
    double* uy = U.component_data(1);
    double* uz = U.component_data(2);
    for (std::size_t f = 0; f < m.n_faces(); ++f) {
        ux[f] = 1.0;  // Ux = 1.0
        uy[f] = 0.0;
        uz[f] = 0.0;
    }
    return U;
}

int main() {
    run_case("interpolate_upwind_scalar", []() {
        Mesh m = make_two_cell_mesh();
        ScalarCellField f(2, "p", "Pa", 1);
        f(0) = 10.0;
        f(1) = 20.0;

        // Create face flux for UPWIND scheme
        auto U = make_constant_velocity(m);
        auto phi = compute_flux(U, m);

        auto face_field = interpolate_cell_to_face(f, m, InterpScheme::UPWIND, &phi);
        EXPECT_TRUE(face_field.size() == 4);
        EXPECT_TRUE(face_field.dimension() == 1);
        // Upwind = owner partout (flux > 0 means outgoing from owner)
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

    run_case("interpolate_limited_requires_gradient_and_flux", []() {
        Mesh m = make_two_cell_mesh();
        ScalarCellField f(2, "p", "Pa", 1);
        f(0) = 10.0;
        f(1) = 30.0;
        EXPECT_THROW(interpolate_cell_to_face(f, m, InterpScheme::LIMITED), std::runtime_error);
    });

    run_case("apply_limiter_bounds_local_extrema", []() {
        EXPECT_NEAR(apply_limiter_tvd(0.0, 100.0, 150.0, LimiterType::NONE), 100.0, 1e-12);
        EXPECT_NEAR(apply_limiter_tvd(100.0, 0.0, -50.0, LimiterType::NONE), 0.0, 1e-12);
        EXPECT_TRUE(apply_limiter_tvd(0.0, 100.0, 50.0, LimiterType::MINMOD) >= 0.0);
        EXPECT_TRUE(apply_limiter_tvd(0.0, 100.0, 50.0, LimiterType::VANLEER) <= 100.0);
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
        // Use LINEAR for empty mesh (no faces, no flux needed)
        auto face_field = interpolate_cell_to_face(f, m, InterpScheme::LINEAR);
        EXPECT_TRUE(face_field.size() == 0);
    });

    run_case("interpolate_owner_out_of_range", []() {
        Mesh m = make_two_cell_mesh();
        ScalarCellField f(2, "p", "Pa", 1);
        m.ownership().set_owner(0, 99);  // owner hors de la portée
        EXPECT_THROW(interpolate_cell_to_face(f, m, InterpScheme::LINEAR),
                     std::runtime_error);
    });

    return run_all();
}
