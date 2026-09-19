// M0.1-T06 — Tests for Mesh topology + validator

#include "cfdx/core/mesh/mesh.h"
#include "test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

// Construit un maillage 2D simple de 2 cellules quadrilatères
// partageant une face interne, avec des faces de frontière.
Mesh make_simple_mesh() {
    Mesh m;

    // 4 points : (0,0), (1,0), (1,1), (0,1)
    m.points().resize(4);
    m.points().set(0, 0.0, 0.0, 0.0);
    m.points().set(1, 1.0, 0.0, 0.0);
    m.points().set(2, 1.0, 1.0, 0.0);
    m.points().set(3, 0.0, 1.0, 0.0);

    // 5 faces :
    //   0 : (0,1,2,3)  bas    (face de frontière, owner=0)
    //   1 : (1,2,3,0)  haut   (face de frontière, owner=1)
    //   2 : (0,3,2,1)  interne (owner=0, neighbour=1)
    //   3 : (0,1)      gauche (face de frontière, owner=0)
    //   4 : (2,3)      droite  (face de frontière, owner=1)
    m.faces().push_face({0, 1, 2, 3});  // face 0
    m.faces().push_face({1, 2, 3, 0});  // face 1
    m.faces().push_face({0, 3, 2, 1});  // face 2 (interne)
    m.faces().push_face({0, 1});        // face 3
    m.faces().push_face({2, 3});        // face 4

    // Owner/neighbour
    m.ownership().resize(5);
    m.ownership().set_owner(0, 0);
    m.ownership().set_neighbour(0, FaceOwnership::BOUNDARY);
    m.ownership().set_owner(1, 1);
    m.ownership().set_neighbour(1, FaceOwnership::BOUNDARY);
    m.ownership().set_owner(2, 0);
    m.ownership().set_neighbour(2, 1);  // interne
    m.ownership().set_owner(3, 0);
    m.ownership().set_neighbour(3, FaceOwnership::BOUNDARY);
    m.ownership().set_owner(4, 1);
    m.ownership().set_neighbour(4, FaceOwnership::BOUNDARY);

    // 2 cellules
    m.cells().push_cell({0, 2, 3});  // cellule 0 : faces 0, 2, 3
    m.cells().push_cell({1, 2, 4});  // cellule 1 : faces 1, 2, 4

    // Patches
    BoundaryPatches bp;
    Patch pin;
    pin.name = "bottom";
    pin.type = PatchType::WALL;
    pin.face_ids = {0};
    bp.add_patch(pin);

    Patch ptop;
    ptop.name = "top";
    ptop.type = PatchType::WALL;
    ptop.face_ids = {1};
    bp.add_patch(ptop);

    Patch pleft;
    pleft.name = "left";
    pleft.type = PatchType::INLET;
    pleft.face_ids = {3};
    bp.add_patch(pleft);

    Patch pright;
    pright.name = "right";
    pright.type = PatchType::OUTLET;
    pright.face_ids = {4};
    bp.add_patch(pright);

    m.boundary() = bp;

    return m;
}

int main() {
    run_case("default_empty", []() {
        Mesh m;
        EXPECT_TRUE(m.n_points() == 0);
        EXPECT_TRUE(m.n_faces() == 0);
        EXPECT_TRUE(m.n_cells() == 0);
        EXPECT_TRUE(m.validate().ok);
    });

    run_case("simple_mesh_stats", []() {
        Mesh m = make_simple_mesh();
        EXPECT_TRUE(m.n_points() == 4);
        EXPECT_TRUE(m.n_faces() == 5);
        EXPECT_TRUE(m.n_cells() == 2);
        EXPECT_TRUE(m.stats().n_patches == 4);
        EXPECT_TRUE(m.stats().n_internal_faces == 1);
        EXPECT_TRUE(m.stats().n_boundary_faces == 4);
    });

    run_case("simple_mesh_validate_ok", []() {
        Mesh m = make_simple_mesh();
        auto result = m.validate();
        if (!result.ok) {
            for (const auto& e : result.errors) {
                std::fprintf(stderr, "  error: %s\n", e.c_str());
            }
        }
        EXPECT_TRUE(result.ok);
    });

    run_case("invalid_points_nan", []() {
        Mesh m = make_simple_mesh();
        m.points().set(0, std::nan(""), 0.0, 0.0);
        EXPECT_FALSE(m.validate().ok);
    });

    run_case("invalid_owner_out_of_range", []() {
        Mesh m = make_simple_mesh();
        m.ownership().set_owner(0, 99);
        EXPECT_FALSE(m.validate().ok);
    });

    run_case("invalid_cell_face_refs", []() {
        Mesh m = make_simple_mesh();
        m.cells().clear();
        m.cells().push_cell({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
        EXPECT_FALSE(m.validate().ok);
    });

    run_case("clear", []() {
        Mesh m = make_simple_mesh();
        m.clear();
        EXPECT_TRUE(m.n_points() == 0);
        EXPECT_TRUE(m.n_faces() == 0);
        EXPECT_TRUE(m.n_cells() == 0);
    });

    run_case("boundary_patch_consistent", []() {
        Mesh m = make_simple_mesh();
        EXPECT_TRUE(m.boundary().is_consistent(m.n_faces()));
    });

    run_case("boundary_patch_overlap", []() {
        Mesh m = make_simple_mesh();
        Patch p;
        p.name = "dup";
        p.type = PatchType::WALL;
        p.face_ids = {0};  // face 0 déjà dans "bottom"
        m.boundary().add_patch(p);
        EXPECT_FALSE(m.boundary().is_consistent(m.n_faces()));
    });

    return run_all();
}