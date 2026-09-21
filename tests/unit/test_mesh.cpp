// M0.3-T01 — Mesh validator tests (topology + geometry + quality + conservation)

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/mesh_validator.h"
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
    run_case("default_empty", []() {
        Mesh m;
        EXPECT_TRUE(m.n_points() == 0);
        EXPECT_TRUE(m.n_faces() == 0);
        EXPECT_TRUE(m.n_cells() == 0);
        EXPECT_TRUE(m.topo_validate().ok);
    });

    run_case("unit_cube_stats", []() {
        Mesh m = make_unit_cube();
        EXPECT_TRUE(m.n_points() == 8);
        EXPECT_TRUE(m.n_faces() == 6);
        EXPECT_TRUE(m.n_cells() == 1);
        EXPECT_TRUE(m.stats().n_patches == 6);
        EXPECT_TRUE(m.stats().n_boundary_faces == 6);
        EXPECT_TRUE(m.stats().n_internal_faces == 0);
    });

    run_case("topology_view_is_geometry_independent", []() {
        Mesh m = make_unit_cube();
        const auto topology = m.topology();
        EXPECT_TRUE(topology.n_cells() == m.n_cells());
        EXPECT_TRUE(topology.n_faces() == m.n_faces());
        EXPECT_TRUE(topology.cells().n_cells() == m.cells().n_cells());
        EXPECT_TRUE(topology.boundary().n_patches() == m.boundary().n_patches());
    });

    run_case("unit_cube_topo_validate_ok", []() {
        Mesh m = make_unit_cube();
        auto result = m.topo_validate();
        if (!result.ok) {
            for (const auto& e : result.errors) {
                std::fprintf(stderr, "  error: %s\n", e.c_str());
            }
        }
        EXPECT_TRUE(result.ok);
    });

    run_case("invalid_points_nan", []() {
        Mesh m = make_unit_cube();
        m.points().set(0, std::nan(""), 0.0, 0.0);
        EXPECT_FALSE(m.topo_validate().ok);
    });

    run_case("invalid_owner_out_of_range", []() {
        Mesh m = make_unit_cube();
        m.ownership().set_owner(0, 99);
        EXPECT_FALSE(m.topo_validate().ok);
    });

    run_case("invalid_cell_face_refs", []() {
        Mesh m = make_unit_cube();
        m.cells().clear();
        m.cells().push_cell({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
        EXPECT_FALSE(m.topo_validate().ok);
    });

    run_case("boundary_patch_consistent", []() {
        Mesh m = make_unit_cube();
        EXPECT_TRUE(m.boundary().is_consistent(m.n_faces()));
    });

    run_case("boundary_patch_overlap", []() {
        Mesh m = make_unit_cube();
        Patch p;
        p.name = "dup";
        p.type = PatchType::WALL;
        p.face_ids = {0};
        m.boundary().add_patch(p);
        EXPECT_FALSE(m.boundary().is_consistent(m.n_faces()));
    });

    run_case("validate_mesh_full_ok", []() {
        Mesh m = make_unit_cube();
        auto report = validate_mesh(m);
        if (!report.ok) {
            for (const auto& e : report.errors) {
                std::fprintf(stderr, "  error: %s\n", e.c_str());
            }
        }
        EXPECT_TRUE(report.ok);
    });

    run_case("validate_mesh_detects_bad_geometry", []() {
        Mesh m = make_unit_cube();
        m.points().set(0, std::nan(""), 0.0, 0.0);
        auto report = validate_mesh(m);
        EXPECT_FALSE(report.ok);
    });

    run_case("validate_mesh_reports_quality", []() {
        Mesh m = make_unit_cube();
        auto report = validate_mesh(m);
        EXPECT_TRUE(report.max_skewness >= 0.0);
        EXPECT_TRUE(report.surface_closure_error >= 0.0);
        // Pour un cube parfait, la fermeture des surfaces est ~0.
        EXPECT_NEAR(report.surface_closure_error, 0.0, 1e-12);
        // Le volume du cube unité est 1.0.
        EXPECT_NEAR(report.max_cell_volume, 1.0, 1e-12);
        EXPECT_NEAR(report.min_cell_volume, 1.0, 1e-12);
    });

    return run_all();
}