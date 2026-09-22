// M0.11 — Tests for GeometryCache

#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/mesh/mesh.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

// Helper: build a unit cube mesh
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

    m.faces().push_face({0, 3, 2, 1});  // bas
    m.faces().push_face({4, 5, 6, 7});  // haut
    m.faces().push_face({0, 1, 5, 4});  // avant
    m.faces().push_face({3, 7, 6, 2});  // arrière
    m.faces().push_face({0, 4, 7, 3});  // gauche
    m.faces().push_face({1, 2, 6, 5});  // droite

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
    run_case("geometry_cache_basic", []() {
        Mesh m = make_unit_cube();
        GeometryCache cache = make_geometry_cache(m);
        
        EXPECT_TRUE(is_valid(cache, m));
        EXPECT_TRUE(cache.valid);
        EXPECT_TRUE(cache.face_centres.size() == 6);
        EXPECT_TRUE(cache.cell_centres.size() == 1);
        EXPECT_TRUE(cache.face_Sf.size() == 6);
        EXPECT_TRUE(cache.cell_volumes.size() == 1);
    });

    run_case("geometry_cache_cube_volume", []() {
        Mesh m = make_unit_cube();
        GeometryCache cache = make_geometry_cache(m);
        
        // Volume should be 1.0
        EXPECT_NEAR(cache.cell_volumes[0], 1.0, 1e-12);
    });

    run_case("geometry_cache_cube_centre", []() {
        Mesh m = make_unit_cube();
        GeometryCache cache = make_geometry_cache(m);
        
        EXPECT_NEAR(cache.cell_centres[0].x, 0.5, 1e-12);
        EXPECT_NEAR(cache.cell_centres[0].y, 0.5, 1e-12);
        EXPECT_NEAR(cache.cell_centres[0].z, 0.5, 1e-12);
    });

    run_case("geometry_cache_face_areas", []() {
        Mesh m = make_unit_cube();
        GeometryCache cache = make_geometry_cache(m);
        
        // Each face area should be 1.0
        for (std::size_t f = 0; f < 6; ++f) {
            EXPECT_NEAR(cache.face_areas[f], 1.0, 1e-12);
        }
    });

    run_case("geometry_cache_face_normals", []() {
        Mesh m = make_unit_cube();
        GeometryCache cache = make_geometry_cache(m);
        
        // Bas (z=0) → normal = (0,0,-1)
        EXPECT_NEAR(cache.face_normals[0].x, 0.0, 1e-12);
        EXPECT_NEAR(cache.face_normals[0].y, 0.0, 1e-12);
        EXPECT_NEAR(cache.face_normals[0].z, -1.0, 1e-12);
        
        // Haut (z=1) → normal = (0,0,+1)
        EXPECT_NEAR(cache.face_normals[1].z, 1.0, 1e-12);
        
        // Avant (y=0) → normal = (0,-1,0)
        EXPECT_NEAR(cache.face_normals[2].y, -1.0, 1e-12);
        
        // Arrière (y=1) → normal = (0,+1,0)
        EXPECT_NEAR(cache.face_normals[3].y, 1.0, 1e-12);
        
        // Gauche (x=0) → normal = (-1,0,0)
        EXPECT_NEAR(cache.face_normals[4].x, -1.0, 1e-12);
        
        // Droite (x=1) → normal = (+1,0,0)
        EXPECT_NEAR(cache.face_normals[5].x, 1.0, 1e-12);
    });

    run_case("geometry_cache_surface_closure", []() {
        Mesh m = make_unit_cube();
        GeometryCache cache = make_geometry_cache(m);
        
        // Surface closure error should be near zero for a valid cube
        EXPECT_NEAR(cache.surface_closure_error[0], 0.0, 1e-10);
        EXPECT_NEAR(cache.global_surface_closure_error, 0.0, 1e-10);
    });

    run_case("geometry_cache_quality_metrics", []() {
        Mesh m = make_unit_cube();
        GeometryCache cache = make_geometry_cache(m);
        
        // For a perfect cube, skewness and non-orthogonality should be 0
        for (std::size_t f = 0; f < 6; ++f) {
            EXPECT_NEAR(cache.face_skewness[f], 0.0, 1e-12);
            EXPECT_NEAR(cache.face_non_orthogonality_deg[f], 0.0, 1e-12);
        }
        EXPECT_NEAR(cache.max_skewness, 0.0, 1e-12);
        EXPECT_NEAR(cache.max_non_orthogonality_deg, 0.0, 1e-12);
    });

    run_case("geometry_cache_delta_coeffs", []() {
        Mesh m = make_unit_cube();
        GeometryCache cache = make_geometry_cache(m);
        
        // For a cube, delta coeffs should all be 0.5 (centre to face distance)
        for (std::size_t f = 0; f < 6; ++f) {
            EXPECT_NEAR(cache.delta_coeffs[f], 0.5, 1e-12);
        }
    });

    run_case("geometry_cache_empty_mesh", []() {
        Mesh m;
        GeometryCache cache = make_geometry_cache(m);
        
        EXPECT_TRUE(is_valid(cache, m));
        EXPECT_TRUE(cache.face_centres.empty());
        EXPECT_TRUE(cache.cell_centres.empty());
    });

    run_case("geometry_cache_internal_face_closure_and_quality", []() {
        Mesh m;
        m.points().resize(12);
        const double p[12][3] = {
            {0,0,0},{1,0,0},{2,0,0},{0,1,0},{1,1,0},{2,1,0},
            {0,0,1},{1,0,1},{2,0,1},{0,1,1},{1,1,1},{2,1,1}};
        for (std::size_t i = 0; i < 12; ++i)
            m.points().set(i, p[i][0], p[i][1], p[i][2]);

        // Two unit cubes sharing x=1. The shared face is deliberately stored
        // with the wrong winding; GeometryCache must canonicalize owner->neighbour.
        m.faces().push_face({0,6,9,3});       // left
        m.faces().push_face({0,1,7,6});       // bottom 0
        m.faces().push_face({3,9,10,4});      // top 0
        m.faces().push_face({0,3,4,1});       // front 0
        m.faces().push_face({6,7,10,9});      // back 0
        m.faces().push_face({7,10,4,1});      // shared, reversed winding
        m.faces().push_face({2,5,11,8});      // right
        m.faces().push_face({1,2,8,7});       // bottom 1
        m.faces().push_face({4,10,11,5});     // top 1
        m.faces().push_face({1,4,5,2});       // front 1
        m.faces().push_face({7,8,11,10});     // back 1

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

        BoundaryPatches bp;
        Patch wall;
        wall.name = "walls";
        wall.type = PatchType::WALL;
        wall.face_ids = {0,1,2,3,4,6,7,8,9,10};
        bp.add_patch(wall);
        m.set_boundary(bp);

        GeometryCache cache = make_geometry_cache(m);
        EXPECT_TRUE(is_valid(cache, m));
        EXPECT_NEAR(cache.cell_volumes[0], 1.0, 1e-12);
        EXPECT_NEAR(cache.cell_volumes[1], 1.0, 1e-12);
        EXPECT_NEAR(cache.surface_closure_error[0], 0.0, 1e-12);
        EXPECT_NEAR(cache.surface_closure_error[1], 0.0, 1e-12);
        EXPECT_NEAR(cache.global_surface_closure_error, 0.0, 1e-12);
        EXPECT_NEAR(cache.face_Sf[5].x, 1.0, 1e-12);
        EXPECT_NEAR(cache.face_Sf[5].y, 0.0, 1e-12);
        EXPECT_NEAR(cache.face_Sf[5].z, 0.0, 1e-12);
        EXPECT_NEAR(cache.face_non_orthogonality_deg[5], 0.0, 1e-12);
        EXPECT_NEAR(cache.face_skewness[5], 0.0, 1e-12);
    });

    return run_all();
}
