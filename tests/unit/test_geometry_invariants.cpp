#include "cfdx/core/geometry/mesh_validator.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "common/test_harness.h"
#include <algorithm>

using namespace cfdx::core;
using namespace cfdx::testing;

static Mesh cube(bool invert_all=false) {
    Mesh m;
    m.points().resize(8);
    m.points().set(0,0,0,0); m.points().set(1,1,0,0);
    m.points().set(2,1,1,0); m.points().set(3,0,1,0);
    m.points().set(4,0,0,1); m.points().set(5,1,0,1);
    m.points().set(6,1,1,1); m.points().set(7,0,1,1);
    const std::vector<std::vector<Index>> faces = {
        {0,3,2,1}, {4,5,6,7}, {0,1,5,4},
        {3,7,6,2}, {0,4,7,3}, {1,2,6,5}};
    for (auto face : faces) {
        if (invert_all) std::reverse(face.begin(), face.end());
        m.faces().push_face(face);
    }
    m.ownership().resize(6);
    for(std::size_t f=0; f<6; ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    Patch wall;
    wall.name = "walls";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,1,2,3,4,5};
    m.boundary().add_patch(wall);
    return m;
}

static Mesh two_cell_x_chain() {
    Mesh m;
    m.points().resize(12);
    const double p[12][3] = {
        {0,0,0},{1,0,0},{2,0,0},{0,1,0},{1,1,0},{2,1,0},
        {0,0,1},{1,0,1},{2,0,1},{0,1,1},{1,1,1},{2,1,1}};
    for (std::size_t i = 0; i < 12; ++i) m.points().set(i, p[i][0], p[i][1], p[i][2]);

    const std::vector<std::vector<Index>> faces = {
        {0,6,9,3}, {0,1,7,6}, {3,9,10,4}, {0,3,4,1}, {6,7,10,9},
        {7,10,4,1}, {2,5,11,8}, {1,2,8,7}, {4,10,11,5}, {1,4,5,2},
        {7,8,11,10}};
    for (const auto& face : faces) m.faces().push_face(face);

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
    run_case("closed_cube_has_finite_positive_geometry", [] {
        auto m=cube();
        auto r=validate_mesh(m);
        EXPECT_TRUE(r.ok);
        EXPECT_NEAR(r.min_cell_volume,1.0,1e-12);
        EXPECT_NEAR(r.max_cell_volume,1.0,1e-12);
        EXPECT_TRUE(r.surface_closure_error < 1e-12);
    });

    run_case("inverted_cell_is_rejected", [] {
        auto m=cube(true);
        auto r=validate_mesh(m);
        EXPECT_TRUE(!r.ok);
    });

    run_case("internal_face_orientation_and_indexing_are_consistent", [] {
        const auto m = two_cell_x_chain();
        const auto g = make_geometry_cache(m);

        EXPECT_TRUE(g.valid);
        EXPECT_NEAR(g.cell_volumes[0], 1.0, 1e-12);
        EXPECT_NEAR(g.cell_volumes[1], 1.0, 1e-12);

        const Vec3 centre_delta = g.cell_centres[1] - g.cell_centres[0];
        EXPECT_TRUE(g.face_Sf[5].dot(centre_delta) > 0.0);
        EXPECT_TRUE(g.face_Sf[5].dot(centre_delta) > 0.999999);

        EXPECT_TRUE(g.surface_closure_error[0] < 1e-12);
        EXPECT_TRUE(g.surface_closure_error[1] < 1e-12);

        EXPECT_TRUE(m.ownership().owner(5) == 0);
        EXPECT_TRUE(m.ownership().neighbour(5) == 1);
    });

    return run_all();
}
