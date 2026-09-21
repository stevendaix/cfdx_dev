#include "cfdx/core/geometry/mesh_validator.h"
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

    return run_all();
}
