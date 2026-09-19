// M0.1-T05 — Tests for BoundaryPatches

#include "cfdx/core/mesh/boundary.h"
#include "test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("default_empty", []() {
        BoundaryPatches bp;
        EXPECT_TRUE(bp.n_patches() == 0);
        EXPECT_TRUE(bp.n_boundary_faces() == 0);
    });

    run_case("add_patch", []() {
        BoundaryPatches bp;
        std::size_t idx = bp.add_patch("inlet", PatchType::INLET);
        EXPECT_TRUE(idx == 0);
        EXPECT_TRUE(bp.n_patches() == 1);
        EXPECT_TRUE(bp.has_patch("inlet"));
        EXPECT_FALSE(bp.has_patch("wall"));
        EXPECT_TRUE(bp.find("inlet") == 0);
        EXPECT_TRUE(bp.find("wall") == 1);
    });

    run_case("duplicate_patch_name", []() {
        BoundaryPatches bp;
        bp.add_patch("inlet", PatchType::INLET);
        EXPECT_THROW(bp.add_patch("inlet", PatchType::OUTLET), std::runtime_error);
    });

    run_case("patch_with_faces", []() {
        BoundaryPatches bp;
        Patch p;
        p.name = "wall";
        p.type = PatchType::WALL;
        p.face_ids = {0, 1, 2, 3};
        bp.add_patch(p);
        EXPECT_TRUE(bp.n_boundary_faces() == 4);
        EXPECT_TRUE(bp.patch(0).size() == 4);
    });

    run_case("consistent", []() {
        BoundaryPatches bp;
        Patch p1;
        p1.name = "inlet"; p1.type = PatchType::INLET;
        p1.face_ids = {0, 1, 2};
        bp.add_patch(p1);

        Patch p2;
        p2.name = "outlet"; p2.type = PatchType::OUTLET;
        p2.face_ids = {3, 4};
        bp.add_patch(p2);

        EXPECT_TRUE(bp.is_consistent(10));
    });

    run_case("inconsistent_overlap", []() {
        BoundaryPatches bp;
        Patch p1;
        p1.name = "inlet"; p1.type = PatchType::INLET;
        p1.face_ids = {0, 1, 2};
        bp.add_patch(p1);

        Patch p2;
        p2.name = "outlet"; p2.type = PatchType::OUTLET;
        p2.face_ids = {2, 3};
        bp.add_patch(p2);

        EXPECT_FALSE(bp.is_consistent(10));
    });

    run_case("inconsistent_out_of_range", []() {
        BoundaryPatches bp;
        Patch p;
        p.name = "inlet"; p.type = PatchType::INLET;
        p.face_ids = {0, 100};
        bp.add_patch(p);
        EXPECT_FALSE(bp.is_consistent(10));
    });

    run_case("remove_patch", []() {
        BoundaryPatches bp;
        bp.add_patch("inlet", PatchType::INLET);
        bp.add_patch("outlet", PatchType::OUTLET);
        bp.remove_patch(0);
        EXPECT_TRUE(bp.n_patches() == 1);
        EXPECT_TRUE(bp.has_patch("outlet"));
    });

    run_case("metadata", []() {
        BoundaryPatches bp;
        Patch p;
        p.name = "inlet";
        p.type = PatchType::INLET;
        p.metadata["velocity"] = "10.0";
        p.metadata["unit"] = "m/s";
        bp.add_patch(p);
        EXPECT_TRUE(bp.patch(0).metadata.at("velocity") == "10.0");
    });

    run_case("patch_type_from_string", []() {
        EXPECT_TRUE(patch_type_from_string("wall") == PatchType::WALL);
        EXPECT_TRUE(patch_type_from_string("inlet") == PatchType::INLET);
        EXPECT_TRUE(patch_type_from_string("outlet") == PatchType::OUTLET);
        EXPECT_TRUE(patch_type_from_string("symmetry") == PatchType::SYMMETRY);
        EXPECT_TRUE(patch_type_from_string("periodic") == PatchType::PERIODIC);
        EXPECT_TRUE(patch_type_from_string("unknown_type") == PatchType::UNKNOWN);
    });

    run_case("index_out_of_range", []() {
        BoundaryPatches bp;
        EXPECT_THROW(bp.patch(5), std::out_of_range);
    });

    return run_all();
}