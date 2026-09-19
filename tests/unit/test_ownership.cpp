// M0.1-T03 — Tests for FaceOwnership

#include "cfdx/core/mesh/ownership.h"
#include "test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("default_empty", []() {
        FaceOwnership fo;
        EXPECT_TRUE(fo.size() == 0);
        EXPECT_TRUE(fo.empty());
        EXPECT_TRUE(fo.n_internal_faces() == 0);
        EXPECT_TRUE(fo.n_boundary_faces() == 0);
    });

    run_case("explicit_size_boundary_default", []() {
        FaceOwnership fo(5);
        EXPECT_TRUE(fo.size() == 5);
        for (std::size_t i = 0; i < 5; ++i) {
            EXPECT_TRUE(fo.owner(i) == 0);
            EXPECT_TRUE(fo.neighbour(i) == FaceOwnership::BOUNDARY);
        }
        EXPECT_TRUE(fo.n_boundary_faces() == 5);
        EXPECT_TRUE(fo.n_internal_faces() == 0);
    });

    run_case("set_owner_neighbour", []() {
        FaceOwnership fo(3);
        fo.set_owner(0, 0);
        fo.set_neighbour(0, 1);
        fo.set_owner(1, 1);
        fo.set_neighbour(1, 2);
        fo.set_owner(2, 2);
        fo.set_neighbour(2, FaceOwnership::BOUNDARY);

        EXPECT_TRUE(fo.owner(0) == 0);
        EXPECT_TRUE(fo.neighbour(0) == 1);
        EXPECT_TRUE(fo.neighbour(2) == FaceOwnership::BOUNDARY);
        EXPECT_TRUE(fo.n_internal_faces() == 2);
        EXPECT_TRUE(fo.n_boundary_faces() == 1);
    });

    run_case("consistent_internal", []() {
        FaceOwnership fo(2);
        fo.set_owner(0, 0);
        fo.set_neighbour(0, 1);
        fo.set_owner(1, 1);
        fo.set_neighbour(1, 0);
        EXPECT_TRUE(fo.is_consistent(2));
    });

    run_case("inconsistent_owner_out_of_range", []() {
        FaceOwnership fo(1);
        fo.set_owner(0, 5);
        fo.set_neighbour(0, FaceOwnership::BOUNDARY);
        EXPECT_FALSE(fo.is_consistent(3));
    });

    run_case("inconsistent_self_connect", []() {
        FaceOwnership fo(1);
        fo.set_owner(0, 0);
        fo.set_neighbour(0, 0);
        EXPECT_FALSE(fo.is_consistent(2));
    });

    run_case("index_out_of_range", []() {
        FaceOwnership fo(1);
        EXPECT_THROW(fo.owner(5), std::out_of_range);
    });

    run_case("bulk_data_access", []() {
        FaceOwnership fo(2);
        fo.set_owner(0, 0);
        fo.set_neighbour(0, 1);
        fo.set_owner(1, 1);
        fo.set_neighbour(1, FaceOwnership::BOUNDARY);

        const auto* od = fo.owner_data();
        const auto* nd = fo.neighbour_data();
        EXPECT_TRUE(od[0] == 0 && od[1] == 1);
        EXPECT_TRUE(nd[0] == 1 && nd[1] == FaceOwnership::BOUNDARY);
    });

    run_case("resize", []() {
        FaceOwnership fo;
        fo.resize(4);
        EXPECT_TRUE(fo.size() == 4);
        EXPECT_TRUE(fo.n_boundary_faces() == 4);
    });

    run_case("clear", []() {
        FaceOwnership fo(3);
        fo.clear();
        EXPECT_TRUE(fo.size() == 0);
        EXPECT_TRUE(fo.empty());
    });

    return run_all();
}