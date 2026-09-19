// M0.1-T04 — Tests for CellConnectivity

#include "cfdx/core/mesh/cell.h"
#include "test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("empty", []() {
        CellConnectivity cc;
        EXPECT_TRUE(cc.n_cells() == 0);
        EXPECT_TRUE(cc.n_face_refs() == 0);
        EXPECT_TRUE(cc.is_consistent());
    });

    run_case("push_single_cell", []() {
        CellConnectivity cc;
        cc.push_cell({0, 1, 2, 3, 4, 5});
        EXPECT_TRUE(cc.n_cells() == 1);
        EXPECT_TRUE(cc.n_face_refs() == 6);
        EXPECT_TRUE(cc.cell_size(0) == 6);
        EXPECT_TRUE(cc.cell_offset(0) == 0);
        EXPECT_TRUE(cc.is_consistent());
    });

    run_case("push_multiple_cells", []() {
        CellConnectivity cc;
        cc.push_cell({0, 1, 2, 3, 4, 5});
        cc.push_cell({6, 7, 8, 9});
        cc.push_cell({10, 11, 12, 13, 14, 15, 16, 17, 18, 19});

        EXPECT_TRUE(cc.n_cells() == 3);
        EXPECT_TRUE(cc.n_face_refs() == 20);
        EXPECT_TRUE(cc.cell_size(0) == 6);
        EXPECT_TRUE(cc.cell_size(1) == 4);
        EXPECT_TRUE(cc.cell_size(2) == 10);
        EXPECT_TRUE(cc.cell_offset(0) == 0);
        EXPECT_TRUE(cc.cell_offset(1) == 6);
        EXPECT_TRUE(cc.cell_offset(2) == 10);
        EXPECT_TRUE(cc.is_consistent());
    });

    run_case("face_ids_valid", []() {
        CellConnectivity cc;
        cc.push_cell({0, 1, 2, 3});
        cc.push_cell({4, 5, 6, 7});
        EXPECT_TRUE(cc.face_ids_valid(10));
        EXPECT_FALSE(cc.face_ids_valid(5));
    });

    run_case("cell_index_out_of_range", []() {
        CellConnectivity cc;
        cc.push_cell({0, 1, 2});
        EXPECT_THROW(cc.cell_size(5), std::out_of_range);
    });

    run_case("reserve", []() {
        CellConnectivity cc;
        cc.reserve(100, 600);
        cc.push_cell({0, 1, 2, 3, 4, 5});
        EXPECT_TRUE(cc.n_cells() == 1);
        EXPECT_TRUE(cc.is_consistent());
    });

    run_case("bulk_data_access", []() {
        CellConnectivity cc;
        cc.push_cell({0, 1, 2, 3});
        cc.push_cell({4, 5, 6});

        const auto* fd = cc.faces_data();
        const auto* od = cc.offsets_data();
        EXPECT_TRUE(fd[0] == 0 && fd[4] == 4);
        EXPECT_TRUE(od[0] == 0 && od[1] == 4 && od[2] == 7);
    });

    run_case("clear", []() {
        CellConnectivity cc;
        cc.push_cell({0, 1, 2});
        cc.clear();
        EXPECT_TRUE(cc.n_cells() == 0);
        EXPECT_TRUE(cc.n_face_refs() == 0);
        EXPECT_TRUE(cc.is_consistent());
    });

    return run_all();
}