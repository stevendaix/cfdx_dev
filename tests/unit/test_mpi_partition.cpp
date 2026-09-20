// M0.11 — Tests for MPI partitioning and halo exchange

#include "cfdx/core/parallel/mesh_partitioner.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::core::parallel;
using namespace cfdx::testing;

int main() {
    run_case("mpi_partition_geometric", []() {
        // Create a simple 2x2x1 mesh (4 cells)
        Mesh m;

        // 9 points in a 2x2 grid
        m.points().add(0.0, 0.0, 0.0);
        m.points().add(1.0, 0.0, 0.0);
        m.points().add(2.0, 0.0, 0.0);
        m.points().add(0.0, 1.0, 0.0);
        m.points().add(1.0, 1.0, 0.0);
        m.points().add(2.0, 1.0, 0.0);
        m.points().add(0.0, 2.0, 0.0);
        m.points().add(1.0, 2.0, 0.0);
        m.points().add(2.0, 2.0, 0.0);

        // 4 cells (2x2 quads)
        m.faces().add({0, 1, 4, 3});
        m.faces().add({1, 4, 3, 0});
        m.faces().add({3, 4, 7, 6});
        m.faces().add({0, 3, 6, 4});

        m.faces().add({1, 2, 5, 4});
        m.faces().add({2, 5, 4, 1});
        m.faces().add({4, 5, 8, 7});

        m.faces().add({3, 4, 7, 6});
        m.faces().add({4, 7, 6, 3});
        m.faces().add({6, 7, 0, 0});
        m.faces().add({3, 6, 0, 0});

        m.faces().add({4, 5, 8, 7});
        m.faces().add({5, 8, 7, 4});
        m.faces().add({7, 8, 0, 0});

        // Set up cell connectivity
        m.cells().add({0, 1, 2, 3});
        m.cells().add({4, 5, 6, 1});
        m.cells().add({7, 8, 9, 3});
        m.cells().add({11, 12, 13, 8});

        // Set up ownership
        m.ownership().owner.resize(14);
        m.ownership().neighbour.resize(14);
        m.ownership().owner[0] = 0;  m.ownership().neighbour[0] = -1;
        m.ownership().owner[1] = 0;  m.ownership().neighbour[1] = 1;
        m.ownership().owner[2] = 0;  m.ownership().neighbour[2] = 2;
        m.ownership().owner[3] = 0;  m.ownership().neighbour[3] = -1;
        m.ownership().owner[4] = 1;  m.ownership().neighbour[4] = -1;
        m.ownership().owner[5] = 1;  m.ownership().neighbour[5] = -1;
        m.ownership().owner[6] = 1;  m.ownership().neighbour[6] = 3;
        m.ownership().owner[1] = 1;  m.ownership().neighbour[1] = 0;
        m.ownership().owner[7] = 2;  m.ownership().neighbour[7] = 0;
        m.ownership().owner[8] = 2;  m.ownership().neighbour[8] = 3;
        m.ownership().owner[9] = 2;  m.ownership().neighbour[9] = -1;
        m.ownership().owner[10] = 2; m.ownership().neighbour[10] = -1;
        m.ownership().owner[11] = 3; m.ownership().neighbour[11] = 1;
        m.ownership().owner[12] = 3; m.ownership().neighbour[12] = -1;
        m.ownership().owner[13] = 3; m.ownership().neighbour[13] = -1;
        m.ownership().owner[8] = 3;  m.ownership().neighbour[8] = 2;
        m.ownership().owner[6] = 3;  m.ownership().neighbour[6] = 1;

        m.boundary().patches.push_back({"bottom", {0, 4}});
        m.boundary().patches.push_back({"right", {5, 12}});
        m.boundary().patches.push_back({"top", {9, 13}});
        m.boundary().patches.push_back({"left", {3, 10}});

        m.topo_validate();

        Partition part = partition_geometric(m, 2);

        EXPECT_EQUAL(part.n_parts, 2);
        EXPECT_EQUAL(part.cell_rank.size(), 4);

        for (int r : part.cell_rank) {
            EXPECT_TRUE(r >= 0 && r < 2);
        }

        for (std::size_t f = 0; f < part.face_owner_rank.size(); ++f) {
            if (m.ownership().owner[f] >= 0) {
                int expected = part.cell_rank[m.ownership().owner[f]];
                EXPECT_EQUAL(part.face_owner_rank[f], expected);
            }
        }

        HaloPlan plan = build_halo_plan(m, part);
        EXPECT_EQUAL(plan.send_faces.size(), 2);
        EXPECT_EQUAL(plan.recv_faces.size(), 2);

        int total_send = 0, total_recv = 0;
        for (int r = 0; r < 2; ++r) {
            total_send += plan.send_faces[r].size();
            total_recv += plan.recv_faces[r].size();
        }
        EXPECT_EQUAL(total_send, total_recv);
    });

    run_case("mpi_halo_exchange_face_field", []() {
        Mesh m;
        m.points().add(0.0, 0.0, 0.0);
        m.points().add(1.0, 0.0, 0.0);
        m.points().add(2.0, 0.0, 0.0);
        m.points().add(0.0, 1.0, 0.0);
        m.points().add(1.0, 1.0, 0.0);
        m.points().add(2.0, 1.0, 0.0);

        m.faces().add({0, 1, 4, 3});
        m.faces().add({1, 4, 3, 0});
        m.faces().add({1, 2, 5, 4});
        m.faces().add({3, 4, 0, 0});
        m.faces().add({4, 5, 0, 0});

        m.cells().add({0, 1, 3});
        m.cells().add({1, 2, 4});

        m.ownership().owner.resize(5);
        m.ownership().neighbour.resize(5);
        m.ownership().owner[0] = 0;  m.ownership().neighbour[0] = -1;
        m.ownership().owner[1] = 0;  m.ownership().neighbour[1] = 1;
        m.ownership().owner[2] = 1;  m.ownership().neighbour[2] = -1;
        m.ownership().owner[3] = 0;  m.ownership().neighbour[3] = -1;
        m.ownership().owner[4] = 1;  m.ownership().neighbour[4] = -1;
        m.ownership().owner[1] = 1;  m.ownership().neighbour[1] = 0;

        m.boundary().patches.push_back({"boundary", {0, 2, 3, 4}});
        m.topo_validate();

        Field<double, Location::FACE> field(1, m.n_faces());
        field.fill(0.0);

        Partition part = partition_geometric(m, 2);
        HaloPlan plan = build_halo_plan(m, part);

        EXPECT_EQUAL(part.n_parts, 2);
        EXPECT_EQUAL(part.cell_rank.size(), 2);
    });

    run_case("mpi_partition_empty_mesh", []() {
        Mesh m;
        Partition part = partition_geometric(m, 2);
        EXPECT_EQUAL(part.n_parts, 2);
        EXPECT_EQUAL(part.cell_rank.size(), 0);
        EXPECT_EQUAL(part.face_owner_rank.size(), 0);
        EXPECT_EQUAL(part.face_ghost_rank.size(), 0);
    });

    run_case("mpi_partition_single_cell", []() {
        Mesh m;
        m.points().add(0.0, 0.0, 0.0);
        m.points().add(1.0, 0.0, 0.0);
        m.points().add(1.0, 1.0, 0.0);
        m.points().add(0.0, 1.0, 0.0);

        m.faces().add({0, 1, 2, 3});
        m.faces().add({1, 2, 0, 0});
        m.faces().add({2, 3, 0, 0});
        m.faces().add({3, 0, 0, 0});

        m.cells().add({0, 1, 2, 3});

        m.ownership().owner.resize(4);
        m.ownership().neighbour.resize(4);
        m.ownership().owner[0] = 0;  m.ownership().neighbour[0] = -1;
        m.ownership().owner[1] = 0;  m.ownership().neighbour[1] = -1;
        m.ownership().owner[2] = 0;  m.ownership().neighbour[2] = -1;
        m.ownership().owner[3] = 0;  m.ownership().neighbour[3] = -1;

        m.boundary().patches.push_back({"all", {0, 1, 2, 3}});
        m.topo_validate();

        Partition part = partition_geometric(m, 2);
        EXPECT_EQUAL(part.cell_rank.size(), 1);
        EXPECT_EQUAL(part.cell_rank[0], 0);
    });

    return run_all();
}
