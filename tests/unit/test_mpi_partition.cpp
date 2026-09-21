// M0.11 — Tests for MPI partitioning and halo exchange

#include "cfdx/core/parallel/mesh_partitioner.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::core::parallel;
using namespace cfdx::testing;

int main() {
    int argc = 0;
    char** argv = nullptr;
    MPI_Init(&argc, &argv);

    run_case("mpi_partition_geometric_basic", []() {
        Mesh m;

        m.points().resize(9);
        m.points().set(0, 0.0, 0.0, 0.0);
        m.points().set(1, 1.0, 0.0, 0.0);
        m.points().set(2, 2.0, 0.0, 0.0);
        m.points().set(3, 0.0, 1.0, 0.0);
        m.points().set(4, 1.0, 1.0, 0.0);
        m.points().set(5, 2.0, 1.0, 0.0);
        m.points().set(6, 0.0, 2.0, 0.0);
        m.points().set(7, 1.0, 2.0, 0.0);
        m.points().set(8, 2.0, 2.0, 0.0);

        m.faces().push_face({0, 3, 2, 1});
        m.faces().push_face({1, 2, 5, 4});
        m.faces().push_face({3, 4, 7, 6});

        m.cells().push_cell({0, 1, 2, 3, 4});
        m.cells().push_cell({1, 2, 5, 4});

        m.ownership().resize(3);
        m.ownership().set_owner(0, 0); m.ownership().set_neighbour(0, FaceOwnership::BOUNDARY);
        m.ownership().set_owner(1, 0); m.ownership().set_neighbour(1, 1);
        m.ownership().set_owner(2, 1); m.ownership().set_neighbour(2, FaceOwnership::BOUNDARY);

        Patch p1;
        p1.name = "bottom";
        p1.type = PatchType::WALL;
        p1.face_ids = {0};
        Patch p2;
        p2.name = "top";
        p2.type = PatchType::WALL;
        p2.face_ids = {2};
        m.boundary().add_patch(p1);
        m.boundary().add_patch(p2);

        m.topo_validate();

        Partition part = partition_geometric(m, 2);

        EXPECT_TRUE(part.n_parts == 2);
        EXPECT_TRUE(part.cell_rank.size() == 2);

        for (const auto& r : part.cell_rank) {
            EXPECT_TRUE(r >= 0 && r < 2);
        }
    });

    run_case("mpi_partition_single_cell", []() {
        Mesh m;
        m.points().resize(4);
        m.points().set(0, 0.0, 0.0, 0.0);
        m.points().set(1, 1.0, 0.0, 0.0);
        m.points().set(2, 1.0, 1.0, 0.0);
        m.points().set(3, 0.0, 1.0, 0.0);

        m.faces().push_face({0, 1, 2, 3});

        m.ownership().resize(1);
        m.ownership().set_owner(0, 0);
        m.ownership().set_neighbour(0, FaceOwnership::BOUNDARY);

        m.cells().push_cell({0, 1, 2, 3});

        m.topo_validate();

        Partition part = partition_geometric(m, 2);
        EXPECT_TRUE(part.cell_rank.size() == 1);
        EXPECT_TRUE(part.cell_rank[0] == 0);
    });

    run_case("mpi_partition_empty_mesh", []() {
        Mesh m;
        Partition part = partition_geometric(m, 2);
        EXPECT_TRUE(part.n_parts == 2);
        EXPECT_TRUE(part.cell_rank.size() == 0);
        EXPECT_TRUE(part.face_owner_rank.size() == 0);
        EXPECT_TRUE(part.face_ghost_rank.size() == 0);
    });

    run_case("mpi_partition_4cells", []() {
        Mesh m;

        m.points().resize(12);
        m.points().set(0, 0.0, 0.0, 0.0);
        m.points().set(1, 1.0, 0.0, 0.0);
        m.points().set(2, 2.0, 0.0, 0.0);
        m.points().set(3, 3.0, 0.0, 0.0);
        m.points().set(4, 0.0, 1.0, 0.0);
        m.points().set(5, 1.0, 1.0, 0.0);
        m.points().set(6, 2.0, 1.0, 0.0);
        m.points().set(7, 3.0, 1.0, 0.0);
        m.points().set(8, 0.0, 0.0, 1.0);
        m.points().set(9, 1.0, 0.0, 1.0);
        m.points().set(10, 2.0, 0.0, 1.0);
        m.points().set(11, 3.0, 0.0, 1.0);

        for (int i = 0; i < 4; ++i) {
            m.faces().push_face({i, i+3, i+4, i+1});
            m.faces().push_face({i+4, i+7, i+8, i+5});
        }

        m.ownership().resize(8);
        for (int i = 0; i < 8; ++i) {
            m.ownership().set_owner(i, i/2);
            m.ownership().set_neighbour(i, FaceOwnership::BOUNDARY);
        }
        m.ownership().set_neighbour(3, 1);
        m.ownership().set_neighbour(7, 3);

        m.cells().push_cell({0, 1, 5, 4});
        m.cells().push_cell({1, 2, 6, 5});
        m.cells().push_cell({2, 3, 7, 6});
        m.cells().push_cell({4, 5, 9, 8});
        m.cells().push_cell({5, 6, 10, 9});

        m.topo_validate();

        Partition part = partition_geometric(m, 4);
        EXPECT_TRUE(part.n_parts == 4);
        EXPECT_TRUE(part.cell_rank.size() == 5);
    });

    const int rc = run_all();
    MPI_Finalize();
    return rc;
}