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

    run_case("mpi_halo_cell_exchange_is_rank_consistent", []() {
        Mesh m;
        m.points().resize(12);
        const double p[12][3] = {
            {0,0,0},{1,0,0},{2,0,0},{0,1,0},{1,1,0},{2,1,0},
            {0,0,1},{1,0,1},{2,0,1},{0,1,1},{1,1,1},{2,1,1}};
        for (std::size_t i=0;i<12;++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);
        m.faces().push_face({0,3,9,6});
        m.faces().push_face({1,4,10,7});
        m.faces().push_face({2,5,11,8});
        m.faces().push_face({0,1,7,6});
        m.faces().push_face({1,2,8,7});
        m.faces().push_face({3,9,10,4});
        m.faces().push_face({4,10,11,5});
        m.faces().push_face({0,6,8,2});
        m.ownership().resize(8);
        for (std::size_t f=0; f<8; ++f) {
            m.ownership().set_owner(f, f == 4 ? 1 : (f >= 6 ? 1 : 0));
            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        }
        m.ownership().set_owner(1, 0);
        m.ownership().set_neighbour(1, 1);
        m.cells().push_cell({0,1,3,5,7});
        m.cells().push_cell({1,2,4,6,7});
        m.topo_validate();

        const int rank = mpi_rank(MPI_COMM_WORLD);
        const Partition part = partition_geometric(m, 2);
        const HaloPlan plan = build_halo_plan(m, part);
        Field<double, Location::CELL> values(m.n_cells(), "phi", "1", 1);
        values.fill(static_cast<double>(rank));
        exchange_halo_cells(values, plan);

        if (mpi_size(MPI_COMM_WORLD) == 2) {
            const int remote = 1 - rank;
            EXPECT_TRUE(part.cell_rank[0] != part.cell_rank[1]);
            const std::size_t remote_cell = (part.cell_rank[0] == remote) ? 0 : 1;
            EXPECT_NEAR(values(remote_cell), static_cast<double>(remote), 1e-14);
            EXPECT_NEAR(values(rank == part.cell_rank[0] ? 0 : 1),
                        static_cast<double>(rank), 1e-14);
        }
    });

    const int rc = run_all();
    MPI_Finalize();
    return rc;
}