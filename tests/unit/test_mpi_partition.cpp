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

        m.cells().push_cell({0, 1});
        m.cells().push_cell({1, 2});

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
        // Four cells represented by a valid 1-D face graph. This test checks
        // partition metadata only; it deliberately does not depend on 3-D
        // geometric volume.
        m.points().resize(10);
        for (std::size_t i=0;i<5;++i) {
            m.points().set(2*i,static_cast<double>(i),0.0,0.0);
            m.points().set(2*i+1,static_cast<double>(i),1.0,0.0);
        }
        for (std::size_t i=0;i<5;++i)
            m.faces().push_face({static_cast<std::uint32_t>(2*i),
                                 static_cast<std::uint32_t>(2*i+1),
                                 static_cast<std::uint32_t>(2*i+1),
                                 static_cast<std::uint32_t>(2*i)});
        m.ownership().resize(5);
        for (std::size_t f=0;f<5;++f) {
            m.ownership().set_owner(f,f==0?0:f-1);
            m.ownership().set_neighbour(f, f<4 ? static_cast<int>(f) + 1 : FaceOwnership::BOUNDARY);
        }
        for (std::size_t c=0;c<4;++c)
            m.cells().push_cell({static_cast<std::uint32_t>(c),
                                 static_cast<std::uint32_t>(c+1)});

        Partition part = partition_geometric(m,4);
        EXPECT_TRUE(part.n_parts==4);
        EXPECT_TRUE(part.cell_rank.size()==4);
        for (const auto r:part.cell_rank) EXPECT_TRUE(r>=0 && r<4);
    });

    return run_all();
}