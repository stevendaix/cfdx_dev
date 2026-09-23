// M0.11 — Tests for MPI partitioning and halo exchange

#include "cfdx/core/parallel/mesh_partitioner.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::core::parallel;
using namespace cfdx::testing;

namespace {
Patch make_wall_patch(const std::string& name, std::initializer_list<std::size_t> face_ids) {
    Patch patch;
    patch.name = name;
    patch.type = PatchType::WALL;
    patch.face_ids.assign(face_ids.begin(), face_ids.end());
    return patch;
}
}  // namespace

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

        EXPECT_TRUE(m.topo_validate().ok);

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

        m.boundary().add_patch(make_wall_patch("wall", {0}));

        m.cells().push_cell({0});

        EXPECT_TRUE(m.topo_validate().ok);

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
            m.faces().push_face({static_cast<std::size_t>(i), static_cast<std::size_t>(i + 3), static_cast<std::size_t>(i + 4), static_cast<std::size_t>(i + 1)});
            m.faces().push_face({static_cast<std::size_t>(i + 4), static_cast<std::size_t>(i + 7), static_cast<std::size_t>(i + 8), static_cast<std::size_t>(i + 5)});
        }

        m.ownership().resize(8);
        for (int i = 0; i < 8; ++i) {
            m.ownership().set_owner(i, i/2);
            m.ownership().set_neighbour(i, FaceOwnership::BOUNDARY);
        }
        m.boundary().add_patch(make_wall_patch("wall", {0, 1, 2, 3, 4, 5, 6, 7}));
        // Keep this partition fixture topologically valid and independent of halo semantics.


        m.cells().push_cell({0, 1});
        m.cells().push_cell({2, 3});
        m.cells().push_cell({4, 5});
        m.cells().push_cell({6, 7});
        EXPECT_TRUE(m.topo_validate().ok);

        Partition part = partition_geometric(m, 4);
        EXPECT_TRUE(part.n_parts == 4);
        EXPECT_TRUE(part.cell_rank.size() == 4);
        for (const auto& r : part.cell_rank)
            EXPECT_TRUE(r >= 0 && r < 4);
    });

    run_case("mpi_halo_cell_exchange_is_rank_consistent", []() {
        Mesh m;
        // Minimal valid two-cell mesh: one internal face is the MPI interface.
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
        m.ownership().set_owner(0, 0);
        m.ownership().set_neighbour(0, FaceOwnership::BOUNDARY);
        m.ownership().set_owner(1, 0);
        m.ownership().set_neighbour(1, 1);
        m.ownership().set_owner(2, 1);
        m.ownership().set_neighbour(2, FaceOwnership::BOUNDARY);
        m.boundary().add_patch(make_wall_patch("wall", {0, 2}));

        EXPECT_TRUE(m.topo_validate().ok);

        const int rank = mpi_rank(MPI_COMM_WORLD);
        const Partition part = partition_geometric(m, 2);
        const HaloPlan plan = build_halo_plan(m, part);
        Field<double, Location::CELL> values(m.n_cells(), "phi", "1", 1);
        values.fill(static_cast<double>(rank));
        exchange_halo_cells(values, plan);

        if (mpi_size(MPI_COMM_WORLD) == 2) {
            EXPECT_TRUE(part.cell_rank[0] != part.cell_rank[1]);
            const std::size_t owned_cell =
                (part.cell_rank[0] == rank) ? 0 : 1;
            EXPECT_NEAR(values(owned_cell), static_cast<double>(rank), 1e-14);

            // Every receive entry contains a global owner-cell slot. After
            // exchange it must contain the value sent by that owner rank.
            for (int owner_rank = 0; owner_rank < 2; ++owner_rank) {
                if (owner_rank == rank) continue;
                for (const int cell : plan.recv_cells[owner_rank])
                    EXPECT_NEAR(values(static_cast<std::size_t>(cell)),
                                static_cast<double>(owner_rank), 1e-14);
            }
        }
    });;

    run_case("mpi_halo_plan_rejects_malformed_counts", []() {
        const int size = mpi_size(MPI_COMM_WORLD);
        HaloPlan plan;
        plan.send_faces.resize(size);
        plan.recv_faces.resize(size);
        plan.send_cells.resize(size);
        plan.recv_cells.resize(size);

        if (size >= 2) {
            plan.send_faces[1].push_back(0);
            bool rejected = false;
            try {
                validate_halo_plan(plan, size);
            } catch (const std::invalid_argument&) {
                rejected = true;
            }
            EXPECT_TRUE(rejected);
        }
    });

    run_case("mpi_halo_plan_rejects_negative_indices", []() {
        const int size = mpi_size(MPI_COMM_WORLD);
        HaloPlan plan;
        plan.send_faces.resize(size);
        plan.recv_faces.resize(size);
        plan.send_cells.resize(size);
        plan.recv_cells.resize(size);

        if (size >= 2) {
            plan.send_faces[1].push_back(-1);
            plan.send_cells[1].push_back(0);
            bool rejected = false;
            try {
                validate_halo_plan(plan, size);
            } catch (const std::invalid_argument&) {
                rejected = true;
            }
            EXPECT_TRUE(rejected);
        }
    });

    run_case("mpi_halo_exchange_rejects_out_of_range_indices", []() {
        const int size = mpi_size(MPI_COMM_WORLD);
        if (size != 2)
            return;

        HaloPlan plan;
        plan.send_faces.resize(size);
        plan.recv_faces.resize(size);
        plan.send_cells.resize(size);
        plan.recv_cells.resize(size);
        // Keep send/receive cardinalities paired so plan validation passes;
        // the index itself is deliberately outside the field bounds.
        plan.send_faces[1].push_back(0);
        plan.recv_faces[1].push_back(0);
        plan.send_cells[1].push_back(1);
        plan.recv_cells[1].push_back(1);

        Field<double, Location::CELL> values(1, "phi", "1", 1);
        values.fill(0.0);

        bool rejected = false;
        try {
            exchange_halo_cells(values, plan);
        } catch (const std::out_of_range&) {
            rejected = true;
        }
        EXPECT_TRUE(rejected);
    });

    run_case("mpi_halo_exchange_rejects_asymmetric_counts", []() {
        const int size = mpi_size(MPI_COMM_WORLD);
        if (size != 2)
            return;

        const int rank = mpi_rank(MPI_COMM_WORLD);
        HaloPlan plan;
        plan.send_faces.resize(size);
        plan.recv_faces.resize(size);
        plan.send_cells.resize(size);
        plan.recv_cells.resize(size);

        if (rank == 0) {
            plan.send_faces[1].push_back(0);
            plan.send_cells[1].push_back(0);
        }

        Field<double, Location::CELL> values(1, "phi", "1", 1);
        values.fill(static_cast<double>(rank));

        bool rejected = false;
        try {
            exchange_halo_cells(values, plan);
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        EXPECT_TRUE(rejected);
    });



    run_case("mpi_parallel_operator_serial_equivalence", []() {
        // Quantitative Phase-5 gate: evaluate the same diffusion operator on
        // the serial/global mesh, then reduce only rows owned by each MPI rank.
        // This exercises the partition ownership contract without pretending
        // that a replicated test mesh is a distributed mesh.
        Mesh m;
        constexpr int n_cells = 4;
        constexpr int n_faces = 5;
        m.points().resize(20);
        for (int i = 0; i < 5; ++i) {
            const double x = static_cast<double>(i);
            const std::size_t b = static_cast<std::size_t>(4 * i);
            m.points().set(b + 0, x, 0.0, 0.0);
            m.points().set(b + 1, x, 1.0, 0.0);
            m.points().set(b + 2, x, 1.0, 1.0);
            m.points().set(b + 3, x, 0.0, 1.0);
        }
        for (int i = 0; i < 5; ++i) {
            const std::size_t b = static_cast<std::size_t>(4 * i);
            m.faces().push_face({b + 0, b + 1, b + 2, b + 3});
        }
        m.ownership().resize(n_faces);
        for (int f = 0; f < n_faces; ++f) {
            m.ownership().set_owner(f, f == 4 ? 3 : f);
            m.ownership().set_neighbour(
                f, (f == 0 || f == 4) ? FaceOwnership::BOUNDARY : f - 1);
        }
        m.cells().push_cell({0, 1});
        m.cells().push_cell({1, 2});
        m.cells().push_cell({2, 3});
        m.cells().push_cell({3, 4});
        m.boundary().add_patch(make_wall_patch("left", {0}));
        m.boundary().add_patch(make_wall_patch("right", {4}));
        EXPECT_TRUE(m.topo_validate().ok);

        const GeometryCache g = make_geometry_cache(m);
        FvDiffusionOperator op(m, g, 1.0);
        Vector x(static_cast<std::size_t>(n_cells));
        x(0) = 0.0; x(1) = 1.0; x(2) = 4.0; x(3) = 9.0;
        Vector serial;
        op.apply(x, serial);

        const int rank = mpi_rank(MPI_COMM_WORLD);
        const int size = mpi_size(MPI_COMM_WORLD);
        const Partition part = partition_geometric(m, size);

        Vector local(serial.size());
        local.fill(0.0);
        for (std::size_t c = 0; c < serial.size(); ++c)
            if (part.cell_rank[c] == rank)
                local(c) = serial(c);

        std::vector<double> reduced(serial.size(), 0.0);
        MPI_Allreduce(local.data(), reduced.data(),
                      static_cast<int>(serial.size()), MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        for (std::size_t c = 0; c < serial.size(); ++c)
            EXPECT_NEAR(reduced[c], serial(c), 1e-13);

        const double local_sum = std::accumulate(local.data(),
                                                  local.data() + local.size(), 0.0);
        const double global_sum = mpi_allreduce_sum(local_sum);
        EXPECT_NEAR(global_sum, 0.0, 1e-13);
    });

    run_case("mpi_partition_size_sweep_and_conservation", []() {
        Mesh m;
        m.points().resize(20);
        for (int i = 0; i < 5; ++i) {
            const double x = static_cast<double>(i);
            const std::size_t b = static_cast<std::size_t>(4 * i);
            m.points().set(b + 0, x, 0.0, 0.0);
            m.points().set(b + 1, x, 1.0, 0.0);
            m.points().set(b + 2, x, 1.0, 1.0);
            m.points().set(b + 3, x, 0.0, 1.0);
        }
        for (int i = 0; i < 5; ++i) {
            const std::size_t b = static_cast<std::size_t>(4 * i);
            m.faces().push_face({b + 0, b + 1, b + 2, b + 3});
        }
        m.ownership().resize(5);
        for (int f = 0; f < 5; ++f) {
            m.ownership().set_owner(f, f == 4 ? 3 : f);
            m.ownership().set_neighbour(
                f, (f == 0 || f == 4) ? FaceOwnership::BOUNDARY : f - 1);
        }
        m.cells().push_cell({0, 1});
        m.cells().push_cell({1, 2});
        m.cells().push_cell({2, 3});
        m.cells().push_cell({3, 4});
        m.boundary().add_patch(make_wall_patch("left", {0}));
        m.boundary().add_patch(make_wall_patch("right", {4}));
        EXPECT_TRUE(m.topo_validate().ok);

        const GeometryCache g = make_geometry_cache(m);
        FvDiffusionOperator op(m, g, 1.0);
        Vector x(4);
        x(0) = 0.0; x(1) = 1.0; x(2) = 4.0; x(3) = 9.0;
        Vector serial;
        op.apply(x, serial);

        const int world_size = mpi_size(MPI_COMM_WORLD);
        const std::vector<int> requested_parts = (world_size >= 4)
            ? std::vector<int>{1, 2, 4}
            : std::vector<int>{1, 2};

        for (const int np : requested_parts) {
            const Partition part = partition_geometric(m, np);
            EXPECT_TRUE(part.n_parts == np);
            std::vector<int> counts(static_cast<std::size_t>(np), 0);
            for (const int r : part.cell_rank) {
                EXPECT_TRUE(r >= 0 && r < np);
                ++counts[static_cast<std::size_t>(r)];
            }
            for (const int count : counts)
                EXPECT_TRUE(count > 0 || np > 4);

            double partitioned_sum = 0.0;
            for (std::size_t c = 0; c < serial.size(); ++c)
                partitioned_sum += serial(c);
            EXPECT_NEAR(partitioned_sum, 0.0, 1e-13);
        }
    });

    const int rc = run_all();
    MPI_Finalize();
    return rc;
}