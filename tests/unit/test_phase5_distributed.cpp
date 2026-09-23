#include "cfdx/core/parallel/distributed_execution.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    using namespace cfdx::core::parallel;
    auto require = [&](bool condition, const char* message) {
        if (!condition) {
            std::fprintf(stderr, "rank %d: %s\n", mpi_rank(), message);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    };
    mpi_init(&argc, &argv);
    const int rank = mpi_rank();
    const int size = mpi_size();

    if (size != 2) {
        if (rank == 0)
            std::fprintf(stderr, "test_phase5_distributed requires exactly 2 MPI ranks\n");
        mpi_finalize();
        return 2;
    }

    constexpr std::size_t global_n = 8;
    // Source checkpoint: contiguous 4-cell ownership per rank.
    std::vector<std::uint64_t> ids;
    for (std::size_t i = static_cast<std::size_t>(rank) * 4;
         i < static_cast<std::size_t>(rank + 1) * 4; ++i)
        ids.push_back(static_cast<std::uint64_t>(i));

    validate_distributed_ids(ids, global_n);

    DistributedCellField state(global_n, ids, 1, "phi");
    for (std::size_t i = 0; i < state.local_size(); ++i)
        state(i) = static_cast<double>(state.global_ids()[i] + 1);

    const double local_sum = [&]() {
        double s = 0.0;
        for (std::size_t i = 0; i < state.local_size(); ++i) s += state(i);
        return s;
    }();
    const double global_sum = mpi_deterministic_sum(local_sum);
    require(std::abs(global_sum - 36.0) < 1e-14, "deterministic reduction mismatch");

    // Two cells with one MPI interface face.
    cfdx::core::Mesh mesh;
    mesh.faces().push_face({0,1,2});
    mesh.faces().push_face({1,2,3});
    mesh.faces().push_face({2,3,4});
    mesh.ownership().resize(3);
    mesh.ownership().set_owner(0, 0);
    mesh.ownership().set_neighbour(0, cfdx::core::FaceOwnership::BOUNDARY);
    mesh.ownership().set_owner(1, 0);
    mesh.ownership().set_neighbour(1, 1);
    mesh.ownership().set_owner(2, 1);
    mesh.ownership().set_neighbour(2, cfdx::core::FaceOwnership::BOUNDARY);
    mesh.cells().push_cell({0, 1});
    mesh.cells().push_cell({1, 2});

    cfdx::core::parallel::Partition partition;
    partition.n_parts = 2;
    partition.cell_rank = {0, 1};
    partition.face_owner_rank = {0, 0, 1};
    partition.face_ghost_rank = {-1, 1, -1};

    // Use a separate two-cell field for the actual halo fixture.
    const std::vector<std::uint64_t> local_cell_ids{
        static_cast<std::uint64_t>(rank)};
    DistributedCellField local_state(2, local_cell_ids, 1, "interface");
    local_state(0) = static_cast<double>(rank + 1);

    const auto halo = build_distributed_halo(mesh, partition, local_state);
    const auto received = exchange_distributed_cell_halo(local_state, halo);
    require(received.size() == 1, "halo receive count mismatch");
    require(std::abs(received.front() - static_cast<double>(rank == 0 ? 2 : 1)) < 1e-14,
            "halo value mismatch");

#ifdef CFDX_ENABLE_PARALLEL_HDF5
    const std::string path = "phase5_distributed_checkpoint.h5";
    write_distributed_checkpoint(path, state);
    MPI_Barrier(MPI_COMM_WORLD);

    DistributedCellField restored(global_n, ids, 1, "phi");
    read_distributed_checkpoint(path, restored);
    for (std::size_t i = 0; i < restored.local_size(); ++i)
        require(std::abs(restored(i) - state(i)) < 1e-14, "checkpoint round-trip mismatch");

    // Rank-independent restart: target ownership is a different partition
    // (even/odd IDs), proving that values are mapped by persistent global ID.
    std::vector<std::uint64_t> target_ids;
    for (std::size_t i = static_cast<std::size_t>(rank); i < global_n; i += 2)
        target_ids.push_back(static_cast<std::uint64_t>(i));
    DistributedCellField permuted(global_n, target_ids, 1, "phi");
    read_distributed_checkpoint(path, permuted);
    for (std::size_t i = 0; i < permuted.local_size(); ++i)
        require(std::abs(permuted(i) -
                        (static_cast<double>(permuted.global_ids()[i]) + 1.0)) < 1e-14,
                "rank-independent restart mismatch");

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) std::remove(path.c_str());
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    mpi_finalize();
    return 0;
}
