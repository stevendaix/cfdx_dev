#include "cfdx/core/parallel/distributed_execution.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

int main(int argc, char** argv) {
    using namespace cfdx::core::parallel;
    mpi_init(&argc, &argv);
    const int rank = mpi_rank();
    const int size = mpi_size();

    if (size != 2) {
        if (rank == 0) std::fprintf(stderr, "test_phase5_distributed requires exactly 2 MPI ranks\n");
        mpi_finalize();
        return 2;
    }

    constexpr std::size_t global_n = 8;
    std::vector<std::uint64_t> ids;
    for (std::size_t i = static_cast<std::size_t>(rank); i < global_n; i += 2)
        ids.push_back(static_cast<std::uint64_t>(i));

    validate_distributed_ids(ids, global_n);

    DistributedCellField state(global_n, ids, 1, "phi");
    for (std::size_t i = 0; i < state.local_size(); ++i)
        state(i) = static_cast<double>(state.global_ids()[i] + 1);

    // Deterministic reduction is exercised on the owned contribution only.
    const double local_sum = [&]() {
        double s = 0.0;
        for (std::size_t i = 0; i < state.local_size(); ++i) s += state(i);
        return s;
    }();
    const double global_sum = mpi_deterministic_sum(local_sum);
    assert(std::abs(global_sum - 36.0) < 1e-14);

    // The distributed halo layer must expose the opposite rank's interface
    // value without requiring a replicated Field.
    cfdx::core::Mesh mesh;
    // Build a minimal two-cell topology: face 0 is internal and connects cells 0/1.
    mesh.points().add({0.0,0.0,0.0});
    mesh.points().add({1.0,0.0,0.0});
    mesh.points().add({0.0,1.0,0.0});
    mesh.points().add({0.0,0.0,1.0});
    mesh.points().add({1.0,1.0,1.0});
    mesh.faces().set_offsets({0,3,6});
    mesh.faces().set_indices({0,1,2, 1,3,4});
    mesh.ownership().set_owner({0,0});
    mesh.ownership().set_neighbour({1,-1});
    mesh.cells().set_offsets({0,1,2});
    mesh.cells().set_faces({0,1});
    cfdx::core::parallel::Partition p;
    p.n_parts = 2;
    p.cell_rank = {0,1};
    p.face_owner_rank = {0,0};
    p.face_ghost_rank = {1,-1};

    auto halo = build_distributed_halo(mesh, p, state);
    const auto received = exchange_distributed_cell_halo(state, halo);
    assert(received.size() == 1);
    assert(std::abs(received.front() - (rank == 0 ? 2.0 : 1.0)) < 1e-14);

#ifdef CFDX_ENABLE_PARALLEL_HDF5
    const std::string path = "phase5_distributed_checkpoint.h5";
    write_distributed_checkpoint(path, state);
    MPI_Barrier(MPI_COMM_WORLD);

    DistributedCellField restored(global_n, ids, 1, "phi");
    read_distributed_checkpoint(path, restored);
    for (std::size_t i = 0; i < restored.local_size(); ++i)
        assert(std::abs(restored(i) - state(i)) < 1e-14);

    // Rank-independent restart: deliberately permute the local order. The
    // reader addresses the file by persistent global ids, not by rank/order.
    std::reverse(ids.begin(), ids.end());
    DistributedCellField permuted(global_n, ids, 1, "phi");
    read_distributed_checkpoint(path, permuted);
    for (std::size_t i = 0; i < permuted.local_size(); ++i)
        assert(std::abs(permuted(i) - (static_cast<double>(permuted.global_ids()[i]) + 1.0)) < 1e-14);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) std::remove(path.c_str());
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    mpi_finalize();
    return 0;
}
