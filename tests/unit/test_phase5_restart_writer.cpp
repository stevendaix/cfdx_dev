#include "cfdx/core/parallel/distributed_execution.h"

#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    using namespace cfdx::core::parallel;
    mpi_init(&argc, &argv);
    const int rank = mpi_rank();
    const int size = mpi_size();

    if (size != 2) {
        if (rank == 0)
            std::fprintf(stderr, "test_phase5_restart_writer requires exactly 2 MPI ranks\n");
        mpi_finalize();
        return 2;
    }

    constexpr std::size_t global_n = 12;
    std::vector<std::uint64_t> ids;
    for (std::size_t i = static_cast<std::size_t>(rank) * 6;
         i < static_cast<std::size_t>(rank + 1) * 6; ++i)
        ids.push_back(static_cast<std::uint64_t>(i));

    validate_distributed_ids(ids, global_n);
    DistributedCellField state(global_n, ids, 3, "phi");
    for (std::size_t i = 0; i < state.local_size(); ++i) {
        const double gid = static_cast<double>(state.global_ids()[i]);
        state(i, 0) = 1.0 + gid;
        state(i, 1) = 100.0 + 2.0 * gid;
        state(i, 2) = -3.0 - gid;
    }

    const std::string path = "phase5_true_n_to_m_checkpoint.h5";
    write_distributed_checkpoint(path, state);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
        std::printf("phase5 N-to-M writer: wrote %zu cells per rank\n", state.local_size());

    mpi_finalize();
    return 0;
}
