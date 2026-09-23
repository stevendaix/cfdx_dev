#include "cfdx/core/parallel/distributed_execution.h"

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    using namespace cfdx::core::parallel;
    mpi_init(&argc, &argv);
    const int rank = mpi_rank();
    const int size = mpi_size();

    if (size != 3) {
        if (rank == 0)
            std::fprintf(stderr, "test_phase5_restart_reader requires exactly 3 MPI ranks\n");
        mpi_finalize();
        return 2;
    }

    constexpr std::size_t global_n = 12;
    std::vector<std::uint64_t> ids;
    for (std::size_t i = static_cast<std::size_t>(rank);
         i < global_n; i += static_cast<std::size_t>(size))
        ids.push_back(static_cast<std::uint64_t>(i));

    validate_distributed_ids(ids, global_n);
    DistributedCellField state(global_n, ids, 3, "phi");

    const std::string path = "phase5_true_n_to_m_checkpoint.h5";
    read_distributed_checkpoint(path, state);

    bool ok = true;
    for (std::size_t i = 0; i < state.local_size(); ++i) {
        const double gid = static_cast<double>(state.global_ids()[i]);
        ok = ok && std::abs(state(i, 0) - (1.0 + gid)) < 1e-14;
        ok = ok && std::abs(state(i, 1) - (100.0 + 2.0 * gid)) < 1e-14;
        ok = ok && std::abs(state(i, 2) - (-3.0 - gid)) < 1e-14;
    }

    const int local_ok = ok ? 1 : 0;
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    if (!global_ok) {
        std::fprintf(stderr, "rank %d: true N-to-M restart value mismatch\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        std::remove(path.c_str());
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
        std::printf("phase5 N-to-M reader: 2 -> 3 rank restart verified\n");

    mpi_finalize();
    return 0;
}
