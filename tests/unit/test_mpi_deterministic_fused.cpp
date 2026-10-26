#include "cfdx/core/linalg/communication_avoiding.h"
#include "cfdx/core/parallel/mpi_utils.h"

#include <cmath>
#include <iostream>

int main(int argc, char** argv) {
    cfdx::core::parallel::mpi_init(&argc, &argv);

    const int rank = cfdx::core::parallel::mpi_rank();
    const int size = cfdx::core::parallel::mpi_size();

    const cfdx::core::ReductionPacket local{
        static_cast<double>(rank + 1),
        static_cast<double>((rank + 1) * (rank + 1)),
        static_cast<double>(2 * rank + 1),
    };
    const auto global = cfdx::core::mpi_deterministic_fused_reduction(local);

    const bool ok =
        std::abs(global.dot - static_cast<double>(size * (size + 1)) / 2.0) <
            1e-14 &&
        std::abs(global.norm2 -
                 static_cast<double>(size * (size + 1) * (2 * size + 1)) /
                     6.0) < 1e-14 &&
        std::abs(global.max_abs - static_cast<double>(2 * size - 1)) < 1e-14;

    if (rank == 0 && !ok) {
        std::cerr << "deterministic fused reduction mismatch\n";
    }

    cfdx::core::parallel::mpi_finalize();
    return ok ? 0 : 1;
}
