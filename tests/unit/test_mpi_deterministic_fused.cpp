#include "cfdx/core/linalg/communication_avoiding.h"
#include "cfdx/core/parallel/mpi_utils.h"

#include <cmath>
#include <iostream>

int main(int argc, char** argv) {
    cfdx::core::parallel::mpi_init(&argc, &argv);

    const int rank = cfdx::core::parallel::mpi_rank();
    const int size = cfdx::core::parallel::mpi_size();
    const double rank_value = static_cast<double>(rank) + 1.0;
    const double size_value = static_cast<double>(size);

    const cfdx::core::ReductionPacket local{
        rank_value,
        rank_value * rank_value,
        2.0 * static_cast<double>(rank) + 1.0,
    };
    const auto global = cfdx::core::mpi_deterministic_fused_reduction(local);

    const double expected_dot = size_value * (size_value + 1.0) / 2.0;
    const double expected_norm2 =
        size_value * (size_value + 1.0) * (2.0 * size_value + 1.0) / 6.0;
    const double expected_max_abs = 2.0 * size_value - 1.0;
    const bool ok =
        std::abs(global.dot - expected_dot) < 1e-14 &&
        std::abs(global.norm2 - expected_norm2) < 1e-14 &&
        std::abs(global.max_abs - expected_max_abs) < 1e-14;

    if (rank == 0 && !ok) {
        std::cerr << "deterministic fused reduction mismatch\n";
    }

    cfdx::core::parallel::mpi_finalize();
    return ok ? 0 : 1;
}
