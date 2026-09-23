#include "cfdx/core/parallel/mpi_utils.h"

#include <cassert>
#include <cmath>

int main(int argc, char** argv) {
    cfdx::core::parallel::mpi_init(&argc, &argv);
    const int rank = cfdx::core::parallel::mpi_rank();
    const int size = cfdx::core::parallel::mpi_size();
    const double expected = static_cast<double>(size * (size + 1)) / 2.0;
    const double value = static_cast<double>(rank + 1);

    const double deterministic = cfdx::core::parallel::mpi_deterministic_sum(value);
    assert(std::abs(deterministic - expected) < 1e-14);

    cfdx::core::parallel::mpi_barrier();
    cfdx::core::parallel::mpi_finalize();
    return 0;
}
