#include "cfdx/core/parallel/mpi_utils.h"

#include <iostream>

int main(int argc, char** argv) {
    cfdx::core::parallel::mpi_init(&argc, &argv);

    const int size = cfdx::core::parallel::mpi_size();
    const bool valid_size = size == 2;
    if (!valid_size && cfdx::core::parallel::mpi_rank() == 0) {
        std::cerr << "MPI runtime control expected exactly two ranks\n";
    }

    cfdx::core::parallel::mpi_finalize();
    return valid_size ? 0 : 1;
}
