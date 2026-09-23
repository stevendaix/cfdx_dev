#include "cfdx/core/linalg/communication_avoiding.h"
#include "cfdx/core/parallel/mpi_utils.h"
#include "common/test_harness.h"

#include <cmath>

using namespace cfdx::core;

int main(int argc, char** argv) {
    parallel::mpi_init(&argc, &argv);

    const int rank = parallel::mpi_rank();
    const int size = parallel::mpi_size();

    Vector a(2), b(2);
    a(0) = static_cast<double>(rank + 1);
    a(1) = 2.0 * static_cast<double>(rank + 1);
    b(0) = 3.0;
    b(1) = -1.0;

    const ReductionPacket local = fused_reduction(a, b);
    const ReductionPacket global = mpi_fused_reduction(local);

    const double rank_sum = static_cast<double>(size * (size + 1)) / 2.0;
    EXPECT_NEAR(global.dot, rank_sum, 1e-12);
    EXPECT_NEAR(global.norm2, 5.0 * (size * (size + 1)) * 0.5, 1e-12);
    EXPECT_NEAR(global.max_abs, 2.0 * static_cast<double>(size), 1e-12);

    const ReductionPacket direct = fused_reduction(a, b, MPI_COMM_WORLD);
    EXPECT_NEAR(direct.dot, rank_sum, 1e-12);
    EXPECT_NEAR(direct.norm2, 5.0 * (size * (size + 1)) * 0.5, 1e-12);
    EXPECT_NEAR(direct.max_abs, 2.0 * static_cast<double>(size), 1e-12);

    parallel::mpi_finalize();
    return run_all();
}
