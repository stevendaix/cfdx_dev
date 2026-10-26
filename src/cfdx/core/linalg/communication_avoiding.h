#pragma once
#include "cfdx/core/linalg/vector.h"
#if defined(CFDX_HAS_MPI) && CFDX_HAS_MPI
#include "cfdx/core/parallel/mpi_utils.h"
#endif
#include <algorithm>
#include <cmath>
#include <cstddef>

namespace cfdx::core {

struct ReductionPacket {
    double dot = 0.0;
    double norm2 = 0.0;
    double max_abs = 0.0;
};

inline ReductionPacket fused_reduction(const Vector& a, const Vector& b) {
    ReductionPacket r;
    const std::size_t n = a.size();
    if (b.size() != n) return r;
    for (std::size_t i = 0; i < n; ++i) {
        r.dot += a(i) * b(i);
        r.norm2 += a(i) * a(i);
        r.max_abs = std::max(r.max_abs, std::abs(a(i)));
    }
    return r;
}

#if defined(CFDX_HAS_MPI) && CFDX_HAS_MPI
// Batch the two additive Krylov reductions into one MPI collective.
// max_abs is reduced separately because it requires MPI_MAX rather than SUM.
// This API is deliberately independent of a solver so distributed Krylov
// implementations can adopt it without coupling linear algebra to physics.
inline ReductionPacket mpi_fused_reduction(const ReductionPacket& local,
                                           MPI_Comm comm = MPI_COMM_WORLD) {
    double local_add[2] = {local.dot, local.norm2};
    double global_add[2] = {0.0, 0.0};
    MPI_Allreduce(local_add, global_add, 2, MPI_DOUBLE, MPI_SUM, comm);

    ReductionPacket global;
    global.dot = global_add[0];
    global.norm2 = global_add[1];
    MPI_Allreduce(&local.max_abs, &global.max_abs, 1, MPI_DOUBLE, MPI_MAX, comm);
    return global;
}

inline ReductionPacket mpi_deterministic_fused_reduction(
    const ReductionPacket& local, MPI_Comm comm = MPI_COMM_WORLD) {
    const int rank = cfdx::core::parallel::mpi_rank(comm);
    const int size = cfdx::core::parallel::mpi_size(comm);

    // Gather all three components together. Rank-ordered accumulation makes
    // the additive reductions reproducible and the max reduction deterministic,
    // while avoiding a separate MPI_MAX operation dispatch on the deterministic
    // path. It also reduces the collective count from three to two.
    const double local_values[3] = {local.dot, local.norm2, local.max_abs};
    std::vector<double> gathered(
        rank == 0 ? static_cast<std::size_t>(size) * 3U : 0U);
    MPI_Gather(local_values, 3, MPI_DOUBLE,
               rank == 0 ? gathered.data() : nullptr, 3, MPI_DOUBLE, 0, comm);

    double global_values[3] = {0.0, 0.0, 0.0};
    if (rank == 0) {
        for (int r = 0; r < size; ++r) {
            const auto offset = static_cast<std::size_t>(r) * 3U;
            global_values[0] += gathered[offset];
            global_values[1] += gathered[offset + 1U];
            global_values[2] = std::max(global_values[2], gathered[offset + 2U]);
        }
    }
    MPI_Bcast(global_values, 3, MPI_DOUBLE, 0, comm);

    ReductionPacket global;
    global.dot = global_values[0];
    global.norm2 = global_values[1];
    global.max_abs = global_values[2];
    return global;
}

inline ReductionPacket fused_reduction(const Vector& a, const Vector& b,
                                       MPI_Comm comm) {
    return mpi_fused_reduction(fused_reduction(a, b), comm);
}

#endif // CFDX_HAS_MPI

} // namespace cfdx::core
