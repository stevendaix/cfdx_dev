#pragma once
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/parallel/mpi_utils.h"
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

inline ReductionPacket fused_reduction(const Vector& a, const Vector& b,
                                       MPI_Comm comm) {
    return mpi_fused_reduction(fused_reduction(a, b), comm);
}

} // namespace cfdx::core
