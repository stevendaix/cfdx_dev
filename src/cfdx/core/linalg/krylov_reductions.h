#pragma once

#include "vector.h"
#include "mixed_precision.h"
#if CFDX_HAS_MPI
#include "cfdx/core/parallel/mpi_utils.h"
#endif
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <limits>

namespace cfdx::core {

#if CFDX_HAS_MPI
using KrylovMPIComm = MPI_Comm;
#else
using KrylovMPIComm = int;
#endif

struct KrylovReductionPolicy {
    bool mpi_enabled = false;
    bool deterministic = false;
    std::size_t local_begin = 0;
    std::size_t local_end = std::numeric_limits<std::size_t>::max();
#if CFDX_HAS_MPI
    KrylovMPIComm comm = MPI_COMM_WORLD;
#else
    KrylovMPIComm comm = 0;
#endif

    std::size_t end_for(std::size_t n) const {
        return std::min(local_end, n);
    }

    bool active(std::size_t n) const {
#if CFDX_HAS_MPI
        return mpi_enabled && cfdx::core::parallel::mpi_size(comm) > 1 &&
               local_begin <= end_for(n);
#else
        (void)n;
        return false;
#endif
    }
};

inline double krylov_sum(double local, const KrylovReductionPolicy& policy) {
#if CFDX_HAS_MPI
    if (!policy.mpi_enabled || cfdx::core::parallel::mpi_size(policy.comm) <= 1)
        return local;
    return policy.deterministic
        ? cfdx::core::parallel::mpi_deterministic_sum(local, policy.comm)
        : cfdx::core::parallel::mpi_allreduce_sum(local, policy.comm);
#else
    (void)policy;
    return local;
#endif
}

inline double krylov_dot(const Vector& a, const Vector& b,
                         SolverPrecision reduction_precision = SolverPrecision::FP64,
                         const KrylovReductionPolicy& policy = {}) {
    const std::size_t n = std::min(a.size(), b.size());
    const std::size_t begin = std::min(policy.local_begin, n);
    const std::size_t end = policy.end_for(n);
    double local = 0.0;
    for (std::size_t i = begin; i < end; ++i) {
        local += reduction_precision == SolverPrecision::FP32
            ? static_cast<double>(static_cast<float>(a(i)) * static_cast<float>(b(i)))
            : a(i) * b(i);
    }
    return krylov_sum(local, policy);
}

inline double krylov_norm2(const Vector& a,
                           SolverPrecision reduction_precision = SolverPrecision::FP64,
                           const KrylovReductionPolicy& policy = {}) {
    const std::size_t begin = std::min(policy.local_begin, a.size());
    const std::size_t end = policy.end_for(a.size());
    double local = 0.0;
    for (std::size_t i = begin; i < end; ++i) {
        local += reduction_precision == SolverPrecision::FP32
            ? static_cast<double>(static_cast<float>(a(i)) * static_cast<float>(a(i)))
            : a(i) * a(i);
    }
    return std::sqrt(std::max(0.0, krylov_sum(local, policy)));
}

} // namespace cfdx::core
