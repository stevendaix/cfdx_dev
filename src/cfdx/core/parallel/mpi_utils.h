// M0.11-T01 — MPI utilities
//
// Basic MPI utilities for CFDX parallel execution.

#pragma once

#include <mpi.h>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>

namespace cfdx {
namespace core {
namespace parallel {

// Initialize MPI if not already initialized
inline void mpi_init(int* argc = nullptr, char*** argv = nullptr) {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
        MPI_Init(argc, argv);
    }
}

// Finalize MPI if we initialized it
inline void mpi_finalize() {
    int initialized = 0;
    int finalized = 0;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    if (initialized && !finalized) {
        MPI_Finalize();
    }
}

// Get communicator rank
inline int mpi_rank(MPI_Comm comm = MPI_COMM_WORLD) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    return rank;
}

// Get communicator size
inline int mpi_size(MPI_Comm comm = MPI_COMM_WORLD) {
    int size = 0;
    MPI_Comm_size(comm, &size);
    return size;
}

// Check if we're in a parallel run
inline bool mpi_is_parallel(MPI_Comm comm = MPI_COMM_WORLD) {
    return mpi_size(comm) > 1;
}

// Barrier synchronization
inline void mpi_barrier(MPI_Comm comm = MPI_COMM_WORLD) {
    MPI_Barrier(comm);
}

// All-reduce sum for scalars
inline double mpi_allreduce_sum(double value, MPI_Comm comm = MPI_COMM_WORLD) {
    double result = 0.0;
    MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_SUM, comm);
    return result;
}

inline int mpi_allreduce_sum(int value, MPI_Comm comm = MPI_COMM_WORLD) {
    int result = 0;
    MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_SUM, comm);
    return result;
}

inline std::int64_t mpi_allreduce_sum(std::int64_t value, MPI_Comm comm = MPI_COMM_WORLD) {
    std::int64_t result = 0;
    MPI_Allreduce(&value, &result, 1, MPI_INT64_T, MPI_SUM, comm);
    return result;
}

// All-reduce max for scalars
inline double mpi_allreduce_max(double value, MPI_Comm comm = MPI_COMM_WORLD) {
    double result = 0.0;
    MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_MAX, comm);
    return result;
}

// All-reduce min for scalars
inline double mpi_allreduce_min(double value, MPI_Comm comm = MPI_COMM_WORLD) {
    double result = 0.0;
    MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_MIN, comm);
    return result;
}

// All-gather for vectors
template <typename T>
inline std::vector<T> mpi_allgather(const std::vector<T>& local_data, MPI_Comm comm = MPI_COMM_WORLD) {
    int rank = mpi_rank(comm);
    int size = mpi_size(comm);
    
    std::vector<int> recv_counts(size);
    std::vector<int> displs(size);
    
    int local_count = static_cast<int>(local_data.size());
    MPI_Allgather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);
    
    displs[0] = 0;
    for (int i = 1; i < size; ++i) {
        displs[i] = displs[i - 1] + recv_counts[i - 1];
    }
    
    int total_count = displs[size - 1] + recv_counts[size - 1];
    std::vector<T> result(total_count);
    
    MPI_Allgatherv(local_data.data(), local_count, MPI_BYTE,
                   result.data(), recv_counts.data(), displs.data(), MPI_BYTE, comm);
    
    return result;
}

// Broadcast from root
template <typename T>
inline void mpi_bcast(T* data, int count, int root = 0, MPI_Comm comm = MPI_COMM_WORLD) {
    MPI_Datatype type;
    if constexpr (std::is_same_v<T, double>) type = MPI_DOUBLE;
    else if constexpr (std::is_same_v<T, int>) type = MPI_INT;
    else if constexpr (std::is_same_v<T, std::int64_t>) type = MPI_INT64_T;
    else if constexpr (std::is_same_v<T, std::uint64_t>) type = MPI_UINT64_T;
    else if constexpr (std::is_same_v<T, float>) type = MPI_FLOAT;
    else if constexpr (std::is_same_v<T, char>) type = MPI_CHAR;
    else {
        throw std::runtime_error("mpi_bcast: unsupported type");
    }
    MPI_Bcast(data, count, type, root, comm);
}

// Broadcast string
inline std::string mpi_bcast_string(const std::string& str, int root = 0, MPI_Comm comm = MPI_COMM_WORLD) {
    int len = static_cast<int>(str.size());
    mpi_bcast(&len, 1, root, comm);
    
    std::string result(len, '\0');
    if (len > 0) {
        mpi_bcast(&result[0], len, root, comm);
    }
    return result;
}

// Reduce to root (sum)
template <typename T>
inline void mpi_reduce_sum(const T* sendbuf, T* recvbuf, int count, int root = 0, MPI_Comm comm = MPI_COMM_WORLD) {
    MPI_Datatype type;
    if constexpr (std::is_same_v<T, double>) type = MPI_DOUBLE;
    else if constexpr (std::is_same_v<T, int>) type = MPI_INT;
    else if constexpr (std::is_same_v<T, std::int64_t>) type = MPI_INT64_T;
    else if constexpr (std::is_same_v<T, float>) type = MPI_FLOAT;
    else {
        throw std::runtime_error("mpi_reduce_sum: unsupported type");
    }
    MPI_Reduce(const_cast<T*>(sendbuf), recvbuf, count, type, MPI_SUM, root, comm);
}

// Print only from rank 0
inline void mpi_print_rank0(const std::string& msg, MPI_Comm comm = MPI_COMM_WORLD) {
    if (mpi_rank(comm) == 0) {
        std::cout << msg << std::flush;
    }
}

inline void mpi_print_rank0(const char* msg, MPI_Comm comm = MPI_COMM_WORLD) {
    if (mpi_rank(comm) == 0) {
        std::cout << msg << std::flush;
    }
}

// Abort with error message
[[noreturn]] inline void mpi_abort(const std::string& msg, int code = 1, MPI_Comm comm = MPI_COMM_WORLD) {
    mpi_print_rank0("MPI ABORT: " + msg + "\n", comm);
    MPI_Abort(comm, code);
    std::exit(code);
}

}  // namespace parallel
}  // namespace core
}  // namespace cfdx
