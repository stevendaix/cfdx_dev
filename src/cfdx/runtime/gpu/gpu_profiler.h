#pragma once

#ifdef CFDX_ENABLE_GPU
#include <cuda_runtime.h>
#include <stdexcept>
#endif

#include <cstddef>

namespace cfdx::runtime::gpu {

// Lightweight CUDA event profiler. Timing is explicit and GPU-side; no CPU
// fallback is used. The measured interval is valid after stop() synchronizes.
class GpuProfiler {
public:
    GpuProfiler() {
#ifdef CFDX_ENABLE_GPU
        if (cudaEventCreate(&start_) != cudaSuccess ||
            cudaEventCreate(&stop_) != cudaSuccess) {
            if (start_) cudaEventDestroy(start_);
            start_ = nullptr;
            if (stop_) cudaEventDestroy(stop_);
            stop_ = nullptr;
            throw std::runtime_error("CUDA event creation failed");
        }
#endif
    }

    ~GpuProfiler() {
#ifdef CFDX_ENABLE_GPU
        if (start_) cudaEventDestroy(start_);
        if (stop_) cudaEventDestroy(stop_);
#endif
    }

    GpuProfiler(const GpuProfiler&) = delete;
    GpuProfiler& operator=(const GpuProfiler&) = delete;

    void start(cudaStream_t stream = nullptr) {
#ifdef CFDX_ENABLE_GPU
        if (cudaEventRecord(start_, stream) != cudaSuccess)
            throw std::runtime_error("CUDA profiler start event failed");
#else
        (void)stream;
        throw std::runtime_error("GPU profiling requested but CUDA is disabled");
#endif
    }

    void stop(cudaStream_t stream = nullptr) {
#ifdef CFDX_ENABLE_GPU
        if (cudaEventRecord(stop_, stream) != cudaSuccess)
            throw std::runtime_error("CUDA profiler stop event failed");
        if (cudaEventSynchronize(stop_) != cudaSuccess)
            throw std::runtime_error("CUDA profiler event synchronization failed");
        float ms = 0.0F;
        if (cudaEventElapsedTime(&ms, start_, stop_) != cudaSuccess)
            throw std::runtime_error("CUDA profiler elapsed-time query failed");
        elapsed_ms_ = static_cast<double>(ms);
#else
        (void)stream;
        throw std::runtime_error("GPU profiling requested but CUDA is disabled");
#endif
    }

    double elapsed_ms() const noexcept { return elapsed_ms_; }

private:
#ifdef CFDX_ENABLE_GPU
    cudaEvent_t start_ = nullptr;
    cudaEvent_t stop_ = nullptr;
#endif
    double elapsed_ms_ = 0.0;
};

} // namespace cfdx::runtime::gpu
