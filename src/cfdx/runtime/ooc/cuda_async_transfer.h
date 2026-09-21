#pragma once

#include "../gpu/cuda_backend.h"
#include <cuda_runtime.h>
#include <cstddef>
#include <stdexcept>

namespace cfdx::runtime::ooc {

class CudaAsyncTransfer {
public:
    explicit CudaAsyncTransfer(cudaStream_t stream) : stream_(stream) {
        if (!stream_) throw std::invalid_argument("CudaAsyncTransfer: null stream");
    }

    void h2d(cfdx::runtime::gpu::CudaDeviceBuffer& dst,
             const void* src, std::size_t bytes) const {
        cfdx::runtime::gpu::async_copy_h2d(dst, src, bytes, stream_);
    }

    void d2h(const cfdx::runtime::gpu::CudaDeviceBuffer& src,
             void* dst, std::size_t bytes) const {
        cfdx::runtime::gpu::async_copy_d2h(src, dst, bytes, stream_);
    }

    void synchronize() const {
        if (cudaStreamSynchronize(stream_) != cudaSuccess)
            throw std::runtime_error("CudaAsyncTransfer synchronization failed");
    }

private:
    cudaStream_t stream_;
};

} // namespace cfdx::runtime::ooc
