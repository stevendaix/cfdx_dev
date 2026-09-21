#pragma once

#include "cuda_device_buffer.h"
#include <cuda_runtime.h>
#include <cstddef>
#include <stdexcept>

namespace cfdx::runtime::gpu {

class CudaStream {
public:
    CudaStream() {
        if (cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking) != cudaSuccess)
            throw std::runtime_error("cudaStreamCreateWithFlags failed");
    }
    ~CudaStream() { if (stream_) cudaStreamDestroy(stream_); }
    CudaStream(const CudaStream&) = delete;
    CudaStream& operator=(const CudaStream&) = delete;
    cudaStream_t get() const noexcept { return stream_; }
    void synchronize() const {
        if (cudaStreamSynchronize(stream_) != cudaSuccess)
            throw std::runtime_error("cudaStreamSynchronize failed");
    }
private:
    cudaStream_t stream_ = nullptr;
};

inline void async_copy_h2d(CudaDeviceBuffer& dst, const void* src,
                           std::size_t bytes, cudaStream_t stream) {
    if (bytes > dst.size_bytes()) throw std::out_of_range("H2D transfer exceeds device buffer");
    if (cudaMemcpyAsync(dst.data(), src, bytes, cudaMemcpyHostToDevice, stream) != cudaSuccess)
        throw std::runtime_error("cudaMemcpyAsync H2D failed");
}

inline void async_copy_d2h(const CudaDeviceBuffer& src, void* dst,
                           std::size_t bytes, cudaStream_t stream) {
    if (bytes > src.size_bytes()) throw std::out_of_range("D2H transfer exceeds device buffer");
    if (cudaMemcpyAsync(dst, src.data(), bytes, cudaMemcpyDeviceToHost, stream) != cudaSuccess)
        throw std::runtime_error("cudaMemcpyAsync D2H failed");
}

} // namespace cfdx::runtime::gpu
