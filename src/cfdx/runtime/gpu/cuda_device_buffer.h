#pragma once

#include "device_buffer.h"
#include <cuda_runtime.h>
#include <stdexcept>
#include <string>

namespace cfdx::runtime::gpu {

class CudaDeviceBuffer final : public DeviceBuffer {
public:
    CudaDeviceBuffer() = default;
    explicit CudaDeviceBuffer(std::size_t bytes) { resize(bytes); }
    ~CudaDeviceBuffer() override { release(); }

    void resize(std::size_t bytes) override {
        if (bytes == size_) return;
        release();
        if (bytes == 0) return;
        if (cudaMalloc(&ptr_, bytes) != cudaSuccess)
            throw std::runtime_error("cudaMalloc failed");
        size_ = bytes;
    }

    void* data() noexcept override { return ptr_; }
    const void* data() const noexcept override { return ptr_; }
    std::size_t size_bytes() const noexcept override { return size_; }

    void copy_from_host(const void* src, std::size_t bytes) override {
        check(bytes);
        if (cudaMemcpy(ptr_, src, bytes, cudaMemcpyHostToDevice) != cudaSuccess)
            throw std::runtime_error("cudaMemcpy H2D failed");
    }

    void copy_to_host(void* dst, std::size_t bytes) const override {
        check(bytes);
        if (cudaMemcpy(dst, ptr_, bytes, cudaMemcpyDeviceToHost) != cudaSuccess)
            throw std::runtime_error("cudaMemcpy D2H failed");
    }

private:
    void check(std::size_t bytes) const {
        if (bytes > size_) throw std::out_of_range("CUDA device buffer too small");
    }
    void release() noexcept {
        if (ptr_) cudaFree(ptr_);
        ptr_ = nullptr; size_ = 0;
    }

    void* ptr_ = nullptr;
    std::size_t size_ = 0;
};

} // namespace cfdx::runtime::gpu
