#pragma once

#include <cuda_runtime.h>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace cfdx::runtime::ooc {

class CudaPinnedBufferPool {
public:
    struct Buffer {
        void* data = nullptr;
        std::size_t bytes = 0;
        bool in_use = false;
    };

    CudaPinnedBufferPool(std::size_t capacity_bytes, std::size_t buffer_bytes)
        : buffer_bytes_(buffer_bytes)
    {
        if (capacity_bytes == 0 || buffer_bytes == 0 || capacity_bytes < buffer_bytes)
            throw std::invalid_argument("CudaPinnedBufferPool: invalid capacity");
        const std::size_t count = capacity_bytes / buffer_bytes;
        buffers_.resize(count);
        for (auto& b : buffers_) {
            if (cudaMallocHost(&b.data, buffer_bytes) != cudaSuccess) {
                release_all();
                throw std::runtime_error("cudaMallocHost failed");
            }
            b.bytes = buffer_bytes;
        }
    }

    ~CudaPinnedBufferPool() { release_all(); }

    CudaPinnedBufferPool(const CudaPinnedBufferPool&) = delete;
    CudaPinnedBufferPool& operator=(const CudaPinnedBufferPool&) = delete;

    static constexpr bool is_pinned() noexcept { return true; }

    std::size_t size() const noexcept { return buffers_.size(); }
    std::size_t buffer_bytes() const noexcept { return buffer_bytes_; }

    Buffer* acquire() {
        for (auto& b : buffers_) if (!b.in_use) {
            b.in_use = true;
            return &b;
        }
        return nullptr;
    }

    void release(Buffer* b) {
        if (!b) return;
        for (auto& owned : buffers_) if (&owned == b) {
            owned.in_use = false;
            return;
        }
        throw std::invalid_argument("CudaPinnedBufferPool: foreign buffer");
    }

private:
    void release_all() noexcept {
        for (auto& b : buffers_) {
            if (b.data) cudaFreeHost(b.data);
            b.data = nullptr;
            b.in_use = false;
        }
    }
    std::size_t buffer_bytes_ = 0;
    std::vector<Buffer> buffers_;
};

} // namespace cfdx::runtime::ooc
