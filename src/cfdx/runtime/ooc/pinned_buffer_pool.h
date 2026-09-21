#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#ifdef CFDX_ENABLE_GPU
#include "cuda_pinned_buffer_pool.h"
#endif

namespace cfdx::runtime::ooc {

// Host staging pool used by CPU-only builds. The name PinnedBufferPool is
// reserved for the CUDA implementation below; std::vector memory is pageable.
class StagingBufferPool {
public:
    struct Buffer {
        std::vector<std::uint8_t> bytes;
        bool in_use = false;
    };

    StagingBufferPool() = default;
    explicit StagingBufferPool(std::size_t capacity_bytes, std::size_t buffer_bytes)
        : buffer_bytes_(buffer_bytes)
    {
        if (capacity_bytes == 0 || buffer_bytes == 0)
            throw std::invalid_argument("StagingBufferPool: capacities must be > 0");
        const std::size_t count = capacity_bytes / buffer_bytes;
        if (count == 0)
            throw std::invalid_argument("StagingBufferPool: buffer larger than capacity");
        buffers_.resize(count);
        for (auto& b : buffers_) b.bytes.resize(buffer_bytes);
    }

    static constexpr bool is_pinned() noexcept { return false; }

    Buffer* acquire() {
        for (auto& b : buffers_) {
            if (!b.in_use) {
                b.in_use = true;
                return &b;
            }
        }
        return nullptr;
    }

    void release(Buffer* b) {
        if (!b) return;
        for (auto& owned : buffers_) {
            if (&owned == b) {
                owned.in_use = false;
                return;
            }
        }
        throw std::invalid_argument("StagingBufferPool: foreign buffer");
    }

    std::size_t size() const noexcept { return buffers_.size(); }
    std::size_t buffer_bytes() const noexcept { return buffer_bytes_; }

private:
    std::size_t buffer_bytes_ = 0;
    std::vector<Buffer> buffers_;
};

#ifdef CFDX_ENABLE_GPU
using PinnedBufferPool = CudaPinnedBufferPool;
#else
// Explicitly pageable CPU emulation: callers cannot accidentally assume that
// asynchronous CUDA transfers are backed by pinned host memory.
using PinnedBufferPool = StagingBufferPool;
#endif

} // namespace cfdx::runtime::ooc
