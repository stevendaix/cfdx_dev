#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::runtime::ooc {

class PinnedBufferPool {
public:
    struct Buffer {
        std::vector<std::uint8_t> bytes;
        bool in_use = false;
    };

    PinnedBufferPool() = default;
    explicit PinnedBufferPool(std::size_t capacity_bytes, std::size_t buffer_bytes)
        : capacity_(capacity_bytes), buffer_bytes_(buffer_bytes)
    {
        if (capacity_bytes == 0 || buffer_bytes == 0)
            throw std::invalid_argument("PinnedBufferPool: capacities must be > 0");
        const std::size_t count = capacity_bytes / buffer_bytes;
        if (count == 0) throw std::invalid_argument("PinnedBufferPool: buffer larger than capacity");
        buffers_.resize(count);
        for (auto& b : buffers_) b.bytes.resize(buffer_bytes);
    }

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
        throw std::invalid_argument("PinnedBufferPool: buffer does not belong to pool");
    }

    std::size_t size() const noexcept { return buffers_.size(); }
    std::size_t buffer_bytes() const noexcept { return buffer_bytes_; }

private:
    std::size_t capacity_ = 0;
    std::size_t buffer_bytes_ = 0;
    std::vector<Buffer> buffers_;
};

} // namespace cfdx::runtime::ooc
