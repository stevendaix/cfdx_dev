#pragma once

#include "cfdx/runtime/gpu/cuda_backend.h"
#include <array>
#include <cstddef>
#include <stdexcept>

namespace cfdx::runtime::gpu {

// Two-slot CUDA pipeline: each slot owns a device buffer and a non-blocking
// stream. A slot must be synchronized before it is reused.
class CudaDoubleBuffer {
public:
    explicit CudaDoubleBuffer(std::size_t bytes) : bytes_(bytes) {
        if (bytes == 0) throw std::invalid_argument("CudaDoubleBuffer: bytes must be > 0");
        for (auto& slot : slots_) slot.buffer.resize(bytes);
    }

    CudaDoubleBuffer(const CudaDoubleBuffer&) = delete;
    CudaDoubleBuffer& operator=(const CudaDoubleBuffer&) = delete;

    std::size_t bytes() const noexcept { return bytes_; }
    static constexpr std::size_t slot_count() noexcept { return 2; }

    CudaDeviceBuffer& buffer(std::size_t slot) {
        return slots_.at(slot & 1U).buffer;
    }
    CudaStream& stream(std::size_t slot) {
        return slots_.at(slot & 1U).stream;
    }

    void synchronize(std::size_t slot) { stream(slot).synchronize(); }

private:
    struct Slot {
        CudaDeviceBuffer buffer;
        CudaStream stream;
    };
    std::array<Slot, 2> slots_;
    std::size_t bytes_;
};

} // namespace cfdx::runtime::gpu
