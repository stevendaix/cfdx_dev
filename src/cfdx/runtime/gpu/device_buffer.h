#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::runtime::gpu {

// Backend-neutral device buffer contract. The implementation is provided by
// the CUDA backend when CFDX_ENABLE_GPU is enabled.
class DeviceBuffer {
public:
    DeviceBuffer() = default;
    explicit DeviceBuffer(std::size_t bytes) { resize(bytes); }
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;
    DeviceBuffer(DeviceBuffer&&) noexcept = default;
    DeviceBuffer& operator=(DeviceBuffer&&) noexcept = default;
    virtual ~DeviceBuffer() = default;

    virtual void resize(std::size_t bytes) = 0;
    virtual void* data() noexcept = 0;
    virtual const void* data() const noexcept = 0;
    virtual std::size_t size_bytes() const noexcept = 0;
    virtual void copy_from_host(const void* src, std::size_t bytes) = 0;
    virtual void copy_to_host(void* dst, std::size_t bytes) const = 0;
};

class HostEmulatedDeviceBuffer final : public DeviceBuffer {
public:
    HostEmulatedDeviceBuffer() = default;
    explicit HostEmulatedDeviceBuffer(std::size_t bytes) : data_(bytes) {}

    void resize(std::size_t bytes) override { data_.resize(bytes); }
    void* data() noexcept override { return data_.data(); }
    const void* data() const noexcept override { return data_.data(); }
    std::size_t size_bytes() const noexcept override { return data_.size(); }

    void copy_from_host(const void* src, std::size_t bytes) override {
        if (bytes > data_.size()) throw std::out_of_range("device buffer too small");
        const auto* p = static_cast<const std::uint8_t*>(src);
        std::copy(p, p + bytes, data_.begin());
    }

    void copy_to_host(void* dst, std::size_t bytes) const override {
        if (bytes > data_.size()) throw std::out_of_range("device buffer too small");
        const auto* p = data_.data();
        std::copy(p, p + bytes, static_cast<std::uint8_t*>(dst));
    }

private:
    std::vector<std::uint8_t> data_;
};

} // namespace cfdx::runtime::gpu
