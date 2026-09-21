#pragma once

#include <cstddef>
#include <cstring>
#include <future>
#include <utility>
#include <vector>

namespace cfdx::runtime::ooc {

class AsyncTransfer {
public:
    template<class T>
    std::future<void> h2d(const std::vector<T>& host, std::vector<T>& device) const {
        return std::async(std::launch::async, [&host, &device] {
            device = host;
        });
    }

    template<class T>
    std::future<void> d2h(const std::vector<T>& device, std::vector<T>& host) const {
        return std::async(std::launch::async, [&device, &host] {
            host = device;
        });
    }
};

template<class T>
class DoubleBuffer {
public:
    explicit DoubleBuffer(std::size_t size = 0) : buffers_{std::vector<T>(size), std::vector<T>(size)} {}
    std::vector<T>& acquire(std::size_t slot) { return buffers_[slot & 1U]; }
    const std::vector<T>& acquire(std::size_t slot) const { return buffers_[slot & 1U]; }
private:
    std::vector<T> buffers_[2];
};

} // namespace cfdx::runtime::ooc
