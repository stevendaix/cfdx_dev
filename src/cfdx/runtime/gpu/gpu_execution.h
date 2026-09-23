#pragma once

#include "cfdx/runtime/gpu/gpu_kernels.h"
#ifdef CFDX_ENABLE_GPU
#include "cfdx/runtime/execution/execution_policy.h"
#include "cfdx/runtime/gpu/cuda_backend.h"
#include "cfdx/runtime/gpu/gpu_profiler.h"
#include "cfdx/runtime/ooc/cuda_pinned_buffer_pool.h"
#endif

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::runtime::gpu {

#ifdef CFDX_ENABLE_GPU
namespace detail {

inline std::size_t max_transfer_bytes(std::initializer_list<std::size_t> sizes)
{
    return *std::max_element(sizes.begin(), sizes.end());
}

inline std::size_t checked_mul(std::size_t a, std::size_t b)
{
    if (b != 0 && a > std::numeric_limits<std::size_t>::max() / b)
        throw std::overflow_error("GPU staging allocation size overflow");
    return a * b;
}

// A transfer batch owns pinned staging buffers until the stream is synchronized.
// This is important: releasing a pinned buffer before an asynchronous CUDA copy
// completes would be a use-after-free on the DMA path.
class PinnedTransferBatch {
public:
    PinnedTransferBatch(CudaStream& stream, std::size_t buffer_bytes,
                        std::size_t max_buffers)
        : stream_(stream),
          pool_(checked_mul(buffer_bytes, max_buffers), buffer_bytes)
    {
        if (buffer_bytes == 0 || max_buffers == 0)
            throw std::invalid_argument("PinnedTransferBatch requires non-zero capacity");
    }

    PinnedTransferBatch(const PinnedTransferBatch&) = delete;
    PinnedTransferBatch& operator=(const PinnedTransferBatch&) = delete;

    ~PinnedTransferBatch() { release_noexcept(); }

    void h2d(CudaDeviceBuffer& dst, const void* src, std::size_t bytes)
    {
        if (bytes > pool_.buffer_bytes())
            throw std::out_of_range("pinned H2D transfer exceeds staging buffer");
        auto* staging = pool_.acquire();
        if (!staging)
            throw std::runtime_error("pinned H2D staging pool exhausted");
        try {
            std::memcpy(staging->data, src, bytes);
            async_copy_h2d(dst, staging->data, bytes, stream_.get());
            h2d_staging_.push_back(staging);
        } catch (...) {
            pool_.release(staging);
            throw;
        }
    }

    void synchronize_h2d()
    {
        stream_.synchronize();
        release(h2d_staging_);
    }

    void d2h(const CudaDeviceBuffer& src, void* dst, std::size_t bytes)
    {
        if (bytes > pool_.buffer_bytes())
            throw std::out_of_range("pinned D2H transfer exceeds staging buffer");
        auto* staging = pool_.acquire();
        if (!staging)
            throw std::runtime_error("pinned D2H staging pool exhausted");
        try {
            async_copy_d2h(src, staging->data, bytes, stream_.get());
            d2h_transfers_.push_back({staging, dst, bytes});
        } catch (...) {
            pool_.release(staging);
            throw;
        }
    }

    void synchronize_d2h()
    {
        stream_.synchronize();
        for (const auto& transfer : d2h_transfers_)
            std::memcpy(transfer.destination, transfer.buffer->data, transfer.bytes);
        release(d2h_buffers());
        d2h_transfers_.clear();
    }

private:
    struct D2HTransfer {
        CudaPinnedBufferPool::Buffer* buffer;
        void* destination;
        std::size_t bytes;
    };

    std::vector<CudaPinnedBufferPool::Buffer*>& d2h_buffers()
    {
        d2h_staging_cache_.clear();
        for (const auto& transfer : d2h_transfers_)
            d2h_staging_cache_.push_back(transfer.buffer);
        return d2h_staging_cache_;
    }

    void release(std::vector<CudaPinnedBufferPool::Buffer*>& buffers)
    {
        for (auto* buffer : buffers)
            pool_.release(buffer);
        buffers.clear();
    }

    void release_noexcept() noexcept
    {
        try {
            stream_.synchronize();
            release(h2d_staging_);
            for (const auto& transfer : d2h_transfers_)
                pool_.release(transfer.buffer);
            d2h_transfers_.clear();
        } catch (...) {
            // Destruction cannot safely propagate CUDA errors.
        }
    }

    CudaStream& stream_;
    CudaPinnedBufferPool pool_;
    std::vector<CudaPinnedBufferPool::Buffer*> h2d_staging_;
    std::vector<D2HTransfer> d2h_transfers_;
    std::vector<CudaPinnedBufferPool::Buffer*> d2h_staging_cache_;
};

inline void validate_device(int device_id)
{
    int device_count = 0;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess || device_count <= 0)
        throw std::runtime_error("CUDA execution requested but no CUDA device is available");
    if (device_id < 0 || device_id >= device_count)
        throw std::invalid_argument("requested GPU device id is outside the available device range");
    if (cudaSetDevice(device_id) != cudaSuccess)
        throw std::runtime_error("cudaSetDevice failed");
}

} // namespace detail

struct GpuExecutionMetrics {
    double h2d_ms = 0.0;
    double kernel_ms = 0.0;
    double d2h_ms = 0.0;
    std::size_t h2d_bytes = 0;
    std::size_t d2h_bytes = 0;
};

inline void execute_gradient_cuda_with_metrics(
    const std::vector<double>& phi,
    const std::vector<double>& face_sx,
    const std::vector<double>& face_sy,
    const std::vector<double>& face_sz,
    const std::vector<std::uint32_t>& owner,
    const std::vector<std::int64_t>& neighbour,
    const std::vector<double>& volume,
    std::vector<double>& grad_x,
    std::vector<double>& grad_y,
    std::vector<double>& grad_z,
    int device_id,
    GpuExecutionMetrics* metrics)
{
    if (owner.size() != neighbour.size() ||
        owner.size() != face_sx.size() ||
        owner.size() != face_sy.size() ||
        owner.size() != face_sz.size() ||
        phi.empty() || volume.empty() || owner.empty() ||
        volume.size() != phi.size()) {
        throw std::invalid_argument("invalid GPU gradient array sizes");
    }

    detail::validate_device(device_id);

    CudaDeviceBuffer d_phi(phi.size() * sizeof(double));
    CudaDeviceBuffer d_sx(face_sx.size() * sizeof(double));
    CudaDeviceBuffer d_sy(face_sy.size() * sizeof(double));
    CudaDeviceBuffer d_sz(face_sz.size() * sizeof(double));
    CudaDeviceBuffer d_owner(owner.size() * sizeof(std::uint32_t));
    CudaDeviceBuffer d_neighbour(neighbour.size() * sizeof(std::int64_t));
    CudaDeviceBuffer d_volume(volume.size() * sizeof(double));
    CudaDeviceBuffer d_gx(volume.size() * sizeof(double));
    CudaDeviceBuffer d_gy(volume.size() * sizeof(double));
    CudaDeviceBuffer d_gz(volume.size() * sizeof(double));

    const std::size_t max_transfer = detail::max_transfer_bytes({
        d_phi.size_bytes(), d_sx.size_bytes(), d_sy.size_bytes(), d_sz.size_bytes(),
        d_owner.size_bytes(), d_neighbour.size_bytes(), d_volume.size_bytes(),
        d_gx.size_bytes(), d_gy.size_bytes(), d_gz.size_bytes()});

    CudaStream stream;
    detail::PinnedTransferBatch transfers(stream, max_transfer, 10U);
    GpuProfiler profiler;

    const std::size_t h2d_bytes =
        d_phi.size_bytes() + d_sx.size_bytes() + d_sy.size_bytes() +
        d_sz.size_bytes() + d_owner.size_bytes() + d_neighbour.size_bytes() +
        d_volume.size_bytes();

    profiler.start(stream.get());
    transfers.h2d(d_phi, phi.data(), d_phi.size_bytes());
    transfers.h2d(d_sx, face_sx.data(), d_sx.size_bytes());
    transfers.h2d(d_sy, face_sy.data(), d_sy.size_bytes());
    transfers.h2d(d_sz, face_sz.data(), d_sz.size_bytes());
    transfers.h2d(d_owner, owner.data(), d_owner.size_bytes());
    transfers.h2d(d_neighbour, neighbour.data(), d_neighbour.size_bytes());
    transfers.h2d(d_volume, volume.data(), d_volume.size_bytes());
    transfers.synchronize_h2d();
    profiler.stop(stream.get());
    if (metrics) metrics->h2d_ms = profiler.elapsed_ms();

    profiler.start(stream.get());
    gradient_gauss_cuda(
        static_cast<const double*>(d_phi.data()),
        static_cast<const double*>(d_sx.data()),
        static_cast<const double*>(d_sy.data()),
        static_cast<const double*>(d_sz.data()),
        static_cast<const std::uint32_t*>(d_owner.data()),
        static_cast<const std::int64_t*>(d_neighbour.data()),
        static_cast<const double*>(d_volume.data()),
        owner.size(), phi.size(),
        static_cast<double*>(d_gx.data()),
        static_cast<double*>(d_gy.data()),
        static_cast<double*>(d_gz.data()),
        stream.get());
    stream.synchronize();
    profiler.stop(stream.get());
    if (metrics) metrics->kernel_ms = profiler.elapsed_ms();

    grad_x.resize(phi.size());
    grad_y.resize(phi.size());
    grad_z.resize(phi.size());

    const std::size_t d2h_bytes =
        d_gx.size_bytes() + d_gy.size_bytes() + d_gz.size_bytes();

    profiler.start(stream.get());
    transfers.d2h(d_gx, grad_x.data(), d_gx.size_bytes());
    transfers.d2h(d_gy, grad_y.data(), d_gy.size_bytes());
    transfers.d2h(d_gz, grad_z.data(), d_gz.size_bytes());
    transfers.synchronize_d2h();
    profiler.stop(stream.get());

    if (metrics) {
        metrics->d2h_ms = profiler.elapsed_ms();
        metrics->h2d_bytes = h2d_bytes;
        metrics->d2h_bytes = d2h_bytes;
    }
}

inline void execute_gradient_cuda(
    const std::vector<double>& phi,
    const std::vector<double>& face_sx,
    const std::vector<double>& face_sy,
    const std::vector<double>& face_sz,
    const std::vector<std::uint32_t>& owner,
    const std::vector<std::int64_t>& neighbour,
    const std::vector<double>& volume,
    std::vector<double>& grad_x,
    std::vector<double>& grad_y,
    std::vector<double>& grad_z)
{
#ifndef CFDX_ENABLE_GPU
    (void)phi; (void)face_sx; (void)face_sy; (void)face_sz;
    (void)owner; (void)neighbour; (void)volume;
    (void)grad_x; (void)grad_y; (void)grad_z;
    throw std::runtime_error("CUDA execution requested but CFDX was built without CUDA");
#else
    execute_gradient_cuda_with_metrics(
        phi, face_sx, face_sy, face_sz, owner, neighbour, volume,
        grad_x, grad_y, grad_z, 0, nullptr);
#endif
}

inline void execute_divergence_cuda(
    const std::vector<double>& phi_face,
    const std::vector<std::uint32_t>& owner,
    const std::vector<std::int64_t>& neighbour,
    std::size_t n_cells,
    std::vector<double>& div)
{
#ifndef CFDX_ENABLE_GPU
    (void)phi_face; (void)owner; (void)neighbour; (void)n_cells; (void)div;
    throw std::runtime_error("CUDA execution requested but CFDX was built without CUDA");
#else
    if (phi_face.empty() || owner.size() != phi_face.size() ||
        neighbour.size() != phi_face.size() || n_cells == 0)
        throw std::invalid_argument("invalid GPU divergence array sizes");

    detail::validate_device(0);

    CudaDeviceBuffer d_flux(phi_face.size() * sizeof(double));
    CudaDeviceBuffer d_owner(owner.size() * sizeof(std::uint32_t));
    CudaDeviceBuffer d_neighbour(neighbour.size() * sizeof(std::int64_t));
    CudaDeviceBuffer d_div(n_cells * sizeof(double));

    const std::size_t max_transfer = detail::max_transfer_bytes({
        d_flux.size_bytes(), d_owner.size_bytes(), d_neighbour.size_bytes(), d_div.size_bytes()});

    CudaStream stream;
    detail::PinnedTransferBatch transfers(stream, max_transfer, 4U);
    transfers.h2d(d_flux, phi_face.data(), d_flux.size_bytes());
    transfers.h2d(d_owner, owner.data(), d_owner.size_bytes());
    transfers.h2d(d_neighbour, neighbour.data(), d_neighbour.size_bytes());
    transfers.synchronize_h2d();

    divergence_cuda(
        static_cast<const double*>(d_flux.data()),
        static_cast<const std::uint32_t*>(d_owner.data()),
        static_cast<const std::int64_t*>(d_neighbour.data()),
        phi_face.size(), n_cells, static_cast<double*>(d_div.data()), stream.get());
    stream.synchronize();

    div.resize(n_cells);
    transfers.d2h(d_div, div.data(), d_div.size_bytes());
    transfers.synchronize_d2h();
#endif
}

} // namespace cfdx::runtime::gpu
