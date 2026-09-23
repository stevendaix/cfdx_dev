#pragma once

#include "cfdx/runtime/gpu/gpu_kernels.h"
#include "cfdx/runtime/gpu/cuda_backend.h"
#include "cfdx/runtime/gpu/gpu_profiler.h"
#include "cfdx/runtime/ooc/cuda_pinned_buffer_pool.h"

#include <cstddef>
#include <cstring>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace cfdx::runtime::gpu {
#ifdef CFDX_ENABLE_GPU
namespace detail {
inline std::size_t checked_mul(std::size_t a, std::size_t b) {
    if (b != 0 && a > std::numeric_limits<std::size_t>::max() / b)
        throw std::overflow_error("GPU staging allocation size overflow");
    return a * b;
}
class PinnedTransferBatch {
public:
    PinnedTransferBatch(CudaStream& stream, std::size_t bytes, std::size_t count)
        : stream_(stream), pool_(checked_mul(bytes, count), bytes) {
        if (bytes == 0 || count == 0)
            throw std::invalid_argument("PinnedTransferBatch requires non-zero capacity");
    }
    ~PinnedTransferBatch() { release_noexcept(); }
    PinnedTransferBatch(const PinnedTransferBatch&) = delete;
    PinnedTransferBatch& operator=(const PinnedTransferBatch&) = delete;
    void h2d(CudaDeviceBuffer& dst, const void* src, std::size_t bytes) {
        if (bytes > pool_.buffer_bytes()) throw std::out_of_range("pinned H2D transfer exceeds staging buffer");
        auto* b=pool_.acquire(); if(!b) throw std::runtime_error("pinned H2D staging pool exhausted");
        try { std::memcpy(b->data,src,bytes); async_copy_h2d(dst,b->data,bytes,stream_.get()); h2d_.push_back(b); }
        catch (...) { pool_.release(b); throw; }
    }
    void synchronize_h2d() { stream_.synchronize(); for(auto* b:h2d_) pool_.release(b); h2d_.clear(); }
    void d2h(const CudaDeviceBuffer& src, void* dst, std::size_t bytes) {
        if (bytes > pool_.buffer_bytes()) throw std::out_of_range("pinned D2H transfer exceeds staging buffer");
        auto* b=pool_.acquire(); if(!b) throw std::runtime_error("pinned D2H staging pool exhausted");
        try { async_copy_d2h(src,b->data,bytes,stream_.get()); d2h_.push_back({b,dst,bytes}); }
        catch (...) { pool_.release(b); throw; }
    }
    void synchronize_d2h() {
        stream_.synchronize();
        for(const auto& t:d2h_) { std::memcpy(t.dst,t.b->data,t.bytes); pool_.release(t.b); }
        d2h_.clear();
    }
private:
    struct T { ooc::CudaPinnedBufferPool::Buffer* b; void* dst; std::size_t bytes; };
    void release_noexcept() noexcept {
        try { stream_.synchronize(); for(auto* b:h2d_) pool_.release(b); for(auto& t:d2h_) pool_.release(t.b); } catch(...) {}
        h2d_.clear(); d2h_.clear();
    }
    CudaStream& stream_;
    ooc::CudaPinnedBufferPool pool_;
    std::vector<ooc::CudaPinnedBufferPool::Buffer*> h2d_;
    std::vector<T> d2h_;
};
} // namespace detail
struct GpuExecutionMetrics {
    double h2d_ms=0.0, kernel_ms=0.0, d2h_ms=0.0;
    std::size_t h2d_bytes=0, d2h_bytes=0;
};
#endif


// Executes the finite-volume gradient entirely on CUDA memory after the
// host-to-device staging copies. No CPU kernel is used on the GPU path.
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
    if (owner.size() != neighbour.size() ||
        owner.size() != face_sx.size() ||
        owner.size() != face_sy.size() ||
        owner.size() != face_sz.size() ||
        phi.empty() || volume.empty() || owner.empty() ||
        volume.size() != phi.size()) {
        throw std::invalid_argument("invalid GPU gradient array sizes");
    }
    const std::size_t n_cells = phi.size();
    for (std::size_t f = 0; f < owner.size(); ++f) {
        if (static_cast<std::size_t>(owner[f]) >= n_cells ||
            neighbour[f] < -1 ||
            (neighbour[f] >= 0 && static_cast<std::size_t>(neighbour[f]) >= n_cells)) {
            throw std::invalid_argument("invalid GPU gradient face ownership");
        }
    }
    int device_count = 0;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess || device_count <= 0)
        throw std::runtime_error("CUDA execution requested but no CUDA device is available");

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

    CudaStream stream;
    const std::size_t max_transfer = std::max({d_phi.size_bytes(),d_sx.size_bytes(),d_sy.size_bytes(),
        d_sz.size_bytes(),d_owner.size_bytes(),d_neighbour.size_bytes(),d_volume.size_bytes(),
        d_gx.size_bytes(),d_gy.size_bytes(),d_gz.size_bytes()});
    detail::PinnedTransferBatch transfers(stream,max_transfer,10U);
    transfers.h2d(d_phi,phi.data(),d_phi.size_bytes());
    transfers.h2d(d_sx,face_sx.data(),d_sx.size_bytes());
    transfers.h2d(d_sy,face_sy.data(),d_sy.size_bytes());
    transfers.h2d(d_sz,face_sz.data(),d_sz.size_bytes());
    transfers.h2d(d_owner,owner.data(),d_owner.size_bytes());
    transfers.h2d(d_neighbour,neighbour.data(),d_neighbour.size_bytes());
    transfers.h2d(d_volume,volume.data(),d_volume.size_bytes());
    transfers.synchronize_h2d();
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

    grad_x.resize(phi.size());
    grad_y.resize(phi.size());
    grad_z.resize(phi.size());
    transfers.d2h(d_gx,grad_x.data(),d_gx.size_bytes());
    transfers.d2h(d_gy,grad_y.data(),d_gy.size_bytes());
    transfers.d2h(d_gz,grad_z.data(),d_gz.size_bytes());
    transfers.synchronize_d2h();
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
    for (std::size_t f = 0; f < owner.size(); ++f) {
        if (static_cast<std::size_t>(owner[f]) >= n_cells ||
            neighbour[f] < -1 ||
            (neighbour[f] >= 0 && static_cast<std::size_t>(neighbour[f]) >= n_cells)) {
            throw std::invalid_argument("invalid GPU divergence face ownership");
        }
    }
    int device_count = 0;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess || device_count <= 0)
        throw std::runtime_error("CUDA execution requested but no CUDA device is available");

    CudaDeviceBuffer d_flux(phi_face.size() * sizeof(double));
    CudaDeviceBuffer d_owner(owner.size() * sizeof(std::uint32_t));
    CudaDeviceBuffer d_neighbour(neighbour.size() * sizeof(std::int64_t));
    CudaDeviceBuffer d_div(n_cells * sizeof(double));
    d_flux.copy_from_host(phi_face.data(), d_flux.size_bytes());
    d_owner.copy_from_host(owner.data(), d_owner.size_bytes());
    d_neighbour.copy_from_host(neighbour.data(), d_neighbour.size_bytes());

    CudaStream stream;
    divergence_cuda(
        static_cast<const double*>(d_flux.data()),
        static_cast<const std::uint32_t*>(d_owner.data()),
        static_cast<const std::int64_t*>(d_neighbour.data()),
        phi_face.size(), n_cells, static_cast<double*>(d_div.data()), stream.get());
    stream.synchronize();

    div.resize(n_cells);
    d_div.copy_to_host(div.data(), d_div.size_bytes());
#endif
}

} // namespace cfdx::runtime::gpu
