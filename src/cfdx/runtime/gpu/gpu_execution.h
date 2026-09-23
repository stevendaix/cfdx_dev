#pragma once

#include "cfdx/runtime/gpu/gpu_kernels.h"
#ifdef CFDX_ENABLE_GPU
#include "cfdx/runtime/gpu/cuda_backend.h"
#endif

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace cfdx::runtime::gpu {

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
    async_copy_h2d(d_phi, phi.data(), d_phi.size_bytes(), stream.get());
    async_copy_h2d(d_sx, face_sx.data(), d_sx.size_bytes(), stream.get());
    async_copy_h2d(d_sy, face_sy.data(), d_sy.size_bytes(), stream.get());
    async_copy_h2d(d_sz, face_sz.data(), d_sz.size_bytes(), stream.get());
    async_copy_h2d(d_owner, owner.data(), d_owner.size_bytes(), stream.get());
    async_copy_h2d(d_neighbour, neighbour.data(), d_neighbour.size_bytes(), stream.get());
    async_copy_h2d(d_volume, volume.data(), d_volume.size_bytes(), stream.get());
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
    async_copy_d2h(d_gx, grad_x.data(), d_gx.size_bytes(), stream.get());
    async_copy_d2h(d_gy, grad_y.data(), d_gy.size_bytes(), stream.get());
    async_copy_d2h(d_gz, grad_z.data(), d_gz.size_bytes(), stream.get());
    stream.synchronize();
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
    int device_count = 0;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess || device_count <= 0)
        throw std::runtime_error("CUDA execution requested but no CUDA device is available");

    CudaDeviceBuffer d_flux(phi_face.size() * sizeof(double));
    CudaDeviceBuffer d_owner(owner.size() * sizeof(std::uint32_t));
    CudaDeviceBuffer d_neighbour(neighbour.size() * sizeof(std::int64_t));
    CudaDeviceBuffer d_div(n_cells * sizeof(double));
    CudaStream stream;
    async_copy_h2d(d_flux, phi_face.data(), d_flux.size_bytes(), stream.get());
    async_copy_h2d(d_owner, owner.data(), d_owner.size_bytes(), stream.get());
    async_copy_h2d(d_neighbour, neighbour.data(), d_neighbour.size_bytes(), stream.get());
    divergence_cuda(
        static_cast<const double*>(d_flux.data()),
        static_cast<const std::uint32_t*>(d_owner.data()),
        static_cast<const std::int64_t*>(d_neighbour.data()),
        phi_face.size(), n_cells, static_cast<double*>(d_div.data()), stream.get());
    stream.synchronize();

    div.resize(n_cells);
    async_copy_d2h(d_div, div.data(), d_div.size_bytes(), stream.get());
    stream.synchronize();
#endif
}

} // namespace cfdx::runtime::gpu
