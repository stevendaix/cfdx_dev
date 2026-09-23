#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace cfdx::runtime::gpu {

// Backend-neutral kernels. The CPU reference path is deterministic and is
// also used by tests when CUDA is unavailable.
void gradient_gauss_reference(
    const double* phi,
    const double* face_sx,
    const double* face_sy,
    const double* face_sz,
    const std::size_t* owner,
    const std::int64_t* neighbour,
    const double* cell_volume,
    std::size_t n_faces,
    std::size_t n_cells,
    double* grad_x,
    double* grad_y,
    double* grad_z);

void divergence_reference(
    const double* phi_face,
    const double* face_sx,
    const double* face_sy,
    const double* face_sz,
    const std::size_t* owner,
    const std::int64_t* neighbour,
    std::size_t n_faces,
    std::size_t n_cells,
    double* div);

} // namespace cfdx::runtime::gpu

#ifdef CFDX_ENABLE_GPU
#include <cuda_runtime.h>
extern "C" void cfdx_cuda_gradient_gauss(
    const double* phi, const double* sx, const double* sy, const double* sz,
    const std::uint32_t* owner, const std::int64_t* neighbour,
    const double* volume, std::size_t n_faces, std::size_t n_cells,
    double* gx, double* gy, double* gz, cudaStream_t stream);

extern "C" void cfdx_cuda_divergence(
    const double* phi_face, const std::uint32_t* owner,
    const std::int64_t* neighbour, std::size_t n_faces, std::size_t n_cells,
    double* div, cudaStream_t stream);
#endif

namespace cfdx::runtime::gpu {

#ifdef CFDX_ENABLE_GPU
inline void gradient_gauss_cuda(
    const double* phi,
    const double* face_sx, const double* face_sy, const double* face_sz,
    const std::uint32_t* owner, const std::int64_t* neighbour,
    const double* volume, std::size_t n_faces, std::size_t n_cells,
    double* grad_x, double* grad_y, double* grad_z,
    cudaStream_t stream = nullptr)
{
    cfdx_cuda_gradient_gauss(phi, face_sx, face_sy, face_sz, owner, neighbour,
                              volume, n_faces, n_cells, grad_x, grad_y, grad_z, stream);
}

inline void divergence_cuda(
    const double* phi_face,
    const std::uint32_t* owner, const std::int64_t* neighbour,
    std::size_t n_faces, std::size_t n_cells,
    double* div, cudaStream_t stream = nullptr)
{
    cfdx_cuda_divergence(phi_face, owner, neighbour, n_faces, n_cells, div, stream);
}
#endif

} // namespace cfdx::runtime::gpu
