#include <cuda_runtime.h>
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <string>

__global__ void cfdx_normalize_gradient_kernel(
    double* gx, double* gy, double* gz, const double* volume, std::size_t n_cells);

__global__ void cfdx_gradient_gauss_kernel(
    const double* phi,
    const double* sx, const double* sy, const double* sz,
    const std::uint32_t* owner, const std::int64_t* neighbour,
    const double* volume, std::size_t n_faces,
    double* gx, double* gy, double* gz)
{
    const std::size_t f = blockIdx.x * blockDim.x + threadIdx.x;
    if (f >= n_faces) return;
    const std::size_t o = owner[f];
    const std::int64_t n = neighbour[f];
    const double vf = (n >= 0) ? 0.5 * (phi[o] + phi[static_cast<std::size_t>(n)]) : phi[o];
    atomicAdd(&gx[o], vf * sx[f]);
    atomicAdd(&gy[o], vf * sy[f]);
    atomicAdd(&gz[o], vf * sz[f]);
    if (n >= 0) {
        const std::size_t j = static_cast<std::size_t>(n);
        atomicAdd(&gx[j], -vf * sx[f]);
        atomicAdd(&gy[j], -vf * sy[f]);
        atomicAdd(&gz[j], -vf * sz[f]);
    }
}

extern "C" void cfdx_cuda_gradient_gauss(
    const double* phi, const double* sx, const double* sy, const double* sz,
    const std::uint32_t* owner, const std::int64_t* neighbour,
    const double* volume, std::size_t n_faces, std::size_t n_cells,
    double* gx, double* gy, double* gz, cudaStream_t stream)
{
    const int block = 256;
    const int grid = static_cast<int>((n_faces + block - 1) / block);
    cudaMemsetAsync(gx, 0, n_cells * sizeof(double), stream);
    cudaMemsetAsync(gy, 0, n_cells * sizeof(double), stream);
    cudaMemsetAsync(gz, 0, n_cells * sizeof(double), stream);
    cfdx_gradient_gauss_kernel<<<grid, block, 0, stream>>>(
        phi, sx, sy, sz, owner, neighbour, volume, n_faces, gx, gy, gz);
    // Normalize in a separate kernel to keep the accumulation kernel simple.
    const int grid_cells = static_cast<int>((n_cells + block - 1) / block);
    cfdx_normalize_gradient_kernel<<<grid_cells, block, 0, stream>>>(
        gx, gy, gz, volume, n_cells);
    const cudaError_t error = cudaPeekAtLastError();
    if (error != cudaSuccess)
        throw std::runtime_error(
            std::string("CFDX CUDA gradient kernel launch failed: ") +
            cudaGetErrorString(error));
}

__global__ void cfdx_normalize_gradient_kernel(double* gx, double* gy, double* gz, const double* volume, std::size_t n_cells) {\n    const std::size_t c = blockIdx.x * blockDim.x + threadIdx.x;\n    if (c >= n_cells) return;\n    const double inv = volume[c] > 0.0 ? 1.0 / volume[c] : 0.0;\n    gx[c] *= inv; gy[c] *= inv; gz[c] *= inv;\n}\n\n__global__ void cfdx_divergence_kernel(
    const double* phi_face, const std::uint32_t* owner,
    const std::int64_t* neighbour, std::size_t n_faces, double* div)
{
    const std::size_t f = blockIdx.x * blockDim.x + threadIdx.x;
    if (f >= n_faces) return;
    const std::size_t o = owner[f];
    atomicAdd(&div[o], phi_face[f]);
    const std::int64_t n = neighbour[f];
    if (n >= 0) atomicAdd(&div[static_cast<std::size_t>(n)], -phi_face[f]);
}

extern "C" void cfdx_cuda_divergence(
    const double* phi_face, const std::uint32_t* owner,
    const std::int64_t* neighbour, std::size_t n_faces, std::size_t n_cells,
    double* div, cudaStream_t stream)
{
    const int block = 256;
    const int grid = static_cast<int>((n_faces + block - 1) / block);
    cudaMemsetAsync(div, 0, n_cells * sizeof(double), stream);
    cfdx_divergence_kernel<<<grid, block, 0, stream>>>(
        phi_face, owner, neighbour, n_faces, div);
    const cudaError_t error = cudaPeekAtLastError();
    if (error != cudaSuccess)
        throw std::runtime_error(
            std::string("CFDX CUDA divergence kernel launch failed: ") +
            cudaGetErrorString(error));
}
