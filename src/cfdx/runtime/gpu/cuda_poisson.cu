#include "cfdx/runtime/gpu/cuda_poisson.h"

#include <cuda_runtime.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::runtime::gpu {
namespace {

inline void cuda_check(cudaError_t status, const char* what)
{
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(status));
}

class DeviceBuffer {
public:
    DeviceBuffer() = default;
    explicit DeviceBuffer(std::size_t bytes) { resize(bytes); }
    ~DeviceBuffer() { if (ptr_) cudaFree(ptr_); }
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    void resize(std::size_t bytes) {
        if (ptr_) cuda_check(cudaFree(ptr_), "cudaFree");
        ptr_ = nullptr;
        bytes_ = bytes;
        if (bytes_ != 0) cuda_check(cudaMalloc(&ptr_, bytes_), "cudaMalloc");
    }
    void* data() noexcept { return ptr_; }
    const void* data() const noexcept { return ptr_; }
    std::size_t bytes() const noexcept { return bytes_; }

private:
    void* ptr_ = nullptr;
    std::size_t bytes_ = 0;
};

__global__ void csr_spmv(
    std::size_t n,
    const double* values,
    const std::uint32_t* columns,
    const std::uint32_t* row_offsets,
    const double* x,
    double* y)
{
    const std::size_t row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= n) return;
    double sum = 0.0;
    for (std::uint32_t k = row_offsets[row]; k < row_offsets[row + 1]; ++k)
        sum += values[k] * x[columns[k]];
    y[row] = sum;
}

__global__ void axpy_residual(
    std::size_t n, double alpha,
    double* x, const double* p,
    double* r, const double* Ap)
{
    const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    x[i] += alpha * p[i];
    r[i] -= alpha * Ap[i];
}

__global__ void jacobi_apply(
    std::size_t n, const double* r, const double* diagonal, double* z)
{
    const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    z[i] = r[i] / diagonal[i];
}

__global__ void direction_update(
    std::size_t n, double beta, const double* z, double* p)
{
    const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    p[i] = z[i] + beta * p[i];
}

__global__ void dot_kernel(
    std::size_t n, const double* a, const double* b, double* result)
{
    extern __shared__ double shared[];
    const std::size_t tid = threadIdx.x;
    const std::size_t i = blockIdx.x * blockDim.x + tid;
    double value = (i < n) ? a[i] * b[i] : 0.0;
    shared[tid] = value;
    __syncthreads();
    for (unsigned int stride = blockDim.x / 2; stride != 0; stride >>= 1) {
        if (tid < stride) shared[tid] += shared[tid + stride];
        __syncthreads();
    }
    if (tid == 0) atomicAdd(result, shared[0]);
}

double device_dot(
    std::size_t n, const double* a, const double* b,
    double* d_result, unsigned int block_size)
{
    cuda_check(cudaMemset(d_result, 0, sizeof(double)), "cudaMemset(dot)");
    const unsigned int blocks = static_cast<unsigned int>(
        (n + block_size - 1) / block_size);
    dot_kernel<<<blocks, block_size, block_size * sizeof(double)>>>(
        n, a, b, d_result);
    cuda_check(cudaGetLastError(), "dot_kernel launch");
    double host = 0.0;
    cuda_check(cudaMemcpy(&host, d_result, sizeof(double),
                          cudaMemcpyDeviceToHost), "cudaMemcpy(dot)");
    return host;
}

void launch_check(const char* name)
{
    cuda_check(cudaGetLastError(), name);
}

} // namespace

bool cuda_poisson_available()
{
    int count = 0;
    const cudaError_t status = cudaGetDeviceCount(&count);
    if (status != cudaSuccess)
        return false;
    return count > 0;
}

cfdx::core::SolverResult solve_poisson_cuda(
    const cfdx::core::SparseMatrix& matrix,
    const cfdx::core::Vector& rhs,
    cfdx::core::Vector& solution,
    std::size_t max_iterations,
    double tolerance)
{
    using cfdx::core::SolverResult;
    using cfdx::core::SolverStatus;

    SolverResult result;
    if (matrix.n_rows() != matrix.n_cols() ||
        rhs.size() != matrix.n_rows() ||
        solution.size() != matrix.n_cols()) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }
    if (!(tolerance >= 0.0) || !std::isfinite(tolerance)) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }
    if (!cuda_poisson_available()) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }

    const std::size_t n = matrix.n_rows();
    if (n == 0) {
        result.status = SolverStatus::CONVERGED;
        return result;
    }
    if (matrix.nnz() > std::numeric_limits<std::uint32_t>::max() ||
        n > std::numeric_limits<std::uint32_t>::max())
        throw std::overflow_error("solve_poisson_cuda: CSR dimensions exceed CUDA index type");

    const auto* values = matrix.values_data();
    const auto* columns = matrix.columns_data();
    const auto* row_offsets = matrix.row_offsets_data();

    std::vector<double> diagonal(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::uint32_t k = row_offsets[i]; k < row_offsets[i + 1]; ++k) {
            if (columns[k] == static_cast<std::uint32_t>(i)) {
                diagonal[i] = values[k];
                break;
            }
        }
        if (!(diagonal[i] > 0.0) || !std::isfinite(diagonal[i]))
            throw std::invalid_argument("solve_poisson_cuda: matrix requires positive diagonal");
        if (!std::isfinite(rhs(i)) || !std::isfinite(solution(i)))
            throw std::invalid_argument("solve_poisson_cuda: non-finite input");
    }

    DeviceBuffer d_values(matrix.nnz() * sizeof(double));
    DeviceBuffer d_columns(matrix.nnz() * sizeof(std::uint32_t));
    DeviceBuffer d_offsets((n + 1) * sizeof(std::uint32_t));
    DeviceBuffer d_diag(n * sizeof(double));
    DeviceBuffer d_b(n * sizeof(double));
    DeviceBuffer d_x(n * sizeof(double));
    DeviceBuffer d_r(n * sizeof(double));
    DeviceBuffer d_z(n * sizeof(double));
    DeviceBuffer d_p(n * sizeof(double));
    DeviceBuffer d_Ap(n * sizeof(double));
    DeviceBuffer d_dot(sizeof(double));

    cuda_check(cudaMemcpy(d_values.data(), values, d_values.bytes(),
                          cudaMemcpyHostToDevice), "cudaMemcpy(values)");
    cuda_check(cudaMemcpy(d_columns.data(), columns, d_columns.bytes(),
                          cudaMemcpyHostToDevice), "cudaMemcpy(columns)");
    cuda_check(cudaMemcpy(d_offsets.data(), row_offsets, d_offsets.bytes(),
                          cudaMemcpyHostToDevice), "cudaMemcpy(row_offsets)");
    cuda_check(cudaMemcpy(d_diag.data(), diagonal.data(), d_diag.bytes(),
                          cudaMemcpyHostToDevice), "cudaMemcpy(diagonal)");

    std::vector<double> host_rhs(n);
    std::vector<double> host_x(n);
    for (std::size_t i = 0; i < n; ++i) {
        host_rhs[i] = rhs(i);
        host_x[i] = solution(i);
    }
    cuda_check(cudaMemcpy(d_b.data(), host_rhs.data(), d_b.bytes(),
                          cudaMemcpyHostToDevice), "cudaMemcpy(rhs)");
    cuda_check(cudaMemcpy(d_x.data(), host_x.data(), d_x.bytes(),
                          cudaMemcpyHostToDevice), "cudaMemcpy(solution)");

    constexpr unsigned int block = 256;
    const unsigned int blocks =
        static_cast<unsigned int>((n + block - 1) / block);

    csr_spmv<<<blocks, block>>>(n,
        static_cast<const double*>(d_values.data()),
        static_cast<const std::uint32_t*>(d_columns.data()),
        static_cast<const std::uint32_t*>(d_offsets.data()),
        static_cast<const double*>(d_x.data()),
        static_cast<double*>(d_Ap.data()));
    launch_check("csr_spmv initial");

    // r = b - A*x, formed by a small host-side kernel to keep the entire
    // Krylov vector state resident on the GPU.
    cuda_check(cudaMemcpy(d_r.data(), d_b.data(), d_b.bytes(),
                          cudaMemcpyDeviceToDevice), "cudaMemcpy residual init");
    axpy_residual<<<blocks, block>>>(n, -1.0,
        static_cast<double*>(d_r.data()),
        static_cast<const double*>(d_Ap.data()),
        static_cast<double*>(d_r.data()),
        static_cast<const double*>(d_Ap.data()));
    launch_check("initial residual kernel");

    jacobi_apply<<<blocks, block>>>(n,
        static_cast<const double*>(d_r.data()),
        static_cast<const double*>(d_diag.data()),
        static_cast<double*>(d_z.data()));
    direction_update<<<blocks, block>>>(n, 0.0,
        static_cast<const double*>(d_z.data()),
        static_cast<double*>(d_p.data()));
    launch_check("initial preconditioner");

    double b2 = 0.0;
    for (double v : host_rhs) b2 += v * v;
    const double b_norm = std::sqrt(b2);
    const double tol_abs = tolerance * std::max(b_norm, 1e-15);
    double rz = device_dot(n,
        static_cast<const double*>(d_r.data()),
        static_cast<const double*>(d_z.data()),
        static_cast<double*>(d_dot.data()), block);

    result.residual = std::sqrt(std::max(0.0, rz));
    result.residual_relative = b_norm > 0.0 ? result.residual / b_norm : 0.0;
    if (result.residual <= tol_abs) {
        result.status = SolverStatus::CONVERGED;
        cuda_check(cudaMemcpy(host_x.data(), d_x.data(), d_x.bytes(),
                              cudaMemcpyDeviceToHost), "cudaMemcpy(solution final)");
        for (std::size_t i = 0; i < n; ++i) solution(i) = host_x[i];
        return result;
    }

    for (std::size_t iter = 1; iter <= max_iterations; ++iter) {
        csr_spmv<<<blocks, block>>>(n,
            static_cast<const double*>(d_values.data()),
            static_cast<const std::uint32_t*>(d_columns.data()),
            static_cast<const std::uint32_t*>(d_offsets.data()),
            static_cast<const double*>(d_p.data()),
            static_cast<double*>(d_Ap.data()));
        launch_check("csr_spmv iteration");

        const double pAp = device_dot(n,
            static_cast<const double*>(d_p.data()),
            static_cast<const double*>(d_Ap.data()),
            static_cast<double*>(d_dot.data()), block);
        if (!(pAp > 0.0) || !std::isfinite(pAp)) {
            result.status = SolverStatus::NOT_APPLICABLE;
            result.iterations = iter;
            break;
        }

        const double alpha = rz / pAp;
        if (!std::isfinite(alpha)) {
            result.status = SolverStatus::DIVERGED;
            result.iterations = iter;
            break;
        }

        axpy_residual<<<blocks, block>>>(n, alpha,
            static_cast<double*>(d_x.data()),
            static_cast<const double*>(d_p.data()),
            static_cast<double*>(d_r.data()),
            static_cast<const double*>(d_Ap.data()));
        launch_check("CG update");

        const double rr = device_dot(n,
            static_cast<const double*>(d_r.data()),
            static_cast<const double*>(d_r.data()),
            static_cast<double*>(d_dot.data()), block);
        result.residual = std::sqrt(std::max(0.0, rr));
        result.residual_relative = b_norm > 0.0 ? result.residual / b_norm : 0.0;
        result.iterations = iter;

        if (!std::isfinite(result.residual)) {
            result.status = SolverStatus::DIVERGED;
            break;
        }
        if (result.residual <= tol_abs) {
            result.status = SolverStatus::CONVERGED;
            break;
        }

        jacobi_apply<<<blocks, block>>>(n,
            static_cast<const double*>(d_r.data()),
            static_cast<const double*>(d_diag.data()),
            static_cast<double*>(d_z.data()));
        launch_check("Jacobi update");

        const double rz_new = device_dot(n,
            static_cast<const double*>(d_r.data()),
            static_cast<const double*>(d_z.data()),
            static_cast<double*>(d_dot.data()), block);
        if (!(rz != 0.0) || !std::isfinite(rz_new)) {
            result.status = SolverStatus::DIVERGED;
            break;
        }
        const double beta = rz_new / rz;
        direction_update<<<blocks, block>>>(n, beta,
            static_cast<const double*>(d_z.data()),
            static_cast<double*>(d_p.data()));
        launch_check("CG direction update");
        rz = rz_new;
    }

    if (result.status == SolverStatus::NOT_APPLICABLE &&
        result.iterations == max_iterations)
        result.status = SolverStatus::MAX_ITER_REACHED;

    cuda_check(cudaMemcpy(host_x.data(), d_x.data(), d_x.bytes(),
                          cudaMemcpyDeviceToHost), "cudaMemcpy(solution final)");
    for (std::size_t i = 0; i < n; ++i) solution(i) = host_x[i];

    return result;
}

} // namespace cfdx::runtime::gpu
