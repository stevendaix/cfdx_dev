#include "cfdx/runtime/gpu/cuda_poisson.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"

#include <cmath>
#include <cstdio>

using namespace cfdx::core;
using namespace cfdx::runtime::gpu;

static SparseMatrix poisson_matrix()
{
    const std::size_t n = 5;
    SparseMatrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        A.push_back(i, i, 2.0);
        if (i > 0) A.push_back(i, i - 1, -1.0);
        if (i + 1 < n) A.push_back(i, i + 1, -1.0);
    }
    A.finalize();
    return A;
}

int main()
{
    if (!cuda_poisson_available()) {
        std::printf("CUDA Poisson: no CUDA device, validation skipped\n");
        return 0;
    }

    const auto A = poisson_matrix();
    Vector exact(5);
    exact(0) = 0.5; exact(1) = 1.0; exact(2) = 1.5; exact(3) = 1.0; exact(4) = 0.5;
    const auto b_values = A.matvec(exact);
    Vector b(5);
    for (std::size_t i = 0; i < 5; ++i) b(i) = b_values[i];

    Vector x(5, 0.0);
    const auto result = solve_poisson_cuda(A, b, x, 100, 1e-11);
    if (result.status != SolverStatus::CONVERGED ||
        result.residual_relative >= 1e-10) {
        std::fprintf(stderr, "CUDA Poisson did not converge: status=%s rel=%.17g\n",
                     to_string(result.status), result.residual_relative);
        return 1;
    }
    for (std::size_t i = 0; i < 5; ++i)
        if (std::abs(x(i) - exact(i)) > 1e-9) {
            std::fprintf(stderr, "CUDA Poisson solution mismatch at %zu\n", i);
            return 1;
        }

    std::printf("CUDA Poisson: device-resident CG validation PASS\n");
    return 0;
}
