#include "cfdx/runtime/gpu/gpu_execution.h"

#include <cmath>
#include <cstdio>
#include <vector>

int main()
{
#ifndef CFDX_ENABLE_GPU
    std::fprintf(stderr, "N13 CUDA equivalence is BLOCKED: CFDX was built without CUDA\n");
    return 2;
#else
    const std::vector<double> phi{0.0, 1.0};
    const std::vector<double> sx{1.0, 1.0, -1.0, -1.0};
    const std::vector<double> sy{0.0, 0.0, 0.0, 0.0};
    const std::vector<double> sz{0.0, 0.0, 0.0, 0.0};
    const std::vector<std::uint32_t> owner{0, 0, 1, 1};
    const std::vector<std::int64_t> neighbour{1, -1, 0, -1};
    const std::vector<double> volume{1.0, 1.0};

    std::vector<double> gx, gy, gz;
    cfdx::runtime::gpu::execute_gradient_cuda(
        phi, sx, sy, sz, owner, neighbour, volume, gx, gy, gz);

    const std::vector<double> cpu_gx{0.5, -0.5};
    if (gx.size() != cpu_gx.size())
        return 1;

    double max_error = 0.0;
    for (std::size_t i = 0; i < gx.size(); ++i)
        max_error = std::max(max_error, std::abs(gx[i] - cpu_gx[i]));

    if (max_error > 1e-12) {
        std::fprintf(stderr, "N13 CPU/CUDA equivalence failed: max error %.17g\n", max_error);
        return 1;
    }

    std::printf("N13 CPU/CUDA equivalence PASS max_error=%.17g\n", max_error);
    return 0;
#endif
}
