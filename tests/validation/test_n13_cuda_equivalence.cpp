#include "cfdx/runtime/gpu/gpu_execution.h"
#include "cfdx/runtime/gpu/gpu_kernels.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace {

struct ErrorMetrics {
    double l2_abs = 0.0;
    double l2_rel = 0.0;
    double linf_abs = 0.0;
    double linf_rel = 0.0;
};

ErrorMetrics compare_field(const std::vector<double>& reference,
                           const std::vector<double>& candidate)
{
    if (reference.size() != candidate.size())
        throw std::runtime_error("CPU/CUDA result size mismatch");

    long double diff2 = 0.0L;
    long double ref2 = 0.0L;
    double linf_abs = 0.0;
    double linf_rel = 0.0;
    for (std::size_t i = 0; i < reference.size(); ++i) {
        const double d = std::abs(candidate[i] - reference[i]);
        diff2 += static_cast<long double>(d) * d;
        ref2 += static_cast<long double>(reference[i]) * reference[i];
        linf_abs = std::max(linf_abs, d);
        const double scale = std::max(std::abs(reference[i]), 1.0e-30);
        linf_rel = std::max(linf_rel, d / scale);
    }
    const double l2_abs = std::sqrt(static_cast<double>(diff2));
    const double l2_ref = std::sqrt(static_cast<double>(ref2));
    return {l2_abs, l2_abs / std::max(l2_ref, 1.0e-30), linf_abs, linf_rel};
}

bool check(const char* name, const ErrorMetrics& m, double tol)
{
    const bool pass = m.l2_rel <= tol && m.linf_rel <= tol;
    std::printf("N13 CPU/CUDA %s %s L2_abs=%.17g L2_rel=%.17g Linf_abs=%.17g Linf_rel=%.17g tolerance=%.17g\n",
                name, pass ? "PASS" : "FAIL", m.l2_abs, m.l2_rel,
                m.linf_abs, m.linf_rel, tol);
    return pass;
}

} // namespace

int main()
{
#ifndef CFDX_ENABLE_GPU
    std::fprintf(stderr, "N13 CUDA equivalence is BLOCKED: CFDX was built without CUDA\n");
    return 2;
#else
    try {
        // Deliberately non-trivial field and mixed-sign geometry. This is compared
        // against the production CPU reference implementation, not against an
        // expected zero field.
        const std::vector<double> phi{0.25, 1.5, -0.75, 2.25};
        const std::vector<double> sx{1.0, 0.7, -1.2, -0.4, 0.3, -0.8};
        const std::vector<double> sy{0.2, -0.6, 0.9, -0.3, 0.5, -0.7};
        const std::vector<double> sz{-0.4, 0.8, 0.1, 0.6, -0.9, 0.2};
        const std::vector<std::uint32_t> owner_u32{0, 0, 1, 1, 2, 3};
        const std::vector<std::size_t> owner{0, 0, 1, 1, 2, 3};
        const std::vector<std::int64_t> neighbour{1, -1, 2, -1, 3, -1};
        const std::vector<double> volume{1.0, 1.3, 0.8, 1.7};

        std::vector<double> cpu_x(phi.size()), cpu_y(phi.size()), cpu_z(phi.size());
        cfdx::runtime::gpu::gradient_gauss_reference(
            phi.data(), sx.data(), sy.data(), sz.data(),
            owner.data(), neighbour.data(),
            volume.data(), owner.size(), phi.size(),
            cpu_x.data(), cpu_y.data(), cpu_z.data());

        std::vector<double> cuda_x, cuda_y, cuda_z;
        GpuExecutionMetrics metrics;
        execute_gradient_cuda_with_metrics(
            phi, sx, sy, sz, owner_u32, neighbour, volume,
            cuda_x, cuda_y, cuda_z, 0, &metrics);

        const auto mx = compare_field(cpu_x, cuda_x);
        const auto my = compare_field(cpu_y, cuda_y);
        const auto mz = compare_field(cpu_z, cuda_z);
        constexpr double tolerance = 1.0e-12;

        const bool pass = check("gradient_x", mx, tolerance) &&
                          check("gradient_y", my, tolerance) &&
                          check("gradient_z", mz, tolerance);

        std::printf("N13 CUDA diagnostics h2d_ms=%.17g kernel_ms=%.17g d2h_ms=%.17g h2d_bytes=%zu d2h_bytes=%zu\n",
                    metrics.h2d_ms, metrics.kernel_ms, metrics.d2h_ms,
                    metrics.h2d_bytes, metrics.d2h_bytes);
        return pass ? 0 : 1;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "N13 CUDA equivalence FAIL exception: %s\n", e.what());
        return 1;
    }
#endif
}
