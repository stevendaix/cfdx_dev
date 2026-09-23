#include "cfdx/runtime/execution/memory_planner.h"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <stdexcept>
#include <vector>

#ifdef CFDX_ENABLE_GPU
#include <cuda_runtime.h>
#endif

namespace {

struct Result {
    std::size_t cells;
    std::size_t estimated_bytes;
    std::size_t peak_bytes;
    double bytes_per_cell;
    double plan_us;
};

Result benchmark(std::size_t cells, int iterations) {
    using clock = std::chrono::steady_clock;
    cfdx::runtime::MemoryPlanner planner;
    const std::size_t field = cells * sizeof(double);
    const std::size_t geometry = cells * 3U * sizeof(double);
    const std::size_t topology = cells * 8U * sizeof(std::uint32_t);
    const std::size_t vectors = cells * 5U * sizeof(double);
    const std::size_t temporary = cells * 2U * sizeof(double);

    planner.add(cfdx::runtime::MemoryKind::Fields, field, "fields");
    planner.add(cfdx::runtime::MemoryKind::Geometry, geometry, "geometry");
    planner.add(cfdx::runtime::MemoryKind::Mesh, topology, "topology");
    planner.add(cfdx::runtime::MemoryKind::LinearAlgebra, vectors, "solver_vectors");
    planner.add(cfdx::runtime::MemoryKind::Temporary, temporary, "temporary");

    std::size_t estimated = 0;
    std::size_t peak = 0;
    const auto begin = clock::now();
    for (int i = 0; i < iterations; ++i) {
        const auto plan = planner.plan(std::numeric_limits<std::size_t>::max());
        if (plan.overflow || !plan.fits())
            throw std::runtime_error("unexpected memory-planner failure");
        estimated = plan.estimated_bytes;
        peak = plan.peak_bytes;
    }
    const auto end = clock::now();

    const double elapsed_us =
        std::chrono::duration<double, std::micro>(end - begin).count();

    return {cells, estimated, peak,
            static_cast<double>(estimated) / static_cast<double>(cells),
            elapsed_us / static_cast<double>(iterations)};
}

} // namespace

int main() {
    try {
        constexpr int iterations = 200;
        const std::vector<std::size_t> sizes{1000, 10000, 100000};

        std::cout << "CFDX Phase 6.12 memory benchmark\n";
        std::cout << "iterations=" << iterations << "\n";
        std::cout << "cells,estimated_bytes,peak_bytes,bytes_per_cell,planner_us_per_call\n";

        std::size_t previous = 0;
        for (const auto cells : sizes) {
            const auto result = benchmark(cells, iterations);
            if (result.estimated_bytes <= previous)
                throw std::runtime_error("memory estimate is not monotonic");
            previous = result.estimated_bytes;

            std::cout << result.cells << ','
                      << result.estimated_bytes << ','
                      << result.peak_bytes << ','
                      << std::fixed << std::setprecision(3)
                      << result.bytes_per_cell << ','
                      << result.plan_us << '\n';
        }

#ifdef CFDX_ENABLE_GPU
        int devices = 0;
        const auto status = cudaGetDeviceCount(&devices);
        if (status == cudaSuccess && devices > 0) {
            std::size_t free_bytes = 0;
            std::size_t total_bytes = 0;
            if (cudaMemGetInfo(&free_bytes, &total_bytes) == cudaSuccess) {
                std::cout << "gpu_devices=" << devices
                          << ",gpu_free_bytes=" << free_bytes
                          << ",gpu_total_bytes=" << total_bytes << '\n';
            }
        }
#else
        std::cout << "gpu_devices=0 (CUDA disabled)\n";
#endif

        return 0;
    } catch (const std::exception& error) {
        std::cerr << "memory benchmark failed: " << error.what() << '\n';
        return 1;
    }
}
