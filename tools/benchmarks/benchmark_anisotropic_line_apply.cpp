// Optional, reproducible local CPU microbenchmark for the anisotropic line-AMG
// preconditioner. This measures SETUP and APPLICATION, not full CFD convergence.
// Compare versions on identical hardware, compiler flags and matrices.
#include "cfdx/core/linalg/anisotropic_line_amg.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;

namespace {
SparseMatrix stretched_diffusion(std::size_t side) {
    const std::size_t n = side * side;
    SparseMatrix A(n, n);
    constexpr double strong = 1e4;
    constexpr double weak = 1.0;
    for (std::size_t y = 0; y < side; ++y) {
        for (std::size_t x = 0; x < side; ++x) {
            const std::size_t i = y * side + x;
            double diagonal = 1.0;
            if (x > 0) { A.push_back(i, i - 1, -strong); diagonal += strong; }
            if (x + 1 < side) { A.push_back(i, i + 1, -strong); diagonal += strong; }
            if (y > 0) { A.push_back(i, i - side, -weak); diagonal += weak; }
            if (y + 1 < side) { A.push_back(i, i + side, -weak); diagonal += weak; }
            A.push_back(i, i, diagonal);
        }
    }
    A.finalize();
    return A;
}
}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: cfdx_benchmark_anisotropic_line_apply grid_side repeats\n";
        return 2;
    }
    try {
        const std::size_t side = std::stoul(argv[1]);
        const int repeats = std::stoi(argv[2]);
        if (side < 4 || side > 256 || repeats < 1 || repeats > 10000)
            throw std::invalid_argument("grid_side must be 4..256 and repeats 1..10000");
        const SparseMatrix A = stretched_diffusion(side);
        AnisotropicLineAMGPreconditioner pc;
        const auto prepare_start = std::chrono::steady_clock::now();
        if (!pc.setup(A)) throw std::runtime_error("anisotropic Line-AMG setup rejected");
        const auto prepare_end = std::chrono::steady_clock::now();

        Vector rhs(A.n_rows(), 0.0), out(A.n_rows(), 0.0);
        for (std::size_t i = 0; i < rhs.size(); ++i)
            rhs(i) = 1.0 + std::sin(0.01 * static_cast<double>(i));
        if (!pc.apply(rhs, out)) throw std::runtime_error("warmup application failed");

        std::vector<double> trial_ms;
        trial_ms.reserve(4);
        for (int trial = 0; trial < 4; ++trial) {
            const auto begin = std::chrono::steady_clock::now();
            for (int i = 0; i < repeats; ++i)
                if (!pc.apply(rhs, out))
                    throw std::runtime_error("preconditioner application failed");
            const auto end = std::chrono::steady_clock::now();
            trial_ms.push_back(std::chrono::duration<double, std::milli>(end - begin).count());
        }
        std::sort(trial_ms.begin(), trial_ms.end());
        Vector Ax(A.n_rows(), 0.0);
        A.matvec_into(out, Ax);
        double residual_sq = 0.0, rhs_sq = 0.0, checksum = 0.0;
        for (std::size_t i = 0; i < rhs.size(); ++i) {
            const double ri = rhs(i) - Ax(i);
            residual_sq += ri * ri;
            rhs_sq += rhs(i) * rhs(i);
            checksum += out(i);
        }
        if (!(rhs_sq > 0.0) || !std::isfinite(residual_sq) ||
            !std::isfinite(checksum))
            throw std::runtime_error("non-finite numerical evidence");
        std::cout << std::setprecision(16)
                  << "n=" << A.n_rows() << " repeats=" << repeats
                  << " setup_ms=" << std::chrono::duration<double, std::milli>(
                      prepare_end - prepare_start).count()
                  << " median_apply_ms=" << (trial_ms[1] + trial_ms[2]) / (2 * repeats)
                  << " residual=" << std::sqrt(residual_sq / rhs_sq)
                  << " checksum=" << checksum << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "LINE_AMG_BENCH_FAIL: " << error.what() << '\n';
        return 1;
    }
}
