// Standalone diagnostic (NOT a CI speed gate):
//   c++ -O3 -std=c++17 -Isrc benchmarks/linear/bench_csr_coalescing.cpp -o bench_csr
//   ./bench_csr [grid_side=192] [matvec_repetitions=150]
// Run both revisions with identical compiler/CPU and repeat to compare setup
// and matvec costs. Never use this as proof of physical CFD convergence.
#include "cfdx/core/linalg/sparse_matrix.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using cfdx::core::SparseMatrix;
using cfdx::core::Vector;

int main(int argc, char** argv) {
    const std::size_t side = argc > 1 ? std::stoul(argv[1]) : 192;
    const std::size_t repeats = argc > 2 ? std::stoul(argv[2]) : 150;
    if (side < 4 || side > 512 || repeats < 1 || repeats > 10000)
        throw std::invalid_argument("invalid CSR benchmark dimensions");
    const std::size_t cells = side * side;
    SparseMatrix A(cells, cells);
    const auto before = std::chrono::steady_clock::now();
    for (std::size_t y = 0; y < side; ++y) {
        for (std::size_t x = 0; x < side; ++x) {
            const std::size_t i = y * side + x;
            double diag = 1.0;
            const auto neighbor = [&](std::size_t j, double a) {
                // Split one physical face contribution over two COO entries.
                A.push_back(i, j, -a / 2);
                A.push_back(i, j, -a / 2);
                diag += a;
            };
            if (x > 0) neighbor(i - 1, 1000.0);
            if (x + 1 < side) neighbor(i + 1, 1000.0);
            if (y > 0) neighbor(i - side, 1.0);
            if (y + 1 < side) neighbor(i + side, 1.0);
            for (int p = 0; p < 4; ++p) A.push_back(i, i, diag / 4);
        }
    }
    const auto staged = A.nnz();
    A.finalize();
    if (!A.is_consistent()) throw std::runtime_error("bad assembled CSR");
    const auto after = std::chrono::steady_clock::now();
    Vector x(cells, 0.0), y(cells, 0.0);
    for (std::size_t i = 0; i < cells; ++i)
        x(i) = std::sin(0.005 * static_cast<double>(i)) +
               0.2 * std::cos(0.121 * static_cast<double>(i));
    double checksum = 0.0;
    const auto solve_start = std::chrono::steady_clock::now();
    for (std::size_t rep = 0; rep < repeats; ++rep) {
        A.matvec_into(x, y);
        checksum += y((rep * 31) % cells);
    }
    const auto solve_stop = std::chrono::steady_clock::now();
    const auto ms = [](auto dt) {
        return std::chrono::duration<double, std::milli>(dt).count();
    };
    std::cout << std::setprecision(17)
              << "cells=" << cells
              << " staged_nnz=" << staged
              << " csr_nnz=" << A.nnz()
              << " assembly_ms=" << ms(after - before)
              << " matvec_us=" << 1e3 * ms(solve_stop - solve_start) / repeats
              << " checksum=" << checksum << '\n';
}
