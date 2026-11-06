// Standalone reproducible Galerkin numeric-refresh microbenchmark.
// Build against V5 or V6 headers with identical flags and input size.
// The matrix is a stretched-grid, anisotropic five-point diffusion stencil.
#include "cfdx/core/linalg/galerkin_csr_product.h"
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace cfdx::core;
int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: benchmark_galerkin_workspace nside refresh_count\n";
        return 2;
    }
    try {
        const std::size_t side = std::stoul(argv[1]);
        const int refreshes = std::stoi(argv[2]);
        if (side < 4 || side > 256 || refreshes < 1 || refreshes > 200)
            throw std::invalid_argument("size or refresh count out of range");
        const std::size_t n = side * side, coarse_n = n / 2;
        SparseMatrix A(n, n);
        for (std::size_t y = 0; y < side; ++y)
            for (std::size_t x = 0; x < side; ++x) {
                const std::size_t i = y*side + x;
                const double ax = 1e3, ay = 1.0;
                double diagonal = 1.0;
                if (x > 0) { A.push_back(i, i-1, -ax); diagonal += ax; }
                if (x+1 < side) { A.push_back(i, i+1, -ax); diagonal += ax; }
                if (y > 0) { A.push_back(i, i-side, -ay); diagonal += ay; }
                if (y+1 < side) { A.push_back(i, i+side, -ay); diagonal += ay; }
                A.push_back(i, i, diagonal);
            }
        A.finalize();
        GalerkinCSRProduct::TransferRows P(n);
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t ci = i/2;
            if (ci + 1 < coarse_n)
                P[i] = {{ci,0.75}, {ci+1,0.25}};
            else P[i] = {{ci,1.0}};
        }
        const auto t0 = std::chrono::steady_clock::now();
        GalerkinCSRProduct plan;
        if (!plan.symbolic(A,P,coarse_n)) throw std::runtime_error("symbolic rejected");
        SparseMatrix coarse;
        if (!plan.numeric(A,P,coarse)) throw std::runtime_error("initial numeric rejected");
        const double symbolic_setup_ms = std::chrono::duration<double,std::milli>(
            std::chrono::steady_clock::now()-t0).count();
        const auto t1 = std::chrono::steady_clock::now();
        for (int t = 0; t < refreshes; ++t) {
            for (std::size_t k = 0; k < A.nnz(); ++k)
                A.values_data()[k] *= (t&1) ? (1.0/1.0001) : 1.0001;
            if (!plan.numeric_into_existing(A,P,coarse))
                throw std::runtime_error("numeric refresh rejected");
        }
        const double refresh_ms = std::chrono::duration<double,std::milli>(
            std::chrono::steady_clock::now()-t1).count();
        // Checksum is robust to the two versions' different temp allocation
        // policies but catches numerical path changes across repeated updates.
        long double checksum = 0.0L;
        for (std::size_t k=0;k<coarse.nnz();++k)
            checksum += static_cast<long double>(coarse.values_data()[k]);
        if (!std::isfinite(static_cast<double>(checksum)))
            throw std::runtime_error("non-finite Galerkin matrix");
        std::cout << std::setprecision(17)
                  << "n=" << n << " coarse_n=" << coarse_n
                  << " refreshes=" << refreshes
                  << " setup_ms=" << symbolic_setup_ms
                  << " refresh_total_ms=" << refresh_ms
                  << " refresh_mean_ms=" << refresh_ms/refreshes
                  << " coarse_nnz=" << coarse.nnz()
                  << " checksum=" << static_cast<double>(checksum) << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "GALERKIN_BENCH_FAIL: " << error.what() << '\n';
        return 1;
    }
}
