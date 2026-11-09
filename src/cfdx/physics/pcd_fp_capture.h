#pragma once

// Byte-exact capture of the ORIGINAL pressure-convection/diffusion operator
// and the exact input to an inner PCD Fp solve. This is NOT a physics oracle.
// CFP01 is a versioned, little-endian IEEE-754/uint32 CSR replay fixture.
// Capture is opt-in, bounded and never enabled by default in production.

#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"

#include <array>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

// Bound diagnostic data by PROCESS rather than by nonlinear/Schur solve.
// Otherwise a physical N9 run can silently overwrite the first Couette
// evidence with later Poiseuille/Ghia systems and perform unbounded disk IO.
// Ordinary and difficult samples have independent three-file budgets.
inline std::size_t reserve_pcd_fp_capture_slot(bool difficult) {
    static std::atomic<std::size_t> first_samples{0};
    static std::atomic<std::size_t> difficult_samples{0};
    auto& counter = difficult ? difficult_samples : first_samples;
    const auto prior = counter.fetch_add(1, std::memory_order_relaxed);
    return prior < 3 ? prior + 1 : 0;
}

struct PcdFpCapture {
    cfdx::core::SparseMatrix matrix;
    cfdx::core::Vector rhs;
    cfdx::core::Vector initial;
};

inline void assert_pcd_capture_platform() {
    const std::uint32_t one = 1;
    if (sizeof(double) != 8 || sizeof(std::uint32_t) != 4 ||
        !std::numeric_limits<double>::is_iec559 ||
        *reinterpret_cast<const unsigned char*>(&one) != 1)
        throw std::runtime_error("CFP01 requires little-endian IEEE754 double and uint32 CSR");
}

template<typename T>
inline void pcd_capture_write(std::ofstream& out, const T* ptr, std::size_t count) {
    if (count && (!ptr || count > static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max()) / sizeof(T)))
        throw std::runtime_error("CFP01 invalid payload size");
    out.write(reinterpret_cast<const char*>(ptr), static_cast<std::streamsize>(count * sizeof(T)));
    if (!out) throw std::runtime_error("CFP01 write failed");
}

template<typename T>
inline void pcd_capture_read(std::ifstream& in, T* ptr, std::size_t count) {
    if (count && (!ptr || count > static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max()) / sizeof(T)))
        throw std::runtime_error("CFP01 invalid input size");
    in.read(reinterpret_cast<char*>(ptr), static_cast<std::streamsize>(count * sizeof(T)));
    if (!in) throw std::runtime_error("CFP01 truncated input");
}

inline void write_pcd_fp_capture(const std::filesystem::path& path,
    const cfdx::core::SparseMatrix& A,
    const cfdx::core::Vector& rhs,
    const cfdx::core::Vector& initial) {
    assert_pcd_capture_platform();
    const auto n = A.n_rows(), nnz = A.nnz();
    if (n == 0 || n != A.n_cols() || !A.is_consistent() ||
        n > std::numeric_limits<std::uint32_t>::max() ||
        nnz > std::numeric_limits<std::uint32_t>::max() ||
        rhs.size() != n || initial.size() != n)
        throw std::invalid_argument("CFP01 requires valid square CSR and matching vectors");
    for (std::size_t i = 0; i < n; ++i)
        if (!std::isfinite(rhs(i)) || !std::isfinite(initial(i)))
            throw std::invalid_argument("CFP01 nonfinite input vectors");
    if (!path.parent_path().empty())
        std::filesystem::create_directories(path.parent_path());
    const auto tmp = path.string() + ".tmp";
    try {
        std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
        if (!out) throw std::runtime_error("CFP01 cannot open temporary file");
        constexpr char magic[8] = {'C', 'F', 'D', 'X', 'F', 'P', '0', '1'};
        const std::uint64_t counts[2] = {
            static_cast<std::uint64_t>(n), static_cast<std::uint64_t>(nnz)};
        pcd_capture_write(out, magic, 8);
        pcd_capture_write(out, counts, 2);
        pcd_capture_write(out, A.row_offsets_data(), n + 1);
        pcd_capture_write(out, A.columns_data(), nnz);
        pcd_capture_write(out, A.values_data(), nnz);
        pcd_capture_write(out, rhs.data(), n);
        pcd_capture_write(out, initial.data(), n);
        out.close();
        if (!out) throw std::runtime_error("CFP01 close failed");
        std::filesystem::rename(tmp, path);
    } catch (...) {
        std::error_code ec;
        std::filesystem::remove(tmp, ec);
        throw;
    }
}

inline PcdFpCapture read_pcd_fp_capture(const std::filesystem::path& path) {
    assert_pcd_capture_platform();
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("CFP01 capture not found: " + path.string());
    char magic[8]{};
    pcd_capture_read(in, magic, 8);
    if (std::memcmp(magic, "CFDXFP01", 8))
        throw std::runtime_error("CFP01 invalid magic/version");
    std::uint64_t counts[2]{};
    pcd_capture_read(in, counts, 2);
    const auto n = static_cast<std::size_t>(counts[0]);
    const auto nnz = static_cast<std::size_t>(counts[1]);
    if (n == 0 || n > 1000000 || nnz > 20000000 || nnz < n ||
        n > std::numeric_limits<std::uint32_t>::max() ||
        nnz > std::numeric_limits<std::uint32_t>::max())
        throw std::runtime_error("CFP01 unreasonable CSR dimensions");
    using Index = cfdx::core::SparseMatrix::Index;
    std::vector<Index> row(n + 1), col(nnz);
    std::vector<double> val(nnz);
    cfdx::core::Vector rhs(n), initial(n);
    pcd_capture_read(in, row.data(), row.size());
    pcd_capture_read(in, col.data(), col.size());
    pcd_capture_read(in, val.data(), val.size());
    pcd_capture_read(in, rhs.data(), n);
    pcd_capture_read(in, initial.data(), n);
    if (in.peek() != std::char_traits<char>::eof())
        throw std::runtime_error("CFP01 extra bytes at end of capture");
    for (std::size_t i = 0; i < n; ++i)
        if (!std::isfinite(rhs(i)) || !std::isfinite(initial(i)))
            throw std::runtime_error("CFP01 nonfinite captured vector");
    PcdFpCapture result;
    result.matrix = cfdx::core::SparseMatrix::from_sorted_csr(
        n, n, std::move(row), std::move(col), std::move(val));
    result.rhs = std::move(rhs);
    result.initial = std::move(initial);
    return result;
}

} // namespace cfdx::physics
