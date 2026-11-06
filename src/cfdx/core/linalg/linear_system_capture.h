#pragma once

// Deterministic opt-in single-rank CSR capture. The file stores the ORIGINAL
// operator, RHS and initial guess, not an AMG/scaled intermediate matrix.
// Round-trip FP64 bit patterns and a non-cryptographic corruption checksum.
// This is a linear-system replay tool, NOT a physical CFD qualification.

#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/linear_solver_models.h"

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

#if defined(CFDX_HAS_MPI) && CFDX_HAS_MPI
#include <mpi.h>
#endif

namespace cfdx::core {

// Strict utility for CLI/validation programs. The solver controls independently
// enforce max_iterations and require a path for an active capture.
inline std::vector<std::size_t> parse_linear_capture_iterations(
    const std::string& csv, std::size_t max_count = 32) {
    if (csv.empty() || max_count == 0)
        throw std::invalid_argument("capture iterations: missing selection");
    std::vector<std::size_t> values;
    std::istringstream input(csv);
    std::string token;
    while (std::getline(input, token, ',')) {
        if (token.empty() || values.size() >= max_count ||
            !std::all_of(token.begin(), token.end(),
                         [](unsigned char ch) { return ch >= '0' && ch <= '9'; }))
            throw std::invalid_argument("capture iterations: malformed list");
        std::size_t parsed = 0;
        const unsigned long long number = std::stoull(token, &parsed);
        if (parsed != token.size() || number == 0 ||
            number > std::numeric_limits<std::size_t>::max())
            throw std::invalid_argument("capture iterations: invalid value");
        values.push_back(static_cast<std::size_t>(number));
    }
    if (values.empty() || csv.back() == ',')
        throw std::invalid_argument("capture iterations: empty selection");
    auto sorted = values;
    std::sort(sorted.begin(), sorted.end());
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
        throw std::invalid_argument("capture iterations: duplicated selection");
    return values;
}

struct LinearSystemCapture {
    SparseMatrix matrix;
    Vector rhs;
    Vector initial_solution;
    LinearProblemKind problem = LinearProblemKind::General;
    NullSpaceModel null_space = NullSpaceModel::None;
    LinearScalingModel scaling = LinearScalingModel::None;
    std::size_t max_iterations = 0;
    double tolerance = 0.0;
    std::uint64_t checksum = 0;
};

namespace capture_detail {
constexpr std::uint64_t offset = 14695981039346656037ULL;
constexpr std::uint64_t prime = 1099511628211ULL;

inline std::uint64_t bits(double x) {
    static_assert(sizeof(double) == sizeof(std::uint64_t) &&
                  std::numeric_limits<double>::is_iec559, "requires IEEE FP64");
    std::uint64_t out = 0;
    std::memcpy(&out, &x, sizeof(out));
    return out;
}
inline double from_bits(std::uint64_t value) {
    double out = 0.0;
    std::memcpy(&out, &value, sizeof(out));
    return out;
}
inline void add(std::uint64_t& hash, std::uint64_t word) {
    // Explicit byte order makes the checksum platform independent for IEEE754.
    for (int i = 0; i < 8; ++i) {
        hash ^= (word >> (8 * i)) & 0xffU;
        hash *= prime;
    }
}
inline void require(bool valid, const char* message) {
    if (!valid) throw std::invalid_argument(message);
}
inline bool valid_problem(int n) {
    return n >= static_cast<int>(LinearProblemKind::General) &&
           n <= static_cast<int>(LinearProblemKind::CoupledPressureVelocity);
}
inline bool valid_scaling(int n) {
    return n >= static_cast<int>(LinearScalingModel::None) &&
           n <= static_cast<int>(LinearScalingModel::SymmetricDiagonal);
}
inline bool valid_nullspace(int n) {
    return n == static_cast<int>(NullSpaceModel::None) ||
           n == static_cast<int>(NullSpaceModel::Constant);
}
inline void verify(const LinearSystemCapture& capture) {
    const auto& a = capture.matrix;
    require(a.n_rows() != 0 && a.n_rows() == a.n_cols() && a.is_consistent(),
            "linear capture: invalid square finalized CSR");
    // A writer must not create captures that the bounded reader cannot replay.
    require(a.n_rows() <= 2000000 && a.nnz() <= 100000000,
            "linear capture: dimensions exceed replay format limits");
    require(valid_problem(static_cast<int>(capture.problem)) &&
            valid_nullspace(static_cast<int>(capture.null_space)) &&
            valid_scaling(static_cast<int>(capture.scaling)),
            "linear capture: invalid solver metadata");
    require(a.n_rows() == capture.rhs.size() &&
            a.n_cols() == capture.initial_solution.size(),
            "linear capture: mismatched vector size");
    require(capture.max_iterations != 0 && capture.max_iterations <= 100000000 &&
            capture.tolerance > 0.0 && std::isfinite(capture.tolerance),
            "linear capture: invalid iteration controls");
    for (std::size_t i = 0; i < capture.rhs.size(); ++i)
        require(std::isfinite(capture.rhs(i)) &&
                std::isfinite(capture.initial_solution(i)),
                "linear capture: nonfinite vector");
}
inline std::uint64_t fingerprint(const LinearSystemCapture& c) {
    verify(c);
    std::uint64_t h = offset;
    add(h, c.matrix.n_rows());
    add(h, c.matrix.n_cols());
    add(h, c.matrix.nnz());
    add(h, static_cast<std::uint64_t>(c.problem));
    add(h, static_cast<std::uint64_t>(c.null_space));
    add(h, static_cast<std::uint64_t>(c.scaling));
    add(h, c.max_iterations);
    add(h, bits(c.tolerance));
    for (std::size_t i = 0; i <= c.matrix.n_rows(); ++i)
        add(h, c.matrix.row_offsets_data()[i]);
    for (std::size_t k = 0; k < c.matrix.nnz(); ++k) {
        add(h, c.matrix.columns_data()[k]);
        add(h, bits(c.matrix.values_data()[k]));
    }
    for (std::size_t i = 0; i < c.rhs.size(); ++i) {
        add(h, bits(c.rhs(i)));
        add(h, bits(c.initial_solution(i)));
    }
    return h;
}
inline std::uint64_t read_hex(std::istream& in) {
    std::uint64_t word = 0;
    if (!(in >> std::hex >> word >> std::dec))
        throw std::runtime_error("linear capture: truncated or invalid hexadecimal number");
    return word;
}
}  // namespace capture_detail

inline void write_linear_system_capture(
    const std::string& path, const SparseMatrix& matrix,
    const Vector& rhs, const Vector& initial_solution,
    LinearProblemKind problem, NullSpaceModel null_space,
    LinearScalingModel scaling, std::size_t max_iterations, double tolerance) {
    if (path.empty()) throw std::invalid_argument("linear capture: empty path");
#if defined(CFDX_HAS_MPI) && CFDX_HAS_MPI
    // These captures contain a complete SERIAL CSR; a distributed local block
    // cannot be replayed as a global physical operator.  Avoid rank races.
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (initialized) {
        int size = 1;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        if (size != 1)
            throw std::runtime_error(
                "linear capture: distributed MPI snapshots require an explicit gather; refusing");
    }
#endif
    LinearSystemCapture c{matrix, rhs, initial_solution, problem, null_space,
                          scaling, max_iterations, tolerance, 0};
    const auto fingerprint = capture_detail::fingerprint(c);
    const std::string temporary = path + ".tmp";
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("linear capture: cannot create temporary file");
    try {
        out << "CFDX_LINEAR_REPLAY_V1\n"
            << static_cast<int>(problem) << ' '
            << static_cast<int>(null_space) << ' '
            << static_cast<int>(scaling) << '\n'
            << matrix.n_rows() << ' ' << matrix.n_cols() << ' '
            << matrix.nnz() << ' ' << max_iterations << '\n';
        out << std::hex << capture_detail::bits(tolerance) << '\n';
        for (std::size_t i = 0; i <= matrix.n_rows(); ++i)
            out << std::hex << matrix.row_offsets_data()[i] << '\n';
        for (std::size_t k = 0; k < matrix.nnz(); ++k)
            out << std::hex << matrix.columns_data()[k] << ' '
                << capture_detail::bits(matrix.values_data()[k]) << '\n';
        for (std::size_t i = 0; i < rhs.size(); ++i)
            out << std::hex << capture_detail::bits(rhs(i)) << ' '
                << capture_detail::bits(initial_solution(i)) << '\n';
        out << std::hex << fingerprint << '\n';
        out.flush();
        if (!out) throw std::runtime_error("linear capture: disk write failed");
        out.close();
        if (!out) throw std::runtime_error("linear capture: disk close failed");
#if defined(_WIN32)
        // Windows rename() refuses an existing destination. Capture retries
        // must replace the previous complete snapshot, never leave a partial one.
        if (!MoveFileExA(temporary.c_str(), path.c_str(),
                         MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
            throw std::runtime_error("linear capture: atomic replacement failed");
#else
        if (std::rename(temporary.c_str(), path.c_str()) != 0)
            throw std::runtime_error("linear capture: atomic rename failed");
#endif
    } catch (...) {
        out.close();
        std::remove(temporary.c_str());
        throw;
    }
}

inline LinearSystemCapture read_linear_system_capture(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("linear capture: cannot open file");
    std::string header;
    if (!(in >> header) || header != "CFDX_LINEAR_REPLAY_V1")
        throw std::runtime_error("linear capture: unsupported format");
    int problem = -1, nullspace = -1, scaling = -1;
    std::size_t rows = 0, columns = 0, nnz = 0, iterations = 0;
    if (!(in >> problem >> nullspace >> scaling >> rows >> columns >> nnz >> iterations))
        throw std::runtime_error("linear capture: invalid header");
    capture_detail::require(capture_detail::valid_problem(problem) &&
        capture_detail::valid_nullspace(nullspace) &&
        capture_detail::valid_scaling(scaling), "linear capture: invalid solver metadata");
    // Prevent accidental/hostile enormous allocations before reading the CSR.
    capture_detail::require(rows > 0 && rows == columns && rows <= 2000000 &&
        nnz <= 100000000 && nnz <= std::numeric_limits<SparseMatrix::Index>::max(),
        "linear capture: dimensions exceed supported limits");
    const double tolerance = capture_detail::from_bits(capture_detail::read_hex(in));
    std::vector<SparseMatrix::Index> offsets(rows + 1);
    std::vector<SparseMatrix::Index> indices(nnz);
    std::vector<double> values(nnz);
    for (auto& v : offsets) {
        const auto word = capture_detail::read_hex(in);
        if (word > std::numeric_limits<SparseMatrix::Index>::max())
            throw std::runtime_error("linear capture: overflowing CSR offset");
        v = static_cast<SparseMatrix::Index>(word);
    }
    for (std::size_t k = 0; k < nnz; ++k) {
        const auto index = capture_detail::read_hex(in);
        if (index >= columns) throw std::runtime_error("linear capture: invalid CSR column");
        indices[k] = static_cast<SparseMatrix::Index>(index);
        values[k] = capture_detail::from_bits(capture_detail::read_hex(in));
    }
    Vector rhs(rows, 0.0), guess(rows, 0.0);
    for (std::size_t i = 0; i < rows; ++i) {
        rhs(i) = capture_detail::from_bits(capture_detail::read_hex(in));
        guess(i) = capture_detail::from_bits(capture_detail::read_hex(in));
    }
    const auto expected = capture_detail::read_hex(in);
    std::string trailing;
    if (in >> trailing) throw std::runtime_error("linear capture: unexpected trailing data");
    LinearSystemCapture c{
        SparseMatrix::from_sorted_csr(rows, columns, std::move(offsets),
            std::move(indices), std::move(values)),
        std::move(rhs), std::move(guess),
        static_cast<LinearProblemKind>(problem),
        static_cast<NullSpaceModel>(nullspace),
        static_cast<LinearScalingModel>(scaling), iterations, tolerance, expected};
    if (capture_detail::fingerprint(c) != expected)
        throw std::runtime_error("linear capture: checksum mismatch");
    return c;
}

} // namespace cfdx::core
