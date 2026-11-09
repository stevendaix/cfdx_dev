#pragma once
// Bounded, optional, binary full-cell observation stream. NOT a convergence gate.
// All multibyte values are explicitly little-endian; floats are IEEE-754.
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {
class CellFieldMapWriter {
public:
    static constexpr std::uint64_t MAX_FILE_BYTES = 512ull * 1024 * 1024;
    CellFieldMapWriter(const std::string& filename, const std::string& mesh_source,
                       std::size_t n_cells, std::size_t n_faces, double reference_tolerance)
      : cells_(n_cells) {
        if (filename.empty() || mesh_source.empty() || mesh_source.size() > 4096 ||
            cells_ == 0 || cells_ > 5'000'000 || n_faces == 0 ||
            !std::isfinite(reference_tolerance) || reference_tolerance <= 0 ||
            std::filesystem::is_symlink(filename))
            throw std::invalid_argument("invalid full-cell map configuration");
        static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559,
                      "cell map requires IEEE-754 binary64");
        out_.open(filename, std::ios::binary | std::ios::trunc);
        if (!out_) throw std::runtime_error("cannot create cell map stream");
        raw("CFDXMAP1", 8);
        u64(static_cast<std::uint64_t>(cells_));
        u64(static_cast<std::uint64_t>(n_faces));
        u32(static_cast<std::uint32_t>(mesh_source.size()));
        f64(reference_tolerance);
        raw(mesh_source.data(), mesh_source.size());
        bytes_ = 36ull + mesh_source.size();
        flush();
    }

    void write_frame(std::size_t iteration, bool comparable,
                     const std::vector<double>& momentum,
                     const std::vector<double>& mass,
                     const std::vector<std::uint16_t>& streak) {
        if (iteration <= last_iteration_ || momentum.size() != cells_ ||
            mass.size() != cells_ || streak.size() != cells_)
            throw std::invalid_argument("cell map frame has inconsistent topology or iteration");
        const std::uint64_t frame_bytes = 4ull + 8ull + 1ull + cells_ * 19ull + 4ull;
        if (bytes_ + frame_bytes + 12 > MAX_FILE_BYTES)
            throw std::runtime_error("full-cell map exceeds 512 MiB I/O budget; increase sampling period");
        raw("FRME", 4);
        u64(static_cast<std::uint64_t>(iteration));
        u8(comparable ? 1 : 0);
        for (std::size_t c = 0; c < cells_; ++c) {
            const bool valid = std::isfinite(momentum[c]) && momentum[c] >= 0.0 &&
                std::isfinite(mass[c]) && mass[c] >= 0.0;
            f64(valid ? momentum[c] : std::numeric_limits<double>::quiet_NaN());
            f64(valid ? mass[c] : std::numeric_limits<double>::quiet_NaN());
            u16(valid ? streak[c] : 0);
            u8(valid ? 1 : 0);
        }
        raw("DONE", 4); // reader ignores a live, uncommitted frame
        bytes_ += frame_bytes;
        last_iteration_ = iteration;
        flush();
    }

    void finish(std::size_t accepted_iterations) {
        if (accepted_iterations < last_iteration_)
            throw std::logic_error("full-cell map end precedes its last frame");
        if (bytes_ + 12 > MAX_FILE_BYTES) throw std::runtime_error("full-cell map size budget exceeded");
        raw("CFND", 4);
        u64(static_cast<std::uint64_t>(accepted_iterations));
        bytes_ += 12;
        flush();
    }
private:
    void raw(const char* data, std::size_t bytes) { out_.write(data, static_cast<std::streamsize>(bytes)); }
    void u8(std::uint8_t value) { const char byte = static_cast<char>(value); raw(&byte, 1); }
    void u16(std::uint16_t value) { for (unsigned i = 0; i < 2; ++i) u8(static_cast<std::uint8_t>(value >> (8*i))); }
    void u32(std::uint32_t value) { for (unsigned i = 0; i < 4; ++i) u8(static_cast<std::uint8_t>(value >> (8*i))); }
    void u64(std::uint64_t value) { for (unsigned i = 0; i < 8; ++i) u8(static_cast<std::uint8_t>(value >> (8*i))); }
    void f64(double value) { std::uint64_t bits = 0; std::memcpy(&bits, &value, 8); u64(bits); }
    void flush() { out_.flush(); if (!out_) throw std::runtime_error("failed to write full-cell map"); }
    std::ofstream out_;
    std::size_t cells_;
    std::size_t last_iteration_ = 0;
    std::uint64_t bytes_ = 0;
};
} // namespace cfdx::physics
