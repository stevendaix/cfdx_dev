#pragma once
// Fail-closed same-mesh NACA checkpoint for an accepted nonlinear U,p,phi step.
// Separate immutable payloads and a single atomically replaced manifest form
// a publication barrier. A crash cannot replace the last complete checkpoint
// with a half-new U,p or face-flux file. Not a physical-convergence certificate.
#include "naca_qualification_state_io.h"
#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <system_error>

namespace cfdx::validation {
struct NacaPairedCheckpoint {
    std::size_t accepted_iteration = 0;
    double momentum_reference = 0.0;
    double nonlinear_monitor_reference = 0.0;
    std::vector<NacaStateRow> cells;
    std::vector<double> mass_flux;
};

inline std::string naca_file_fnv64(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("NACA checkpoint: cannot hash file: " + path.string());
    std::uint64_t h = UINT64_C(14695981039346656037);
    char buf[32768];
    while (in) {
        in.read(buf, sizeof(buf));
        for (std::streamsize k=0;k<in.gcount();++k) {
            h ^= static_cast<unsigned char>(buf[k]);
            h *= UINT64_C(1099511628211);
        }
    }
    if (!in.eof()) throw std::runtime_error("NACA checkpoint: hash input I/O failure");
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << h;
    return out.str();
}

inline bool naca_checkpoint_safe_name(const std::string& name) {
    if (name.empty() || name == "." || name == ".." || name.find("..") != std::string::npos)
        return false;
    return std::all_of(name.begin(), name.end(), [](char ch) {
        return (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') ||
               (ch >= '0' && ch <= '9') || ch == '.' || ch == '_' || ch == '-';
    });
}

inline void naca_write_paired_checkpoint_atomic(
    const std::string& manifest_path, std::size_t accepted_iteration,
    const std::vector<NacaStateRow>& cells,
    const cfdx::core::Mesh& mesh, const cfdx::physics::FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    double momentum_reference, double nonlinear_monitor_reference) {
    namespace fs=std::filesystem;
    if (manifest_path.empty() || accepted_iteration == 0 ||
        cells.size() != mesh.n_cells() || mass_flux.size() != mesh.n_faces() ||
        !std::isfinite(momentum_reference) || !(momentum_reference > 0.0) ||
        !std::isfinite(nonlinear_monitor_reference) ||
        !(nonlinear_monitor_reference > 0.0))
        throw std::invalid_argument("NACA paired checkpoint: invalid shapes/iteration/path");
    const fs::path manifest(manifest_path);
    static std::atomic<std::uint64_t> generation{0};
    const auto nonce = std::to_string(
        std::chrono::steady_clock::now().time_since_epoch().count()) + "_" +
        std::to_string(generation.fetch_add(1));
    const auto prefix = "naca_step_" + std::to_string(accepted_iteration) + "_" + nonce;
    const std::string state_name = prefix + ".state.csv";
    const std::string flux_name = prefix + ".flux.csv";
    const fs::path parent = manifest.parent_path();
    const fs::path state_file = parent / state_name;
    const fs::path flux_file = parent / flux_name;
    // Immutable names: interrupted writes only leave unreachable orphan files.
    naca_write_state_csv(state_file.string(), cells);
    naca_write_face_flux_csv(flux_file.string(), mesh, geometry, mass_flux);
    const auto state_hash=naca_file_fnv64(state_file);
    const auto flux_hash=naca_file_fnv64(flux_file);
    const fs::path temporary=manifest_path + ".incomplete";
    {
        std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
        if (!out) throw std::runtime_error("NACA checkpoint: cannot create manifest");
        out << std::setprecision(std::numeric_limits<double>::max_digits10)
            << "CFDX_NACA_PAIRED_V2\n"
            << accepted_iteration << ',' << mesh.n_cells() << ',' << mesh.n_faces()
            << ',' << momentum_reference << ',' << nonlinear_monitor_reference
            << ',' << state_name << ',' << flux_name << ',' << state_hash
            << ',' << flux_hash << '\n';
        out.flush();
        if (!out) throw std::runtime_error("NACA checkpoint: cannot write manifest");
    }
    if (std::rename(temporary.string().c_str(), manifest.string().c_str()) != 0)
        throw std::runtime_error("NACA checkpoint: cannot atomically publish manifest");
}

inline NacaPairedCheckpoint naca_read_paired_checkpoint(
    const std::string& manifest_path, const cfdx::core::Mesh& mesh,
    const cfdx::physics::FvGeometry& geometry) {
    namespace fs=std::filesystem;
    std::ifstream in(manifest_path, std::ios::binary);
    if (!in) throw std::runtime_error("NACA paired checkpoint: manifest missing");
    std::string magic, line, extra;
    if (!std::getline(in,magic) || magic != "CFDX_NACA_PAIRED_V2" ||
        !std::getline(in,line) || std::getline(in,extra))
        throw std::invalid_argument("NACA paired checkpoint: invalid manifest schema");
    const auto parts=naca_csv_fields(line);
    if (parts.size()!=9) throw std::invalid_argument("NACA checkpoint: invalid manifest fields");
    const auto strict_number=[](const std::string& value) {
        std::size_t consumed=0;
        const auto parsed=std::stoull(value, &consumed);
        if (consumed != value.size()) throw std::invalid_argument("NACA checkpoint: invalid count");
        return parsed;
    };
    const auto iteration=strict_number(parts[0]);
    if (iteration == 0 || strict_number(parts[1]) != mesh.n_cells() ||
        strict_number(parts[2]) != mesh.n_faces())
        throw std::invalid_argument("NACA paired checkpoint: mesh size/iteration mismatch");
    const double momentum_reference=naca_csv_finite_double(parts[3]);
    const double nonlinear_reference=naca_csv_finite_double(parts[4]);
    if (!(momentum_reference > 0.0) || !(nonlinear_reference > 0.0))
        throw std::invalid_argument("NACA checkpoint: invalid normalization baseline");
    if (!naca_checkpoint_safe_name(parts[5]) || !naca_checkpoint_safe_name(parts[6]))
        throw std::invalid_argument("NACA paired checkpoint: invalid payload filename");
    const auto parent=fs::path(manifest_path).parent_path();
    const auto state_path=parent / parts[5];
    const auto flux_path=parent / parts[6];
    if (naca_file_fnv64(state_path) != parts[7] ||
        naca_file_fnv64(flux_path) != parts[8])
        throw std::invalid_argument("NACA paired checkpoint: payload checksum mismatch");
    NacaPairedCheckpoint result;
    result.accepted_iteration=static_cast<std::size_t>(iteration);
    result.momentum_reference=momentum_reference;
    result.nonlinear_monitor_reference=nonlinear_reference;
    result.cells=naca_read_state_csv(state_path.string(),geometry.cell_centres);
    result.mass_flux=naca_read_face_flux_csv(flux_path.string(),mesh,geometry);
    return result;
}
} // namespace cfdx::validation
