#pragma once

// Strict, optional CSV exchange for NACA validation initial guesses.
// These are initial guesses only, never substitutes for a converged CFD solve.
// The import validates every cell index and center before mutating the fields.

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/physics/finite_volume_transport.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::validation {
struct NacaStateRow {
    cfdx::core::Vec3 center;
    double ux = 0.0;
    double uy = 0.0;
    double pressure = 0.0;
};

inline std::vector<std::string> naca_csv_fields(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> out;
    std::string token;
    while (std::getline(stream, token, ',')) out.push_back(token);
    return out;
}

inline double naca_csv_finite_double(const std::string& token) {
    std::size_t end = 0;
    const double value = std::stod(token, &end);
    if (end != token.size() || !std::isfinite(value))
        throw std::invalid_argument("NACA initial state has invalid number");
    return value;
}

inline std::vector<NacaStateRow> naca_read_state_csv(
    const std::string& path, const std::vector<cfdx::core::Vec3>& expected_centers) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open NACA initial field: " + path);
    std::string line;
    if (!std::getline(input, line))
        throw std::invalid_argument("NACA initial field is empty");
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line != "cell,x,y,z,ux,uy,p")
        throw std::invalid_argument("NACA initial field has wrong CSV schema");
    std::vector<NacaStateRow> values;
    values.reserve(expected_centers.size());
    for (std::size_t i = 0; i < expected_centers.size(); ++i) {
        if (!std::getline(input, line))
            throw std::invalid_argument("NACA initial field is incomplete");
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto parts = naca_csv_fields(line);
        if (parts.size() != 7 || parts[0] != std::to_string(i))
            throw std::invalid_argument("NACA initial field index/order mismatch");
        const cfdx::core::Vec3 center{
            naca_csv_finite_double(parts[1]),
            naca_csv_finite_double(parts[2]),
            naca_csv_finite_double(parts[3])};
        const auto& expected = expected_centers[i];
        const double scale = std::max({1.0, expected.mag(), center.mag()});
        if ((center - expected).mag() > 2e-10 * scale)
            throw std::invalid_argument("NACA initial field cell center mismatch");
        values.push_back({center, naca_csv_finite_double(parts[4]),
                                  naca_csv_finite_double(parts[5]),
                                  naca_csv_finite_double(parts[6])});
    }
    if (std::getline(input, line))
        throw std::invalid_argument("NACA initial field contains extra cells");
    return values;
}

inline void naca_write_state_csv(const std::string& path,
                                  const std::vector<NacaStateRow>& values) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot create NACA field CSV: " + path);
    out << "cell,x,y,z,ux,uy,p\n"
        << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i=0;i<values.size();++i) {
        const auto& r=values[i];
        if (!std::isfinite(r.center.x) || !std::isfinite(r.center.y) ||
            !std::isfinite(r.center.z) || !std::isfinite(r.ux) ||
            !std::isfinite(r.uy) || !std::isfinite(r.pressure))
            throw std::invalid_argument("cannot export non-finite NACA field");
        out << i << ',' << r.center.x << ',' << r.center.y << ','
            << r.center.z << ',' << r.ux << ',' << r.uy << ','
            << r.pressure << '\n';
    }
    if (!out.good()) throw std::runtime_error("failed to write NACA field CSV");
}

// Replace only after the entire snapshot has been serialized and closed.
// An interrupted run preserves the preceding checkpoint.
inline void naca_write_state_csv_atomic(
    const std::string& path, const std::vector<NacaStateRow>& values) {
    const std::string temporary = path + ".incomplete";
    naca_write_state_csv(temporary, values);
    if (std::rename(temporary.c_str(), path.c_str()) != 0)
        throw std::runtime_error("NACA checkpoint atomic rename failed");
}

inline void naca_write_centers_csv(
    const std::string& path, const std::vector<cfdx::core::Vec3>& centers) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot create NACA centers CSV");
    out << "cell,x,y,z\n"
        << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i=0;i<centers.size();++i) {
        const auto& x=centers[i];
        if (!std::isfinite(x.x) || !std::isfinite(x.y) || !std::isfinite(x.z))
            throw std::invalid_argument("invalid NACA cell center");
        out << i << ',' << x.x << ',' << x.y << ',' << x.z << '\n';
    }
    if (!out.good()) throw std::runtime_error("failed to write NACA centers");
}
// Same-mesh restart: the conservative mass flux at the accepted iterate
// cannot in general be reconstructed from cell-centred U and p alone.
// Bind every row to both connectivity and geometry so a face sequence from
// another mesh or a reversed orientation cannot silently be accepted.
inline void naca_write_face_flux_csv(
    const std::string& path, const cfdx::core::Mesh& mesh,
    const cfdx::physics::FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& flux) {
    const auto nf = mesh.n_faces();
    if (flux.size() != nf || flux.dimension() != 1 ||
        geometry.face_centres.size() != nf || geometry.face_area_vectors.size() != nf)
        throw std::invalid_argument("NACA face-flux output shape mismatch");
    std::ofstream out(path, std::ios::trunc);
    if (!out) throw std::runtime_error("cannot open NACA converged face-flux CSV");
    out << "face,owner,neighbour,cx,cy,cz,sx,sy,sz,mass_flux\n"
        << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (std::size_t f=0; f<nf; ++f) {
        const auto& cf=geometry.face_centres[f];
        const auto& sf=geometry.face_area_vectors[f];
        if (!std::isfinite(cf.x) || !std::isfinite(cf.y) || !std::isfinite(cf.z) ||
            !std::isfinite(sf.x) || !std::isfinite(sf.y) || !std::isfinite(sf.z) ||
            !std::isfinite(flux(f)))
            throw std::invalid_argument("NACA face-flux export nonfinite");
        out << f << ',' << mesh.ownership().owner(f) << ','
            << mesh.ownership().neighbour(f) << ','
            << cf.x << ',' << cf.y << ',' << cf.z << ','
            << sf.x << ',' << sf.y << ',' << sf.z << ',' << flux(f) << '\n';
    }
    out.flush();
    if (!out) throw std::runtime_error("NACA converged face-flux CSV write failure");
}

inline std::vector<double> naca_read_face_flux_csv(
    const std::string& path, const cfdx::core::Mesh& mesh,
    const cfdx::physics::FvGeometry& geometry) {
    const auto nf=mesh.n_faces();
    if (geometry.face_centres.size() != nf || geometry.face_area_vectors.size() != nf)
        throw std::invalid_argument("NACA restart face geometry size mismatch");
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open NACA initial face flux: " + path);
    std::string line;
    if (!std::getline(in,line)) throw std::invalid_argument("NACA face-flux CSV empty");
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line != "face,owner,neighbour,cx,cy,cz,sx,sy,sz,mass_flux")
        throw std::invalid_argument("NACA face-flux CSV schema mismatch");
    std::vector<double> values;
    values.reserve(nf);
    for (std::size_t f=0; f<nf; ++f) {
        if (!std::getline(in,line)) throw std::invalid_argument("NACA face-flux CSV incomplete");
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto parts=naca_csv_fields(line);
        if (parts.size() != 10 || parts[0] != std::to_string(f) ||
            parts[1] != std::to_string(mesh.ownership().owner(f)) ||
            parts[2] != std::to_string(mesh.ownership().neighbour(f)))
            throw std::invalid_argument("NACA face-flux CSV connectivity mismatch");
        const auto& cf=geometry.face_centres[f];
        const auto& sf=geometry.face_area_vectors[f];
        const double reference[]{cf.x,cf.y,cf.z,sf.x,sf.y,sf.z};
        for (std::size_t j=0; j<6; ++j) {
            const double value=naca_csv_finite_double(parts[j+3]);
            if (std::abs(value-reference[j]) >
                2e-10*std::max({1.0,std::abs(value),std::abs(reference[j])}))
                throw std::invalid_argument("NACA face-flux CSV geometry mismatch");
        }
        values.push_back(naca_csv_finite_double(parts[9]));
    }
    if (std::getline(in,line))
        throw std::invalid_argument("NACA face-flux CSV extra rows");
    return values;
}

} // namespace cfdx::validation
