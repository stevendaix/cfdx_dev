#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace cfdx::io {

// DAT format versions.
//   1 — velocity and pressure only.
//   2 — adds an optional T / k / second turbulence section.
//   3 — adds the mesh identity that a restart reader validates against.
enum : int {
    DatRestartVersionFields = 2,
    DatRestartVersionIdentity = 3,
};

// Geometry identity of the mesh a checkpoint belongs to.
//
// A cell count cannot distinguish two different meshes. Validating a restart
// by cell count alone let a checkpoint written on one mesh load its field
// values into any other mesh of the same size, which is a silent wrong-state
// restart rather than a reported failure.
//
// The identity describes geometry: the point/face/cell counts and an FNV-1a
// checksum over the point coordinates. It does not describe the cell ordering,
// so a checkpoint whose cells are permuted over identical geometry is not
// distinguished here. That case is a different problem and belongs to the
// persistent cell-ID remapping used by the distributed checkpoint path.
struct DatRestartIdentity {
    std::size_t n_points = 0;
    std::size_t n_faces = 0;
    std::size_t n_cells = 0;
    std::uint64_t geometry_checksum = 0;

    bool operator==(const DatRestartIdentity& other) const {
        return n_points == other.n_points && n_faces == other.n_faces &&
               n_cells == other.n_cells &&
               geometry_checksum == other.geometry_checksum;
    }
    bool operator!=(const DatRestartIdentity& other) const { return !(*this == other); }
};

inline std::uint64_t dat_restart_geometry_checksum(const cfdx::core::Mesh& mesh)
{
    std::uint64_t hash = 1469598103934665603ULL;
    const auto mix = [&hash](std::uint64_t word) {
        for (int byte = 0; byte < 8; ++byte) {
            hash ^= static_cast<std::uint8_t>((word >> (8 * byte)) & 0xFFULL);
            hash *= 1099511628211ULL;
        }
    };
    mix(static_cast<std::uint64_t>(mesh.n_points()));
    mix(static_cast<std::uint64_t>(mesh.n_faces()));
    mix(static_cast<std::uint64_t>(mesh.n_cells()));
    for (std::size_t i = 0; i < mesh.n_points(); ++i) {
        const double coordinates[3] = {
            mesh.points().x(i), mesh.points().y(i), mesh.points().z(i)};
        for (const double coordinate : coordinates) {
            if (!std::isfinite(coordinate))
                throw std::invalid_argument(
                    "dat_restart_identity: non-finite point coordinate");
            // -0.0 and 0.0 are the same point, so they must checksum alike.
            const double normalised = (coordinate == 0.0) ? 0.0 : coordinate;
            std::uint64_t bits = 0;
            std::memcpy(&bits, &normalised, sizeof(bits));
            mix(bits);
        }
    }
    return hash;
}

inline DatRestartIdentity dat_restart_identity(const cfdx::core::Mesh& mesh)
{
    return DatRestartIdentity{
        mesh.n_points(), mesh.n_faces(), mesh.n_cells(),
        dat_restart_geometry_checksum(mesh)};
}

inline std::string describe_dat_restart_identity(const DatRestartIdentity& id)
{
    return "points=" + std::to_string(id.n_points) +
           " faces=" + std::to_string(id.n_faces) +
           " cells=" + std::to_string(id.n_cells) +
           " geometry_checksum=" + std::to_string(id.geometry_checksum);
}

// Settle the mesh identity of a checkpoint before any field value is read.
//
// A checkpoint written before this contract carries no identity at all. It is
// rejected rather than accepted on the strength of its cell count, because a
// cell count is exactly the check that cannot tell two meshes apart.
inline void validate_dat_restart_identity(
    const cfdx::core::Mesh& mesh, int version, std::istream& in,
    const char* reader = "read_dat_restart")
{
    if (version < DatRestartVersionIdentity)
        throw std::runtime_error(
            std::string(reader) +
            ": checkpoint predates the mesh identity contract (DAT version " +
            std::to_string(version) + "); regenerate it with a current writer");
    std::string key;
    DatRestartIdentity recorded;
    if (!(in >> key) || key != "identity")
        throw std::runtime_error(std::string(reader) + ": missing mesh identity");
    if (!(in >> recorded.n_points >> recorded.n_faces >> recorded.n_cells >>
          recorded.geometry_checksum))
        throw std::runtime_error(std::string(reader) + ": malformed mesh identity");
    const auto expected = dat_restart_identity(mesh);
    if (recorded != expected)
        throw std::runtime_error(
            std::string(reader) + ": checkpoint belongs to a different mesh (" +
            describe_dat_restart_identity(recorded) + ", expected " +
            describe_dat_restart_identity(expected) + ")");
}

struct DatRestartState {
    std::size_t iteration = 0;
    double time = 0.0;
};

struct DatRestartFields {
    cfdx::core::Field<double, cfdx::core::Location::CELL>* temperature = nullptr;
    cfdx::core::Field<double, cfdx::core::Location::CELL>* k = nullptr;
    cfdx::core::Field<double, cfdx::core::Location::CELL>* second_turbulence = nullptr;
};

inline void validate_optional_restart_field(
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* field,
    const char* name)
{
    if (!field) return;
    if (field->dimension() != 1 || field->size() != mesh.n_cells())
        throw std::invalid_argument(std::string("DAT restart: invalid ") + name + " field");
}

inline void write_dat_restart_fields(
    const std::string& path,
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    const DatRestartFields& fields,
    std::size_t iteration = 0,
    double time = 0.0)
{
    if (U.dimension() != 3 || U.size() != mesh.n_cells())
        throw std::invalid_argument("write_dat_restart_fields: invalid velocity field");
    if (p.dimension() != 1 || p.size() != mesh.n_cells())
        throw std::invalid_argument("write_dat_restart_fields: invalid pressure field");
    if (!std::isfinite(time))
        throw std::invalid_argument("write_dat_restart_fields: non-finite time");
    validate_optional_restart_field(mesh, fields.temperature, "temperature");
    validate_optional_restart_field(mesh, fields.k, "k");
    validate_optional_restart_field(mesh, fields.second_turbulence, "second turbulence");

    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        if (!std::isfinite(U.component_data(0)[c]) || !std::isfinite(U.component_data(1)[c]) ||
            !std::isfinite(U.component_data(2)[c]) || !std::isfinite(p(c)))
            throw std::invalid_argument("write_dat_restart_fields: non-finite field value");
    }
    // Every optional field is checked here, before the output file is opened.
    // Validating them inside the write loop would leave the temporary file
    // behind on a non-finite value, because the stream is still open there.
    const auto reject_non_finite = [](const cfdx::core::Field<double,cfdx::core::Location::CELL>* field,
                                      const char* name) {
        if (!field) return;
        for (std::size_t c = 0; c < field->size(); ++c)
            if (!std::isfinite((*field)(c)))
                throw std::invalid_argument(
                    std::string("write_dat_restart_fields: non-finite ") + name);
    };
    reject_non_finite(fields.temperature, "temperature");
    reject_non_finite(fields.k, "k");
    reject_non_finite(fields.second_turbulence, "second turbulence");

    const auto identity = dat_restart_identity(mesh);
    const std::string tmp_path = path + ".tmp";
    std::ofstream out(tmp_path);
    if (!out) throw std::runtime_error("write_dat_restart_fields: cannot open " + path);
    out.precision(std::numeric_limits<double>::max_digits10);
    out << "CFDX-DAT " << DatRestartVersionIdentity << "\n";
    out << "cells " << mesh.n_cells() << "\n";
    out << "iteration " << iteration << "\n";
    out << "time " << time << "\n";
    out << "identity " << identity.n_points << ' ' << identity.n_faces << ' '
        << identity.n_cells << ' ' << identity.geometry_checksum << "\n";
    out << "field U 3\n";
    for (std::size_t c=0;c<mesh.n_cells();++c)
        out << U.component_data(0)[c] << ' ' << U.component_data(1)[c] << ' ' << U.component_data(2)[c] << "\n";
    out << "field p 1\n";
    for (std::size_t c=0;c<mesh.n_cells();++c) out << p(c) << "\n";
    out << "optional_fields\n";
    const auto write_field = [&](const char* name, const cfdx::core::Field<double,cfdx::core::Location::CELL>* field) {
        if (!field) return;
        out << "field " << name << " 1\n";
        for (std::size_t c=0;c<mesh.n_cells();++c) out << (*field)(c) << "\n";
    };
    write_field("T", fields.temperature);
    write_field("k", fields.k);
    write_field("second_turbulence", fields.second_turbulence);
    if (!out) {
        out.close();
        std::remove(tmp_path.c_str());
        throw std::runtime_error("write_dat_restart_fields: write failed for " + path);
    }
    out.close();
    if (std::rename(tmp_path.c_str(), path.c_str()) != 0) {
        std::remove(tmp_path.c_str());
        throw std::runtime_error("write_dat_restart_fields: cannot replace " + path);
    }
}
inline DatRestartState read_dat_restart_fields(
    const std::string& path,
    const cfdx::core::Mesh& mesh,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    DatRestartFields fields)
{
    if (U.dimension() != 3 || U.size() != mesh.n_cells())
        throw std::invalid_argument("read_dat_restart_fields: invalid velocity field");
    if (p.dimension() != 1 || p.size() != mesh.n_cells())
        throw std::invalid_argument("read_dat_restart_fields: invalid pressure field");
    validate_optional_restart_field(mesh, fields.temperature, "temperature");
    validate_optional_restart_field(mesh, fields.k, "k");
    validate_optional_restart_field(mesh, fields.second_turbulence, "second turbulence");
    std::ifstream in(path);
    if (!in) throw std::runtime_error("read_dat_restart_fields: cannot open " + path);

    std::string magic; int version = 0;
    if (!(in >> magic >> version) || magic != "CFDX-DAT" ||
        (version != 1 && version != DatRestartVersionFields &&
         version != DatRestartVersionIdentity))
        throw std::runtime_error("read_dat_restart_fields: unsupported DAT format");
    std::size_t cells = 0; std::string key;
    if (!(in >> key >> cells) || key != "cells" || cells != mesh.n_cells())
        throw std::runtime_error("read_dat_restart_fields: mesh cell count mismatch");
    DatRestartState state;
    if (!(in >> key >> state.iteration) || key != "iteration")
        throw std::runtime_error("read_dat_restart_fields: invalid iteration");
    if (!(in >> key >> state.time) || key != "time" || !std::isfinite(state.time))
        throw std::runtime_error("read_dat_restart_fields: invalid time");
    validate_dat_restart_identity(mesh, version, in);

    auto read_scalar_field = [&](const char* expected,
                                 cfdx::core::Field<double, cfdx::core::Location::CELL>* field) {
        if (!(in >> key)) throw std::runtime_error("read_dat_restart_fields: missing field");
        std::size_t dimension = 0;
        if (!(in >> key >> dimension) || key != expected || dimension != 1)
            throw std::runtime_error(std::string("read_dat_restart_fields: invalid ") + expected + " field");
        if (!field) {
            for (std::size_t c = 0; c < mesh.n_cells(); ++c) { double ignored; if (!(in >> ignored)) throw std::runtime_error("read_dat_restart_fields: invalid field values"); }
            return;
        }
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            double value = 0.0;
            if (!(in >> value) || !std::isfinite(value)) throw std::runtime_error("read_dat_restart_fields: invalid field value");
            (*field)(c) = value;
        }
    };

    if (!(in >> key) || key != "field") throw std::runtime_error("read_dat_restart_fields: missing U field");
    std::size_t dimension = 0;
    if (!(in >> key >> dimension) || key != "U" || dimension != 3)
        throw std::runtime_error("read_dat_restart_fields: invalid U field");
    for (std::size_t c=0;c<mesh.n_cells();++c) {
        double ux,uy,uz;
        if (!(in>>ux>>uy>>uz) || !std::isfinite(ux)||!std::isfinite(uy)||!std::isfinite(uz))
            throw std::runtime_error("read_dat_restart_fields: invalid U values");
        U.component_data(0)[c]=ux; U.component_data(1)[c]=uy; U.component_data(2)[c]=uz;
    }
    if (!(in >> key) || key != "field") throw std::runtime_error("read_dat_restart_fields: missing p field");
    if (!(in >> key >> dimension) || key != "p" || dimension != 1)
        throw std::runtime_error("read_dat_restart_fields: invalid p field");
    for (std::size_t c=0;c<mesh.n_cells();++c) {
        double value;
        if (!(in>>value) || !std::isfinite(value)) throw std::runtime_error("read_dat_restart_fields: invalid p values");
        p(c)=value;
    }

    // Every version that survives the identity gate is >= 3 and always carries
    // the optional-fields section, so the former v1/v2 tails are unreachable.
    if (!(in >> key) || key != "optional_fields")
        throw std::runtime_error("read_dat_restart_fields: missing optional fields section");
    bool got_t=false, got_k=false, got_second=false;
    while (in >> key) {
        if (key != "field") throw std::runtime_error("read_dat_restart_fields: invalid optional section");
        std::string name; if (!(in >> name >> dimension) || dimension != 1)
            throw std::runtime_error("read_dat_restart_fields: invalid optional field header");
        auto* target = name == "T" ? fields.temperature :
                       name == "k" ? fields.k :
                       name == "second_turbulence" ? fields.second_turbulence : nullptr;
        if (name == "T") got_t=true;
        else if (name == "k") got_k=true;
        else if (name == "second_turbulence") got_second=true;
        else throw std::runtime_error("read_dat_restart_fields: unknown optional field " + name);
        if (!target) {
            for (std::size_t i=0;i<mesh.n_cells();++i) { double ignored; if (!(in>>ignored)) throw std::runtime_error("read_dat_restart_fields: invalid optional values"); }
        } else {
            for (std::size_t i=0;i<mesh.n_cells();++i) {
                double value; if (!(in>>value) || !std::isfinite(value)) throw std::runtime_error("read_dat_restart_fields: invalid optional value");
                (*target)(i)=value;
            }
        }
    }
    if (fields.temperature && !got_t) throw std::runtime_error("read_dat_restart_fields: requested T but DAT has none");
    if (fields.k && !got_k) throw std::runtime_error("read_dat_restart_fields: requested k but DAT has none");
    if (fields.second_turbulence && !got_second) throw std::runtime_error("read_dat_restart_fields: requested turbulence field but DAT has none");
    return state;
}



// Velocity/pressure checkpoint. This is the same on-disk format as the
// `_fields` variant with no optional fields requested, so it delegates instead
// of maintaining a second writer. Two independent writers previously produced
// two different layouts under the same version number, which made a checkpoint
// written by one unreadable by the other.
inline void write_dat_restart(
    const std::string& path,
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    std::size_t iteration = 0,
    double time = 0.0)
{
    write_dat_restart_fields(path, mesh, U, p, DatRestartFields{}, iteration, time);
}

// Velocity/pressure restart, discarding any optional fields present in the
// checkpoint. Requesting an optional field that the caller does not ask for is
// how a reader declines state it cannot consume.
inline DatRestartState read_dat_restart(
    const std::string& path,
    const cfdx::core::Mesh& mesh,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p)
{
    return read_dat_restart_fields(path, mesh, U, p, DatRestartFields{});
}

} // namespace cfdx::io
