#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"

#include <cmath>
#include <cstddef>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

namespace cfdx::io {

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

    std::ofstream out(path);
    if (!out) throw std::runtime_error("write_dat_restart_fields: cannot open " + path);
    out.precision(std::numeric_limits<double>::max_digits10);
    out << "CFDX-DAT 2\n";
    out << "cells " << mesh.n_cells() << "\n";
    out << "iteration " << iteration << "\n";
    out << "time " << time << "\n";
    out << "field U 3\n";
    for (std::size_t c=0;c<mesh.n_cells();++c)
        out << U.component_data(0)[c] << ' ' << U.component_data(1)[c] << ' ' << U.component_data(2)[c] << "\n";
    out << "field p 1\n";
    for (std::size_t c=0;c<mesh.n_cells();++c) out << p(c) << "\n";
    out << "optional_fields\n";
    const auto write_field = [&](const char* name, const cfdx::core::Field<double,cfdx::core::Location::CELL>* field) {
        if (!field) return;
        out << "field " << name << " 1\n";
        for (std::size_t c=0;c<mesh.n_cells();++c) {
            if (!std::isfinite((*field)(c)))
                throw std::invalid_argument(std::string("write_dat_restart_fields: non-finite ") + name);
            out << (*field)(c) << "\n";
        }
    };
    write_field("T", fields.temperature);
    write_field("k", fields.k);
    write_field("second_turbulence", fields.second_turbulence);
}
inline DatRestartState read_dat_restart_fields(
    const std::string& path,
    const cfdx::core::Mesh& mesh,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    DatRestartFields fields)
{
    validate_optional_restart_field(mesh, fields.temperature, "temperature");
    validate_optional_restart_field(mesh, fields.k, "k");
    validate_optional_restart_field(mesh, fields.second_turbulence, "second turbulence");
    std::ifstream in(path);
    if (!in) throw std::runtime_error("read_dat_restart_fields: cannot open " + path);

    std::string magic; int version = 0;
    if (!(in >> magic >> version) || magic != "CFDX-DAT" || (version != 1 && version != 2))
        throw std::runtime_error("read_dat_restart_fields: unsupported DAT format");
    std::size_t cells = 0; std::string key;
    if (!(in >> key >> cells) || key != "cells" || cells != mesh.n_cells())
        throw std::runtime_error("read_dat_restart_fields: mesh cell count mismatch");
    DatRestartState state;
    if (!(in >> key >> state.iteration) || key != "iteration")
        throw std::runtime_error("read_dat_restart_fields: invalid iteration");
    if (!(in >> key >> state.time) || key != "time" || !std::isfinite(state.time))
        throw std::runtime_error("read_dat_restart_fields: invalid time");

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

    if (version == 1) {
        if (in >> key) throw std::runtime_error("read_dat_restart_fields: unexpected trailing data");
        return state;
    }
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


inline void write_dat_restart(
    const std::string& path,
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    std::size_t iteration = 0,
    double time = 0.0)
{
    using namespace cfdx::core;
    if (U.dimension() != 3 || U.size() != mesh.n_cells())
        throw std::invalid_argument("write_dat_restart: invalid velocity field");
    if (p.dimension() != 1 || p.size() != mesh.n_cells())
        throw std::invalid_argument("write_dat_restart: invalid pressure field");
    if (!std::isfinite(time))
        throw std::invalid_argument("write_dat_restart: non-finite time");

    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("write_dat_restart: cannot open " + path);

    out.precision(std::numeric_limits<double>::max_digits10);
    out << "CFDX-DAT 1\n";
    out << "cells " << mesh.n_cells() << "\n";
    out << "iteration " << iteration << "\n";
    out << "time " << time << "\n";
    out << "field U 3\n";
    for (std::size_t c = 0; c < mesh.n_cells(); ++c)
        out << U.component_data(0)[c] << ' '
            << U.component_data(1)[c] << ' '
            << U.component_data(2)[c] << "\n";
    out << "field p 1\n";
    for (std::size_t c = 0; c < mesh.n_cells(); ++c)
        out << p(c) << "\n";
    if (!out)
        throw std::runtime_error("write_dat_restart: write failed for " + path);
}

inline DatRestartState read_dat_restart(
    const std::string& path,
    const cfdx::core::Mesh& mesh,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p)
{
    using namespace cfdx::core;
    if (U.dimension() != 3 || U.size() != mesh.n_cells())
        throw std::invalid_argument("read_dat_restart: invalid velocity field");
    if (p.dimension() != 1 || p.size() != mesh.n_cells())
        throw std::invalid_argument("read_dat_restart: invalid pressure field");

    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("read_dat_restart: cannot open " + path);

    auto expect = [&in](const char* key) {
        std::string actual;
        if (!(in >> actual) || actual != key)
            throw std::runtime_error(
                std::string("read_dat_restart: expected '") + key + "'");
    };

    std::string magic;
    int version = 0;
    if (!(in >> magic >> version) || magic != "CFDX-DAT" || version != 1)
        throw std::runtime_error("read_dat_restart: unsupported DAT format");

    std::size_t cells = 0;
    expect("cells");
    if (!(in >> cells) || cells != mesh.n_cells())
        throw std::runtime_error("read_dat_restart: mesh cell count mismatch");

    DatRestartState state;
    expect("iteration");
    if (!(in >> state.iteration))
        throw std::runtime_error("read_dat_restart: invalid iteration");

    expect("time");
    if (!(in >> state.time) || !std::isfinite(state.time))
        throw std::runtime_error("read_dat_restart: invalid time");

    std::string field_name;
    std::size_t dimension = 0;
    expect("field");
    if (!(in >> field_name >> dimension) || field_name != "U" || dimension != 3)
        throw std::runtime_error("read_dat_restart: invalid U field header");

    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        double ux = 0.0, uy = 0.0, uz = 0.0;
        if (!(in >> ux >> uy >> uz) ||
            !std::isfinite(ux) || !std::isfinite(uy) || !std::isfinite(uz))
            throw std::runtime_error("read_dat_restart: invalid U values");
        U.component_data(0)[c] = ux;
        U.component_data(1)[c] = uy;
        U.component_data(2)[c] = uz;
    }

    expect("field");
    if (!(in >> field_name >> dimension) || field_name != "p" || dimension != 1)
        throw std::runtime_error("read_dat_restart: invalid p field header");

    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        double value = 0.0;
        if (!(in >> value) || !std::isfinite(value))
            throw std::runtime_error("read_dat_restart: invalid p values");
        p(c) = value;
    }

    std::string trailing;
    if (in >> trailing)
        throw std::runtime_error("read_dat_restart: unexpected trailing data");

    return state;
}

} // namespace cfdx::io
