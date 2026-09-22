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
    out << "CFDX-DAT 1\\n";
    out << "cells " << mesh.n_cells() << "\\n";
    out << "iteration " << iteration << "\\n";
    out << "time " << time << "\\n";
    out << "field U 3\\n";
    for (std::size_t c = 0; c < mesh.n_cells(); ++c)
        out << U.component_data(0)[c] << ' '
            << U.component_data(1)[c] << ' '
            << U.component_data(2)[c] << "\\n";
    out << "field p 1\\n";
    for (std::size_t c = 0; c < mesh.n_cells(); ++c)
        out << p(c) << "\\n";
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
