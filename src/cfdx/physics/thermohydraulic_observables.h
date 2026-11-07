#pragma once

// Raw, solver-owned evidence for thermohydraulic benchmarks. This does not
// perform reference comparisons or claim qualification. All velocities and
// pressures come from the accepted CFDX state; mdot is the authoritative
// Rhie--Chow face flux, never reconstructed from interpolated velocities.

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/sst_production_coupling.h"
#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct ThermohydraulicObservableSummary {
    std::size_t cells = 0;
    std::size_t faces = 0;
    std::size_t boundary_faces = 0;
    std::size_t no_slip_wall_faces = 0;
    double boundary_net_mass_flux_kg_s = 0.0;
    double boundary_absolute_mass_flux_kg_s = 0.0;
    double global_mass_imbalance_relative = 0.0;
    double maximum_cell_mass_imbalance_relative = 0.0;
};

inline void thermohydraulic_check_finite(double value, const char* label) {
    if (!std::isfinite(value))
        throw std::invalid_argument(std::string("thermohydraulic export: non-finite ") + label);
}

inline std::string thermohydraulic_csv_string(const std::string& value) {
    std::string escaped = "\"";
    for (const char c : value) {
        if (c == '"') escaped += "\"\"";
        else if (c == '\n' || c == '\r')
            throw std::invalid_argument("thermohydraulic export: invalid newline in patch name");
        else escaped += c;
    }
    return escaped + '"';
}

// Install a complete CSV only after the stream has been checked. A benchmark
// cannot consume half-written output after an interrupt or I/O failure.
template <typename Writer>
inline void thermohydraulic_atomic_csv(
    const std::filesystem::path& file, Writer writer)
{
    const auto temporary = file.string() + ".tmp";
    try {
        std::ofstream stream(temporary, std::ios::trunc);
        if (!stream)
            throw std::runtime_error("thermohydraulic export: cannot open " + temporary);
        stream << std::setprecision(std::numeric_limits<double>::max_digits10);
        writer(stream);
        stream.close();
        if (!stream)
            throw std::runtime_error("thermohydraulic export: incomplete write to " + temporary);
        std::filesystem::rename(temporary, file);
    } catch (...) {
        std::error_code ec;
        std::filesystem::remove(temporary, ec);
        throw;
    }
}

inline ThermohydraulicObservableSummary write_thermohydraulic_observables(
    const std::filesystem::path& directory,
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    const VelocityBoundaryConditions& velocity_boundary_conditions,
    double density,
    double molecular_kinematic_viscosity,
    const SSTProductionState* sst = nullptr,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* temperature = nullptr)
{
    using cfdx::core::PatchType;
    const auto cells = mesh.n_cells();
    const auto faces = mesh.n_faces();
    if (cells == 0 || faces == 0 || velocity.size() != cells ||
        velocity.dimension() != 3 || pressure.size() != cells ||
        pressure.dimension() != 1 || mass_flux.size() != faces ||
        mass_flux.dimension() != 1 ||
        !(density > 0.0) || !std::isfinite(density) ||
        !(molecular_kinematic_viscosity > 0.0) ||
        !std::isfinite(molecular_kinematic_viscosity) ||
        (temperature && (temperature->size() != cells || temperature->dimension() != 1)))
        throw std::invalid_argument("thermohydraulic export: invalid state, shape or material");
    if (sst) validate_sst_production_state(*sst, cells);

    const auto geometry = build_fv_geometry(mesh);
    std::vector<double> per_cell_net(cells, 0.0);
    std::vector<double> per_cell_abs(cells, 0.0);
    ThermohydraulicObservableSummary summary;
    summary.cells = cells;
    summary.faces = faces;

    for (std::size_t f = 0; f < faces; ++f) {
        const double flux = mass_flux(f);
        thermohydraulic_check_finite(flux, "mass flux");
        const auto owner = static_cast<std::size_t>(mesh.ownership().owner(f));
        if (owner >= cells)
            throw std::invalid_argument("thermohydraulic export: invalid face owner");
        per_cell_net[owner] += flux;
        per_cell_abs[owner] += std::abs(flux);
        const auto neighbour = mesh.ownership().neighbour(f);
        if (neighbour >= 0) {
            if (static_cast<std::size_t>(neighbour) >= cells)
                throw std::invalid_argument("thermohydraulic export: invalid neighbour");
            per_cell_net[static_cast<std::size_t>(neighbour)] -= flux;
            per_cell_abs[static_cast<std::size_t>(neighbour)] += std::abs(flux);
        } else {
            ++summary.boundary_faces;
            summary.boundary_net_mass_flux_kg_s += flux;
            summary.boundary_absolute_mass_flux_kg_s += std::abs(flux);
        }
    }
    summary.global_mass_imbalance_relative =
        summary.boundary_absolute_mass_flux_kg_s > 0.0
            ? std::abs(summary.boundary_net_mass_flux_kg_s) /
                (0.5 * summary.boundary_absolute_mass_flux_kg_s)
            : 0.0;
    for (std::size_t c = 0; c < cells; ++c) {
        if (per_cell_abs[c] > 0.0)
            summary.maximum_cell_mass_imbalance_relative = std::max(
                summary.maximum_cell_mass_imbalance_relative,
                std::abs(per_cell_net[c]) / per_cell_abs[c]);
    }
    std::filesystem::create_directories(directory);

    thermohydraulic_atomic_csv(directory / "cell_fields.csv", [&](std::ostream& out) {
        out << "cell_id,x_m,y_m,z_m,volume_m3,u_x_m_s,u_y_m_s,u_z_m_s,"
               "p_solver_Pa,p_static_Pa,k_m2_s2,omega_1_s,nu_t_m2_s,"
               "temperature_K,wall_distance_m,cell_net_mass_flux_kg_s\n";
        for (std::size_t c = 0; c < cells; ++c) {
            for (std::size_t d = 0; d < 3; ++d)
                thermohydraulic_check_finite(velocity(c, d), "velocity");
            thermohydraulic_check_finite(pressure(c), "pressure");
            const auto& x = geometry.cell_centres[c];
            out << c << ',' << x.x << ',' << x.y << ',' << x.z << ','
                << geometry.cell_volumes[c] << ',' << velocity(c, 0) << ','
                << velocity(c, 1) << ',' << velocity(c, 2) << ','
                << pressure(c) << ','
                // The SST carrier solves for p*=p_static + 2 rho k/3.
                << (pressure(c) - (sst ? (2.0/3.0)*density*sst->k(c) : 0.0)) << ',';
            if (sst) out << sst->k(c) << ',' << sst->omega(c) << ','
                         << sst->nu_t(c);
            else out << ",,";
            out << ',';
            if (temperature) {
                thermohydraulic_check_finite((*temperature)(c), "temperature");
                out << (*temperature)(c);
            }
            out << ',';
            if (sst) out << sst->wall_distance(c);
            out << ',' << per_cell_net[c] << '\n';
        }
    });

    thermohydraulic_atomic_csv(directory / "face_fluxes.csv", [&](std::ostream& out) {
        out << "face_id,owner,neighbour,patch,x_m,y_m,z_m,"
               "area_x_m2,area_y_m2,area_z_m2,area_m2,"
               "mass_flux_kg_s,wall_normal_distance_m,"
               "wall_shear_x_Pa,wall_shear_y_Pa,wall_shear_z_Pa,"
               "wall_shear_method,wall_y_plus_estimate,owner_nu_t_over_nu\n";
        for (std::size_t f = 0; f < faces; ++f) {
            const auto owner = static_cast<std::size_t>(mesh.ownership().owner(f));
            const auto neighbour = mesh.ownership().neighbour(f);
            const auto& x = geometry.face_centres[f];
            const auto& S = geometry.face_area_vectors[f];
            const double area = S.mag();
            if (!(area > 0.0) || !std::isfinite(area))
                throw std::invalid_argument("thermohydraulic export: invalid face area");
            const std::size_t patch_id = geometry.face_patch[f];
            const bool boundary = neighbour < 0;
            if (boundary && patch_id >= mesh.boundary().n_patches())
                throw std::invalid_argument("thermohydraulic export: unpatched boundary face");
            const std::string patch_name = boundary ?
                mesh.boundary().patch(patch_id).name : "";
            out << f << ',' << owner << ',' << neighbour << ','
                << thermohydraulic_csv_string(patch_name) << ','
                << x.x << ',' << x.y << ',' << x.z << ','
                << S.x << ',' << S.y << ',' << S.z << ',' << area
                << ',' << mass_flux(f) << ',';
            if (boundary && mesh.boundary().patch(patch_id).type == PatchType::WALL) {
                const auto it = velocity_boundary_conditions.find(patch_name);
                if (it != velocity_boundary_conditions.end() &&
                    it->second.type == VelocityBoundaryCondition::Type::FIXED_VALUE) {
                    // First-order, one-sided shear. This is explicitly a
                    // diagnostic: mesh convergence and higher-order wall
                    // reconstruction are still required for Cf qualification.
                    const double nx = S.x / area, ny = S.y / area, nz = S.z / area;
                    const auto& centre = geometry.cell_centres[owner];
                    const double distance = -((centre.x-x.x)*nx +
                                              (centre.y-x.y)*ny +
                                              (centre.z-x.z)*nz);
                    if (!(distance > 0.0) || !std::isfinite(distance))
                        throw std::invalid_argument(
                            "thermohydraulic export: invalid wall-normal distance");
                    const auto& wall = it->second.value;
                    const double ux = velocity(owner,0) - wall.x;
                    const double uy = velocity(owner,1) - wall.y;
                    const double uz = velocity(owner,2) - wall.z;
                    const double dot = ux*nx + uy*ny + uz*nz;
                    const double multiplier = density*molecular_kinematic_viscosity/distance;
                    out << distance << ',' << multiplier*(ux - dot*nx)
                        << ',' << multiplier*(uy - dot*ny)
                        << ',' << multiplier*(uz - dot*nz)
                        << ",first_order_one_sided_diagnostic";
                    // y+ is a *diagnostic estimate* from the one-sided
                    // molecular shear, not a resolved-wall qualification.
                    const double tau_tangent = multiplier *
                        std::sqrt((ux - dot*nx)*(ux - dot*nx) +
                                  (uy - dot*ny)*(uy - dot*ny) +
                                  (uz - dot*nz)*(uz - dot*nz));
                    const double u_tau = std::sqrt(tau_tangent / density);
                    out << ',' << distance * u_tau / molecular_kinematic_viscosity;
                    out << ',';
                    if (sst)
                        out << sst->nu_t(owner)/molecular_kinematic_viscosity;
                    ++summary.no_slip_wall_faces;
                } else {
                    out << ",,,,not_available_non_dirichlet_wall,,";
                }
            } else {
                out << ",,,,not_a_no_slip_wall,,";
            }
            out << '\n';
        }
    });

    thermohydraulic_atomic_csv(directory / "mass_conservation.json", [&](std::ostream& out) {
        out << "{\n  \"format\": \"CFDX-THERMOHYDRAULIC-CONSERVATION\",\n"
               "  \"schema_version\": 1,\n"
               "  \"physical_qualification\": \"NOT_QUALIFIED\",\n"
            << "  \"cells\": " << summary.cells << ",\n"
            << "  \"faces\": " << summary.faces << ",\n"
            << "  \"boundary_faces\": " << summary.boundary_faces << ",\n"
            << "  \"no_slip_wall_faces\": " << summary.no_slip_wall_faces << ",\n"
            << "  \"boundary_net_mass_flux_kg_s\": "
            << summary.boundary_net_mass_flux_kg_s << ",\n"
            << "  \"boundary_absolute_mass_flux_kg_s\": "
            << summary.boundary_absolute_mass_flux_kg_s << ",\n"
            << "  \"global_mass_imbalance_relative\": "
            << summary.global_mass_imbalance_relative << ",\n"
            << "  \"maximum_cell_mass_imbalance_relative\": "
            << summary.maximum_cell_mass_imbalance_relative << "\n}\n";
    });
    return summary;
}

}  // namespace cfdx::physics
