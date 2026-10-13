// cfdx_thermal_solver: production solid-conduction thermal path.
//
// Stage 1 of the multiphysics solver-level orchestration capability
// (issue #885); qualification case T1 (issue #880). It executes a steady
// multi-material conduction solve with volumetric power sources, contact
// resistances at material interfaces and WALL_THERMAL boundary conditions,
// then persists machine-readable artifacts:
//   execution.json           (CFDX-EXECUTION v1, same process semantics as
//                             the production solver: 0/1/2)
//   thermal_convergence.json (linear solve status and residual)
//   thermal_summary.json     (temperature range, per-zone stats, energy
//                             balance, per-patch boundary flux)
//   result_0.vtu             (temperature field)
//   probe CSV                (--probe samples, nearest cell centre)
//
// Case conventions (documented, fail-closed):
//   * Only CFDX case files (.cfdx.h5) are accepted: the thermal path needs
//     the declared case setup. Mesh-only files are rejected.
//   * case_setup.energy_model must be "conjugate" or "solid_conduction".
//   * One material = one thermal zone. Each material must declare a positive
//     thermal_conductivity. With more than one material, every material
//     must carry extra_properties.zone_id, and the case file must carry a
//     per-cell scalar field "material_id" (or "zone_id") mapping cells to
//     zones. A single-material case may omit both.
//   * Optional per-material extra_properties:
//       volumetric_power    W/m^3 (uniform within the zone, default 0)
//       contact_resistance m^2.K/W at interfaces with other zones
//                           (default 0; both sides must agree if declared)
//   * Boundary conditions: a patch with value_type WALL_THERMAL and a
//     non-zero wall_temperature is a fixed-temperature wall. Otherwise a
//     raw_params entry "heat_flux" (W/m^2, positive into the domain) turns
//     it into a fixed-gradient wall (the patch must touch a single zone).
//     Every other patch is adiabatic. At least one fixed-temperature patch
//     is required, otherwise the conduction system is singular.
//   * Initial temperature: initial_condition.temperature when declared,
//     otherwise ref_temperature. The linear solve result does not depend
//     on this initial guess.
#include "cfdx/physics/thermal_zones.h"
#include "cfdx/io/hdf5/case_hdf5_io.h"
#include "cfdx/io/runtime/execution_summary.h"
#include "cfdx/io/vtu/vtu_writer.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::io;
using namespace cfdx::physics;

namespace {

struct ThermalProbe {
    std::string name;
    Vec3 location;
};

struct Options {
    std::string mesh;
    std::filesystem::path output_dir;
    std::vector<ThermalProbe> probes;
    std::filesystem::path probe_csv;
};

ThermalProbe parse_probe(const std::string& spec)
{
    const auto first = spec.find(':');
    if (first == std::string::npos)
        throw std::invalid_argument("--probe expects <name>:<x>,<y>,<z>:T");
    const auto second = spec.find(':', first + 1);
    if (second == std::string::npos)
        throw std::invalid_argument("--probe expects <name>:<x>,<y>,<z>:T");
    const std::string name = spec.substr(0, first);
    if (name.empty())
        throw std::invalid_argument("--probe name must not be empty");

    std::vector<double> coordinates;
    const std::string point = spec.substr(first + 1, second - first - 1);
    std::size_t begin = 0;
    while (begin <= point.size()) {
        const auto comma = point.find(',', begin);
        const std::string token = point.substr(
            begin, comma == std::string::npos ? std::string::npos : comma - begin);
        try {
            coordinates.push_back(std::stod(token));
        } catch (const std::exception&) {
            throw std::invalid_argument("--probe coordinates must be numeric: " + token);
        }
        if (comma == std::string::npos) break;
        begin = comma + 1;
    }
    if (coordinates.size() != 3)
        throw std::invalid_argument("--probe expects exactly three coordinates");

    const std::string field_name = spec.substr(second + 1);
    if (field_name != "T" && field_name != "temperature")
        throw std::invalid_argument(
            "--probe field must be T (temperature): " + field_name);

    ThermalProbe probe;
    probe.name = name;
    probe.location = Vec3{coordinates[0], coordinates[1], coordinates[2]};
    return probe;
}

Options parse(int argc, char** argv)
{
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto value = [&](const char* name) -> std::string {
            if (i + 1 >= argc)
                throw std::invalid_argument(std::string("missing value for ") + name);
            return argv[++i];
        };
        if (arg == "--mesh") o.mesh = value("--mesh");
        else if (arg == "--output-dir") o.output_dir = value("--output-dir");
        else if (arg == "--probe") o.probes.push_back(parse_probe(value("--probe")));
        else if (arg == "--probe-csv") o.probe_csv = value("--probe-csv");
        else if (arg == "--help") {
            std::cout << "cfdx_thermal_solver --mesh CASE.cfdx.h5 --output-dir DIR\n"
                         "  [--probe <name>:<x>,<y>,<z>:T]... [--probe-csv PATH]\n"
                         "  Steady multi-material solid conduction with volumetric\n"
                         "  power, contact resistance and WALL_THERMAL boundaries.\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + arg);
        }
    }
    if (o.mesh.empty() || o.output_dir.empty())
        throw std::invalid_argument("--mesh and --output-dir are required");
    if (o.mesh.size() < 8 || o.mesh.substr(o.mesh.size() - 8) != ".cfdx.h5")
        throw std::invalid_argument(
            "thermal solver requires a CFDX case file (.cfdx.h5): the thermal "
            "path consumes the declared case setup, not a bare mesh");
    if (!o.probe_csv.empty() && o.probes.empty())
        throw std::invalid_argument("--probe-csv requires at least one --probe");
    for (std::size_t i = 0; i < o.probes.size(); ++i)
        for (std::size_t j = i + 1; j < o.probes.size(); ++j)
            if (o.probes[i].name == o.probes[j].name)
                throw std::invalid_argument("duplicate probe name: " + o.probes[i].name);
    std::filesystem::create_directories(o.output_dir);
    return o;
}

double parse_finite_number(const std::string& text, const std::string& context)
{
    double value = 0.0;
    try {
        std::size_t consumed = 0;
        value = std::stod(text, &consumed);
        if (consumed != text.size())
            throw std::invalid_argument("trailing characters");
    } catch (const std::exception&) {
        throw std::invalid_argument(
            context + " must be a finite number, got '" + text + "'");
    }
    if (!std::isfinite(value))
        throw std::invalid_argument(context + " must be a finite number");
    return value;
}

int parse_integer(const std::string& text, const std::string& context)
{
    const double value = parse_finite_number(text, context);
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1e-9 || std::abs(rounded) > 2.0e9)
        throw std::invalid_argument(context + " must be an integer, got '" + text + "'");
    return static_cast<int>(rounded);
}

// One declared material = one thermal zone.
std::vector<ThermalZoneSpec> build_zone_specs(const CaseSetup& setup)
{
    if (setup.materials.empty())
        throw std::invalid_argument(
            "thermal solver: case declares no material; a material with a "
            "positive thermal_conductivity is required");
    std::vector<ThermalZoneSpec> specs;
    specs.reserve(setup.materials.size());
    for (const auto& material : setup.materials) {
        ThermalZoneSpec spec;
        spec.name = material.name;
        spec.solid = true;
        spec.conductivity = material.thermal_conductivity;
        if (!(spec.conductivity > 0.0) || !std::isfinite(spec.conductivity))
            throw std::invalid_argument(
                "thermal solver: material '" + material.name +
                "' must declare a positive finite thermal_conductivity");

        const auto zone_it = material.extra_properties.find("zone_id");
        if (zone_it != material.extra_properties.end())
            spec.zone_id = parse_integer(
                zone_it->second,
                "material '" + material.name + "' extra_properties.zone_id");
        else if (setup.materials.size() > 1)
            throw std::invalid_argument(
                "thermal solver: material '" + material.name +
                "' must declare extra_properties.zone_id when the case "
                "declares several materials");
        else
            spec.zone_id = 0;

        const auto power_it = material.extra_properties.find("volumetric_power");
        if (power_it != material.extra_properties.end())
            spec.volumetric_power = parse_finite_number(
                power_it->second,
                "material '" + material.name + "' extra_properties.volumetric_power");

        const auto resistance_it =
            material.extra_properties.find("contact_resistance");
        if (resistance_it != material.extra_properties.end()) {
            spec.contact_resistance = parse_finite_number(
                resistance_it->second,
                "material '" + material.name + "' extra_properties.contact_resistance");
            if (spec.contact_resistance < 0.0)
                throw std::invalid_argument(
                    "material '" + material.name +
                    "' extra_properties.contact_resistance must be non-negative");
        }
        specs.push_back(spec);
    }
    return specs;
}

// Cell -> zone map from the case file's scalar fields.
std::vector<int> resolve_cell_zones(
    const Mesh& mesh,
    const std::string& case_path,
    const std::vector<ThermalZoneSpec>& specs)
{
    const std::size_t nc = mesh.n_cells();
    std::vector<std::pair<std::string, ScalarCellField>> fields;
    if (!read_scalar_fields_hdf5(case_path, fields))
        fields.clear();

    const ScalarCellField* zone_field = nullptr;
    std::string zone_field_name;
    for (const auto& entry : fields) {
        const std::string& name = entry.first;
        const ScalarCellField& field = entry.second;
        if (name == "material_id" || name == "zone_id") {
            if (zone_field)
                throw std::invalid_argument(
                    "thermal solver: case carries both material_id and "
                    "zone_id cell fields; declare exactly one");
            zone_field = &field;
            zone_field_name = name;
        }
    }

    if (!zone_field) {
        if (specs.size() == 1)
            return std::vector<int>(nc, specs.front().zone_id);
        throw std::invalid_argument(
            "thermal solver: case carries no material_id/zone_id cell field; "
            "a multi-material case must map cells onto zones");
    }
    if (zone_field->size() != nc)
        throw std::invalid_argument(
            "thermal solver: cell field '" + zone_field_name +
            "' does not match the mesh cell count");

    std::vector<int> cell_zone(nc, 0);
    for (std::size_t c = 0; c < nc; ++c) {
        const double value = (*zone_field)(c);
        const double rounded = std::round(value);
        if (std::abs(value - rounded) > 1e-6 || std::abs(rounded) > 2.0e9)
            throw std::invalid_argument(
                "thermal solver: cell field '" + zone_field_name +
                "' must contain integer zone ids (cell " + std::to_string(c) +
                " has " + std::to_string(value) + ")");
        cell_zone[c] = static_cast<int>(rounded);
    }
    return cell_zone;
}

// Thermal boundary conditions from the declared case (fail-closed).
ScalarBoundaryConditions build_thermal_boundary_conditions(
    const Mesh& mesh,
    const FvGeometry& geometry,
    const CaseSetup& setup,
    const ThermalZoneModel& zones)
{
    ScalarBoundaryConditions bcs;
    std::size_t fixed_temperature_patches = 0;
    for (std::size_t i = 0; i < mesh.boundary().n_patches(); ++i) {
        const auto& patch = mesh.boundary().patch(i);
        const auto* bc = setup.find_boundary(patch.name);
        if (bc == nullptr)
            throw std::invalid_argument(
                "thermal solver: case has no boundary condition for mesh patch '" +
                patch.name + "'");

        ScalarBoundaryCondition thermal;  // default: adiabatic
        if (bc->value_type == BCValueType::WALL_THERMAL &&
            bc->wall_temperature != 0.0) {
            thermal = {ScalarBoundaryType::FIXED_VALUE, bc->wall_temperature, 0.0};
            ++fixed_temperature_patches;
        } else {
            const auto flux_it = bc->raw_params.find("heat_flux");
            if (flux_it != bc->raw_params.end()) {
                const double heat_flux = parse_finite_number(
                    flux_it->second, "patch '" + patch.name + "' heat_flux");
                // The fixed gradient g = q''/k depends on the local zone
                // conductivity, so a heat-flux patch must touch exactly one
                // declared zone.
                bool has_face = false;
                bool mono_zone = true;
                int patch_zone = 0;
                for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                    if (geometry.face_patch[f] != i) continue;
                    if (mesh.ownership().neighbour(f) >= 0) continue;
                    const std::size_t owner = mesh.ownership().owner(f);
                    if (!has_face) {
                        patch_zone = zones.cell_zone[owner];
                        has_face = true;
                    } else if (zones.cell_zone[owner] != patch_zone) {
                        mono_zone = false;
                    }
                }
                if (!has_face || !mono_zone)
                    throw std::invalid_argument(
                        "thermal solver: heat_flux on patch '" + patch.name +
                        "' requires boundary faces from a single material zone");
                double zone_conductivity = 0.0;
                bool found = false;
                for (const auto& spec : zones.specs) {
                    if (spec.zone_id == patch_zone) {
                        zone_conductivity = spec.conductivity;
                        found = true;
                        break;
                    }
                }
                if (!found)
                    throw std::invalid_argument(
                        "thermal solver: patch '" + patch.name +
                        "' references undeclared zone " +
                        std::to_string(patch_zone));
                // Positive gradient = heat flowing into the domain.
                thermal = {ScalarBoundaryType::FIXED_GRADIENT, 0.0,
                           heat_flux / zone_conductivity};
            }
        }
        bcs[patch.name] = thermal;
    }
    if (fixed_temperature_patches == 0)
        throw std::invalid_argument(
            "thermal solver: case declares no fixed wall temperature; steady "
            "conduction requires at least one Dirichlet patch (otherwise the "
            "system is singular)");
    return bcs;
}

std::size_t nearest_cell(const FvGeometry& geometry, const Vec3& location)
{
    std::size_t best = 0;
    double best_distance = std::numeric_limits<double>::infinity();
    for (std::size_t c = 0; c < geometry.cell_centres.size(); ++c) {
        const double distance = (geometry.cell_centres[c] - location).mag();
        if (distance < best_distance) {
            best_distance = distance;
            best = c;
        }
    }
    return best;
}

void write_text_file_atomic(const std::filesystem::path& target,
                            const std::string& content)
{
    const auto temporary = target.string() + ".tmp";
    {
        std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
        if (!out)
            throw std::runtime_error("cannot open temporary file: " + temporary);
        out << content;
        out.close();
        if (!out)
            throw std::runtime_error("failed writing: " + temporary);
    }
    std::error_code ec;
    std::filesystem::rename(temporary, target, ec);
    if (ec)
        throw std::runtime_error(
            "failed to install " + target.string() + ": " + ec.message());
}

std::string thermal_convergence_json(const ThermalConductionResult& result)
{
    std::ostringstream out;
    out << std::setprecision(17);
    out << "{\n"
        << "  \"format\": \"CFDX-THERMAL-CONVERGENCE\",\n"
        << "  \"schema_version\": 1,\n"
        << "  \"converged\": " << (result.converged ? "true" : "false") << ",\n"
        << "  \"status\": " << static_cast<int>(result.status) << ",\n"
        << "  \"iterations\": " << result.iterations << ",\n"
        << "  \"residual_inf\": " << result.residual_inf << "\n"
        << "}\n";
    return out.str();
}

std::string thermal_summary_json(
    const Mesh& mesh,
    const ThermalZoneModel& zones,
    const ThermalConductionResult& result)
{
    std::ostringstream out;
    out << std::setprecision(17);
    out << "{\n"
        << "  \"format\": \"CFDX-THERMAL-SUMMARY\",\n"
        << "  \"schema_version\": 1,\n"
        << "  \"n_cells\": " << mesh.n_cells() << ",\n"
        << "  \"n_zones\": " << zones.diagnostics.n_zones << ",\n"
        << "  \"min_temperature\": " << result.min_temperature << ",\n"
        << "  \"max_temperature\": " << result.max_temperature << ",\n"
        << "  \"total_source_power\": " << result.total_source_power << ",\n"
        << "  \"boundary_flux_out\": " << result.boundary_flux_out << ",\n"
        << "  \"energy_imbalance_relative\": "
        << result.energy_imbalance_relative << ",\n"
        << "  \"zones\": [\n";
    for (std::size_t i = 0; i < zones.specs.size(); ++i) {
        const auto& spec = zones.specs[i];
        const auto cells_it = zones.diagnostics.cells_per_zone.find(spec.name);
        const std::size_t cells =
            cells_it == zones.diagnostics.cells_per_zone.end()
                ? 0 : cells_it->second;
        out << "    {\"name\": \"" << spec.name
            << "\", \"zone_id\": " << spec.zone_id
            << ", \"conductivity\": " << spec.conductivity
            << ", \"volumetric_power\": " << spec.volumetric_power
            << ", \"contact_resistance\": " << spec.contact_resistance
            << ", \"cells\": " << cells << "}"
            << (i + 1 < zones.specs.size() ? "," : "") << "\n";
    }
    out << "  ],\n"
        << "  \"patch_flux_out\": {";
    bool first = true;
    for (const auto& entry : result.patch_flux_out) {
        if (!first) out << ", ";
        first = false;
        out << "\"" << entry.first << "\": " << entry.second;
    }
    out << "}\n"
        << "}\n";
    return out.str();
}

}  // namespace

int main(int argc, char** argv)
{
    std::filesystem::path execution_output_dir;
    std::optional<ThermalConductionResult> execution_result;
    try {
        const Options options = parse(argc, argv);
        execution_output_dir = options.output_dir;

        Mesh mesh;
        CaseSetup case_setup;
        SourceInfo case_source;
        GapAnalysis case_gap;
        if (!read_case_cfdx_h5(options.mesh, mesh, case_source, case_setup, case_gap))
            throw std::runtime_error("thermal solver: CFDX case load failed");
        if (mesh.n_cells() == 0)
            throw std::runtime_error("thermal solver: case carries no cells");

        // Energy-model gate: the thermal path is only meaningful for a case
        // that actually declares a thermal problem.
        if (case_setup.energy_model != "conjugate" &&
            case_setup.energy_model != "solid_conduction")
            throw std::invalid_argument(
                "thermal solver: case energy_model is '" +
                case_setup.energy_model + "'; the thermal path requires "
                "'conjugate' or 'solid_conduction'");

        const auto geometry = build_fv_geometry(mesh);
        const auto specs = build_zone_specs(case_setup);
        const auto cell_zone = resolve_cell_zones(mesh, options.mesh, specs);
        const auto zones = build_thermal_zones(mesh, cell_zone, specs);
        const auto bcs = build_thermal_boundary_conditions(
            mesh, geometry, case_setup, zones);

        double initial_temperature = case_setup.ref_temperature;
        if (case_setup.has_initial_condition &&
            case_setup.initial_condition.temperature > 0.0)
            initial_temperature = case_setup.initial_condition.temperature;

        Vector temperature(mesh.n_cells(), initial_temperature);
        const ScalarSolveControls controls{5000, 1e-12, 1.0};
        const auto result = solve_thermal_conduction(
            mesh, geometry, zones, bcs, temperature, controls);
        execution_result = result;

        const auto execution = options.output_dir / "execution.json";
        const int solver_exit_code = execution_exit_code_for_solver(result.converged);
        // Persist the solver outcome before any secondary artifact.
        write_execution_summary(execution.string(), ExecutionSummary{
            solver_exit_code,
            result.converged,
            result.iterations,
            static_cast<int>(result.status),
            result.converged
                ? std::string("linear conduction solve converged")
                : std::string("linear conduction solve did not converge"),
            false,
            false});

        const auto convergence = options.output_dir / "thermal_convergence.json";
        write_text_file_atomic(convergence, thermal_convergence_json(result));

        const auto summary = options.output_dir / "thermal_summary.json";
        write_text_file_atomic(summary, thermal_summary_json(mesh, zones, result));

        ScalarCellField temperature_field(
            mesh.n_cells(), "T", "K", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            temperature_field(c) = temperature(c);
        const auto vtu = options.output_dir / "result_0.vtu";
        VtuWriter writer;
        if (!writer.write(vtu.string(), mesh, {{"temperature", temperature_field}},
                          {}, {}, 0.0, 0, true))
            throw std::runtime_error("thermal solver: VTU export failed");

        write_execution_summary(execution.string(), ExecutionSummary{
            solver_exit_code,
            result.converged,
            result.iterations,
            static_cast<int>(result.status),
            result.converged
                ? std::string("linear conduction solve converged")
                : std::string("linear conduction solve did not converge"),
            false,
            std::filesystem::is_regular_file(convergence)});

        std::cout << "Converged " << (result.converged ? "YES" : "NO")
                  << " Iterations " << result.iterations
                  << " Residual " << result.residual_inf << "\n"
                  << "Temperature min " << result.min_temperature
                  << " max " << result.max_temperature << " K\n"
                  << "Energy source " << result.total_source_power
                  << " W boundary_out " << result.boundary_flux_out
                  << " W imbalance " << result.energy_imbalance_relative << "\n"
                  << "Convergence history " << convergence.string() << "\n"
                  << "Thermal summary " << summary.string() << "\n"
                  << "Field export " << vtu.string() << "\n"
                  << "Execution summary " << execution.string() << "\n";

        if (!options.probe_csv.empty()) {
            std::ostringstream csv;
            csv << "name,x,y,z,cell_index,temperature\n";
            for (const auto& probe : options.probes) {
                const std::size_t cell = nearest_cell(geometry, probe.location);
                csv << probe.name << "," << std::setprecision(17)
                    << probe.location.x << "," << probe.location.y << ","
                    << probe.location.z << "," << cell << ","
                    << temperature(cell) << "\n";
            }
            write_text_file_atomic(options.probe_csv, csv.str());
            std::cout << "Probe CSV " << options.probe_csv.string()
                      << " (" << options.probes.size() << " probes)\n";
        }

        return solver_exit_code;
    } catch (const std::exception& exc) {
        std::cerr << "CFDX thermal solver error: " << exc.what() << "\n";
        if (!execution_output_dir.empty()) {
            try {
                const auto execution = execution_output_dir / "execution.json";
                const bool has_convergence = std::filesystem::is_regular_file(
                    execution_output_dir / "thermal_convergence.json");
                write_execution_summary(execution.string(), ExecutionSummary{
                    kExecutionExitError,
                    execution_result.has_value() ? execution_result->converged : false,
                    execution_result.has_value() ? execution_result->iterations : 0,
                    execution_result.has_value()
                        ? static_cast<int>(execution_result->status)
                        : 0,
                    exc.what(),
                    false,
                    has_convergence});
            } catch (const std::exception& summary_exc) {
                std::cerr << "CFDX thermal solver: failed to persist execution "
                             "summary: " << summary_exc.what() << "\n";
            }
        }
        return kExecutionExitError;
    }
}
