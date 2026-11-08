
#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/physics/production_boundary_contract.h"
#include "cfdx/physics/continuation.h"
#include "cfdx/physics/production_energy_stage.h"
#include "cfdx/physics/production_turbulence_contract.h"
#include "cfdx/physics/production_sst_thermal_coupling.h"
#include "cfdx/physics/production_radiation_contract.h"
#include "cfdx/physics/production_radiation_s2s_stage.h"
#include "cfdx/physics/production_buoyancy_stage.h"
#include "cfdx/io/gmsh/gmsh_importer.h"
#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/io/hdf5/case_hdf5_io.h"
#include "cfdx/io/restart/dat_restart.h"
#include "cfdx/io/probe/probe_csv.h"
#include "cfdx/io/runtime/convergence_history.h"
#include "cfdx/io/runtime/execution_summary.h"
#include "cfdx/io/vtu/vtu_writer.h"
#include "cfdx/physics/production_time_contract.h"
#include "cfdx/physics/wall_quadratic_forces.h"
#include "cfdx/physics/thermohydraulic_observables.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <csignal>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#if defined(CFDX_HAS_HDF5) && CFDX_HAS_HDF5
#include <hdf5.h>
#endif

using namespace cfdx::core;
using namespace cfdx::io;
using namespace cfdx::physics;

namespace {
// A CFDX case with a generic .h5 name MUST NOT silently take the mesh-only
// loader: that route omits its physical BCs and can falsely converge at U=0.
// Inspect only the root metadata. Ordinary HDF5 meshes stay supported.
inline bool hdf5_contains_case_setup(const std::string& path)
{
#if defined(CFDX_HAS_HDF5) && CFDX_HAS_HDF5
    hid_t file = -1;
    H5E_BEGIN_TRY { file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT); }
    H5E_END_TRY;
    if (file < 0) return false; // Mesh importer gives the canonical error.
    const htri_t attribute = H5Aexists(file, "case_setup_json");
    (void)H5Fclose(file);
    if (attribute < 0)
        throw std::runtime_error("production solver: cannot inspect HDF5 case metadata");
    return attribute > 0;
#else
    (void)path;
    return false;
#endif
}

volatile std::sig_atomic_t g_stop_requested = 0;

void handle_stop_signal(int) noexcept
{
    g_stop_requested = 1;
}

struct Options {
    std::string mesh;
    std::filesystem::path output_dir;
    std::filesystem::path restart;
    std::size_t iterations = 2;
    bool adaptive_convergence = false;
    bool initial_velocity_field = false;
    bool initial_pressure_field = false;
    bool initial_sst_fields = false;
    bool initial_face_mass_flux = false;
    bool experimental_momentum_first_continuation = false;
    bool wall_forces_quadratic = false;
    // Explicit opt-in only; this path remains experimental and unqualified.
    bool experimental_sst = false;
    bool benchmark_export = false;
    // Opt-in exact pressure CSR capture for production BFS and other cases.
    std::string pressure_replay_prefix;
    std::vector<std::size_t> pressure_replay_iterations;
    std::vector<IncompressiblePointProbe> probes;
    std::filesystem::path probe_csv;
    bool probe_stream = false;
};

std::string encode_probe_name_hex(const std::string& name)
{
    constexpr char digits[] = "0123456789abcdef";
    std::string encoded;
    encoded.reserve(name.size() * 2);
    for (const char raw_byte : name) {
        const auto byte = static_cast<unsigned char>(raw_byte);
        encoded.push_back(digits[byte >> 4]);
        encoded.push_back(digits[byte & 0x0f]);
    }
    return encoded;
}


// Parses `--probe <name>:<x>,<y>,<z>:<field>`.
//
// A probe has to be declarable from a run, otherwise sampling exists only inside
// the C++ API and no run can produce a probe artifact. The field is the probe's
// enum name, so an unsupported one is rejected at parse time rather than
// silently producing a column of the wrong quantity.
IncompressiblePointProbe parse_probe(const std::string& spec)
{
    const auto first = spec.find(':');
    if (first == std::string::npos)
        throw std::invalid_argument("--probe expects <name>:<x>,<y>,<z>:<field>");
    const auto second = spec.find(':', first + 1);
    if (second == std::string::npos)
        throw std::invalid_argument("--probe expects <name>:<x>,<y>,<z>:<field>");
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
    IncompressibleProbeField field;
    if (field_name == "p" || field_name == "pressure")
        field = IncompressibleProbeField::PRESSURE;
    else if (field_name == "u_x")
        field = IncompressibleProbeField::U_X;
    else if (field_name == "u_y")
        field = IncompressibleProbeField::U_Y;
    else if (field_name == "u_z")
        field = IncompressibleProbeField::U_Z;
    else if (field_name == "u_mag" || field_name == "u_magnitude")
        field = IncompressibleProbeField::U_MAGNITUDE;
    else
        throw std::invalid_argument(
            "--probe field must be one of p, u_x, u_y, u_z, u_mag: " + field_name);

    IncompressiblePointProbe probe;
    probe.name = name;
    probe.location = Vec3{coordinates[0], coordinates[1], coordinates[2]};
    probe.field = field;
    return probe;
}


struct LoadedCaseNumerics {
    bool from_case = false;
};

void apply_explicit_case_numerics(
    const cfdx::io::CaseSetup& setup,
    IncompressibleSolverControls& controls)
{
    if (!setup.has_explicit_numerics || !setup.numerical_report.valid())
        throw std::invalid_argument(
            "production solver requires a valid explicit numerics.selection block");

    bool pressure_velocity_selected = false;
    bool convection_selected = false;
    bool gradient_selected = false;
    bool schur_selected = false;
    for (const auto& selection : setup.numerical_report.resolved) {
        if (selection.family == NumericalMethodFamily::PressureVelocity) {
            pressure_velocity_selected = true;
            if (selection.method_id == "pressure_velocity.simple")
                controls.algorithm = PressureVelocityAlgorithm::SIMPLE;
            else if (selection.method_id == "pressure_velocity.simplec")
                controls.algorithm = PressureVelocityAlgorithm::SIMPLEC;
            else if (selection.method_id == "pressure_velocity.piso")
                controls.algorithm = PressureVelocityAlgorithm::PISO;
            else if (selection.method_id == "pressure_velocity.pimple")
                controls.algorithm = PressureVelocityAlgorithm::PIMPLE;
            else if (selection.method_id == "pressure_velocity.fractional_step")
                controls.algorithm = PressureVelocityAlgorithm::FRACTIONAL_STEP;
            else if (selection.method_id == "pressure_velocity.coupled")
                controls.algorithm = PressureVelocityAlgorithm::COUPLED;
            else
                throw std::invalid_argument(
                    "unsupported resolved pressure-velocity selection: " +
                    selection.method_id);
        } else if (selection.family == NumericalMethodFamily::Gradient) {
            gradient_selected = true;
            if (selection.method_id != "gradient.gauss_cell")
                throw std::invalid_argument(
                    "production solver does not yet consume resolved gradient selection: " +
                    selection.method_id);
        } else if (selection.family == NumericalMethodFamily::Convection) {
            convection_selected = true;
            if (selection.method_id == "convection.upwind")
                controls.convection_scheme = ConvectionScheme::UPWIND;
            else if (selection.method_id == "convection.second_order_upwind")
                controls.convection_scheme = ConvectionScheme::SECOND_ORDER_UPWIND;
            else if (selection.method_id.rfind("convection.tvd.", 0) == 0)
                controls.convection_scheme = ConvectionScheme::TVD;
            else
                throw std::invalid_argument(
                    "unsupported resolved convection selection: " +
                    selection.method_id);
        } else if (selection.family == NumericalMethodFamily::LinearSolver) {
            if (selection.method_id == "linear.cg") {
                controls.momentum_linear_solver.krylov = KrylovModel::CG;
                controls.pressure_linear_solver.krylov = KrylovModel::CG;
                controls.coupled_linear_solver.krylov = KrylovModel::CG;
            } else if (selection.method_id == "linear.bicgstab") {
                controls.momentum_linear_solver.krylov = KrylovModel::BiCGStab;
                controls.pressure_linear_solver.krylov = KrylovModel::BiCGStab;
                controls.coupled_linear_solver.krylov = KrylovModel::BiCGStab;
            } else if (selection.method_id == "linear.gmres") {
                controls.momentum_linear_solver.krylov = KrylovModel::GMRES;
                controls.pressure_linear_solver.krylov = KrylovModel::GMRES;
                controls.coupled_linear_solver.krylov = KrylovModel::GMRES;
            } else if (selection.method_id == "linear.fgmres") {
                controls.momentum_linear_solver.krylov = KrylovModel::FGMRES;
                controls.pressure_linear_solver.krylov = KrylovModel::FGMRES;
                controls.coupled_linear_solver.krylov = KrylovModel::FGMRES;
            } else {
                throw std::invalid_argument(
                    "unsupported resolved linear-solver selection: " +
                    selection.method_id);
            }
        } else if (selection.family == NumericalMethodFamily::Initialization) {
            if (selection.method_id == "initialization.provided") {
                controls.initialization.mode = InitializationMode::Provided;
            } else if (selection.method_id == "initialization.uniform") {
                controls.initialization.mode = InitializationMode::Uniform;
            } else if (selection.method_id == "initialization.restart") {
                controls.initialization.mode = InitializationMode::Restart;
            } else {
                throw std::invalid_argument(
                    "unsupported resolved initialization selection: " +
                    selection.method_id);
            }
        } else if (selection.family == NumericalMethodFamily::Temporal) {
            throw std::invalid_argument(
                "steady production solver does not consume a temporal selection: " +
                selection.method_id);
        } else if (selection.family == NumericalMethodFamily::Schur) {
            schur_selected = true;
            if (selection.method_id.rfind("schur.", 0) != 0)
                throw std::invalid_argument(
                    "unsupported resolved Schur selection: " + selection.method_id);
            CoupledSchurModel model = CoupledSchurModel::BlockLocal;
            if (!parse_coupled_schur_model(selection.method_id.substr(6), model))
                throw std::invalid_argument(
                    "unsupported resolved Schur selection: " + selection.method_id);
            controls.coupling.schur_model = model;
        } else if (selection.family == NumericalMethodFamily::Preconditioner) {
            PreconditionerModel model = PreconditionerModel::Auto;
            if (selection.method_id == "preconditioner.native_amg")
                model = PreconditionerModel::NativeAMG;
            else if (selection.method_id == "preconditioner.smoothed_aggregation_amg")
                model = PreconditionerModel::SmoothedAggregationAMG;
            else if (selection.method_id == "preconditioner.native_fieldsplit")
                model = PreconditionerModel::NativeFieldSplit;
            else if (selection.method_id == "preconditioner.coupled_block_schur")
                model = PreconditionerModel::CoupledBlockSchur;
            else if (selection.method_id == "preconditioner.pcd")
                model = PreconditionerModel::PCD;
            else
                throw std::invalid_argument(
                    "unsupported resolved preconditioner selection: " +
                    selection.method_id);
            controls.momentum_linear_solver.preconditioner = model;
            controls.pressure_linear_solver.preconditioner = model;
            controls.coupled_linear_solver.preconditioner = model;
        }
    }

    if (!pressure_velocity_selected)
        throw std::invalid_argument(
            "case numerics.selection does not select pressure_velocity");
    if (!convection_selected)
        throw std::invalid_argument(
            "case numerics.selection does not select convection");
    if (!gradient_selected)
        throw std::invalid_argument(
            "case numerics.selection does not select a consumed gradient method");
    if (controls.algorithm == PressureVelocityAlgorithm::COUPLED && !schur_selected)
        throw std::invalid_argument(
            "coupled production solver requires an explicit Schur numerical selection");
}

// Writes a text artifact atomically: the payload goes to a temporary
// sibling first and is then renamed, so a reader never observes a
// half-written artifact.
void write_text_file_atomic(
    const std::filesystem::path& path, const std::string& content)
{
    auto temporary = path;
    temporary += ".tmp";
    {
        std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
        if (!out)
            throw std::runtime_error(
                "production solver: cannot open " + temporary.string());
        out << content;
        if (!out)
            throw std::runtime_error(
                "production solver: cannot write " + temporary.string());
    }
    std::filesystem::rename(temporary, path);
}

// JSON string escaping: patch names come from the case/mesh and are not
// guaranteed to be JSON-safe. Quotes, backslashes and control characters
// are escaped; everything else passes through as UTF-8. Without this, a
// patch name like a"b would produce an artifact no parser can read.
std::string json_escape(const std::string& text)
{
    static const char* hex = "0123456789abcdef";
    std::string out;
    out.reserve(text.size() + 2);
    for (const char c : text) {
        switch (c) {
        case '"': out += "\\\""; break;
        case '\\': out += "\\\\"; break;
        case '\b': out += "\\b"; break;
        case '\f': out += "\\f"; break;
        case '\n': out += "\\n"; break;
        case '\r': out += "\\r"; break;
        case '\t': out += "\\t"; break;
        default: {
            const auto uc = static_cast<unsigned char>(c);
            if (uc < 0x20) {
                out += "\\u00";
                out += hex[(uc >> 4) & 0xF];
                out += hex[uc & 0xF];
            } else {
                out += c;
            }
            break;
        }
        }
    }
    return out;
}

// JSON number: a non-finite double must never serialize as NaN/Infinity,
// which are not valid JSON numbers. It becomes null so the artifact stays
// machine-readable in every failure mode while the explicit convergence
// gates carry the actual failure.
std::string json_number(double value)
{
    if (!std::isfinite(value))
        return "null";
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

// Machine-readable convergence history of the T2 energy stage. The
// independent conservation residual is reported per iteration so the
// artifact carries the audit, not only the assembled residual. The
// blocking gates (linear, global balance, local balance) and their
// normalization scale are explicit: a reader can reconstruct why a run
// was declared converged or not without re-running the audit.
std::string energy_convergence_json(
    const ProductionEnergyStageResult& energy)
{
    std::ostringstream out;
    out << "{\n";
    out << "  \"converged\": " << (energy.converged ? "true" : "false")
        << ",\n";
    out << "  \"linear_solver_converged\": "
        << (energy.linear_solver_converged ? "true" : "false") << ",\n";
    out << "  \"linear_solver_status\": " << energy.linear_solver_status << ",\n";
    out << "  \"linear_residual_relative\": "
        << json_number(energy.linear_residual_relative) << ",\n";
    out << "  \"global_energy_balance_passed\": "
        << (energy.global_energy_balance_passed ? "true" : "false") << ",\n";
    out << "  \"local_energy_balance_passed\": "
        << (energy.local_energy_balance_passed ? "true" : "false") << ",\n";
    out << "  \"temperature_admissibility_passed\": "
        << (energy.temperature_admissibility_passed ? "true" : "false") << ",\n";
    out << "  \"nonpositive_temperature\": " << energy.nonpositive_temperature
        << ",\n";
    out << "  \"iterations\": " << energy.iterations << ",\n";
    out << "  \"independent_conservation_residual\": "
        << json_number(energy.independent_conservation_residual) << ",\n";
    out << "  \"independent_conservation_normalized\": "
        << json_number(energy.independent_conservation_normalized) << ",\n";
    out << "  \"audit_power_scale\": "
        << json_number(energy.audit_power_scale) << ",\n";
    out << "  \"global_balance_normalized\": "
        << json_number(energy.global_balance_normalized) << ",\n";
    out << "  \"max_cell_residual\": "
        << json_number(energy.max_cell_residual) << ",\n";
    out << "  \"worst_cell\": " << energy.worst_cell << ",\n";
    out << "  \"local_balance_normalized\": "
        << json_number(energy.local_balance_normalized) << ",\n";
    out << "  \"max_cell_residual_normalized\": "
        << json_number(energy.max_cell_residual_normalized) << ",\n";
    out << "  \"worst_cell_normalized\": "
        << energy.worst_cell_normalized << ",\n";
    out << "  \"nonfinite_faces\": " << energy.nonfinite_faces << ",\n";
    out << "  \"nonfinite_cells\": " << energy.nonfinite_cells << ",\n";
    out << "  \"history\": [\n";
    for (std::size_t i = 0; i < energy.history.size(); ++i) {
        const auto& entry = energy.history[i];
        out << "    {\"iteration\": " << entry.iteration
            << ", \"residual\": " << json_number(entry.residual)
            << ", \"backward_error\": "
            << json_number(entry.backward_error)
            << ", \"linear_residual_relative\": "
            << json_number(entry.linear_residual_relative)
            << ", \"linear_solver_status\": " << entry.linear_solver_status
            << ", \"linear_solver_converged\": "
            << (entry.linear_solver_converged ? "true" : "false")
            << ", \"independent_conservation_normalized\": "
            << json_number(entry.independent_conservation_normalized)
            << "}";
        out << (i + 1 == energy.history.size() ? "\n" : ",\n");
    }
    out << "  ]\n";
    out << "}\n";
    return out.str();
}

// Machine-readable physical summary of the T2 energy stage: temperature
// extrema, source and boundary power, and per-patch fluxes (positive
// leaving the domain).
std::string energy_summary_json(const ProductionEnergyStageResult& energy)
{
    std::ostringstream out;
    out << "{\n";
    out << "  \"min_temperature\": "
        << json_number(energy.min_temperature) << ",\n";
    out << "  \"max_temperature\": "
        << json_number(energy.max_temperature) << ",\n";
    out << "  \"total_source_power\": "
        << json_number(energy.total_source_power) << ",\n";
    out << "  \"boundary_flux_out\": "
        << json_number(energy.boundary_flux_out) << ",\n";
    out << "  \"nonfinite_temperature\": " << energy.nonfinite_temperature
        << ",\n";
    out << "  \"nonpositive_temperature\": " << energy.nonpositive_temperature
        << ",\n";
    out << "  \"temperature_admissibility_passed\": "
        << (energy.temperature_admissibility_passed ? "true" : "false") << ",\n";
    out << "  \"volumetric_power\": "
        << json_number(energy.volumetric_power) << ",\n";
    out << "  \"patches\": [\n";
    std::size_t index = 0;
    for (const auto& entry : energy.patch_flux_out) {
        out << "    {\"name\": \"" << json_escape(entry.first)
            << "\", \"flux_out\": " << json_number(entry.second) << "}";
        out << (++index == energy.patch_flux_out.size() ? "\n" : ",\n");
    }
    out << "  ]\n";
    out << "}\n";
    return out.str();
}

// JSON summary of the external S2S radiation + energy iteration (issue
// #885 step 3): the coupled fixed-point verdict, its separate blocking
// criteria (unrelaxed radiative fixed-point residual, wall-temperature
// stationarity), the per-iteration convergence history, and the
// per-surface radiosity data of the final radiative response.
std::string radiation_s2s_json(const ProductionS2SOuterLoopResult& stage)
{
    std::ostringstream out;
    out << "{\n";
    out << "  \"converged\": " << (stage.converged ? 1 : 0) << ",\n";
    out << "  \"fixed_point_flux_passed\": "
        << (stage.fixed_point_flux_passed ? 1 : 0) << ",\n";
    out << "  \"fixed_point_wall_temperature_passed\": "
        << (stage.fixed_point_wall_temperature_passed ? 1 : 0) << ",\n";
    out << "  \"outer_iterations\": " << stage.outer_iterations << ",\n";
    out << "  \"max_flux_fixed_point_normalized\": "
        << json_number(stage.max_flux_fixed_point_normalized) << ",\n";
    out << "  \"max_wall_temperature_change_normalized\": "
        << json_number(stage.max_wall_temperature_change_normalized)
        << ",\n";
    out << "  \"max_wall_temperature_lag\": "
        << json_number(stage.max_wall_temperature_lag) << ",\n";
    out << "  \"outer_history\": [\n";
    for (std::size_t i = 0; i < stage.history.size(); ++i) {
        const auto& entry = stage.history[i];
        out << "    {\"iteration\": " << entry.iteration
            << ", \"energy_converged\": "
            << (entry.energy_converged ? 1 : 0)
            << ", \"radiation_converged\": "
            << (entry.radiation_converged ? 1 : 0)
            << ", \"max_flux_fixed_point_normalized\": "
            << json_number(entry.max_flux_fixed_point_normalized)
            << ", \"max_wall_temperature_change_normalized\": "
            << json_number(entry.max_wall_temperature_change_normalized)
            << "}";
        out << (i + 1 == stage.history.size() ? "\n" : ",\n");
    }
    out << "  ],\n";
    const ProductionS2SRadiationResult& radiation = stage.radiation;
    out << "  \"radiation\": {\n";
    out << "    \"balance_passed\": "
        << (radiation.balance_passed ? 1 : 0) << ",\n";
    out << "    \"surface_identity_passed\": "
        << (radiation.surface_identity_passed ? 1 : 0) << ",\n";
    out << "    \"total_net_power\": "
        << json_number(radiation.total_net_power) << ",\n";
    out << "    \"balance_scale\": "
        << json_number(radiation.balance_scale) << ",\n";
    out << "    \"balance_normalized\": "
        << json_number(radiation.balance_normalized) << ",\n";
    out << "    \"max_surface_identity_normalized\": "
        << json_number(radiation.max_surface_identity_normalized) << ",\n";
    out << "    \"surface_count\": " << radiation.surface_count << ",\n";
    out << "    \"surfaces\": [\n";
    for (std::size_t index = 0; index < radiation.surface_names.size();
         ++index) {
        const std::string& name = radiation.surface_names[index];
        out << "      {\"name\": \"" << json_escape(name)
            << "\", \"area\": " << json_number(radiation.areas[index])
            << ", \"emissivity\": "
            << json_number(radiation.emissivities[index])
            << ", \"wall_temperature\": "
            << json_number(radiation.wall_temperatures[index])
            << ", \"net_flux\": "
            << json_number(radiation.net_flux[index])
            << ", \"net_power\": "
            << json_number(radiation.net_power[index])
            << ", \"injected_heat_flux\": "
            << json_number(radiation.injected_heat_flux.at(name))
            << ", \"applied_heat_flux\": "
            << json_number(radiation.applied_heat_flux.at(name)) << "}";
        out << (index + 1 == radiation.surface_names.size() ? "\n"
                                                            : ",\n");
    }
    out << "    ]\n";
    out << "  }\n";
    out << "}\n";
    return out.str();
}

void write_quadratic_wall_forces(
    const Mesh& mesh,
    const Field<double, Location::CELL>& velocity,
    const Field<double, Location::CELL>& pressure,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    double dynamic_viscosity,
    const std::filesystem::path& output_dir)
{
    if (!(dynamic_viscosity > 0.0) || !std::isfinite(dynamic_viscosity))
        throw std::invalid_argument(
            "quadratic wall forces require finite positive dynamic viscosity");

    const auto geometry = build_fv_geometry(mesh);
    std::vector<double> cell_pressure(mesh.n_cells());
    std::vector<Vec3> cell_velocity(mesh.n_cells());
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        cell_pressure[cell] = pressure(cell);
        cell_velocity[cell] = Vec3{
            velocity(cell, 0), velocity(cell, 1), velocity(cell, 2)};
    }

    std::vector<forces::SurfaceSample> surface;
    std::vector<std::string> wall_patch_names;
    double total_area = 0.0;
    for (std::size_t patch_id = 0;
         patch_id < mesh.boundary().n_patches(); ++patch_id) {
        const auto& patch = mesh.boundary().patch(patch_id);
        if (patch.type != PatchType::WALL)
            continue;

        const auto velocity_it = velocity_bcs.find(patch.name);
        if (velocity_it == velocity_bcs.end() ||
            velocity_it->second.type != VelocityBoundaryCondition::Type::FIXED_VALUE)
            throw std::invalid_argument(
                "quadratic wall forces require prescribed FIXED_VALUE velocity on WALL patch '" +
                patch.name + "'");
        const auto pressure_it = pressure_bcs.find(patch.name);
        if (pressure_it == pressure_bcs.end())
            throw std::invalid_argument(
                "quadratic wall forces have no pressure boundary contract for WALL patch '" +
                patch.name + "'");

        std::vector<double> wall_pressure(mesh.n_faces(), 0.0);
        std::vector<Vec3> wall_velocity(mesh.n_faces());
        for (const auto face : patch.face_ids) {
            if (face >= mesh.n_faces())
                throw std::invalid_argument(
                    "quadratic wall forces: patch face index is out of range");
            const auto owner = mesh.ownership().owner(face);
            if (owner >= mesh.n_cells())
                throw std::invalid_argument(
                    "quadratic wall forces: invalid owner cell");
            if (pressure_it->second.type == ScalarBoundaryType::FIXED_VALUE) {
                wall_pressure[face] = pressure_it->second.value;
            } else if (pressure_it->second.type == ScalarBoundaryType::ZERO_GRADIENT) {
                wall_pressure[face] = pressure(owner);
            } else {
                throw std::invalid_argument(
                    "quadratic wall forces do not support FIXED_GRADIENT pressure on WALL patch '" +
                    patch.name + "'");
            }
            wall_velocity[face] = velocity_it->second.value;
        }

        auto patch_surface = reconstruct_quadratic_wall_surface(
            mesh, geometry, patch_id, cell_pressure, cell_velocity,
            wall_pressure, wall_velocity, dynamic_viscosity, 4);
        for (const auto& sample : patch_surface)
            total_area += sample.area;
        surface.insert(surface.end(), patch_surface.begin(), patch_surface.end());
        wall_patch_names.push_back(patch.name);
    }

    if (surface.empty() || !(total_area > 0.0) || !std::isfinite(total_area))
        throw std::invalid_argument(
            "quadratic wall forces found no finite, positive-area WALL surface");

    // The artifact reports dimensional forces/moments only. A user-defined
    // rho_inf/U_inf/reference area is required before publishing Cd/Cl/Cm.
    forces::Reference reference;
    reference.area = total_area;
    reference.length = 1.0;
    reference.pressure_ref = 0.0;
    const auto integrated = forces::integrate(surface, reference);

    std::ostringstream out;
    out << "{\n"
        << "  \"format\": \"CFDX-WALL-FORCES\",\n"
        << "  \"schema_version\": 1,\n"
        << "  \"method\": \"quadratic_wls_3d\",\n"
        << "  \"normal_convention\": \"body_to_fluid\",\n"
        << "  \"pressure_reference_pa\": 0,\n"
        << "  \"moment_center_m\": [0, 0, 0],\n"
        << "  \"force_definition\": \"fluid_on_body = -p*n_body + tau*n_body\",\n"
        << "  \"dynamic_viscosity_pa_s\": " << json_number(dynamic_viscosity) << ",\n"
        << "  \"surface_face_count\": " << surface.size() << ",\n"
        << "  \"surface_area_m2\": " << json_number(total_area) << ",\n"
        << "  \"wall_patches\": [";
    for (std::size_t i = 0; i < wall_patch_names.size(); ++i) {
        if (i != 0) out << ", ";
        out << "\"" << json_escape(wall_patch_names[i]) << "\"";
    }
    out << "],\n"
        << "  \"pressure_force_n\": ["
        << json_number(integrated.pressure_force.x) << ", "
        << json_number(integrated.pressure_force.y) << ", "
        << json_number(integrated.pressure_force.z) << "],\n"
        << "  \"viscous_force_n\": ["
        << json_number(integrated.viscous_force.x) << ", "
        << json_number(integrated.viscous_force.y) << ", "
        << json_number(integrated.viscous_force.z) << "],\n"
        << "  \"total_force_n\": ["
        << json_number(integrated.total_force.x) << ", "
        << json_number(integrated.total_force.y) << ", "
        << json_number(integrated.total_force.z) << "],\n"
        << "  \"total_moment_nm\": ["
        << json_number(integrated.total_moment.x) << ", "
        << json_number(integrated.total_moment.y) << ", "
        << json_number(integrated.total_moment.z) << "],\n"
        << "  \"nondimensional_coefficients\": null\n"
        << "}\n";
    write_text_file_atomic(output_dir / "wall_forces.json", out.str());
    std::cout << "Quadratic 3-D wall forces: " << surface.size()
              << " faces, area " << total_area
              << " m2; wrote wall_forces.json\n";
}

Options parse(int argc, char** argv)
{
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto value = [&](const char* name) -> std::string {
            if (i + 1 >= argc) throw std::invalid_argument(std::string("missing value for ") + name);
            return argv[++i];
        };
        if (arg == "--mesh") o.mesh = value("--mesh");
        else if (arg == "--output-dir") o.output_dir = value("--output-dir");
        else if (arg == "--restart") o.restart = value("--restart");
        else if (arg == "--iterations") o.iterations = std::stoull(value("--iterations"));
        else if (arg == "--adaptive-convergence") o.adaptive_convergence = true;
        else if (arg == "--initial-velocity-field") o.initial_velocity_field = true;
        else if (arg == "--initial-pressure-field") o.initial_pressure_field = true;
        else if (arg == "--initial-sst-fields") o.initial_sst_fields = true;
        else if (arg == "--initial-face-mass-flux") o.initial_face_mass_flux = true;
        else if (arg == "--experimental-momentum-first-continuation")
            o.experimental_momentum_first_continuation = true;
        else if (arg == "--wall-forces-quadratic") o.wall_forces_quadratic = true;
        else if (arg == "--experimental-sst") o.experimental_sst = true;
        else if (arg == "--benchmark-export") o.benchmark_export = true;
        else if (arg == "--pressure-replay-prefix")
            o.pressure_replay_prefix = value("--pressure-replay-prefix");
        else if (arg == "--pressure-replay-iterations")
            o.pressure_replay_iterations = cfdx::core::parse_linear_capture_iterations(
                value("--pressure-replay-iterations"));
        else if (arg == "--probe") o.probes.push_back(parse_probe(value("--probe")));
        else if (arg == "--probe-csv") o.probe_csv = value("--probe-csv");
        else if (arg == "--probe-stream") o.probe_stream = true;
        else if (arg == "--help") {
            std::cout << "cfdx_production_solver --mesh PATH --output-dir DIR "
                         "[--restart DAT] [--iterations N] [--adaptive-convergence]\n"
                         "  [--initial-velocity-field] (initialize U from /fields/vector/U)\n"
                         "  [--initial-pressure-field] (initialize p from /fields/scalar/p; requires imported U)\n"
                         "  [--initial-sst-fields] (initial k and omega from /fields/scalar)\n"
                         "  [--initial-face-mass-flux] (projected /fields/face/mdot; requires initial velocity)\n"
                         "  [--experimental-momentum-first-continuation] (momentum-only then strict SST; NOT qualified)\n"
                         "  [--wall-forces-quadratic] (opt-in 3-D quadratic WLS wall-force artifact; rank-deficient stencils fail)\n"
                         "  [--experimental-sst] (explicit experimental k-omega SST opt-in; not qualified)\n"
                         "  [--benchmark-export] (raw cell/face CSV, wall shear diagnostic and mass balance)\n"
                         "  [--pressure-replay-prefix PATH --pressure-replay-iterations 1,20]\n"
                         "    (diagnostic post-gauge pressure CSR, RHS, initial guess; opt-in)\n"
                         "  [--probe <name>:<x>,<y>,<z>:<field>]... [--probe-csv PATH]\n"
                         "  probe field: p | u_x | u_y | u_z | u_mag\n";
            std::cout << "  --probe-stream emits committed samples as CFDX_PROBE_SAMPLE_V1 stdout records\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + arg);
        }
    }
    if (o.mesh.empty() || o.output_dir.empty() || o.iterations == 0)
        throw std::invalid_argument("--mesh, --output-dir and positive --iterations are required");
    // Asking for an export with nothing to sample would otherwise produce an
    // empty file, or worse a file with a header and no data.
    if (!o.probe_csv.empty() && o.probes.empty())
        throw std::invalid_argument("--probe-csv requires at least one --probe");
    if (o.probe_stream && o.probes.empty())
        throw std::invalid_argument("--probe-stream requires at least one --probe");
    if (o.pressure_replay_prefix.empty() != o.pressure_replay_iterations.empty())
        throw std::invalid_argument(
            "--pressure-replay-prefix and --pressure-replay-iterations must be supplied together");
    for (std::size_t i = 0; i < o.probes.size(); ++i)
        for (std::size_t j = i + 1; j < o.probes.size(); ++j)
            if (o.probes[i].name == o.probes[j].name)
                throw std::invalid_argument("duplicate probe name: " + o.probes[i].name);
    std::filesystem::create_directories(o.output_dir);
    return o;
}
}

int main(int argc, char** argv)
{
    std::signal(SIGINT, handle_stop_signal);
    std::signal(SIGTERM, handle_stop_signal);
    std::filesystem::path execution_output_dir;
    std::optional<IncompressibleSolveResult> execution_result;
    try {
        const Options options = parse(argc, argv);
        execution_output_dir = options.output_dir;
        Mesh mesh;
        const std::string mesh_path = options.mesh;
        CaseSetup case_setup;
        SourceInfo case_source;
        GapAnalysis case_gap;
        // A CFDX case is distinct from a generic mesh HDF5 file.  The
        // suffix is 8 characters (".cfdx.h5"); the previous 9-character
        // check could never match a valid case path and silently routed
        // .cfdx.h5 files through the mesh-only HDF5 reader.  That bypassed
        // case numerics, boundary conditions and initialization metadata.
        const bool is_case_hdf5 =
            mesh_path.size() >= 8 &&
            mesh_path.substr(mesh_path.size() - 8) == ".cfdx.h5";
        bool ok = false;
        if (is_case_hdf5) {
            ok = read_case_cfdx_h5(mesh_path, mesh, case_source, case_setup, case_gap);
            if (!ok)
                throw std::runtime_error("production solver: CFDX case load failed");
        } else if (mesh_path.size() >= 3 &&
                   mesh_path.substr(mesh_path.size() - 3) == ".h5") {
            if (hdf5_contains_case_setup(mesh_path))
                throw std::invalid_argument(
                    "production solver: HDF5 contains case_setup_json but file "
                    "does not end in .cfdx.h5; rename it to load the physical case "
                    "(case boundaries/initialization must not be ignored)");
            ok = read_mesh_hdf5(mesh_path, mesh);
        } else {
            ok = gmsh::import_gmsh_mesh(mesh_path, mesh);
        }
        if (!ok || mesh.n_cells() == 0)
            throw std::runtime_error("production solver: mesh import failed");
        if (options.experimental_sst && !is_case_hdf5)
            throw std::invalid_argument(
                "production solver: --experimental-sst requires a .cfdx.h5 case");

        // The declared energy model is validated fail-closed as soon as the
        // case is loaded, before any expensive flow computation: an unknown
        // or unsupported declaration must fail the run up front instead of
        // being discovered after the carrier solve, even when the energy
        // stage itself would not run. The buoyant model additionally
        // requires its coupling declarations (gravity vector, thermal
        // expansion coefficient, reference temperature) to be valid at
        // load time; the flow-energy buoyancy outer loop below is the only
        // solver path that may run it.
        if (is_case_hdf5) {
            validate_production_steady_time_selection(case_setup, "production solver");
            validate_production_turbulence_selection(
                case_setup, options.experimental_sst);
            validate_production_energy_model(case_setup);
            validate_production_radiation_selection(case_setup);
            validate_production_buoyancy_model(case_setup);
            // Supported-by-this-executable gate: the schema accepts more
            // energy models than this incompressible frozen-flow solver
            // implements. "solid_conduction" needs the dedicated thermal
            // solver and "compressible" needs a compressible solver;
            // running this executable on them would silently compute a
            // different physics than the case declares (an incompressible
            // isothermal flow), so they are rejected with an explicit
            // routing message instead of being accepted and ignored.
            if (case_setup.energy_model == "solid_conduction")
                throw std::invalid_argument(
                    "production solver: energy_model 'solid_conduction' is "
                    "not solved by this executable; use the dedicated "
                    "thermal solver");
            if (case_setup.energy_model == "compressible")
                throw std::invalid_argument(
                    "production solver: energy_model 'compressible' is not "
                    "solved by this executable (incompressible carrier "
                    "only); use a compressible solver for this case");
        }

        if (options.initial_pressure_field && (!options.initial_velocity_field ||
            !options.restart.empty() || !is_case_hdf5))
            throw std::invalid_argument(
                "--initial-pressure-field requires --initial-velocity-field, HDF5, and no DAT restart");
        if (options.initial_sst_fields && (!options.experimental_sst ||
            !options.restart.empty() || !is_case_hdf5))
            throw std::invalid_argument(
                "--initial-sst-fields requires --experimental-sst, HDF5, and no DAT restart");
        if (options.initial_face_mass_flux && (!options.initial_velocity_field ||
            !options.restart.empty() || !is_case_hdf5))
            throw std::invalid_argument(
                "--initial-face-mass-flux requires --initial-velocity-field, HDF5, and no DAT restart");
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);

        if (options.initial_velocity_field) {
            if (!is_case_hdf5)
                throw std::invalid_argument(
                    "production solver: --initial-velocity-field requires a .cfdx.h5 case");
            std::vector<std::pair<std::string, Vec3CellField>> imported_vectors;
            if (!read_vector_fields_hdf5(mesh_path, imported_vectors))
                throw std::runtime_error(
                    "production solver: failed to read vector fields from case HDF5");
            const Vec3CellField* imported_u = nullptr;
            for (const auto& entry : imported_vectors) {
                if (entry.first == "U") {
                    if (imported_u != nullptr)
                        throw std::invalid_argument(
                            "production solver: duplicate imported U fields");
                    imported_u = &entry.second;
                }
            }
            if (imported_u == nullptr)
                throw std::invalid_argument(
                    "production solver: --initial-velocity-field requested but field 'U' is absent");
            if (imported_u->size() != mesh.n_cells())
                throw std::invalid_argument(
                    "production solver: imported U cell count does not match mesh");
            for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
                const auto value = (*imported_u)(cell);
                if (!std::isfinite(value.x) || !std::isfinite(value.y) ||
                    !std::isfinite(value.z))
                    throw std::invalid_argument(
                        "production solver: imported U contains non-finite values");
                U.set(cell, value.x, value.y, value.z);
            }
            std::cout << "Imported cell-centred velocity field U from case HDF5; "
                         "pressure initialization remains case-declared.\n";

            // Preserve machine-readable evidence of the exact imported field.
            // This artifact is written before solver setup can reject unrelated
            // case contracts, making import failures distinguishable from later
            // production-solver qualification failures.
            std::ostringstream imported_velocity_json;
            imported_velocity_json << "{\n  \"field\": \"U\",\n"
                                   << "  \"location\": \"cell\",\n"
                                   << "  \"cell_count\": " << U.size() << ",\n"
                                   << "  \"values\": [";
            for (std::size_t cell = 0; cell < U.size(); ++cell) {
                if (cell != 0)
                    imported_velocity_json << ",";
                imported_velocity_json << "\n    [" << json_number(U(cell, 0))
                                       << ", " << json_number(U(cell, 1))
                                       << ", " << json_number(U(cell, 2)) << "]";
            }
            if (U.size() != 0)
                imported_velocity_json << "\n  ";
            imported_velocity_json << "]\n}\n";
            write_text_file_atomic(
                options.output_dir / "imported_velocity_field.json",
                imported_velocity_json.str());
        }

        // The initial state is declared, not inferred. When the case carries an
        // explicit initial_condition block it is mapped onto the fields through
        // the uniform strategy; otherwise the caller-provided state is kept.
        // A restart, when requested, replaces it below.
        InitializationControls initialization;
        if (is_case_hdf5 && case_setup.has_initial_condition) {
            const auto& ic = case_setup.initial_condition;
            initialization.mode = options.initial_velocity_field
                ? InitializationMode::Provided : InitializationMode::Uniform;
            if (options.initial_velocity_field)
                p.fill(ic.pressure);
            if (!ic.velocity_vector.empty()) {
                if (ic.velocity_vector.size() != 3)
                    throw std::runtime_error(
                        "production solver: initial_condition.velocity_vector must have 3 components");
                initialization.uniform_velocity = cfdx::core::Vec3{
                    ic.velocity_vector[0], ic.velocity_vector[1], ic.velocity_vector[2]};
            } else {
                initialization.uniform_velocity =
                    cfdx::core::Vec3{ic.velocity, ic.velocity, ic.velocity};
            }
            initialization.uniform_pressure = ic.pressure;
        }
        if (options.initial_velocity_field && !case_setup.has_initial_condition)
            initialization.mode = InitializationMode::Provided;

        // A nonuniform pressure seed is useful for high-aspect-ratio BFS
        // meshes whose first segregated pressure correction is very large.
        // Import it only on explicit request, after uniform initialization
        // has been established. The solver still reassembles its own FV
        // equations and uses every original residual/conservation gate.
        if (options.initial_pressure_field) {
            std::vector<std::pair<std::string, ScalarCellField>> scalars;
            if (!read_scalar_fields_hdf5(mesh_path, scalars))
                throw std::invalid_argument(
                    "initial pressure fields cannot be read from case HDF5");
            const ScalarCellField* imported_p = nullptr;
            for (const auto& [name, field] : scalars) {
                if (name != "p") continue;
                if (imported_p)
                    throw std::invalid_argument("duplicate imported pressure field");
                imported_p = &field;
            }
            if (!imported_p)
                throw std::invalid_argument(
                    "--initial-pressure-field requested but scalar field 'p' is absent");
            if (imported_p->size() != mesh.n_cells() || imported_p->dimension() != 1)
                throw std::invalid_argument("imported pressure field shape mismatch");
            double pmin = std::numeric_limits<double>::infinity();
            double pmax = -std::numeric_limits<double>::infinity();
            for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
                const double value = (*imported_p)(cell);
                if (!std::isfinite(value))
                    throw std::invalid_argument("imported pressure field is not finite");
                p(cell) = value;
                pmin = std::min(pmin, value);
                pmax = std::max(pmax, value);
            }
            std::cout << "INITIAL_PRESSURE_FIELD_IMPORTED cells=" << mesh.n_cells()
                      << " min=" << std::setprecision(17) << pmin
                      << " max=" << pmax
                      << " (initial guess only; NOT a convergence certificate)\n";
        }

        if (is_case_hdf5)
            validate_production_boundary_unique_names(case_setup);
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        for (std::size_t i = 0; i < mesh.boundary().n_patches(); ++i) {
            const auto& patch = mesh.boundary().patch(i);
            if (!is_case_hdf5) {
                ubc[patch.name] = {
                    VelocityBoundaryCondition::Type::FIXED_VALUE,
                    {0.0, 0.0, 0.0}};
                pbc[patch.name] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
                continue;
            }

            const auto* bc = case_setup.find_boundary(patch.name);
            if (bc == nullptr)
                throw std::invalid_argument(
                    "production solver: case has no boundary condition for mesh patch '" +
                    patch.name + "'");
            // A physical label must never be silently interpreted as an
            // unrelated generic zero-gradient condition.
            validate_production_incompressible_boundary(*bc, patch.type);

            if (bc->value_type == BCValueType::FIXED ||
                bc->value_type == BCValueType::WALL_NO_SLIP ||
                bc->value_type == BCValueType::WALL_THERMAL) {
                cfdx::core::Vec3 value{0.0, 0.0, 0.0};
                if (!bc->velocity_vector.empty()) {
                    if (bc->velocity_vector.size() != 3)
                        throw std::invalid_argument(
                            "production solver: velocity_vector must have 3 components for patch '" +
                            patch.name + "'");
                    value = cfdx::core::Vec3{
                        bc->velocity_vector[0],
                        bc->velocity_vector[1],
                        bc->velocity_vector[2]};
                } else if (bc->velocity_magnitude != 0.0) {
                    value = cfdx::core::Vec3{bc->velocity_magnitude, 0.0, 0.0};
                }
                ubc[patch.name] = {
                    VelocityBoundaryCondition::Type::FIXED_VALUE, value};
            } else if (bc->value_type == BCValueType::ZERO_GRADIENT ||
                       bc->value_type == BCValueType::WALL_SLIP ||
                       bc->type == BCType::OUTLET ||
                       bc->type == BCType::PRESSURE_OUTLET ||
                       bc->type == BCType::SYMMETRY) {
                ubc[patch.name] = {
                    VelocityBoundaryCondition::Type::ZERO_GRADIENT,
                    {0.0, 0.0, 0.0}};
            } else {
                throw std::invalid_argument(
                    "production solver: unsupported velocity BC '" +
                    patch.name + "'");
            }

            if (bc->value_type == BCValueType::OUTLET_PRESSURE) {
                pbc[patch.name] = {
                    ScalarBoundaryType::FIXED_VALUE, bc->pressure, 0.0};
            } else {
                pbc[patch.name] = {
                    ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
            }
        }

        IncompressibleSolverControls controls;
        controls.convergence.max_iterations = options.iterations;
        if (is_case_hdf5) {
            apply_explicit_case_numerics(case_setup, controls);
            // A parsed relaxation value that is silently ignored prevents
            // a difficult SST separation case from being stabilized by its
            // declared numerical input. Carry the case-level controls into
            // SIMPLE/SIMPLEC/PISO; do not substitute hidden defaults.
            const double requested_alpha_u =
                case_setup.numerics.under_relaxation_momentum;
            const double requested_alpha_p =
                case_setup.numerics.under_relaxation_pressure;
            if (!std::isfinite(requested_alpha_u) ||
                !(requested_alpha_u > 0.0) || requested_alpha_u > 1.0 ||
                !std::isfinite(requested_alpha_p) ||
                !(requested_alpha_p > 0.0) || requested_alpha_p > 1.0)
                throw std::invalid_argument(
                    "production solver: case momentum/pressure relaxation "
                    "must be finite and in (0,1]");
            controls.coupling.alpha_u = requested_alpha_u;
            controls.coupling.alpha_p = requested_alpha_p;
            std::cout << "Case under-relaxation alpha_u=" << requested_alpha_u
                      << " alpha_p=" << requested_alpha_p << '\n';
            std::cout << "Resolved numerical selections:\n"
                      << cfdx::core::format_numerics_report(case_setup.numerical_report);
            if (case_setup.numerics.max_iterations > 0)
                controls.convergence.max_iterations =
                    std::min<std::size_t>(
                        controls.convergence.max_iterations,
                        static_cast<std::size_t>(case_setup.numerics.max_iterations));
            // Case-level uniform body force. The solver consumes it per unit
            // mass (it multiplies by density internally), so the case declares
            // an acceleration vector in m/s^2. An absent block leaves the
            // solver default untouched.
            if (!case_setup.body_force.empty()) {
                if (case_setup.body_force.size() != 3)
                    throw std::invalid_argument(
                        "production solver: case body_force must have 3 components");
                controls.body_force = cfdx::core::Vec3{
                    case_setup.body_force[0],
                    case_setup.body_force[1],
                    case_setup.body_force[2]};
            }
            // Material properties: a single material declaring ANY fluid
            // property must declare BOTH a finite positive density and a
            // finite positive dynamic viscosity (nu = mu / rho). A partial
            // declaration is rejected instead of silently letting the
            // solver defaults substitute for the declared fluid: the flow
            // must be solved with the material the case declares, or the
            // velocity field and mass flux are computed with a different
            // physics than the case advertises. A material block that
            // declares no fluid property at all (e.g. thermal-only) keeps
            // the solver defaults, and more than one declared material is
            // rejected: the production path must not silently pick a
            // winner.
            if (case_setup.materials.size() == 1) {
                const auto& material = case_setup.materials.front();
                const bool declares_fluid_properties =
                    material.density != 0.0 ||
                    material.dynamic_viscosity != 0.0;
                if (declares_fluid_properties) {
                    if (!std::isfinite(material.density) ||
                        !(material.density > 0.0) ||
                        !std::isfinite(material.dynamic_viscosity) ||
                        !(material.dynamic_viscosity > 0.0))
                        throw std::invalid_argument(
                            "production solver: the case material declares "
                            "fluid properties; density and "
                            "dynamic_viscosity must then both be finite and "
                            "positive instead of silently falling back to "
                            "solver defaults");
                    controls.density = material.density;
                    controls.kinematic_viscosity =
                        material.dynamic_viscosity / material.density;
                }
            } else if (case_setup.materials.size() > 1) {
                throw std::invalid_argument(
                    "production solver: case declares more than one material; "
                    "material selection is not supported");
            }
        } else {
            controls.algorithm = PressureVelocityAlgorithm::SIMPLE;
            controls.convection_scheme = ConvectionScheme::UPWIND;
        }
        controls.convergence.continuity_tolerance = 1e-12;
        controls.linear_max_iterations = 1000;
        controls.linear_tolerance = 1e-11;
        // Optional algebraic equilibration for ill-scaled imported meshes.
        // The linear dispatcher reconstructs the physical solution and checks
        // its *unscaled true residual* at the original tolerance. This does
        // not relax any of the nonlinear/physical qualification criteria.
        if (is_case_hdf5) {
            // Auxiliary pressure-space BCs for the PCD Schur model only.
            // A preconditioner experiment, never a modification of physical p/U BCs.
            const auto pcd_inflow = case_setup.numerics.raw_settings.find(
                "pcd_pressure_inflow_dirichlet_experimental");
            if (pcd_inflow != case_setup.numerics.raw_settings.end()) {
                if (pcd_inflow->second != "1" ||
                    controls.algorithm != PressureVelocityAlgorithm::COUPLED ||
                    controls.coupling.schur_model != CoupledSchurModel::PCD)
                    throw std::invalid_argument(
                        "pcd_pressure_inflow_dirichlet_experimental=1 "
                        "requires COUPLED with PCD Schur");
                controls.coupling.pcd_pressure_inflow_dirichlet_experimental = true;
                std::cout << "PCD auxiliary pressure-inflow Dirichlet: experimental "
                          << "(physical BCs and convergence gates unchanged)\n";
            }
            const auto schur_guard = case_setup.numerics.raw_settings.find(
                "schur_action_growth_guard_experimental");
            if (schur_guard != case_setup.numerics.raw_settings.end()) {
                if (schur_guard->second != "1" ||
                    controls.algorithm != PressureVelocityAlgorithm::COUPLED)
                    throw std::invalid_argument(
                        "schur_action_growth_guard_experimental=1 requires COUPLED");
                controls.coupling.schur_action_growth_guard_experimental = true;
                std::cout << "COUPLED six-Schur true-matrix action guard: "
                          << "1e10 (experimental; no fallback)\n";
            }
            const auto scaling = case_setup.numerics.raw_settings.find(
                "pressure_linear_scaling");
            if (scaling != case_setup.numerics.raw_settings.end()) {
                if (scaling->second == "none")
                    controls.pressure_linear_solver.scaling = LinearScalingModel::None;
                else if (scaling->second == "row_column")
                    controls.pressure_linear_solver.scaling = LinearScalingModel::RowColumn;
                else if (scaling->second == "symmetric_diagonal")
                    controls.pressure_linear_solver.scaling = LinearScalingModel::SymmetricDiagonal;
                else
                    throw std::invalid_argument(
                        "production solver: pressure_linear_scaling must be "
                        "none, row_column, or symmetric_diagonal");
                std::cout << "Pressure algebraic scaling: " << scaling->second
                          << " (physical true-residual gate retained)\n";
            }
            // Pressure and momentum are different linear operators. Keep the
            // case's global momentum selection, but allow an explicit SPD
            // pressure-only Krylov/AMG selection on stretched meshes.
            const auto pressure_krylov = case_setup.numerics.raw_settings.find(
                "pressure_krylov");
            if (pressure_krylov != case_setup.numerics.raw_settings.end()) {
                const auto& v = pressure_krylov->second;
                if (v == "auto") controls.pressure_linear_solver.krylov = KrylovModel::Auto;
                else if (v == "cg") controls.pressure_linear_solver.krylov = KrylovModel::CG;
                else if (v == "gmres") controls.pressure_linear_solver.krylov = KrylovModel::GMRES;
                else if (v == "fgmres") controls.pressure_linear_solver.krylov = KrylovModel::FGMRES;
                else if (v == "bicgstab") controls.pressure_linear_solver.krylov = KrylovModel::BiCGStab;
                else throw std::invalid_argument("invalid pressure_krylov; expected auto, cg, gmres, fgmres or bicgstab");
                std::cout << "Pressure Krylov: " << v << "\n";
            }
            const auto pressure_pc = case_setup.numerics.raw_settings.find(
                "pressure_preconditioner");
            if (pressure_pc != case_setup.numerics.raw_settings.end()) {
                const auto& v = pressure_pc->second;
                if (v == "auto") controls.pressure_linear_solver.preconditioner = PreconditionerModel::Auto;
                else if (v == "none") controls.pressure_linear_solver.preconditioner = PreconditionerModel::None;
                else if (v == "jacobi") controls.pressure_linear_solver.preconditioner = PreconditionerModel::Jacobi;
                else if (v == "ilu0") controls.pressure_linear_solver.preconditioner = PreconditionerModel::ILU0;
                else if (v == "anisotropic_line_amg") controls.pressure_linear_solver.preconditioner = PreconditionerModel::AnisotropicLineAMG;
                else if (v == "native_amg") controls.pressure_linear_solver.preconditioner = PreconditionerModel::NativeAMG;
                else if (v == "smoothed_aggregation_amg") controls.pressure_linear_solver.preconditioner = PreconditionerModel::SmoothedAggregationAMG;
                else throw std::invalid_argument("invalid pressure_preconditioner; expected auto, none, jacobi, ilu0, native_amg, anisotropic_line_amg or smoothed_aggregation_amg");
                std::cout << "Pressure preconditioner: " << v << "\n";
            }
            const auto interpolation = case_setup.numerics.raw_settings.find(
                "pressure_face_interpolation");
            if (interpolation != case_setup.numerics.raw_settings.end()) {
                if (interpolation->second == "arithmetic")
                    controls.pressure_projected_face_interpolation = false;
                else if (interpolation->second == "projected")
                    controls.pressure_projected_face_interpolation = true;
                else throw std::invalid_argument(
                    "invalid pressure_face_interpolation; expected arithmetic or projected");
                std::cout << "Pressure face interpolation: "
                          << interpolation->second << " (SIMPLE-only experimental)\n";
            }
            const auto pressure_tol = case_setup.numerics.raw_settings.find(
                "pressure_linear_tolerance");
            if (pressure_tol != case_setup.numerics.raw_settings.end()) {
                std::size_t consumed = 0;
                double value = 0.0;
                try { value = std::stod(pressure_tol->second, &consumed); }
                catch (const std::exception&) { throw std::invalid_argument("invalid pressure_linear_tolerance"); }
                if (consumed != pressure_tol->second.size() ||
                    !(value > 0.0) || !std::isfinite(value))
                    throw std::invalid_argument("invalid pressure_linear_tolerance");
                controls.pressure_linear_tolerance = value;
                std::cout << "Pressure true-residual target: " << value
                          << " (physical conservation gate unchanged)\n";
            }
            const auto refine_setting = case_setup.numerics.raw_settings.find(
                "pressure_iterative_refinement_steps");
            if (refine_setting != case_setup.numerics.raw_settings.end()) {
                const auto& value = refine_setting->second;
                if (value.size() != 1 || value[0] < '0' || value[0] > '8')
                    throw std::invalid_argument(
                        "pressure_iterative_refinement_steps must be an integer from 0 to 8");
                controls.pressure_linear_solver.residual_refinement_steps =
                    static_cast<std::size_t>(value[0] - '0');
                std::cout << "Pressure original-system refinement steps: "
                          << controls.pressure_linear_solver.residual_refinement_steps
                          << " (opt-in; original FP64 tolerance unchanged)\n";
            }
            const auto pressure_budget = case_setup.numerics.raw_settings.find(
                "pressure_linear_max_iterations");
            if (pressure_budget != case_setup.numerics.raw_settings.end()) {
                std::size_t consumed = 0;
                unsigned long value = 0;
                try { value = std::stoul(pressure_budget->second, &consumed, 10); }
                catch (const std::exception&) { throw std::invalid_argument("invalid pressure_linear_max_iterations"); }
                if (consumed != pressure_budget->second.size() || value < 1 || value > 100000)
                    throw std::invalid_argument("invalid pressure_linear_max_iterations");
                controls.pressure_linear_max_iterations = value;
                std::cout << "Pressure max iterations: " << value << "\n";
            }
            const auto restart = case_setup.numerics.raw_settings.find(
                "pressure_gmres_restart");
            if (restart != case_setup.numerics.raw_settings.end()) {
                std::size_t consumed = 0;
                int value = 0;
                try {
                    value = std::stoi(restart->second, &consumed, 10);
                } catch (const std::exception&) {
                    throw std::invalid_argument(
                        "production solver: pressure_gmres_restart must be an integer 1..2048");
                }
                if (consumed != restart->second.size() || value < 1 || value > 2048)
                    throw std::invalid_argument(
                        "production solver: pressure_gmres_restart must be an integer 1..2048");
                controls.pressure_linear_solver.gmres_restart = value;
                std::cout << "Pressure GMRES restart: " << value
                          << " (linear tolerance unchanged)\n";
            }
        }
        controls.acceleration.adaptive_linear_tolerance =
            options.adaptive_convergence;
        controls.acceleration.adaptive_pressure_correctors =
            options.adaptive_convergence;
        // N10 diagnostic mode: expose the solver's existing per-iteration
        // convergence metrics without changing any convergence criterion.
        // This is intentionally diagnostic-only; the production result still
        // requires the authoritative convergence gates below.
        controls.diagnostics.iteration_trace = true;
        controls.diagnostics.iteration_trace_frequency = 1;
        if (is_case_hdf5) {
            const auto& raw=case_setup.numerics.raw_settings;
            const auto flag=[&](const char* name) {
                const auto it=raw.find(name);
                if (it==raw.end()) return false;
                if (it->second!="true" && it->second!="false")
                    throw std::invalid_argument(std::string(name)+" must be true or false");
                return it->second=="true";
            };
            const auto number=[&](const char* name, double fallback) {
                const auto it=raw.find(name);
                if (it==raw.end()) return fallback;
                std::size_t used=0;
                double result=0.0;
                try { result=std::stod(it->second,&used); }
                catch (const std::exception&) {
                    throw std::invalid_argument(std::string(name)+" must be a number");
                }
                if (used!=it->second.size() || !std::isfinite(result))
                    throw std::invalid_argument(std::string(name)+" must be finite");
                return result;
            };
            controls.diagnostics.flux_patch_audit=flag("v17_flux_patch_audit");
            controls.pressure_compensated_assembly=flag("pressure_compensated_assembly");
            controls.pressure_compensated_refinement_diagnostic=flag("pressure_compensated_refinement_diagnostic");
            controls.v18_full_flux_pressure_oracle=flag("v18_full_flux_pressure_oracle");
            controls.v18_sparse_flux_pressure_jacobian=flag("v18_sparse_flux_pressure_jacobian");
            controls.v19_full_pressure_flux_step=flag("v19_full_pressure_flux_step");
            // V21 opt-in for the existing hysteretic relaxation controller.
            // It uses momentum/flux coupling, not projected continuity, for
            // alpha_p whenever V19 full-flux SIMPLE is active.
            controls.adaptive_relaxation.enabled=flag("v21_adaptive_relaxation");
            if (controls.adaptive_relaxation.enabled &&
                !controls.v19_full_pressure_flux_step)
                throw std::invalid_argument(
                    "V21 coupling adaptation requires V19 full pressure flux step");
            if (controls.adaptive_relaxation.enabled)
                std::cout << "V21 momentum/flux adaptive relaxation: EXPERIMENTAL\n";
            if (controls.v19_full_pressure_flux_step)
                std::cout << "V19 full pressure flux step, relaxed physical pressure: "
                          << "EXPERIMENTAL, not qualified\n";
            if (controls.v18_full_flux_pressure_oracle)
                std::cout << "V18 reference pressure Schur Jacobian: ENABLED "
                          << "(small-mesh diagnostic only, never qualified)\n";
            controls.pressure_continuity_line_search.enabled=
                flag("v17_pressure_line_search");
            if (controls.pressure_continuity_line_search.enabled)
                controls.pressure_continuity_line_search.max_trials=12;
            auto& pseudo=controls.pseudo_momentum;
            pseudo.enabled=flag("v17_pseudo_momentum");
            pseudo.initial_cfl=number("v17_pseudo_initial_cfl",pseudo.initial_cfl);
            pseudo.minimum_cfl=number("v17_pseudo_min_cfl",pseudo.minimum_cfl);
            pseudo.maximum_cfl=number("v17_pseudo_max_cfl",pseudo.maximum_cfl);
            pseudo.growth=number("v17_pseudo_growth",pseudo.growth);
            pseudo.shrink=number("v17_pseudo_shrink",pseudo.shrink);
            validate_v17_pseudo_momentum_controls(pseudo);
            if (pseudo.enabled)
                std::cout << "V17 experimental steady pseudo-momentum: initial CFL="
                          << pseudo.initial_cfl << " range=[" << pseudo.minimum_cfl
                          << "," << pseudo.maximum_cfl << "] (not a physical transient)\n";
            if (controls.diagnostics.flux_patch_audit)
                std::cout << "V17 authoritative Rhie-Chow patch flux audit enabled\n";
        }
        controls.pressure_reference_cell = 0;
        controls.pressure_reference_value = 0.0;
        controls.probes = options.probes;
        if (options.probe_stream) {
            controls.probe_callback = [&options](
                const IncompressibleProbeSample& sample) {
                const auto probe = std::find_if(
                    options.probes.begin(), options.probes.end(),
                    [&sample](const IncompressiblePointProbe& candidate) {
                        return candidate.name == sample.name;
                    });
                if (probe == options.probes.end())
                    throw std::logic_error(
                        "probe stream received a sample without a configured probe");
                std::cout << "CFDX_PROBE_SAMPLE_V1\t"
                          << encode_probe_name_hex(sample.name) << '\t'
                          << probe_field_cli(probe->field) << '\t'
                          << probe_field_unit(probe->field) << '\t'
                          << sample.iteration << '\t'
                          << format_probe_csv_value(sample.time) << '\t'
                          << format_probe_csv_value(sample.value) << '\n'
                          << std::flush;
            };
        }

        // The declared initialization strategy is applied by the solver, before
        // the first residual is formed. An explicit --restart request therefore
        // overrides the case-derived uniform state.
        controls.initialization = initialization;
        if (!options.restart.empty())
            controls.initialization.mode = InitializationMode::Restart;
        if (options.initial_face_mass_flux) {
#if defined(CFDX_HAS_HDF5) && CFDX_HAS_HDF5
            const auto fail = [](const char* why) -> void {
                throw std::invalid_argument(
                    std::string("hybrid /fields/face/mdot: ") + why);
            };
            hid_t file=H5Fopen(mesh_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            if (file<0) fail("cannot open HDF5");
            hid_t ds=H5Dopen2(file,"/fields/face/mdot",H5P_DEFAULT);
            if (ds<0) { H5Fclose(file); fail("missing dataset"); }
            hid_t space=H5Dget_space(ds);
            hsize_t dims[1]={0};
            if (space<0 || H5Sget_simple_extent_ndims(space)!=1 ||
                H5Sget_simple_extent_dims(space,dims,nullptr)<0 ||
                dims[0]!=mesh.n_faces()) {
                if (space>=0) H5Sclose(space);
                H5Dclose(ds); H5Fclose(file);
                fail("invalid face field shape");
            }
            controls.initial_face_mass_flux.resize(mesh.n_faces());
            const auto read_result=H5Dread(ds,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,
                                        H5P_DEFAULT,controls.initial_face_mass_flux.data());
            H5Sclose(space); H5Dclose(ds); H5Fclose(file);
            if (read_result<0) fail("cannot read dataset");
#else
            throw std::invalid_argument("hybrid face flux requires HDF5-enabled CFDX");
#endif
        }

        if (!options.pressure_replay_prefix.empty()) {
            controls.pressure_replay_prefix = options.pressure_replay_prefix;
            controls.pressure_replay_iterations = options.pressure_replay_iterations;
            const std::filesystem::path capture_path(options.pressure_replay_prefix);
            if (!capture_path.parent_path().empty())
                std::filesystem::create_directories(capture_path.parent_path());
            std::cout << "Pressure CSR replay capture enabled: "
                      << options.pressure_replay_prefix << " ("
                      << options.pressure_replay_iterations.size()
                      << " selected iterations; retries can overwrite a capture)\n";
        }

        // Own all SST runtime objects in main for the complete solver call.
        // The default production selector stays fail-closed; this branch can
        // only be entered through the explicit --experimental-sst opt-in.
        std::optional<SSTProductionState> experimental_sst_state;
        SSTProductionRuntime experimental_sst_runtime;
        DatRestartFields solver_restart_fields{};
        // The restored Rhie-Chow face mass flux is part of the iteration
        // state: it is read from the DAT into this caller-owned field and
        // consumed by the solver after its own initialization transaction.
        std::optional<Field<double, Location::FACE>> restart_face_flux;
        if (!options.restart.empty() &&
            (options.experimental_sst || controls.v19_full_pressure_flux_step)) {
            // DAT v4 checkpoints created before the optional face_fields
            // section remain valid for established non-SST restart paths.
            // Experimental SST and the V19 split pressure/flux policy need
            // the authoritative mdot to replay the same momentum predictor.
            // Do not silently reconstruct mdot from U and relaxed p here.
            restart_face_flux.emplace(mesh.n_faces(), "mdot", "kg/s", 1);
            solver_restart_fields.face_mass_flux = &*restart_face_flux;
        }
        if (options.experimental_sst) {
            if (production_s2s_radiation_requested(case_setup))
                throw std::invalid_argument(
                    "production solver: experimental SST-radiation coupling is "
                    "not implemented; cannot assert a physical energy balance");
            if (production_energy_requested(case_setup))
                (void)require_production_sst_turbulent_prandtl(case_setup);

            TurbulenceTransportControls sst_controls;
            sst_controls.model = TurbulenceModel::SST;
            sst_controls.density = controls.density;
            sst_controls.molecular_viscosity = controls.kinematic_viscosity;
            const auto damping = case_setup.numerics.raw_settings.find(
                "sst_picard_damping");
            if (damping != case_setup.numerics.raw_settings.end()) {
                std::size_t used = 0;
                double value=0.0;
                try { value=std::stod(damping->second, &used); }
                catch (const std::exception&) {
                    throw std::invalid_argument("invalid sst_picard_damping");
                }
                if (used != damping->second.size() ||
                    !std::isfinite(value) || !(value>0.0 && value<=1.0))
                    throw std::invalid_argument("sst_picard_damping must be in (0,1]");
                sst_controls.sst_picard_damping = value;
                std::cout << "SST Picard damping: " << value
                          << " (equation and convergence gates unchanged)\n";
            }
            const auto relative_gate = case_setup.numerics.raw_settings.find(
                "sst_scale_aware_convergence");
            if (relative_gate != case_setup.numerics.raw_settings.end()) {
                if (relative_gate->second != "true" && relative_gate->second != "false")
                    throw std::invalid_argument(
                        "sst_scale_aware_convergence must be true or false");
                sst_controls.sst_scale_aware_convergence =
                    relative_gate->second == "true";
                std::cout << "SST scale-aware convergence: "
                          << (sst_controls.sst_scale_aware_convergence ? "enabled" : "disabled")
                          << " (true equation residual, linear and physical gates retained)\n";
            }
            const auto sst_newton_setting = case_setup.numerics.raw_settings.find(
                "sst_omega_newton_sink");
            if (sst_newton_setting != case_setup.numerics.raw_settings.end()) {
                if (sst_newton_setting->second != "true" &&
                    sst_newton_setting->second != "false")
                    throw std::invalid_argument("sst_omega_newton_sink must be true or false");
                sst_controls.sst_omega_newton_sink = sst_newton_setting->second == "true";
                std::cout << "SST omega Newton sink: "
                          << (sst_controls.sst_omega_newton_sink ? "enabled" : "disabled") << "\n";
            }
            const auto sst_gs_setting = case_setup.numerics.raw_settings.find(
                "sst_sequential_omega_first");
            if (sst_gs_setting != case_setup.numerics.raw_settings.end()) {
                if (sst_gs_setting->second != "true" &&
                    sst_gs_setting->second != "false")
                    throw std::invalid_argument("sst_sequential_omega_first must be true or false");
                sst_controls.sst_sequential_omega_first = sst_gs_setting->second == "true";
                std::cout << "SST omega-first Picard: "
                          << (sst_controls.sst_sequential_omega_first ? "enabled" : "disabled") << "\n";
            }
            const auto sst_wall_setting = case_setup.numerics.raw_settings.find(
                "sst_wall_omega_factor");
            if (sst_wall_setting != case_setup.numerics.raw_settings.end()) {
                if (sst_wall_setting->second != "6" &&
                    sst_wall_setting->second != "60")
                    throw std::invalid_argument("sst_wall_omega_factor must be 6 or 60");
                sst_controls.sst_wall_omega_factor = std::stod(sst_wall_setting->second);
                std::cout << "SST wall omega face factor: "
                          << sst_controls.sst_wall_omega_factor
                          << " (not an OpenFOAM wall function)\n";
            }
            const auto sst_linear_relaxation = case_setup.numerics.raw_settings.find(
                "sst_transport_linear_relaxation");
            if (sst_linear_relaxation != case_setup.numerics.raw_settings.end()) {
                std::size_t parsed=0;
                double alpha=0.0;
                try { alpha=std::stod(sst_linear_relaxation->second,&parsed); }
                catch (const std::exception&) {
                    throw std::invalid_argument("invalid sst_transport_linear_relaxation");
                }
                if (parsed!=sst_linear_relaxation->second.size() ||
                    !std::isfinite(alpha) || alpha<=0.0 || alpha>1.0)
                    throw std::invalid_argument("sst_transport_linear_relaxation must be in (0,1]");
                sst_controls.sst_linear_relaxation=alpha;
            }
            std::cout << "SST scalar post-solve relaxation: "
                      << sst_controls.sst_linear_relaxation << " (default 0.7)\n";
            for (const auto& setting : {std::pair<const char*,bool*>{
                     "sst_compensated_residual_gate", &sst_controls.sst_compensated_residual_gate},
                     {"sst_nonlinear_line_search", &sst_controls.sst_nonlinear_line_search}}) {
                const auto found=case_setup.numerics.raw_settings.find(setting.first);
                if (found != case_setup.numerics.raw_settings.end()) {
                    if (found->second != "true" && found->second != "false")
                        throw std::invalid_argument(std::string(setting.first)+" must be true or false");
                    *setting.second = found->second == "true";
                    std::cout << "Experimental " << setting.first << "=" << found->second
                              << " (global physical gates unchanged)\n";
                }
            }
            const auto molecular_wall_diffusion = case_setup.numerics.raw_settings.find(
                "sst_molecular_wall_diffusion");
            if (molecular_wall_diffusion != case_setup.numerics.raw_settings.end()) {
                if (molecular_wall_diffusion->second != "true" &&
                    molecular_wall_diffusion->second != "false")
                    throw std::invalid_argument("sst_molecular_wall_diffusion must be true or false");
                sst_controls.sst_molecular_wall_diffusion =
                    molecular_wall_diffusion->second == "true";
            }
            std::cout << "SST molecular boundary diffusion: "
                      << (sst_controls.sst_molecular_wall_diffusion ? "experimental" : "legacy")
                      << " (only fixed-value wall faces; no wall function)\n";
            validate_turbulence_controls(sst_controls);

            const auto boundary_contract =
                make_sst_runtime_boundary_contract(case_setup);
            std::vector<std::string> wall_patch_names;
            for (std::size_t patch_index = 0;
                 patch_index < mesh.boundary().n_patches(); ++patch_index) {
                const auto& patch = mesh.boundary().patch(patch_index);
                const auto* bc = case_setup.find_boundary(patch.name);
                if (bc == nullptr)
                    throw std::invalid_argument(
                        "production solver: SST runtime has no case boundary metadata for patch '" +
                        patch.name + "'");
                const bool declared_wall =
                    bc->type == BCType::WALL ||
                    bc->value_type == BCValueType::WALL_NO_SLIP ||
                    bc->value_type == BCValueType::WALL_SLIP ||
                    bc->value_type == BCValueType::WALL_THERMAL;
                if (declared_wall) {
                    if (patch.type != PatchType::WALL)
                        throw std::invalid_argument(
                            "production solver: SST case declares wall patch '" +
                            patch.name + "' but mesh patch is not WALL");
                    wall_patch_names.push_back(patch.name);
                }
            }
            if (wall_patch_names.empty())
                throw std::invalid_argument(
                    "production solver: experimental SST requires at least one "
                    "explicit wall patch for wall-distance reconstruction");

            std::map<std::string, std::array<double, 3>> fixed_velocity_values;
            for (const auto& [name, bc] : ubc) {
                if (bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                    fixed_velocity_values.emplace(
                        name, std::array<double, 3>{bc.value.x, bc.value.y, bc.value.z});
            }

            const FvGeometry sst_geometry = build_fv_geometry(mesh);
            if (!options.restart.empty()) {
                experimental_sst_state.emplace(mesh.n_cells());
                const auto wall_distance = compute_mesh_cell_wall_distance(
                    mesh, sst_geometry, wall_patch_names);
                // Stage and validate U, p, k, omega and derived SST fields
                // before changing live fields. The solver re-reads the same
                // checkpoint through its normal initialization transaction.
                (void)restore_sst_production_state_from_dat(
                    options.restart.string(), mesh, U, p,
                    *experimental_sst_state, wall_distance.distance,
                    sst_controls, fixed_velocity_values);
            } else {
                if (!case_setup.has_initial_condition ||
                    case_setup.initial_condition.scalar_fields.find("k") ==
                        case_setup.initial_condition.scalar_fields.end() ||
                    case_setup.initial_condition.scalar_fields.find("omega") ==
                        case_setup.initial_condition.scalar_fields.end())
                    throw std::invalid_argument(
                        "production solver: experimental SST requires explicit "
                        "initial_condition.scalar_fields.k and .omega, or a DAT restart");
                experimental_sst_state.emplace(
                    initialize_sst_production_state_from_mesh(
                        mesh, case_setup, sst_controls, U,
                        fixed_velocity_values, wall_patch_names));
                if (options.initial_sst_fields) {
                    std::vector<std::pair<std::string,ScalarCellField>> scalars;
                    if (!read_scalar_fields_hdf5(mesh_path,scalars))
                        throw std::invalid_argument("hybrid SST fields cannot be read from HDF5");
                    bool found_k=false, found_omega=false;
                    for (const auto& [name,field]:scalars) {
                        if (name!="k" && name!="omega") continue;
                        if (field.size()!=mesh.n_cells() || field.dimension()!=1)
                            throw std::invalid_argument("hybrid SST scalar shape mismatch");
                        bool& found=name=="k" ? found_k : found_omega;
                        if (found)
                            throw std::invalid_argument("duplicate hybrid SST scalar");
                        found=true;
                        for (std::size_t c=0;c<mesh.n_cells();++c) {
                            const double value=field(c);
                            if (!std::isfinite(value) ||
                                !(value>=(name=="k" ? sst_controls.k_min : sst_controls.omega_min)))
                                throw std::invalid_argument("hybrid SST field contains an invalid or subfloor value");
                            if (name=="k") experimental_sst_state->k(c)=value;
                            else experimental_sst_state->omega(c)=value;
                        }
                    }
                    if (!found_k || !found_omega)
                        throw std::invalid_argument("hybrid SST requires k and omega datasets");
                    refresh_sst_production_derived_state_from_mesh(
                        *experimental_sst_state,mesh,sst_controls,U,
                        fixed_velocity_values,wall_patch_names);
                    std::cout << "HYBRID_INITIAL_SST_FIELDS imported nonuniform k and omega"
                                 " (derived nut refreshed; NOT qualified)\n";
                }
            }

            experimental_sst_runtime.state = &*experimental_sst_state;
            experimental_sst_runtime.transport = sst_controls;
            experimental_sst_runtime.k_boundary_conditions = boundary_contract.k;
            experimental_sst_runtime.omega_boundary_conditions = boundary_contract.omega;
            experimental_sst_runtime.wall_patch_names = wall_patch_names;
            experimental_sst_runtime.wall_omega_sst_patches =
                boundary_contract.wall_omega_sst_patches;
            // Larger SST inner budgets may be needed for high-Re separated
            // flows. This setting changes only the work budget, never the
            // fixed 1e-8 acceptance tolerance or conservation gates.
            std::size_t sst_inner_budget = 100;
            const auto budget_it = case_setup.numerics.raw_settings.find(
                "sst_transport_max_iterations");
            if (budget_it != case_setup.numerics.raw_settings.end()) {
                const auto& declared = budget_it->second;
                std::size_t end = 0;
                unsigned long parsed = 0;
                try {
                    parsed = std::stoul(declared, &end, 10);
                } catch (const std::exception&) {
                    throw std::invalid_argument(
                        "SST transport budget must be a decimal integer in [10, 2000]");
                }
                if (end != declared.size() || parsed < 10 || parsed > 2000)
                    throw std::invalid_argument(
                        "SST transport budget must be a decimal integer in [10, 2000]");
                sst_inner_budget = static_cast<std::size_t>(parsed);
            }
            experimental_sst_runtime.max_transport_iterations = sst_inner_budget;
            experimental_sst_runtime.transport_tolerance = 1e-8;
            const auto provisional = case_setup.numerics.raw_settings.find(
                "sst_outer_provisional");
            if (provisional != case_setup.numerics.raw_settings.end()) {
                if (provisional->second == "enabled")
                    experimental_sst_runtime.allow_bounded_provisional_transport = true;
                else if (provisional->second != "disabled")
                    throw std::invalid_argument(
                        "sst_outer_provisional must be enabled or disabled");
                std::cout << "SST bounded provisional outer coupling: "
                          << provisional->second
                          << " (final SST and conservation gates UNCHANGED)\n";
            }
            std::cout << "SST transport max inner iterations "
                      << sst_inner_budget << " (unchanged tolerance 1e-8)\n";
            if (!boundary_contract.wall_omega_sst_patches.empty()) {
                experimental_sst_runtime.omega_face_values =
                    make_sst_wall_omega_face_values(
                        mesh, sst_geometry,
                        boundary_contract.wall_omega_sst_patches, sst_controls);
                require_supported_sst_runtime_boundaries(
                    boundary_contract, mesh, sst_geometry,
                    &experimental_sst_runtime.omega_face_values);
            } else {
                require_supported_sst_runtime_boundaries(
                    boundary_contract, mesh, sst_geometry, nullptr);
            }
            controls.sst_runtime = &experimental_sst_runtime;
            solver_restart_fields.k = &experimental_sst_state->k;
            solver_restart_fields.omega = &experimental_sst_state->omega;
            std::cerr
                << "WARNING: experimental SST runtime explicitly enabled. "
                   "This path is not physically qualified and does not enable "
                   "the general RANS selector.\n";
        }

        std::cout << "Initialization strategy "
                  << cfdx::physics::to_string(controls.initialization.mode) << "\n";

        controls.initialized_state_callback =
            [&](const Mesh& initialized_mesh,
                const Field<double, Location::CELL>& initialized_u,
                const Field<double, Location::CELL>& initialized_p) {
                (void)initialized_p;
                std::ostringstream initialized_json;
                initialized_json << "{\n  \"field\": \"U\",\n"
                                << "  \"location\": \"cell\",\n"
                                << "  \"stage\": \"after_initialization_before_first_residual\",\n"
                                << "  \"cell_count\": " << initialized_mesh.n_cells() << ",\n"
                                << "  \"values\": [";
                for (std::size_t cell = 0; cell < initialized_u.size(); ++cell) {
                    if (cell != 0) initialized_json << ",";
                    initialized_json << "\n    ["
                                    << json_number(initialized_u(cell, 0)) << ", "
                                    << json_number(initialized_u(cell, 1)) << ", "
                                    << json_number(initialized_u(cell, 2)) << "]";
                }
                if (initialized_u.size() != 0) initialized_json << "\n  ";
                initialized_json << "]\n}\n";
                write_text_file_atomic(
                    options.output_dir / "initialized_velocity_field.json",
                    initialized_json.str());
            };

        controls.iteration_output_callback =
            [&](std::size_t iteration, double time, const Mesh& state_mesh,
                const Field<double, Location::CELL>& state_u,
                const Field<double, Location::CELL>& state_p) {
                ScalarCellField pressure(state_p);
                ScalarCellField ux(state_u.size(), "u_x", "m/s");
                ScalarCellField uy(state_u.size(), "u_y", "m/s");
                ScalarCellField uz(state_u.size(), "u_z", "m/s");
                for (std::size_t c = 0; c < state_u.size(); ++c) {
                    ux(c) = state_u(c, 0);
                    uy(c) = state_u(c, 1);
                    uz(c) = state_u(c, 2);
                }
                const auto path = options.output_dir /
                    ("result_" + std::to_string(iteration) + ".vtu");
                // V21: diagnostics for the actual accepted coupled SST state.
                // Omitting these cell fields made it impossible to assess
                // upstream boundary layers, nut, k and omega in the same
                // snapshot as the velocity. No extra arrays for laminar cases.
                std::map<std::string, ScalarCellField> cell_fields{
                    {"p", pressure}, {"u_x", ux}, {"u_y", uy}, {"u_z", uz}};
                if (experimental_sst_state) {
                    validate_sst_production_state(*experimental_sst_state,
                                                  state_mesh.n_cells());
                    cell_fields.emplace("k", ScalarCellField(experimental_sst_state->k));
                    cell_fields.emplace("omega", ScalarCellField(experimental_sst_state->omega));
                    cell_fields.emplace("nu_t", ScalarCellField(experimental_sst_state->nu_t));
                    cell_fields.emplace("wall_distance", ScalarCellField(experimental_sst_state->wall_distance));
                    cell_fields.emplace("sst_strain_rate", ScalarCellField(experimental_sst_state->strain_rate));
                }
                VtuWriter writer;
                if (!writer.write(path.string(), state_mesh,
                                  cell_fields, {}, {}, time, iteration, true))
                    return false;
                std::cout << "Iteration " << iteration
                          << " Time = " << time << "\n";
                return g_stop_requested == 0;
            };

        const bool energy_requested = is_case_hdf5 &&
            production_energy_requested(case_setup);
        Field<double, Location::CELL> temperature(
            mesh.n_cells(), "temperature", "K", 1);
        const bool thermal_restart = energy_requested && !options.restart.empty();
        if (thermal_restart) {
            // Validate the complete thermal checkpoint before running the
            // carrier. The carrier's own initialization subsequently consumes
            // U/p and its iteration metadata; no live fields are changed here.
            auto candidate_velocity = U;
            auto candidate_pressure = p;
            DatRestartFields fields;
            fields.temperature = &temperature;
            read_dat_restart_fields(options.restart.string(), mesh,
                candidate_velocity, candidate_pressure, fields);
            for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
                if (!std::isfinite(temperature(cell)) || !(temperature(cell) > 0.0))
                    throw std::invalid_argument(
                        "production solver: restarted temperature must be finite and positive (kelvin)");
        }

        IncompressibleSolveResult result;
        bool buoyancy_ran = false;
        ProductionBuoyancyOuterLoopResult buoyancy_result;
        if (is_case_hdf5 && production_buoyancy_requested(case_setup)) {
            // Coupled flow-energy buoyancy iteration (issue #885 step 4):
            // each outer pass rebuilds the cellwise Boussinesq buoyancy
            // source from the relaxed temperature, re-solves the carrier
            // (warm-started from the previous pass), then solves the
            // energy response - with the energy-S2S fixed point when
            // radiation is declared - on the new authoritative mass
            // flux. The blocking verdict is the coupled fixed point,
            // never an intermediate Picard pass.
            buoyancy_result = run_production_buoyant_outer_loop(
                mesh, ubc, pbc, controls, case_setup, U, p,
                options.restart.string(), thermal_restart ? &temperature : nullptr,
                options.experimental_sst ? &experimental_sst_runtime : nullptr);
            result = buoyancy_result.flow;
            buoyancy_ran = true;
            std::cout << "Buoyant flow-energy outer loop "
                      << (buoyancy_result.converged ? "converged"
                                                    : "NOT converged")
                      << " outer iterations "
                      << buoyancy_result.outer_iterations
                      << " fixed point temperature "
                      << buoyancy_result.max_temperature_change_normalized
                      << " velocity "
                      << buoyancy_result.max_velocity_change_normalized
                      << " pressure "
                      << buoyancy_result.max_pressure_change_normalized
                      << "\n";
        } else if (options.experimental_momentum_first_continuation) {
            if (!options.experimental_sst || controls.sst_runtime == nullptr)
                throw std::invalid_argument(
                    "--experimental-momentum-first-continuation requires --experimental-sst");
            if (!options.restart.empty())
                throw std::invalid_argument(
                    "--experimental-momentum-first-continuation does not accept a DAT restart");
            ContinuationControls schedule;
            schedule.enabled = true;
            schedule.momentum_only_before_final = true;
            schedule.upwind_before_final = true;
            schedule.initial_step = 0.5;
            schedule.maximum_step = 0.5;
            schedule.minimum_step = 0.0625;
            const auto staged = solve_steady_incompressible_continuation(
                mesh, U, p, ubc, pbc, controls, schedule);
            for (const auto& stage : staged.stages)
                std::cout << "V25_CONTINUATION_STAGE target=" << stage.parameter
                          << " attempts=" << stage.attempts
                          << " converged=" << (stage.converged ? 1 : 0)
                          << " full_physics=" <<
                              (continuation_final_physical_stage(stage.parameter) ? 1 : 0)
                          << "\n";
            result = staged.final_solver_result;
            // A converged intermediate momentum-only stage must never
            // masquerade as a converged final SST result.
            if (!staged.converged || staged.final_parameter != 1.0)
                result.converged = false;
        } else {
            result = solve_steady_incompressible(
                mesh, U, p, ubc, pbc, controls, options.restart.string(),
                solver_restart_fields);
        }
        execution_result = result;
        const auto execution = options.output_dir / "execution.json";
        const auto flow_exit_code = execution_exit_code_for_solver(result.converged);
        // Persist the solver outcome before any secondary artifact is written.
        // This guarantees that a restart/convergence/probe write failure still
        // leaves an execution artifact which the catch block can finalize with
        // the actual process error code (2).
        write_execution_summary(execution.string(), ExecutionSummary{
            flow_exit_code,
            result.converged,
            result.iterations,
            static_cast<int>(result.convergence_status),
            result.convergence_reason,
            false,
            false});
        if (options.wall_forces_quadratic) {
            if (!result.converged)
                throw std::runtime_error(
                    "--wall-forces-quadratic requires a converged flow solution");
            write_quadratic_wall_forces(
                mesh, U, p, ubc, pbc,
                controls.density * controls.kinematic_viscosity,
                options.output_dir);
        }
        if (!options.restart.empty())
            std::cout << "Restart checkpoint " << options.restart.string() << "\n";

        // T2 production energy stage (capability #885 stage 2, issues
        // #880 / #118): once the frozen flow has converged, couple it
        // one-way into the steady energy equation on the solver's
        // authoritative conservative face mass flux. The stage runs only
        // for a converged carrier flow: an energy field carried by a
        // diverged flow is not evidence, and the global exit code below
        // is flow AND energy. The energy model was already validated
        // fail-closed right after case load.
        bool energy_ran = false;
        ProductionEnergyStageResult energy_result;
        ProductionS2SOuterLoopResult s2s_outer_result;
        bool s2s_ran = false;
        if (buoyancy_ran) {
            // The buoyant outer loop already produced the coupled energy
            // response (and, when declared, the energy-S2S fixed point)
            // on the final carrier solve: no frozen-flow pass is run
            // again. The artifacts below reuse the loop's recorded
            // responses; they are gated on the loop having completed at
            // least one full outer pass, exactly like the frozen-flow
            // path.
            energy_ran = buoyancy_result.energy_ran;
            s2s_ran = buoyancy_result.s2s_ran;
            energy_result = buoyancy_result.energy;
            if (s2s_ran)
                s2s_outer_result = buoyancy_result.s2s;
            temperature = buoyancy_result.temperature;
        } else if (energy_requested && result.converged) {
            const auto energy_geometry = build_fv_geometry(mesh);
            // Initial temperature: the declared initial condition, else
            // the case reference temperature. Zero means "not declared"
            // (stage-1 convention) and falls back to ref_temperature.
            double initial_temperature = thermal_restart
                ? temperature(0) : case_setup.ref_temperature;
            if (!thermal_restart && case_setup.has_initial_condition &&
                case_setup.initial_condition.temperature != 0.0)
                initial_temperature =
                    case_setup.initial_condition.temperature;
            // Absolute temperature contract: a non-finite or non-positive
            // value cannot seed the energy stage (whose kernel also rejects
            // it), so fail the run here with the solver-level message.
            if (!std::isfinite(initial_temperature) ||
                !(initial_temperature > 0.0))
                throw std::invalid_argument(
                    "production solver: initial temperature must be finite "
                    "and positive (kelvin)");
            if (!thermal_restart)
                for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
                    temperature(cell) = initial_temperature;
            if (production_s2s_radiation_requested(case_setup)) {
                // External energy-radiation iteration (issue #885 step
                // 3): energy pass on the frozen flow, S2S radiosity from
                // the resulting wall temperatures, net radiative fluxes
                // re-injected with controlled under-relaxation, iterated
                // until the unrelaxed radiative fixed point and the
                // wall-temperature stationarity both hold. The verdict
                // is the coupled fixed point on the frozen flow, never
                // an intermediate Picard pass.
                s2s_outer_result = run_production_s2s_energy_outer_loop(
                    mesh, energy_geometry, case_setup,
                    result.authoritative_mass_flux, temperature);
                energy_result = s2s_outer_result.energy;
                energy_ran = true;
                s2s_ran = true;
            } else {
                // Couple the final *accepted* SST viscosity to energy.
                // Keep this vector alive until the assembled energy solve
                // and its independent face-flux audit have finished.
                std::vector<double> lambda_eff;
                if (experimental_sst_state)
                    lambda_eff = build_production_sst_thermal_conductivity(
                        case_setup, *experimental_sst_state, mesh.n_cells(),
                        require_production_sst_turbulent_prandtl(case_setup));
                energy_result = run_production_energy_stage(
                    mesh, energy_geometry, case_setup,
                    result.authoritative_mass_flux, temperature, {},
                    lambda_eff.empty() ? nullptr : &lambda_eff);
                energy_ran = true;
            }
        } else if (energy_requested) {
            std::cout
                << "Energy stage skipped: carrier flow did not converge\n";
        }
        // Energy artifacts are written from the recorded responses of
        // whichever path ran (frozen-flow stage or buoyant outer loop),
        // gated on a completed energy response exactly like the
        // frozen-flow path.
        if (energy_ran && temperature.size() == mesh.n_cells()) {
            write_text_file_atomic(
                options.output_dir / "energy_convergence.json",
                energy_convergence_json(energy_result));
            write_text_file_atomic(
                options.output_dir / "energy_summary.json",
                energy_summary_json(energy_result));
            VtuWriter energy_writer;
            const auto energy_vtu =
                options.output_dir / "energy_result_0.vtu";
            // The VTU is written through a temporary sibling and renamed,
            // matching the JSON artifacts: a reader never observes a
            // half-written energy result, whatever the interruption.
            const auto energy_vtu_tmp =
                options.output_dir / "energy_result_0.vtu.tmp";
            if (!energy_writer.write(energy_vtu_tmp.string(), mesh,
                                     {{"temperature", temperature}},
                                     {}, {}, 0.0, 0, true))
                throw std::runtime_error(
                    "production solver: energy result VTU write failed");
            std::filesystem::rename(energy_vtu_tmp, energy_vtu);
            std::cout << "Energy stage "
                      << (energy_result.converged ? "converged"
                                                   : "NOT converged")
                      << " iterations " << energy_result.iterations
                      // Blocking-gate diagnostics, not the legacy
                      // normalized audit: global_balance_normalized and
                      // local_balance_normalized (gross-scale) plus
                      // max_cell_residual_normalized (per-cell scale) are
                      // the criteria converged actually requires.
                      << " conservation global="
                      << energy_result.global_balance_normalized
                      << " local=" << energy_result.local_balance_normalized
                      << " local_cell="
                      << energy_result.max_cell_residual_normalized
                      << "\n";
            std::cout << "Energy gates linear="
                      << (energy_result.linear_solver_converged ? 1 : 0)
                      << " global_balance="
                      << (energy_result.global_energy_balance_passed ? 1 : 0)
                      << " local_balance="
                      << (energy_result.local_energy_balance_passed ? 1 : 0)
                      << " temperature_admissibility="
                      << (energy_result.temperature_admissibility_passed ? 1 : 0)
                      << " nonpositive_temperature="
                      << energy_result.nonpositive_temperature
                      << " audit_scale=" << energy_result.audit_power_scale
                      << "\n";
            std::cout << "Energy temperatures ["
                      << energy_result.min_temperature << ", "
                      << energy_result.max_temperature << "] K\n";
            std::cout << "Energy convergence history "
                      << (options.output_dir / "energy_convergence.json")
                             .string()
                      << "\n";
            std::cout << "Energy summary "
                      << (options.output_dir / "energy_summary.json")
                             .string()
                      << "\n";
            std::cout << "Energy result " << energy_vtu.string() << "\n";
            if (s2s_ran) {
                write_text_file_atomic(
                    options.output_dir / "radiation_s2s.json",
                    radiation_s2s_json(s2s_outer_result));
                std::cout << "S2S radiation outer loop "
                          << (s2s_outer_result.converged
                                  ? "converged"
                                  : "NOT converged")
                          << " outer iterations "
                          << s2s_outer_result.outer_iterations
                          << " surfaces "
                          << s2s_outer_result.radiation.surface_count
                          << " flux fixed point "
                          << s2s_outer_result
                                 .max_flux_fixed_point_normalized
                          << " wall temperature stationarity "
                          << s2s_outer_result
                                 .max_wall_temperature_change_normalized
                          << " balance_normalized "
                          << s2s_outer_result.radiation.balance_normalized
                          << " surface_identity_normalized "
                          << s2s_outer_result.radiation
                                 .max_surface_identity_normalized
                          << "\n";
                std::cout << "S2S last wall temperature lag "
                          << s2s_outer_result.max_wall_temperature_lag
                          << " K\n";
                std::cout << "S2S radiation summary "
                          << (options.output_dir / "radiation_s2s.json")
                                 .string()
                          << "\n";
            }
        }
        // Native benchmark fields are emitted only from the accepted final
        // solver state, never from the intermediate p-only iteration VTU.
        // A failed flow remains marked unqualified in the manifest; the
        // exports are still useful for diagnosing why it failed.
        if (options.benchmark_export) {
            const auto evidence = write_thermohydraulic_observables(
                options.output_dir / "benchmark_observables",
                mesh, U, p, result.authoritative_mass_flux, ubc,
                controls.density, controls.kinematic_viscosity,
                experimental_sst_state ? &*experimental_sst_state : nullptr,
                energy_ran ? &temperature : nullptr);
            std::cout << "Benchmark raw evidence: " << evidence.cells
                      << " cells, global face-flux mass imbalance="
                      << evidence.global_mass_imbalance_relative << "\n";
        }
        const bool overall_converged =
            result.converged &&
            (!buoyancy_ran || buoyancy_result.converged) &&
            (!energy_requested ||
             (energy_ran && energy_result.converged &&
              (!s2s_ran || s2s_outer_result.converged)));
        const auto solver_exit_code =
            execution_exit_code_for_solver(overall_converged);
        auto overall_status = result.convergence_status;
        std::string overall_reason = result.convergence_reason;
        if (result.converged && !overall_converged) {
            overall_status = cfdx::core::ConvergenceStatus::DIVERGED;
            if (energy_ran && !energy_result.temperature_admissibility_passed)
                overall_reason = "energy temperature admissibility gate failed";
            else if (energy_ran && !energy_result.converged)
                overall_reason = "energy stage linear or conservation gate failed";
            else if (s2s_ran && !s2s_outer_result.converged)
                overall_reason = "energy-radiation fixed point did not converge";
            else if (buoyancy_ran && !buoyancy_result.converged)
                overall_reason = "flow-energy buoyancy fixed point did not converge";
            else
                overall_reason = "requested energy stage did not complete";
        }

        const auto dat = options.output_dir / "restart.dat";
        DatRestartFields checkpoint_fields;
        if (energy_ran)
            checkpoint_fields.temperature = &temperature;
        if (experimental_sst_state) {
            // The experimental opt-in owns k and omega for the whole run, so
            // the final checkpoint persists them under their stable DAT v4
            // named fields. Without this wiring a persisted SST case can
            // initialize but never restart: the staged DAT reader requires
            // the named k and omega fields whenever the active runtime
            // requests them.
            checkpoint_fields.k = &experimental_sst_state->k;
            checkpoint_fields.omega = &experimental_sst_state->omega;
        }
        // Persist the authoritative face mass flux so the next run replays
        // the same iteration sequence instead of re-deriving a different
        // initial flux from the restored U and p.
        if (result.authoritative_mass_flux.dimension() == 1 &&
            result.authoritative_mass_flux.size() == mesh.n_faces())
            checkpoint_fields.face_mass_flux = &result.authoritative_mass_flux;
        write_dat_restart_fields(dat.string(), mesh, U, p, checkpoint_fields,
                                 result.iterations, 0.0);
        const auto convergence = options.output_dir / "convergence.json";
        write_convergence_history(convergence.string(), result);
        write_execution_summary(execution.string(), ExecutionSummary{
            solver_exit_code,
            overall_converged,
            result.iterations,
            static_cast<int>(overall_status),
            overall_reason,
            std::filesystem::is_regular_file(dat),
            std::filesystem::is_regular_file(convergence)});
        std::cout << "Checkpoint " << dat.string() << "\n";
        std::cout << "Convergence history " << convergence.string() << "\n";
        std::cout << "Execution summary " << execution.string() << "\n";

        if (!options.probe_csv.empty()) {
            write_probe_csv(options.probe_csv.string(), result.probe_samples, options.probes);
            std::cout << "Probe CSV " << options.probe_csv.string() << " ("
                      << result.probe_samples.size() << " samples)\n";
        }

        std::cout << "Converged " << (overall_converged ? "YES" : "NO")
                  << " Iterations " << result.iterations << "\n";
        return solver_exit_code;
    } catch (const std::exception& exc) {
        std::cerr << "CFDX production solver error: " << exc.what() << "\n";
        if (!execution_output_dir.empty()) {
            try {
                const auto execution = execution_output_dir / "execution.json";
                const bool has_restart =
                    std::filesystem::is_regular_file(execution_output_dir / "restart.dat");
                const bool has_convergence =
                    std::filesystem::is_regular_file(execution_output_dir / "convergence.json");
                write_execution_summary(execution.string(), ExecutionSummary{
                    kExecutionExitError,
                    false,
                    execution_result.has_value() ? execution_result->iterations : 0,
                    static_cast<int>(cfdx::core::ConvergenceStatus::DIVERGED),
                    exc.what(),
                    has_restart,
                    has_convergence});
            } catch (const std::exception& summary_exc) {
                std::cerr << "CFDX production solver: failed to persist execution summary: "
                          << summary_exc.what() << "\n";
            }
        }
        return kExecutionExitError;
    }
}
