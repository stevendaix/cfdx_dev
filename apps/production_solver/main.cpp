
#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/physics/coupled_matrix_capture_v16.h"
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
#include <chrono>
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
    std::string high_low_checkpoint_prefix;
    // Opt-in exact pressure CSR capture for production BFS and other cases.
    std::string pressure_replay_prefix;
    std::vector<std::size_t> pressure_replay_iterations;
    std::string high_low_checkpoint_prefix;
    bool validate_numerics_only = false;
    // Opt-in exact pressure CSR capture for production BFS and other cases.
    std::vector<IncompressiblePointProbe> probes;
    std::filesystem::path probe_csv;
    bool probe_stream = false;
    bool cell_hotspots = false;
    bool cell_map_full = false;
    std::size_t cell_hotspots_every = 5;
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
            else if (selection.method_id == "preconditioner.mgr")
                model = PreconditionerModel::MGR;
            else if (selection.method_id == "preconditioner.braess_sarazin")
                model = PreconditionerModel::BraessSarazin;
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
        else if (arg == "--high-low-checkpoint-prefix")
            o.high_low_checkpoint_prefix = value("--high-low-checkpoint-prefix");
        else if (arg == "--pressure-replay-prefix")
            o.pressure_replay_prefix = value("--pressure-replay-prefix");
        else if (arg == "--pressure-replay-iterations")
            o.pressure_replay_iterations = cfdx::core::parse_linear_capture_iterations(
                value("--pressure-replay-iterations"));
        else if (arg == "--high-low-checkpoint-prefix")
            o.high_low_checkpoint_prefix = value("--high-low-checkpoint-prefix");
        else if (arg == "--validate-numerics-only") o.validate_numerics_only = true;
        else if (arg == "--probe") o.probes.push_back(parse_probe(value("--probe")));
        else if (arg == "--probe-csv") o.probe_csv = value("--probe-csv");
        else if (arg == "--probe-stream") o.probe_stream = true;
        else if (arg == "--cell-hotspots") o.cell_hotspots = true;
        else if (arg == "--cell-map-full") { o.cell_map_full = true; o.cell_hotspots = true; }
        else if (arg == "--cell-hotspots-every") {
            const auto raw = value("--cell-hotspots-every");
            std::size_t used = 0;
            try { o.cell_hotspots_every = std::stoull(raw, &used); }
            catch (const std::exception&) { throw std::invalid_argument("invalid --cell-hotspots-every"); }
            if (used != raw.size() || o.cell_hotspots_every < 1 || o.cell_hotspots_every > 100000)
                throw std::invalid_argument("--cell-hotspots-every must be in [1,100000]");
        }
        else if (arg == "--help") {
            std::cout << "cfdx_production_solver --mesh PATH --output-dir DIR "
                         "[--restart DAT] [--iterations N] [--adaptive-convergence]\n"
                         "  [--initial-velocity-field] (initialize U from /fields/vector/U)\n"
                         "  [--initial-pressure-field] (initialize p from /fields/scalar/p; requires imported U)\n"
                         "  [--initial-pressure-field] (initialize p from /fields/scalar/p; requires imported U)\n"
                         "  [--initial-sst-fields] (initial k and omega from /fields/scalar)\n"
                         "  [--initial-face-mass-flux] (projected /fields/face/mdot; requires initial velocity)\n"
                         "  [--experimental-momentum-first-continuation] (momentum-only then strict SST; NOT qualified)\n"
                         "  [--wall-forces-quadratic] (opt-in 3-D quadratic WLS wall-force artifact; rank-deficient stencils fail)\n"
                         "  [--experimental-sst] (explicit experimental k-omega SST opt-in; not qualified)\n"
                         "  [--benchmark-export] (raw cell/face CSV, wall shear diagnostic and mass balance)\n"
                         "  [--high-low-checkpoint-prefix PATH] (DAT6 per committed high/low iteration)\n"
                         "  [--pressure-replay-prefix PATH --pressure-replay-iterations 1,20]\n"
                         "    (diagnostic post-gauge pressure CSR, RHS, initial guess; opt-in)\n"
                         "  [--high-low-checkpoint-prefix PATH] (DAT6 per committed high/low iteration)\n"
                         "  [--validate-numerics-only] (validate and print resolved case methods; do not solve)\n"
                         "  [--pressure-replay-prefix PATH --pressure-replay-iterations 1,20]\n"
                         "    (diagnostic post-gauge pressure CSR, RHS, initial guess; opt-in)\n"
                         "  [--probe <name>:<x>,<y>,<z>:<field>]... [--probe-csv PATH]\n"
                         "  probe field: p | u_x | u_y | u_z | u_mag\n"
                         "  [--cell-hotspots] [--cell-hotspots-every N] [--cell-map-full] (full map <=512 MiB)\n";
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
    const auto process_started = std::chrono::steady_clock::now();
    const auto elapsed_wall_seconds = [&]() -> double {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now() - process_started).count();
    };
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
        if (is_case_hdf5)
            resolve_production_mesh_patch_types(mesh, case_setup);
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
        if (options.validate_numerics_only) {
            if (!is_case_hdf5)
                throw std::invalid_argument(
                    "--validate-numerics-only requires a .cfdx.h5 case");
            IncompressibleSolverControls selection_controls;
            apply_explicit_case_numerics(case_setup, selection_controls);
            std::cout << "Resolved numerical selections:\\n"
                      << cfdx::core::format_numerics_report(case_setup.numerical_report);
            return 0;
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
            // A physical label must never be silently interpreted as an
            // unrelated generic zero-gradient condition.
            validate_production_incompressible_boundary(*bc, patch.type);

            // Explicit per-field conditions override the legacy shared value_type.
            // Never infer the pressure boundary from the velocity mode.
            if (const auto* explicit_u = bc->find_field("velocity")) {
                if (explicit_u->type == "FIXED_VALUE") {
                    ubc[patch.name] = {VelocityBoundaryCondition::Type::FIXED_VALUE,
                        {explicit_u->value.at(0), explicit_u->value.at(1), explicit_u->value.at(2)}};
                } else if (explicit_u->type == "ZERO_GRADIENT") {
                    ubc[patch.name] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,
                        {0.0, 0.0, 0.0}};
                } else throw std::invalid_argument("unsupported native velocity field condition");
            } else if (bc->value_type == BCValueType::FIXED ||
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

            if (const auto* explicit_p = bc->find_field("pressure")) {
                if (explicit_p->type == "FIXED_VALUE")
                    pbc[patch.name] = {ScalarBoundaryType::FIXED_VALUE, explicit_p->value.at(0), 0.0};
                else if (explicit_p->type == "FIXED_GRADIENT")
                    pbc[patch.name] = {ScalarBoundaryType::FIXED_GRADIENT, 0.0, explicit_p->gradient};
                else if (explicit_p->type == "ZERO_GRADIENT")
                    pbc[patch.name] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
                else throw std::invalid_argument("unsupported native pressure field condition");
            } else if (bc->value_type == BCValueType::OUTLET_PRESSURE) {
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
