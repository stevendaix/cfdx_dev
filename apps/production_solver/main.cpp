
#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/physics/production_energy_stage.h"
#include "cfdx/physics/production_turbulence_contract.h"
#include "cfdx/physics/production_radiation_contract.h"
#include "cfdx/physics/production_radiation_s2s_stage.h"
#include "cfdx/io/gmsh/gmsh_importer.h"
#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/io/hdf5/case_hdf5_io.h"
#include "cfdx/io/restart/dat_restart.h"
#include "cfdx/io/probe/probe_csv.h"
#include "cfdx/io/runtime/convergence_history.h"
#include "cfdx/io/runtime/execution_summary.h"
#include "cfdx/io/vtu/vtu_writer.h"

#include <algorithm>
#include <cmath>
#include <csignal>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace cfdx::core;
using namespace cfdx::io;
using namespace cfdx::physics;

namespace {
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
    std::vector<IncompressiblePointProbe> probes;
    std::filesystem::path probe_csv;
};


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
    out << "  \"global_energy_balance_passed\": "
        << (energy.global_energy_balance_passed ? "true" : "false") << ",\n";
    out << "  \"local_energy_balance_passed\": "
        << (energy.local_energy_balance_passed ? "true" : "false") << ",\n";
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
            << json_number(radiation.injected_heat_flux.at(name)) << "}";
        out << (index + 1 == radiation.surface_names.size() ? "\n"
                                                            : ",\n");
    }
    out << "    ]\n";
    out << "  }\n";
    out << "}\n";
    return out.str();
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
        else if (arg == "--probe") o.probes.push_back(parse_probe(value("--probe")));
        else if (arg == "--probe-csv") o.probe_csv = value("--probe-csv");
        else if (arg == "--help") {
            std::cout << "cfdx_production_solver --mesh PATH --output-dir DIR "
                         "[--restart DAT] [--iterations N] [--adaptive-convergence]\n"
                         "  [--initial-velocity-field] (initialize U from /fields/vector/U)\n"
                         "  [--probe <name>:<x>,<y>,<z>:<field>]... [--probe-csv PATH]\n"
                         "  probe field: p | u_x | u_y | u_z | u_mag\n";
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
            ok = read_mesh_hdf5(mesh_path, mesh);
        } else {
            ok = gmsh::import_gmsh_mesh(mesh_path, mesh);
        }
        if (!ok || mesh.n_cells() == 0)
            throw std::runtime_error("production solver: mesh import failed");

        // The declared energy model is validated fail-closed as soon as the
        // case is loaded, before any expensive flow computation: an unknown
        // declaration (e.g. "buoyant") must fail the run up front instead of
        // being discovered after the carrier solve, even when the energy
        // stage itself would not run.
        if (is_case_hdf5) {
            validate_production_turbulence_selection(case_setup);
            validate_production_energy_model(case_setup);
            validate_production_radiation_selection(case_setup);
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
        controls.pressure_reference_cell = 0;
        controls.pressure_reference_value = 0.0;
        controls.probes = options.probes;

        // The declared initialization strategy is applied by the solver, before
        // the first residual is formed. An explicit --restart request therefore
        // overrides the case-derived uniform state.
        controls.initialization = initialization;
        if (!options.restart.empty())
            controls.initialization.mode = InitializationMode::Restart;
        std::cout << "Initialization strategy "
                  << cfdx::physics::to_string(controls.initialization.mode) << "\n";

        controls.iteration_output_callback =
            [&](std::size_t iteration, double time, const Mesh& state_mesh,
                const Field<double, Location::CELL>&,
                const Field<double, Location::CELL>& state_p) {
                ScalarCellField pressure(state_p);
                const auto path = options.output_dir /
                    ("result_" + std::to_string(iteration) + ".vtu");
                VtuWriter writer;
                if (!writer.write(path.string(), state_mesh, {{"p", pressure}},
                                  {}, {}, time, iteration, true))
                    return false;
                std::cout << "Iteration " << iteration
                          << " Time = " << time << "\n";
                return g_stop_requested == 0;
            };

        const auto result = solve_steady_incompressible(
            mesh, U, p, ubc, pbc, controls, options.restart.string());
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
        bool energy_requested = false;
        bool energy_ran = false;
        ProductionEnergyStageResult energy_result;
        ProductionS2SOuterLoopResult s2s_outer_result;
        bool s2s_ran = false;
        if (is_case_hdf5)
            energy_requested = production_energy_requested(case_setup);
        if (energy_requested && result.converged) {
            const auto energy_geometry = build_fv_geometry(mesh);
            // Initial temperature: the declared initial condition, else
            // the case reference temperature. Zero means "not declared"
            // (stage-1 convention) and falls back to ref_temperature.
            double initial_temperature = case_setup.ref_temperature;
            if (case_setup.has_initial_condition &&
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
            Field<double, Location::CELL> temperature(
                mesh.n_cells(), "temperature", "K", 1);
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
                energy_result = run_production_energy_stage(
                    mesh, energy_geometry, case_setup,
                    result.authoritative_mass_flux, temperature);
                energy_ran = true;
            }

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

        } else if (energy_requested) {
            std::cout
                << "Energy stage skipped: carrier flow did not converge\n";
        }
        const bool overall_converged =
            result.converged &&
            (!energy_requested ||
             (energy_ran && energy_result.converged &&
              (!s2s_ran || s2s_outer_result.converged)));
        const auto solver_exit_code =
            execution_exit_code_for_solver(overall_converged);

        const auto dat = options.output_dir / "restart.dat";
        write_dat_restart(dat.string(), mesh, U, p, result.iterations, 0.0);
        const auto convergence = options.output_dir / "convergence.json";
        write_convergence_history(convergence.string(), result);
        write_execution_summary(execution.string(), ExecutionSummary{
            solver_exit_code,
            overall_converged,
            result.iterations,
            static_cast<int>(result.convergence_status),
            result.convergence_reason,
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
                    execution_result.has_value() ? execution_result->converged : false,
                    execution_result.has_value() ? execution_result->iterations : 0,
                    execution_result.has_value()
                        ? static_cast<int>(execution_result->convergence_status)
                        : 0,
                    execution_result.has_value()
                        ? execution_result->convergence_reason
                        : exc.what(),
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

