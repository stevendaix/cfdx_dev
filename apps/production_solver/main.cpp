#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/io/gmsh/gmsh_importer.h"
#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/io/hdf5/case_hdf5_io.h"
#include "cfdx/io/restart/dat_restart.h"
#include "cfdx/io/probe/probe_csv.h"
#include "cfdx/io/runtime/convergence_history.h"
#include "cfdx/io/runtime/execution_summary.h"
#include "cfdx/io/vtu/vtu_writer.h"

#include <algorithm>
#include <csignal>
#include <filesystem>
#include <iostream>
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
        else if (arg == "--probe") o.probes.push_back(parse_probe(value("--probe")));
        else if (arg == "--probe-csv") o.probe_csv = value("--probe-csv");
        else if (arg == "--help") {
            std::cout << "cfdx_production_solver --mesh PATH --output-dir DIR "
                         "[--restart DAT] [--iterations N] [--adaptive-convergence]\n"
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
    try {
        const Options options = parse(argc, argv);
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

        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);

        // The initial state is declared, not inferred. When the case carries an
        // explicit initial_condition block it is mapped onto the fields through
        // the uniform strategy; otherwise the caller-provided state is kept.
        // A restart, when requested, replaces it below.
        InitializationControls initialization;
        if (is_case_hdf5 && case_setup.has_initial_condition) {
            const auto& ic = case_setup.initial_condition;
            initialization.mode = InitializationMode::Uniform;
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
        if (!options.restart.empty())
            std::cout << "Restart checkpoint " << options.restart.string() << "\n";

        const auto dat = options.output_dir / "restart.dat";
        write_dat_restart(dat.string(), mesh, U, p, result.iterations, 0.0);
        const auto convergence = options.output_dir / "convergence.json";
        write_convergence_history(convergence.string(), result);
        const auto execution = options.output_dir / "execution.json";
        write_execution_summary(execution.string(), ExecutionSummary{
            result.converged ? 0 : 1,
            result.converged,
            result.iterations,
            static_cast<int>(result.convergence_status),
            result.convergence_reason,
            std::filesystem::exists(dat),
            std::filesystem::exists(convergence)});
        std::cout << "Checkpoint " << dat.string() << "\n";
        std::cout << "Convergence history " << convergence.string() << "\n";
        std::cout << "Execution summary " << execution.string() << "\n";

        if (!options.probe_csv.empty()) {
            write_probe_csv(options.probe_csv.string(), result.probe_samples, options.probes);
            std::cout << "Probe CSV " << options.probe_csv.string() << " ("
                      << result.probe_samples.size() << " samples)\n";
        }

        std::cout << "Converged " << (result.converged ? "YES" : "NO")
                  << " Iterations " << result.iterations << "\n";
        return result.converged ? 0 : 1;
    } catch (const std::exception& exc) {
        std::cerr << "CFDX production solver error: " << exc.what() << "\n";
        return 2;
    }
}
