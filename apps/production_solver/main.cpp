#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/io/gmsh/gmsh_importer.h"
#include "cfdx/io/restart/dat_restart.h"
#include "cfdx/io/vtu/vtu_writer.h"

#include <filesystem>
#include <iostream>
#include <string>

using namespace cfdx::core;
using namespace cfdx::io;
using namespace cfdx::physics;

namespace {
struct Options {
    std::string mesh;
    std::filesystem::path output_dir;
    std::filesystem::path restart;
    std::size_t iterations = 2;
    bool adaptive_convergence = false;
};

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
        else if (arg == "--help") {
            std::cout << "cfdx_production_solver --mesh PATH --output-dir DIR "
                         "[--restart DAT] [--iterations N] [--adaptive-convergence]\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + arg);
        }
    }
    if (o.mesh.empty() || o.output_dir.empty() || o.iterations == 0)
        throw std::invalid_argument("--mesh, --output-dir and positive --iterations are required");
    std::filesystem::create_directories(o.output_dir);
    return o;
}
}

int main(int argc, char** argv)
{
    try {
        const Options options = parse(argc, argv);
        Mesh mesh;
        if (!gmsh::import_gmsh_mesh(options.mesh, mesh) || mesh.n_cells() == 0)
            throw std::runtime_error("production solver: mesh import failed");

        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        U.fill(0.0);
        p.fill(0.0);

        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        for (std::size_t i = 0; i < mesh.boundary().n_patches(); ++i) {
            const auto& patch = mesh.boundary().patch(i);
            ubc[patch.name] = {
                VelocityBoundaryCondition::Type::FIXED_VALUE,
                {0.0, 0.0, 0.0}};
            pbc[patch.name] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        }

        IncompressibleSolverControls controls;
        controls.algorithm = PressureVelocityAlgorithm::SIMPLE;
        controls.convergence.max_iterations = options.iterations;
        controls.convergence.continuity_tolerance = 1e-12;
        controls.linear_max_iterations = 1000;
        controls.linear_tolerance = 1e-11;
        controls.acceleration.adaptive_linear_tolerance =
            options.adaptive_convergence;
        controls.acceleration.adaptive_pressure_correctors =
            options.adaptive_convergence;
        controls.pressure_reference_cell = 0;
        controls.pressure_reference_value = 0.0;

        if (!options.restart.empty()) {
            const auto restart_state = read_dat_restart(
                options.restart.string(), mesh, U, p);
            std::cout << "Restart iteration=" << restart_state.iteration
                      << " time=" << restart_state.time << "\n";
        }

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
                return true;
            };

        const auto result = solve_steady_incompressible(
            mesh, U, p, ubc, pbc, controls, options.restart.string());

        const auto dat = options.output_dir / "restart.dat";
        write_dat_restart(dat.string(), mesh, U, p, result.iterations, 0.0);
        std::cout << "Checkpoint " << dat.string() << "\n";
        std::cout << "Converged " << (result.converged ? "YES" : "NO")
                  << " Iterations " << result.iterations << "\n";
        return result.converged ? 0 : 1;
    } catch (const std::exception& exc) {
        std::cerr << "CFDX production solver error: " << exc.what() << "\n";
        return 2;
    }
}
