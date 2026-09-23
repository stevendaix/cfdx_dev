#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/io/restart/dat_restart.h"
#include "cfdx/io/vtu/vtu_writer.h"
#include "cfdx/physics/steady_incompressible_solver.h"

#include <chrono>
#include <filesystem>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <thread>

using cfdx::core::Field;
using cfdx::core::Location;
using cfdx::core::Mesh;
using cfdx::core::ScalarCellField;
using cfdx::physics::IncompressibleSolverControls;
using cfdx::physics::ScalarBoundaryConditions;
using cfdx::physics::ScalarBoundaryType;
using cfdx::physics::VelocityBoundaryConditions;
using cfdx::physics::VelocityBoundaryCondition;

namespace {
std::string step_name(std::size_t iteration, const char* extension)
{
    const std::string number = std::to_string(iteration);
    const std::string padding(number.size() < 4 ? 4 - number.size() : 0, '0');
    return "step_" + padding + number + extension;
}

void write_outputs(
    const std::filesystem::path& output_dir,
    const Mesh& mesh,
    const Field<double, Location::CELL>& U,
    const Field<double, Location::CELL>& p,
    std::size_t iteration,
    double time)
{
    std::filesystem::create_directories(output_dir);
    cfdx::io::write_dat_restart(
        (output_dir / step_name(iteration, ".dat")).string(),
        mesh, U, p, iteration, time);

    ScalarCellField ux(mesh.n_cells(), "Ux", "m/s", 1);
    ScalarCellField uy(mesh.n_cells(), "Uy", "m/s", 1);
    ScalarCellField uz(mesh.n_cells(), "Uz", "m/s", 1);
    ScalarCellField pressure = p;
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        ux(c) = U.component_data(0)[c];
        uy(c) = U.component_data(1)[c];
        uz(c) = U.component_data(2)[c];
    }

    std::map<std::string, ScalarCellField> fields{
        {"Ux", ux}, {"Uy", uy}, {"Uz", uz}, {"p", pressure}
    };
    if (!cfdx::io::VtuWriter{}.write(
            (output_dir / step_name(iteration, ".vtu")).string(),
            mesh, fields, {}, {}, time, iteration, true)) {
        throw std::runtime_error("failed to write VTU output");
    }
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "usage: cfdx_application_e2e_driver <case.h5> <output_dir> "
                     "[--restart <dat>] [--iterations <n>] [--delay-ms <n>]\n";
        return 2;
    }

    try {
        const std::filesystem::path case_path = argv[1];
        const std::filesystem::path output_dir = argv[2];
        std::filesystem::path restart_path;
        int delay_ms = 0;
        std::size_t iterations = 2;

        for (int i = 3; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--restart") {
                if (++i >= argc) throw std::invalid_argument("--restart requires a path");
                restart_path = argv[i];
            } else if (arg == "--iterations") {
                if (++i >= argc) throw std::invalid_argument("--iterations requires an integer");
                iterations = static_cast<std::size_t>(std::stoull(argv[i]));
                if (iterations == 0) throw std::invalid_argument("--iterations must be positive");
            } else if (arg == "--delay-ms") {
                if (++i >= argc) throw std::invalid_argument("--delay-ms requires an integer");
                delay_ms = std::stoi(argv[i]);
                if (delay_ms < 0) throw std::invalid_argument("--delay-ms must be non-negative");
            } else {
                throw std::invalid_argument("unknown argument: " + arg);
            }
        }

        if (!std::filesystem::is_regular_file(case_path)) {
            throw std::invalid_argument("case file does not exist");
        }

        Mesh mesh;
        if (!cfdx::io::read_mesh_hdf5(case_path.string(), mesh)) {
            throw std::runtime_error("native HDF5 mesh import failed");
        }
        if (mesh.n_cells() != 1 || mesh.n_faces() != 6) {
            throw std::runtime_error("E2E driver expects the one-cell unit-cube acceptance mesh");
        }

        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        if (restart_path.empty()) {
            U.set(0, 0.25, -0.15, 0.05);
            p(0) = 37.5;
        } else {
            U.fill(0.0);
            p.fill(0.0);
        }

        VelocityBoundaryConditions velocity_bcs;
        velocity_bcs["wall"] = {
            VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}
        };
        ScalarBoundaryConditions pressure_bcs;
        pressure_bcs["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        IncompressibleSolverControls controls;
        controls.algorithm = cfdx::physics::PressureVelocityAlgorithm::SIMPLE;
        controls.convergence.max_iterations = iterations;
        controls.convergence.relative_tolerance = 1e-30;
        controls.convergence.continuity_tolerance = 1e-30;
        controls.convergence.velocity_tolerance = 1e-30;
        controls.convergence.pressure_tolerance = 1e-30;
        controls.linear_tolerance = 1e-12;
        controls.pressure_reference_cell = 0;
        controls.pressure_reference_value = 0.0;
        controls.iteration_output_callback =
            [&](std::size_t iteration, double time, const Mesh& callback_mesh,
                const Field<double, Location::CELL>& callback_U,
                const Field<double, Location::CELL>& callback_p) {
                write_outputs(output_dir, callback_mesh, callback_U, callback_p, iteration, time);
                std::cout << "Iteration " << iteration << " Time = " << time << std::endl;
                if (delay_ms > 0) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(delay_ms));
                }
                return true;
            };

        const auto result = cfdx::physics::solve_steady_incompressible(
            mesh, U, p, velocity_bcs, pressure_bcs, controls, restart_path.string());
        if (result.history.empty()) {
            throw std::runtime_error("solver produced no nonlinear iteration");
        }

        std::cout << "CFDX_E2E iterations=" << result.iterations
                  << " converged=" << (result.converged ? 1 : 0)
                  << " cells=" << mesh.n_cells() << std::endl;
        return 0;
    } catch (const std::exception& exc) {
        std::cerr << "CFDX_E2E_ERROR " << exc.what() << std::endl;
        return 1;
    }
}
