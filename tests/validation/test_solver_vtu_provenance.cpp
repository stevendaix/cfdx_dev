#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/io/vtu/vtu_writer.h"
#include "common/test_harness.h"

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::io;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_unit_cube()
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i = 0; i < 8; ++i)
        m.points().set(i, p[i][0], p[i][1], p[i][2]);

    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});

    m.ownership().resize(6);
    for (std::size_t f = 0; f < 6; ++f) {
        m.ownership().set_owner(f, 0);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }

    m.cells().push_cell({0,1,2,3,4,5});

    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,1,2,3,4,5};
    m.boundary().add_patch(wall);
    return m;
}

int main()
{
    run_case("solver_generated_vtu_carries_authoritative_iteration_and_time", [] {
        const Mesh mesh = make_unit_cube();

        Field<double, Location::CELL> U(1, "U", "m/s", 3);
        Field<double, Location::CELL> p(1, "p", "Pa", 1);
        U.fill(0.0);
        p.fill(0.0);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {
            VelocityBoundaryCondition::Type::FIXED_VALUE,
            {0.0, 0.0, 0.0}
        };

        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        IncompressibleSolverControls controls;
        controls.algorithm = PressureVelocityAlgorithm::SIMPLE;
        controls.convergence.max_iterations = 2;
        controls.convergence.continuity_tolerance = 1e-12;
        controls.linear_tolerance = 1e-12;
        controls.pressure_reference_cell = 0;
        controls.pressure_reference_value = 0.0;

        std::vector<std::filesystem::path> outputs;
        controls.iteration_output_callback =
            [&](std::size_t iteration, double time, const Mesh& state_mesh,
                const Field<double, Location::CELL>& state_U,
                const Field<double, Location::CELL>& state_p) {
                const auto path = std::filesystem::temp_directory_path() /
                    ("cfdx_solver_vtu_provenance_" + std::to_string(iteration) + ".vtu");

                ScalarCellField pressure(state_p);
                VtuWriter writer;
                EXPECT_TRUE(writer.write(
                    path.string(), state_mesh, {{"p", pressure}}, {}, {}, time, iteration, true));
                outputs.push_back(path);
                return true;
            };

        const auto result =
            solve_steady_incompressible(mesh, U, p, ubc, pbc, controls);

        EXPECT_TRUE(result.iterations >= std::size_t{1});
        EXPECT_TRUE(!outputs.empty());

        const auto final_path = outputs.back();
        std::ifstream in(final_path);
        const std::string xml(
            (std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

        EXPECT_TRUE(xml.find("Name=\"physical_time\"") != std::string::npos);
        EXPECT_TRUE(xml.find("Name=\"iteration\"") != std::string::npos);
        EXPECT_TRUE(xml.find("0.000000000000e+00") != std::string::npos);
        EXPECT_TRUE(
            xml.find(std::to_string(result.iterations)) != std::string::npos);

        for (const auto& path : outputs)
            std::filesystem::remove(path);
    });

    return run_all();
}
