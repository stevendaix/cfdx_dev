#include "cfdx/core/geometry/mesh_validator.h"
#include "cfdx/io/mesh/mesh_importer.h"
#include "cfdx/io/hdf5/case_hdf5_io.h"
#include "cfdx/io/cfdx_io/case_schema.h"

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace fs = std::filesystem;

using cfdx::core::Mesh;
using cfdx::core::NumericalMethodFamily;
using cfdx::io::BCType;
using cfdx::io::BCValueType;
using cfdx::io::BoundarySpec;
using cfdx::io::CaseSetup;
using cfdx::io::GapAnalysis;
using cfdx::io::SourceInfo;

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "usage: test_n10_fixture_production_case <source> <output.cfdx.h5>\n";
        return 2;
    }

    const fs::path source_path = fs::absolute(argv[1]);
    const fs::path output_path = fs::absolute(argv[2]);

    try {
        Mesh mesh;
        if (!cfdx::io::mesh::import_mesh(source_path.string(), mesh))
            throw std::runtime_error("production C++ importer rejected fixture");
        if (mesh.n_cells() == 0 || mesh.n_faces() == 0)
            throw std::runtime_error("imported fixture has no cells/faces");

        const auto topology = mesh.topo_validate();
        if (!topology.ok)
            throw std::runtime_error("imported fixture failed topology validation");

        const auto quality = cfdx::core::validate_mesh(mesh);
        if (!quality.ok)
            throw std::runtime_error("imported fixture failed geometry validation");

        CaseSetup setup;
        setup.source.solver = "meshio";
        setup.source.version = "pinned-public-fixture";
        setup.source.case_path = source_path.string();
        setup.source.case_name = source_path.stem().string();
        setup.source.format = "vtk";
        setup.mesh_info.n_vertices = mesh.n_points();
        setup.mesh_info.n_faces = mesh.n_faces();
        setup.mesh_info.n_cells = mesh.n_cells();
        setup.mesh_info.n_boundary_faces = mesh.stats().n_boundary_faces;
        setup.mesh_info.n_internal_faces = mesh.stats().n_internal_faces;
        setup.mesh_info.n_patches = mesh.stats().n_patches;
        setup.mesh_info.dimension = 3;
        setup.mesh_info.mesh_type = "unstructured";
        setup.physics_model = "incompressible_laminar";
        setup.turbulence_model = "laminar";
        setup.energy_model = "isothermal";
        setup.transient = false;
        setup.solver_mode = "steady";
        setup.units = "SI";
        setup.numerics.max_iterations = 500;

        // This case is deliberately an execution fixture, not a physical
        // validation case. The public meshio VTK fixture does not define CFD
        // semantics for its patch ordering, so the execution case must not
        // invent a "moving lid" from patch index. Homogeneous stationary
        // no-slip walls with a zero initial state provide a deterministic,
        // well-posed smoke case while still exercising the complete
        // production import -> case serialization -> solver path.
        for (std::size_t i = 0; i < mesh.boundary().n_patches(); ++i) {
            const auto& patch = mesh.boundary().patch(i);
            BoundarySpec bc;
            bc.patch_name = patch.name;
            bc.type = BCType::WALL;
            bc.value_type = BCValueType::WALL_NO_SLIP;
            bc.velocity_vector = {0.0, 0.0, 0.0};
            bc.source_zone_name = patch.name;
            bc.source_type_name = "CFDX N10 production execution fixture";
            setup.boundary_conditions.push_back(std::move(bc));
        }
        if (setup.boundary_conditions.empty())
            throw std::runtime_error("imported fixture has no boundary patches");

        std::cout << "N10 boundary patches:";
        for (std::size_t i = 0; i < mesh.boundary().n_patches(); ++i) {
            const auto& imported_patch = mesh.boundary().patch(i);
            std::cout << " [" << i << "]" << imported_patch.name
                      << " faces=" << imported_patch.face_ids.size();
        }
        std::cout << "\n";

        setup.has_initial_condition = true;
        setup.initial_condition.velocity_vector = {0.0, 0.0, 0.0};
        setup.initial_condition.pressure = 0.0;

        // Explicit N1 selections: these are resolved by the canonical HDF5
        // writer and consumed by the production solver. No solver fallback is
        // introduced by this fixture path.
        using F = NumericalMethodFamily;
        setup.numerical_config.entries = {
            {F::Gradient, "numerics.gradient.gauss"},
            {F::Convection, "numerics.convection.upwind"},
            {F::PressureVelocity, "pressure_velocity.simple"},
        };
        setup.numerical_config.required_families = {
            F::Gradient,
            F::Convection,
            F::PressureVelocity,
        };
        setup.has_explicit_numerics = true;

        SourceInfo source;
        source.solver = setup.source.solver;
        source.version = setup.source.version;
        source.case_path = setup.source.case_path;
        source.case_name = setup.source.case_name;
        source.format = setup.source.format;

        GapAnalysis gap;
        gap.supported("mesh", "vtk", "Imported through the production C++ mesh importer");
        gap.supported("boundary_condition", "homogeneous_wall",
                      "Deterministic N10 numerical execution fixture");
        gap.supported("numerics", "explicit_selection",
                      "Gradient, convection and pressure-velocity are explicitly selected");

        fs::create_directories(output_path.parent_path());
        if (!cfdx::io::write_case_cfdx_h5(
                output_path.string(), mesh, source, setup, gap))
            throw std::runtime_error("failed to write CFDX HDF5 production case");

        std::cout << "N10 production case created: " << output_path << "\n";
        std::cout << "  cells=" << mesh.n_cells()
                  << " faces=" << mesh.n_faces()
                  << " patches=" << mesh.boundary().n_patches() << "\n";
        return 0;
    } catch (const std::exception& exc) {
        std::cerr << "N10 production case creation FAILED: " << exc.what() << "\n";
        return 1;
    }
}
