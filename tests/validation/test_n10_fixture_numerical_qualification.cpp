#include "cfdx/core/field/field.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/geometry/mesh_validator.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/io/mesh/mesh_importer.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using cfdx::core::Field;
using cfdx::core::Location;
using cfdx::core::Mesh;
using cfdx::core::Vec3;

namespace fs = std::filesystem;

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

Field<double, Location::CELL> make_field(
    const Mesh& mesh,
    const cfdx::core::GeometryCache& geometry,
    bool linear)
{
    Field<double, Location::CELL> field(mesh.n_cells(), "phi", "1", 1);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        const auto& p = geometry.cell_centres[c];
        field(c) = linear ? (2.0 * p.x - 3.0 * p.y + 0.5 * p.z) : 7.0;
    }
    return field;
}

struct Error {
    double linf{0.0};
    std::size_t count{0};
    std::size_t excluded{0};
};

Error gradient_error(
    const Mesh& mesh,
    const Field<double, Location::CELL>& gradient,
    bool linear,
    bool interior_only = false)
{
    const Vec3 exact{2.0, -3.0, 0.5};
    Error error;
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        if (mesh.cells().cell_size(c) == 0) continue;

        bool has_boundary_face = false;
        if (interior_only) {
            const auto offset = mesh.cells().cell_offset(c);
            const auto size = mesh.cells().cell_size(c);
            for (std::size_t k = 0; k < size; ++k) {
                const auto face = mesh.cells().faces_data()[offset + k];
                if (mesh.ownership().neighbour(face) < 0) {
                    has_boundary_face = true;
                    break;
                }
            }
            if (has_boundary_face) {
                ++error.excluded;
                continue;
            }
        }

        const Vec3 got{gradient(c, 0), gradient(c, 1), gradient(c, 2)};
        const Vec3 want = linear ? exact : Vec3{0.0, 0.0, 0.0};
        const double err = (got - want).mag();
        require(std::isfinite(err), "gradient error is non-finite");
        error.linf = std::max(error.linf, err);
        ++error.count;
    }
    return error;
}

fs::path normalize_vtk_fixture(const fs::path& fixture)
{
    if (fs::is_regular_file(fixture) && fixture.extension() == ".vtk") {
        return fixture;
    }

    if (fs::is_regular_file(fixture)) {
        const fs::path staging =
            fs::temp_directory_path() / "cfdx_n10_numerical_vtk_fixture.vtk";
        fs::copy_file(fixture, staging, fs::copy_options::overwrite_existing);
        return staging;
    }

    require(fs::is_directory(fixture),
            "N10.11 fixture path is neither a VTK file nor a fixture directory");

    const fs::path acquired = fixture / "06_unstructured.vtk";
    require(fs::is_regular_file(acquired),
            "N10.11 VTK fixture directory does not contain 06_unstructured.vtk");

    const fs::path staging =
        fs::temp_directory_path() / "cfdx_n10_numerical_vtk_fixture.vtk";
    fs::copy_file(acquired, staging, fs::copy_options::overwrite_existing);
    return staging;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "usage: test_n10_fixture_numerical_qualification <fixture> <report.json>\n";
        return 2;
    }

    const fs::path fixture = fs::absolute(argv[1]);
    const fs::path report_path = fs::absolute(argv[2]);
    fs::path staged_fixture;

    try {
        staged_fixture = normalize_vtk_fixture(fixture);

        Mesh mesh;
        require(cfdx::io::mesh::import_mesh(staged_fixture.string(), mesh),
                "production importer rejected numerical qualification fixture");
        require(mesh.n_cells() > 0 && mesh.n_faces() > 0,
                "qualification fixture is empty");

        const auto topology = mesh.topo_validate();
        require(topology.ok, "qualification fixture failed topology validation");

        const auto quality = cfdx::core::validate_mesh(mesh);
        require(quality.ok, "qualification fixture failed geometry validation");

        const auto geometry = cfdx::core::make_geometry_cache(mesh);

        const auto constant = make_field(mesh, geometry, false);
        const auto constant_gradient =
            cfdx::core::compute_gradient_least_squares(constant, mesh);
        const auto constant_error =
            gradient_error(mesh, constant_gradient, false);
        require(constant_error.count == mesh.n_cells(),
                "constant-field qualification did not cover all cells");
        require(constant_error.linf <= 1e-12,
                "constant field must have zero least-squares gradient");

        const auto linear = make_field(mesh, geometry, true);
        const auto linear_gradient =
            cfdx::core::compute_gradient_least_squares(linear, mesh);
        const auto linear_error =
            gradient_error(mesh, linear_gradient, true, true);
        require(linear_error.count > 0,
                "linear-field qualification found no interior cells with a full face-neighbour stencil");
        require(linear_error.linf <= 1e-9,
                "linear field must be exact for least-squares reconstruction");

        fs::create_directories(report_path.parent_path());
        std::ofstream report(report_path);
        require(report.good(), "unable to create qualification report");
        report << std::setprecision(17)
               << "{\n"
               << "  \"campaign\": \"N10.11 public fixture numerical qualification\",\n"
               << "  \"fixture\": \"meshio-vtk-unstructured\",\n"
               << "  \"status\": \"PASS\",\n"
               << "  \"qualification\": \"linear_reconstruction_on_imported_3d_mesh\",\n"
               << "  \"method\": \"least_squares\",\n"
               << "  \"cells_total\": " << mesh.n_cells() << ",\n"
               << "  \"cells_checked\": " << linear_error.count << ",\n"
               << "  \"boundary_cells_excluded\": " << linear_error.excluded << ",\n"
               << "  \"constant_gradient_linf\": " << constant_error.linf << ",\n"
               << "  \"linear_gradient_linf\": " << linear_error.linf << ",\n"
               << "  \"physical_validation\": \"NOT_CLAIMED\",\n"
               << "  \"observed_order\": \"NOT_CLAIMED\",\n"
               << "  \"numerical_tolerances_changed\": false\n"
               << "}\n";
        report.close();

        if (staged_fixture != fixture) fs::remove(staged_fixture);

        std::cout << "N10.11 fixture numerical qualification: PASS\n"
                  << "fixture=meshio-vtk-unstructured cells=" << linear_error.count
                  << " constant_gradient_Linf=" << constant_error.linf
                  << " linear_gradient_Linf=" << linear_error.linf << "\n";
        return 0;
    } catch (const std::exception& exc) {
        if (!staged_fixture.empty() && staged_fixture != fixture) fs::remove(staged_fixture);
        std::cerr << "N10.11 fixture numerical qualification FAILED: "
                  << exc.what() << "\n";
        return 1;
    }
}