#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/solvers/scalar_diffusion.h"
#include "cfdx/io/gmsh/gmsh_importer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr double PI = 3.14159265358979323846;

double exact_solution(double x, double y)
{
    return std::sin(PI * x) * std::sin(PI * y);
}

double exact_source(double x, double y)
{
    // -laplacian(phi) for phi = sin(pi*x) sin(pi*y).
    return 2.0 * PI * PI * exact_solution(x, y);
}

} // namespace

int main(int argc, char** argv)
{
    std::string mesh_path = "/tmp/mms_meshes/mms_square_32x32.msh";
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--mesh" && i + 1 < argc)
            mesh_path = argv[++i];
    }

    try {
        cfdx::core::Mesh mesh;
        if (!cfdx::io::gmsh::import_gmsh_mesh(mesh_path, mesh)) {
            std::cerr << "Failed to load mesh: " << mesh_path << "\n";
            return 1;
        }

        const cfdx::core::GeometryCache geometry =
            cfdx::core::make_geometry_cache(mesh);

        std::vector<double> source(mesh.n_cells(), 0.0);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            const auto& x = geometry.cell_centres[c];
            source[c] = exact_source(x.x, x.y);
        }

        auto boundary =
            cfdx::core::PoissonBoundaryCondition::dirichlet(mesh.n_faces());
        for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
            if (mesh.ownership().neighbour(f) < 0) {
                const auto& x = geometry.face_centres[f];
                boundary.face_values[f] = exact_solution(x.x, x.y);
            }
        }

        const cfdx::core::ScalarDiffusionConfig config{
            .diffusivity = 1.0,
            .max_iterations = 5000,
            .tolerance = 1e-11};

        const auto result = cfdx::core::solve_poisson_dirichlet(
            mesh, boundary, source, config);

        double max_error = 0.0;
        double l2_error_squared = 0.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            const auto& x = geometry.cell_centres[c];
            const double error =
                result.solution(c) - exact_solution(x.x, x.y);
            max_error = std::max(max_error, std::abs(error));
            l2_error_squared +=
                error * error * geometry.cell_volumes[c];
        }

        const double l2_error = std::sqrt(l2_error_squared);
        std::cout << "Poisson MMS (canonical Dirichlet path)\n";
        std::cout << "  cells      = " << mesh.n_cells() << "\n";
        std::cout << "  L_inf      = " << max_error << "\n";
        std::cout << "  L_2        = " << l2_error << "\n";
        std::cout << "  iterations = " << result.linear_result.iterations << "\n";
        std::cout << "  residual    = " << result.linear_result.residual << "\n";

        return result.linear_result.converged() ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "Poisson MMS exception: " << e.what() << "\n";
        return 1;
    }
}
