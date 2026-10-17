#include "cfdx/core/numerics/gradient.h"
#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

double linear_field(double y) { return y; }
double quadratic_field(double y) { return y * (1.0 - y); }
double quadratic_derivative(double y) { return 1.0 - 2.0 * y; }

struct Errors {
    double l2 = 0.0;
    double linf = 0.0;
    std::size_t count = 0;
};

Errors gradient_error(const Field<double, Location::CELL>& gradient,
                      const FvGeometry& geometry,
                      std::size_t nx,
                      std::size_t ny,
                      bool quadratic)
{
    Errors e;
    double sum2 = 0.0;
    // Exclude inlet/outlet-adjacent cells: gauss_gradient_with_boundary accepts
    // patch-wise constants, whereas the exact inlet/outlet data vary with y.
    for (std::size_t c = 0; c < geometry.cell_centres.size(); ++c) {
        const auto& x = geometry.cell_centres[c];
        const auto i = static_cast<std::size_t>(
            std::llround(x.x * static_cast<double>(nx) / 8.0));
        if (i == 0 || i + 1 >= nx) continue;
        const double exact_y = quadratic ? quadratic_derivative(x.y) : 1.0;
        const double ex = gradient.component_data(0)[c];
        const double ey = gradient.component_data(1)[c];
        const double ez = gradient.component_data(2)[c];
        const double err = std::sqrt(ex * ex + (ey - exact_y) * (ey - exact_y) + ez * ez);
        sum2 += err * err;
        e.linf = std::max(e.linf, err);
        ++e.count;
    }
    if (e.count == 0) throw std::runtime_error("no interior channel cells in gradient audit");
    e.l2 = std::sqrt(sum2 / static_cast<double>(e.count));
    (void)ny;
    return e;
}

void require_finite(const Field<double, Location::CELL>& gradient, const char* name)
{
    for (std::size_t c = 0; c < gradient.size(); ++c) {
        for (std::size_t d = 0; d < 3; ++d) {
            if (!std::isfinite(gradient.component_data(d)[c]))
                throw std::runtime_error(std::string(name) + ": non-finite gradient");
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    try {
        if (argc != 2)
            throw std::invalid_argument("usage: test_wall_shear_gradient_audit channel_n16.h5");

        constexpr std::size_t nx = 16;
        constexpr std::size_t ny = 16;
        Mesh mesh;
        if (!cfdx::io::read_mesh_hdf5(argv[1], mesh))
            throw std::runtime_error(std::string("failed to read channel mesh: ") + argv[1]);
        const auto topo = mesh.topo_validate();
        if (!topo.ok) throw std::runtime_error("channel mesh topology is invalid");

        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::CELL> linear(mesh.n_cells(), "phi_linear", "m", 1);
        Field<double, Location::CELL> quadratic(mesh.n_cells(), "u_quadratic", "m/s", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            const double y = geometry.cell_centres[c].y;
            linear(c) = linear_field(y);
            quadratic(c) = quadratic_field(y);
        }

        // Face-indexed WLS boundary conditions: exact Dirichlet values at each
        // boundary face; front/back are EMPTY and therefore use zero-gradient ghosts.
        std::vector<BoundaryGradientCondition> linear_bc(mesh.n_faces());
        std::vector<BoundaryGradientCondition> quadratic_bc(mesh.n_faces());
        for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
            if (mesh.ownership().neighbour(f) >= 0) continue;
            const auto patch_id = geometry.face_patch[f];
            if (patch_id < mesh.boundary().n_patches() &&
                mesh.boundary().patch(patch_id).type == PatchType::EMPTY) {
                linear_bc[f] = {BoundaryGradientConditionType::ZERO_GRADIENT, 0.0};
                quadratic_bc[f] = {BoundaryGradientConditionType::ZERO_GRADIENT, 0.0};
            } else {
                const double y = geometry.face_centres[f].y;
                linear_bc[f] = {BoundaryGradientConditionType::DIRICHLET, linear_field(y)};
                quadratic_bc[f] = {BoundaryGradientConditionType::DIRICHLET, quadratic_field(y)};
            }
        }

        const auto linear_wls = compute_gradient_weighted_least_squares(
            linear, mesh, GradientWeighting::INVERSE_DISTANCE_SQUARED,
            std::numeric_limits<double>::infinity(),
            BoundaryGradientPolicy::ZERO_GRADIENT_GHOST, &linear_bc);
        const auto quadratic_wls = compute_gradient_weighted_least_squares(
            quadratic, mesh, GradientWeighting::INVERSE_DISTANCE_SQUARED,
            std::numeric_limits<double>::infinity(),
            BoundaryGradientPolicy::ZERO_GRADIENT_GHOST, &quadratic_bc);

        ScalarBoundaryConditions quadratic_gg_bc;
        quadratic_gg_bc["bottom"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
        quadratic_gg_bc["top"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
        quadratic_gg_bc["front"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        quadratic_gg_bc["back"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        quadratic_gg_bc["inlet"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        quadratic_gg_bc["outlet"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        const auto quadratic_gauss =
            gauss_gradient_with_boundary(quadratic, mesh, geometry, quadratic_gg_bc);

        require_finite(linear_wls, "linear WLS");
        require_finite(quadratic_wls, "quadratic WLS");
        require_finite(quadratic_gauss, "quadratic Green-Gauss");

        const auto linear_error = gradient_error(linear_wls, geometry, nx, ny, false);
        const auto quadratic_wls_error = gradient_error(quadratic_wls, geometry, nx, ny, true);
        const auto quadratic_gg_error = gradient_error(quadratic_gauss, geometry, nx, ny, true);

        std::cout << "WALL_SHEAR_GRADIENT_AUDIT"
                  << " linear_wls_l2=" << linear_error.l2
                  << " linear_wls_linf=" << linear_error.linf
                  << " quadratic_wls_l2=" << quadratic_wls_error.l2
                  << " quadratic_wls_linf=" << quadratic_wls_error.linf
                  << " quadratic_gauss_l2=" << quadratic_gg_error.l2
                  << " quadratic_gauss_linf=" << quadratic_gg_error.linf
                  << " cells=" << linear_error.count << '\n';

        // The affine-field WLS result is a consistency gate. Quadratic results
        // are diagnostic only until their expected order is qualified on the
        // production NACA mesh family; no method is selected from this test alone.
        if (linear_error.linf > 1e-9)
            throw std::runtime_error("WLS failed affine-field exactness on channel mesh");

        std::cout << "WALL_SHEAR_GRADIENT_AUDIT: PASS (quadratic metrics diagnostic)\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "WALL_SHEAR_GRADIENT_AUDIT: FAIL: " << e.what() << '\n';
        return 1;
    }
}
