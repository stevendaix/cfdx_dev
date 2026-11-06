#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/physics/wall_quadratic_forces.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

constexpr double kMu = 0.7;

double exact_pressure(const Vec3& x)
{
    return 2.0 + 0.3*x.x - 0.2*x.y + 0.4*x.z
        + 0.05*x.x*x.x + 0.02*x.x*x.y - 0.03*x.x*x.z
        + 0.04*x.y*x.y + 0.01*x.y*x.z + 0.06*x.z*x.z;
}

Vec3 exact_velocity(const Vec3& x)
{
    return {
        1.0 + 0.2*x.x + 0.3*x.y - 0.1*x.z
            + 0.04*x.x*x.x + 0.02*x.y*x.z,
        -0.5 + 0.1*x.x - 0.2*x.y + 0.4*x.z
            + 0.03*x.y*x.y - 0.01*x.x*x.z,
        0.7 - 0.3*x.x + 0.2*x.y + 0.1*x.z
            + 0.02*x.z*x.z + 0.015*x.x*x.y};
}

std::array<Vec3, 3> exact_velocity_gradient(const Vec3& x)
{
    return {{
        {0.2 + 0.08*x.x, 0.3 + 0.02*x.z, -0.1 + 0.02*x.y},
        {0.1 - 0.01*x.z, -0.2 + 0.06*x.y, 0.4 - 0.01*x.x},
        {-0.3 + 0.015*x.y, 0.2 + 0.015*x.x, 0.1 + 0.04*x.z}}};
}

Vec3 exact_viscous_traction(const std::array<Vec3, 3>& g,
                            const Vec3& body_normal)
{
    return {
        kMu*(2.0*g[0].x*body_normal.x
             + (g[0].y + g[1].x)*body_normal.y
             + (g[0].z + g[2].x)*body_normal.z),
        kMu*((g[1].x + g[0].y)*body_normal.x
             + 2.0*g[1].y*body_normal.y
             + (g[1].z + g[2].y)*body_normal.z),
        kMu*((g[2].x + g[0].z)*body_normal.x
             + (g[2].y + g[1].z)*body_normal.y
             + 2.0*g[2].z*body_normal.z)};
}

double error(const Vec3& a, const Vec3& b)
{
    return (a - b).mag();
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2)
        throw std::invalid_argument(
            "usage: test_naca3d_wall_quadratic_mesh mesh.h5");

    Mesh mesh;
    if (!cfdx::io::read_mesh_hdf5(argv[1], mesh) || !mesh.topo_validate().ok)
        throw std::runtime_error("invalid extruded NACA0012 mesh");

    const auto geometry = build_fv_geometry(mesh);
    const std::size_t airfoil = mesh.boundary().find("airfoil");
    if (airfoil >= mesh.boundary().n_patches() ||
        mesh.boundary().patch(airfoil).type != PatchType::WALL)
        throw std::runtime_error("NACA airfoil WALL patch is missing");

    std::vector<double> pressure(mesh.n_cells());
    std::vector<Vec3> velocity(mesh.n_cells());
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        pressure[cell] = exact_pressure(geometry.cell_centres[cell]);
        velocity[cell] = exact_velocity(geometry.cell_centres[cell]);
    }

    std::vector<double> wall_pressure(mesh.n_faces(), 0.0);
    std::vector<Vec3> wall_velocity(mesh.n_faces());
    for (const auto face : mesh.boundary().patch(airfoil).face_ids) {
        wall_pressure[face] = exact_pressure(geometry.face_centres[face]);
        wall_velocity[face] = exact_velocity(geometry.face_centres[face]);
    }

    // This is a real curved, extruded NACA O-grid, not a Cartesian stencil.
    // Four spanwise layers provide enough independent z coordinates for a
    // full 3-D quadratic basis; rank-deficient local stencils must fail closed.
    const auto surface = reconstruct_quadratic_wall_surface(
        mesh, geometry, airfoil, pressure, velocity,
        wall_pressure, wall_velocity, kMu, 4);
    const auto& faces = mesh.boundary().patch(airfoil).face_ids;
    if (surface.size() != faces.size() || surface.empty())
        throw std::runtime_error("NACA quadratic wall reconstruction returned wrong face count");

    double max_pressure_error = 0.0;
    double max_traction_error = 0.0;
    double total_area = 0.0;
    for (std::size_t i = 0; i < surface.size(); ++i) {
        const auto& sample = surface[i];
        const auto face = faces[i];
        const Vec3 x = geometry.face_centres[face];
        const Vec3 body_normal{sample.normal.x, sample.normal.y, sample.normal.z};
        max_pressure_error = std::max(
            max_pressure_error, std::abs(sample.pressure - exact_pressure(x)));
        max_traction_error = std::max(
            max_traction_error,
            error(Vec3{sample.viscous_traction.x,
                       sample.viscous_traction.y,
                       sample.viscous_traction.z},
                  exact_viscous_traction(exact_velocity_gradient(x), body_normal)));
        total_area += sample.area;
    }

    std::cout << "NACA3D_WALL_QUADRATIC_MESH faces=" << surface.size()
              << " area=" << total_area
              << " max_pressure_error=" << max_pressure_error
              << " max_traction_error=" << max_traction_error << '\n';

    if (!std::isfinite(max_pressure_error) ||
        !std::isfinite(max_traction_error) ||
        max_pressure_error > 2e-8 || max_traction_error > 2e-8)
        throw std::runtime_error(
            "quadratic wall reconstruction failed on curved NACA O-grid");

    forces::Reference reference;
    reference.area = total_area;
    const auto integrated = forces::integrate(surface, reference);
    if (!std::isfinite(integrated.cd) || !std::isfinite(integrated.cmz))
        throw std::runtime_error("NACA quadratic wall force integration is non-finite");

    std::cout << "NACA3D_WALL_QUADRATIC_MESH PASS"
              << " cd=" << integrated.cd
              << " cmz=" << integrated.cmz << '\n';
    return 0;
}
