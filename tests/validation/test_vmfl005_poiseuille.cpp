#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {

Mesh make_pipe(std::size_t ns, std::size_t nz, double R, double L)
{
    Mesh m;
    const std::size_t p_layer = ns + 1;
    m.points().resize(p_layer * (nz + 1));
    constexpr double pi = 3.1415926535897932384626433832795;
    const double dz = L / static_cast<double>(nz);

    for (std::size_t k = 0; k <= nz; ++k) {
        const std::size_t b = k * p_layer;
        const double z = dz * static_cast<double>(k);
        m.points().set(b, 0.0, 0.0, z);
        for (std::size_t s = 0; s < ns; ++s) {
            const double a = 2.0 * pi * static_cast<double>(s) / static_cast<double>(ns);
            m.points().set(b + 1 + s, R * std::cos(a), R * std::sin(a), z);
        }
    }

    std::vector<std::vector<std::size_t>> end_faces(nz + 1, std::vector<std::size_t>(ns));
    std::vector<std::vector<std::size_t>> radial(nz, std::vector<std::size_t>(ns));
    std::vector<std::vector<std::size_t>> wall(nz, std::vector<std::size_t>(ns));

    // End triangles. k=0 points outward in -z; k=nz outward in +z.
    for (std::size_t k = 0; k <= nz; ++k) {
        const std::size_t b = k * p_layer;
        for (std::size_t s = 0; s < ns; ++s) {
            const std::size_t sn = (s + 1) % ns;
            end_faces[k][s] = m.faces().n_faces();
            if (k == 0)
                m.faces().push_face({b, b + 1 + sn, b + 1 + s});
            else
                m.faces().push_face({b, b + 1 + s, b + 1 + sn});
        }
    }

    for (std::size_t k = 0; k < nz; ++k) {
        const std::size_t b0 = k * p_layer;
        const std::size_t b1 = (k + 1) * p_layer;
        for (std::size_t s = 0; s < ns; ++s) {
            const std::size_t sn = (s + 1) % ns;
            // Radial interface at angle s. The ordering gives the outward
            // normal of the owner sector s; the neighbour sector s-1 receives
            // the opposite orientation through CellGeometry's owner/neighbour
            // convention.
            radial[k][s] = m.faces().n_faces();
            m.faces().push_face({b0, b0 + 1 + s, b1 + 1 + s, b1});

            // Circular wall face, outward in the radial direction.
            wall[k][s] = m.faces().n_faces();
            m.faces().push_face({b0 + 1 + s, b0 + 1 + sn,
                                 b1 + 1 + sn, b1 + 1 + s});
        }
    }

    m.ownership().resize(m.n_faces());
    for (std::size_t s = 0; s < ns; ++s) {
        for (std::size_t k = 0; k <= nz; ++k) {
            const auto f = end_faces[k][s];
            if (k == 0) {
                m.ownership().set_owner(f, 0);
                m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
            } else if (k == nz) {
                m.ownership().set_owner(f, (nz - 1) * ns + s);
                m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
            } else {
                m.ownership().set_owner(f, (k - 1) * ns + s);
                m.ownership().set_neighbour(f, static_cast<int>(k * ns + s));
            }
        }
    }
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t s = 0; s < ns; ++s) {
            const auto r = radial[k][s];
            m.ownership().set_owner(r, k * ns + s);
            m.ownership().set_neighbour(r, static_cast<int>(k * ns + ((s + ns - 1) % ns)));
            m.ownership().set_owner(wall[k][s], k * ns + s);
            m.ownership().set_neighbour(wall[k][s], FaceOwnership::BOUNDARY);
        }
    }

    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t s = 0; s < ns; ++s) {
            const std::size_t sn = (s + 1) % ns;
            m.cells().push_cell({
                end_faces[k][s], end_faces[k + 1][s],
                radial[k][s], wall[k][s], radial[k][sn]
            });
        }
    }

    Patch p;
    p.name = "inlet"; p.type = PatchType::INLET;
    p.face_ids = end_faces[0]; m.boundary().add_patch(p);
    p = {}; p.name = "outlet"; p.type = PatchType::OUTLET;
    p.face_ids = end_faces[nz]; m.boundary().add_patch(p);
    std::vector<std::size_t> wall_faces;
    wall_faces.reserve(ns * nz);
    for (const auto& row : wall) wall_faces.insert(wall_faces.end(), row.begin(), row.end());
    p = {}; p.name = "wall"; p.type = PatchType::WALL;
    p.face_ids = std::move(wall_faces); m.boundary().add_patch(p);

    if (!m.ownership().is_consistent(m.n_cells()))
        throw std::runtime_error("VMFL005: inconsistent face ownership");
    return m;
}

void check_case(std::size_t nz)
{
    constexpr double R = 0.00125;
    constexpr double L = 0.1;
    constexpr double rho = 1.0;
    constexpr double mu = 1.0e-5;
    constexpr double dp = 10.24;
    constexpr std::size_t ns = 32;
    constexpr double pi = 3.1415926535897932384626433832795;

    Mesh mesh = make_pipe(ns, nz, R, L);
    Field<double,Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
    Field<double,Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
    U.fill(0.0); p.fill(0.0);

    VelocityBoundaryConditions ubc;
    ubc["inlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0.0,0.0,0.0}};
    ubc["outlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0.0,0.0,0.0}};
    ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};

    ScalarBoundaryConditions pbc;
    pbc["inlet"] = {ScalarBoundaryType::FIXED_VALUE,dp,0.0};
    pbc["outlet"] = {ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
    pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

    IncompressibleSolverControls c;
    c.algorithm = PressureVelocityAlgorithm::SIMPLE;
    c.convergence.max_iterations = 500;
    c.convergence.relative_tolerance = 1e-8;
    c.convergence.continuity_tolerance = 1e-10;
    c.linear_max_iterations = 5000;
    c.linear_tolerance = 1e-8;
    c.density = rho;
    c.kinematic_viscosity = mu / rho;

    const auto result = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
    if (!result.converged)
        throw std::runtime_error("VMFL005 coupled solver did not converge");

    const auto geometry = build_fv_geometry(mesh);
    double flow_rate = 0.0;
    double volume = 0.0;
    double max_profile_error = 0.0;
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        const double r = std::hypot(geometry.cell_centres[cell].x,
                                    geometry.cell_centres[cell].y);
        const double exact = dp * (R * R - r * r) / (4.0 * mu * L);
        max_profile_error = std::max(max_profile_error, std::abs(U(cell,0) - exact));
        flow_rate += U(cell,0) * geometry.cell_volumes[cell];
        volume += geometry.cell_volumes[cell];
    }
    const double mean_u = flow_rate / volume;
    const double q_exact = pi * std::pow(R,4) * dp / (8.0 * mu * L);
    const double mean_exact = 2.0;
    const double rel_q = std::abs(flow_rate - q_exact) / q_exact;
    const double rel_mean = std::abs(mean_u - mean_exact) / mean_exact;
    const auto& last = result.history.back();

    std::cout << "VMFL005 nz=" << nz
              << " flow_rate=" << flow_rate
              << " q_rel_error=" << rel_q
              << " mean_velocity=" << mean_u
              << " mean_rel_error=" << rel_mean
              << " profile_abs_error=" << max_profile_error
              << " continuity_linf=" << last.continuity_linf
              << " iterations=" << result.iterations << "\n";

    const double max_exact = dp * R * R / (4.0 * mu * L);
    const double rel_profile = max_profile_error / max_exact;
    if (rel_q > 5e-2 || rel_mean > 5e-2 || rel_profile > 5e-2)
        throw std::runtime_error("VMFL005 quantitative mismatch");
    if (last.continuity_linf > 1e-8)
        throw std::runtime_error("VMFL005 continuity error too large");
}

} // namespace

int main()
{
    try {
        check_case(16);
        check_case(32);
        std::cout << "VMFL005_VALIDATION: PASS\\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "VMFL005_VALIDATION: FAIL: " << e.what() << "\n";
        return 1;
    }
}
