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

Mesh make_pipe(std::size_t ns, std::size_t nz, double R, double L, std::size_t nr)
{
    Mesh m;
    const std::size_t layer_points = 1 + nr * ns;
    m.points().resize(layer_points * (nz + 1));
    constexpr double pi = 3.1415926535897932384626433832795;
    const double dz = L / static_cast<double>(nz);

    auto point = [&](std::size_t k, std::size_t j, std::size_t s) -> std::size_t {
        const std::size_t base = k * layer_points;
        if (j == 0) return base;
        return base + 1 + (j - 1) * ns + s;
    };

    for (std::size_t k = 0; k <= nz; ++k) {
        const double z = dz * static_cast<double>(k);
        m.points().set(point(k, 0, 0), 0.0, 0.0, z);
        for (std::size_t j = 1; j <= nr; ++j) {
            const double r = R * static_cast<double>(j) / static_cast<double>(nr);
            for (std::size_t s = 0; s < ns; ++s) {
                const double a = 2.0 * pi * static_cast<double>(s) / static_cast<double>(ns);
                m.points().set(point(k, j, s), r * std::cos(a), r * std::sin(a), z);
            }
        }
    }

    using EndFaces = std::vector<std::vector<std::vector<std::size_t>>>;
    using AxialFaces = std::vector<std::vector<std::vector<std::size_t>>>;
    EndFaces end_faces(nz + 1,
                       std::vector<std::vector<std::size_t>>(nr,
                           std::vector<std::size_t>(ns)));
    AxialFaces angular(nz,
                      std::vector<std::vector<std::size_t>>(nr,
                          std::vector<std::size_t>(ns)));
    AxialFaces radial(nz,
                     std::vector<std::vector<std::size_t>>(nr - 1,
                         std::vector<std::size_t>(ns)));
    std::vector<std::vector<std::size_t>> wall(nz, std::vector<std::size_t>(ns));

    // End faces: bottom outward in -z, top outward in +z.
    for (std::size_t k = 0; k <= nz; ++k) {
        for (std::size_t j = 0; j < nr; ++j) {
            for (std::size_t s = 0; s < ns; ++s) {
                const std::size_t sn = (s + 1) % ns;
                const auto c = point(k, 0, 0);
                const auto a = point(k, j + 1, s);
                const auto b = point(k, j + 1, sn);
                end_faces[k][j][s] = m.faces().n_faces();
                if (j == 0) {
                    if (k == 0)
                        m.faces().push_face({c, b, a});
                    else
                        m.faces().push_face({c, a, b});
                } else if (k == 0) {
                    // Bottom annulus: outward normal is -z.
                    m.faces().push_face({point(k, j, s), point(k, j, sn),
                                         point(k, j + 1, sn), point(k, j + 1, s)});
                } else {
                    // Top annulus: outward normal is +z.
                    m.faces().push_face({point(k, j, s), point(k, j + 1, s),
                                         point(k, j + 1, sn), point(k, j, sn)});
                }
            }
        }
    }

    // Angular faces are radial planes. For sector s the lower face has an
    // outward normal in the -tangential direction; the upper face is shared
    // with sector s+1 and gets the opposite local orientation automatically.
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < nr; ++j) {
            for (std::size_t s = 0; s < ns; ++s) {
                const auto lower = point(k, j == 0 ? 0 : j, s);
                const auto upper = point(k, j + 1, s);
                angular[k][j][s] = m.faces().n_faces();
                m.faces().push_face({lower, upper, point(k + 1, j + 1, s),
                                     point(k + 1, j == 0 ? 0 : j, s)});
            }
        }
    }

    // Radial interfaces at ring j+1. The outer cell is the owner and the
    // face is oriented inward, i.e. outward from that owner cell.
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j + 1 < nr; ++j) {
            const std::size_t ring = j + 1;
            for (std::size_t s = 0; s < ns; ++s) {
                const std::size_t sn = (s + 1) % ns;
                radial[k][j][s] = m.faces().n_faces();
                m.faces().push_face({point(k, ring, s), point(k + 1, ring, s),
                                     point(k + 1, ring, sn), point(k, ring, sn)});
            }
        }

        for (std::size_t s = 0; s < ns; ++s) {
            const std::size_t sn = (s + 1) % ns;
            wall[k][s] = m.faces().n_faces();
            m.faces().push_face({point(k, nr, s), point(k, nr, sn),
                                 point(k + 1, nr, sn), point(k + 1, nr, s)});
        }
    }

    m.ownership().resize(m.n_faces());
    const std::size_t cells_per_plane = nr * ns;
    auto cell_id = [&](std::size_t k, std::size_t j, std::size_t s) {
        return static_cast<CellIndex>(k * cells_per_plane + j * ns + s);
    };

    for (std::size_t k = 0; k <= nz; ++k) {
        for (std::size_t j = 0; j < nr; ++j) {
            for (std::size_t s = 0; s < ns; ++s) {
                const auto f = end_faces[k][j][s];
                if (k == 0) {
                    m.ownership().set_owner(f, cell_id(0, j, s));
                    m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
                } else if (k == nz) {
                    m.ownership().set_owner(f, cell_id(nz - 1, j, s));
                    m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
                } else {
                    m.ownership().set_owner(f, cell_id(k - 1, j, s));
                    m.ownership().set_neighbour(f, static_cast<std::int64_t>(cell_id(k, j, s)));
                }
            }
        }
    }

    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < nr; ++j) {
            for (std::size_t s = 0; s < ns; ++s) {
                const auto lower = angular[k][j][s];
                const auto upper = angular[k][j][(s + 1) % ns];
                m.ownership().set_owner(lower, cell_id(k, j, s));
                m.ownership().set_neighbour(
                    lower, static_cast<std::int64_t>(cell_id(k, j, (s + ns - 1) % ns)));
                m.ownership().set_owner(upper, cell_id(k, j, (s + 1) % ns));
                m.ownership().set_neighbour(
                    upper, static_cast<std::int64_t>(cell_id(k, j, s)));

                if (j + 1 < nr) {
                    const auto r = radial[k][j][s];
                    m.ownership().set_owner(r, cell_id(k, j + 1, s));
                    m.ownership().set_neighbour(
                        r, static_cast<std::int64_t>(cell_id(k, j, s)));
                } else {
                    const auto w = wall[k][s];
                    m.ownership().set_owner(w, cell_id(k, j, s));
                    m.ownership().set_neighbour(w, FaceOwnership::BOUNDARY);
                }
            }
        }
    }

    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < nr; ++j) {
            for (std::size_t s = 0; s < ns; ++s) {
                std::vector<FaceIndex> faces = {
                    end_faces[k][j][s], end_faces[k + 1][j][s],
                    angular[k][j][s], angular[k][j][(s + 1) % ns]
                };
                if (j > 0) faces.push_back(radial[k][j - 1][s]);
                if (j + 1 < nr) faces.push_back(radial[k][j][s]);
                else faces.push_back(wall[k][s]);
                m.cells().push_cell(faces);
            }
        }
    }

    Patch p;
    p.name = "inlet"; p.type = PatchType::INLET;
    for (std::size_t j = 0; j < nr; ++j)
        for (const auto f : end_faces[0][j]) p.face_ids.push_back(f);
    m.boundary().add_patch(p);
    p = {}; p.name = "outlet"; p.type = PatchType::OUTLET;
    for (std::size_t j = 0; j < nr; ++j)
        for (const auto f : end_faces[nz][j]) p.face_ids.push_back(f);
    m.boundary().add_patch(p);
    p = {}; p.name = "wall"; p.type = PatchType::WALL;
    p.face_ids.reserve(ns * nz);
    for (const auto& row : wall)
        p.face_ids.insert(p.face_ids.end(), row.begin(), row.end());
    m.boundary().add_patch(p);

    const auto topology = m.topo_validate();
    if (!topology.ok) {
        throw std::runtime_error("VMFL005 topology invalid: " + topology.errors.front());
    }
    return m;
}
void check_case(std::size_t nz)
{
    constexpr double R = 0.00125;
    constexpr double L = 0.1;
    constexpr double rho = 1.0;
    constexpr double mu = 1.0e-5;
    constexpr double dp = 0.1024;
    constexpr std::size_t ns = 32;
    constexpr std::size_t nr = 8;
    constexpr double pi = 3.1415926535897932384626433832795;

    Mesh mesh = make_pipe(ns, nz, R, L, nr);
    Field<double,Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
    Field<double,Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
    U.fill(0.0); p.fill(0.0);

    VelocityBoundaryConditions ubc;
    ubc["inlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0.0,0.0,0.0}};
    ubc["outlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0.0,0.0,0.0}};
    ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};

    ScalarBoundaryConditions pbc;
    pbc["inlet"] = {ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
    pbc["outlet"] = {ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
    pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

    IncompressibleSolverControls c;
    c.algorithm = PressureVelocityAlgorithm::SIMPLE;
    c.active_velocity_components = {false, false, true};
    c.convergence.max_iterations = 500;
    c.convergence.relative_tolerance = 1e-8;
    c.convergence.continuity_tolerance = 1e-10;
    c.linear_max_iterations = 5000;
    c.linear_tolerance = 1e-8;
    c.density = rho;
    c.kinematic_viscosity = mu / rho;
    c.body_force = {0.0, 0.0, dp / L};
    // The Poiseuille case is laminar and pressure-driven; disable the
    // bounded-convection correction so the validation isolates the viscous
    // pressure-driven solution without a nonlinear convective stabilization
    // term dominating the nearly zero transverse fluxes.
    c.use_bounded_convection = false;

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
        max_profile_error = std::max(max_profile_error, std::abs(U(cell,2) - exact));
        flow_rate += U(cell,2) * geometry.cell_volumes[cell];
        volume += geometry.cell_volumes[cell];
    }
    const double area = volume / L;
    flow_rate /= L;
    const double mean_u = flow_rate / area;
    const double q_exact = pi * std::pow(R,4) * dp / (8.0 * mu * L);
    const double mean_exact = 0.2;
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
