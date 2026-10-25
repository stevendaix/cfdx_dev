#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <unordered_map>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_four_cell_stress_channel()
{
    Mesh mesh;
    constexpr std::size_t nx = 4;
    mesh.points().resize((nx + 1) * 4);
    for (std::size_t i = 0; i <= nx; ++i) {
        const double x = static_cast<double>(i);
        const double xyz[4][3] = {
            {x,0,0}, {x,1,0}, {x,1,1}, {x,0,1}};
        for (std::size_t j = 0; j < 4; ++j)
            mesh.points().set(4*i+j, xyz[j][0], xyz[j][1], xyz[j][2]);
    }

    std::vector<std::size_t> x_faces(nx + 1);
    std::vector<std::array<std::size_t, 4>> side_faces(nx);
    for (std::size_t i = 0; i <= nx; ++i) {
        const std::size_t b = 4*i;
        mesh.faces().push_face({b, b+1, b+2, b+3});
        x_faces[i] = i;
    }
    for (std::size_t i = 0; i < nx; ++i) {
        const std::size_t l = 4*i, r = 4*(i+1);
        const std::size_t first = nx + 1 + 4*i;
        mesh.faces().push_face({l, r, r+3, l+3});       // y = 0
        mesh.faces().push_face({l+1, l+2, r+2, r+1});   // y = 1
        mesh.faces().push_face({l, l+1, r+1, r});       // z = 0
        mesh.faces().push_face({l+3, r+3, r+2, l+2});   // z = 1
        side_faces[i] = {first, first+1, first+2, first+3};
    }

    const std::size_t n_faces = (nx + 1) + 4*nx;
    mesh.ownership().resize(n_faces);
    for (std::size_t face = 0; face < n_faces; ++face)
        mesh.ownership().set_neighbour(face, FaceOwnership::BOUNDARY);
    for (std::size_t i = 0; i <= nx; ++i) {
        mesh.ownership().set_owner(x_faces[i], i == 0 ? 0 : i-1);
        if (i > 0 && i < nx)
            mesh.ownership().set_neighbour(x_faces[i], static_cast<std::int64_t>(i));
    }
    for (std::size_t i = 0; i < nx; ++i) {
        const auto& sf = side_faces[i];
        mesh.ownership().set_owner(sf[0], i);
        mesh.ownership().set_owner(sf[1], i);
        mesh.ownership().set_owner(sf[2], i);
        mesh.ownership().set_owner(sf[3], i);
        mesh.cells().push_cell({
            x_faces[i], x_faces[i+1], sf[0], sf[1], sf[2], sf[3]});
    }

    const auto add_wall_patch = [&](const std::string& name,
                                    std::vector<std::size_t> faces) {
        Patch patch;
        patch.name = name;
        patch.type = PatchType::WALL;
        patch.face_ids = std::move(faces);
        mesh.boundary().add_patch(patch);
    };
    add_wall_patch("x_min", {x_faces.front()});
    add_wall_patch("x_max", {x_faces.back()});
    std::vector<std::size_t> y_min, y_max, z_min, z_max;
    for (const auto& sf : side_faces) {
        y_min.push_back(sf[0]);
        y_max.push_back(sf[1]);
        z_min.push_back(sf[2]);
        z_max.push_back(sf[3]);
    }
    add_wall_patch("y_min", std::move(y_min));
    add_wall_patch("y_max", std::move(y_max));
    add_wall_patch("z_min", std::move(z_min));
    add_wall_patch("z_max", std::move(z_max));
    return mesh;
}


static Mesh make_skewed_stress_cube(std::size_t n, double shear)
{
    Mesh mesh;
    const auto point_id = [n](std::size_t i, std::size_t j, std::size_t k) {
        return (k * (n + 1) + j) * (n + 1) + i;
    };
    mesh.points().resize((n + 1) * (n + 1) * (n + 1));
    for (std::size_t k = 0; k <= n; ++k)
        for (std::size_t j = 0; j <= n; ++j)
            for (std::size_t i = 0; i <= n; ++i) {
                const auto stretched_coordinate = [n](std::size_t index) {
                    const double eta = static_cast<double>(index) / static_cast<double>(n);
                    return eta + 0.20 * eta * eta;
                };
                const double x = stretched_coordinate(i);
                const double y = stretched_coordinate(j);
                const double z = stretched_coordinate(k);
                mesh.points().set(point_id(i, j, k), x + shear * y, y, z);
            }

    struct FaceKey {
        std::size_t i, j, k, axis;
        bool operator==(const FaceKey& other) const {
            return i == other.i && j == other.j && k == other.k && axis == other.axis;
        }
    };
    struct FaceKeyHash {
        std::size_t operator()(const FaceKey& f) const {
            std::size_t h = f.i;
            h = h * 1315423911u + f.j;
            h = h * 1315423911u + f.k;
            h = h * 1315423911u + f.axis;
            return h;
        }
    };
    const std::size_t nc = n * n * n;
    const auto cell_id = [n](std::size_t i, std::size_t j, std::size_t k) {
        return (k * n + j) * n + i;
    };
    std::unordered_map<FaceKey, std::size_t, FaceKeyHash> face_map;
    std::vector<std::vector<std::size_t>> faces_per_cell(nc);
    const auto add_face = [&](FaceKey key, std::initializer_list<std::size_t> vertices) {
        const auto found = face_map.find(key);
        if (found != face_map.end()) return found->second;
        const std::size_t id = mesh.faces().n_faces();
        mesh.faces().push_face(std::vector<FaceIndex>(vertices.begin(), vertices.end()));
        face_map.emplace(key, id);
        return id;
    };
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i) {
                const auto p = [&](std::size_t pi, std::size_t pj, std::size_t pk) {
                    return point_id(pi, pj, pk);
                };
                faces_per_cell[cell_id(i,j,k)] = {
                    add_face({i,j,k,0}, {p(i,j,k),p(i,j,k+1),p(i,j+1,k+1),p(i,j+1,k)}),
                    add_face({i+1,j,k,0}, {p(i+1,j,k),p(i+1,j+1,k),p(i+1,j+1,k+1),p(i+1,j,k+1)}),
                    add_face({i,j,k,1}, {p(i,j,k),p(i+1,j,k),p(i+1,j,k+1),p(i,j,k+1)}),
                    add_face({i,j+1,k,1}, {p(i,j+1,k),p(i,j+1,k+1),p(i+1,j+1,k+1),p(i+1,j+1,k)}),
                    add_face({i,j,k,2}, {p(i,j,k),p(i,j+1,k),p(i+1,j+1,k),p(i+1,j,k)}),
                    add_face({i,j,k+1,2}, {p(i,j,k+1),p(i+1,j,k+1),p(i+1,j+1,k+1),p(i,j+1,k+1)})
                };
            }

    mesh.ownership().resize(mesh.n_faces());
    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    for (const auto& entry : face_map) {
        const FaceKey key = entry.first;
        const std::size_t face = entry.second;
        const bool lower = key.axis == 0 ? key.i == 0 :
                           key.axis == 1 ? key.j == 0 : key.k == 0;
        const bool upper = key.axis == 0 ? key.i == n :
                           key.axis == 1 ? key.j == n : key.k == n;
        const std::size_t oi = key.axis == 0 ? (lower ? 0 : key.i - 1) : key.i;
        const std::size_t oj = key.axis == 1 ? (lower ? 0 : key.j - 1) : key.j;
        const std::size_t ok = key.axis == 2 ? (lower ? 0 : key.k - 1) : key.k;
        const std::size_t owner = cell_id(oi, oj, ok);
        mesh.ownership().set_owner(face, owner);
        if (lower || upper) {
            mesh.ownership().set_neighbour(face, FaceOwnership::BOUNDARY);
            wall.face_ids.push_back(face);
        } else {
            std::size_t ni = oi, nj = oj, nk = ok;
            if (key.axis == 0) ++ni;
            else if (key.axis == 1) ++nj;
            else ++nk;
            mesh.ownership().set_neighbour(face,
                static_cast<std::int64_t>(cell_id(ni, nj, nk)));
        }
    }
    for (std::size_t cell = 0; cell < nc; ++cell)
        mesh.cells().push_cell(faces_per_cell[cell]);
    mesh.boundary().add_patch(wall);
    return mesh;
}

static VelocityBoundaryConditions affine_velocity_bcs(double a)
{
    VelocityBoundaryConditions bcs;
    bcs["x_min"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.5*a, 0.0, 0.0}};
    bcs["x_max"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.5*a, 0.0, 0.0}};
    bcs["y_min"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
    bcs["y_max"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {a, 0.0, 0.0}};
    bcs["z_min"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.5*a, 0.0, 0.0}};
    bcs["z_max"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.5*a, 0.0, 0.0}};
    return bcs;
}

int main()
{
    run_case("variable_viscosity_transpose_stress_affine_face_flux", [] {
        constexpr double a = 2.5;
        const Mesh mesh = make_four_cell_stress_channel();
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        std::vector<double> mu(mesh.n_cells(), 0.0);
        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
            U.set(cell, a*geometry.cell_centres[cell].y, 0.0, 0.0);
            mu[cell] = 1.0 + geometry.cell_centres[cell].x;
        }
        const auto bcs = affine_velocity_bcs(a);

        const auto face_flux = momentum_transpose_diffusion_face_fluxes(
            mesh, geometry, U, bcs, mu, 1);
        // For U=(a*y,0,0), div(U)=0 and the y-momentum transpose-stress
        // flux through an x-normal face is mu(x_face)*a*S_x. Linear mu(x)
        // is reconstructed exactly at the three internal face centres.
        for (std::size_t face = 1; face <= 3; ++face) {
            const double expected =
                (1.0 + geometry.face_centres[face].x) * a *
                geometry.face_area_vectors[face].x;
            EXPECT_NEAR(face_flux[face], expected, 1e-12);
        }

        const auto rhs = momentum_transpose_diffusion_rhs(
            mesh, geometry, U, bcs, mu, 1);
        // Cells 1 and 2 are internal in x. Their y-momentum stress divergence
        // is d/dx(mu*a)=a; side-wall stress flux has S_x=0.
        EXPECT_NEAR(rhs[1], a*geometry.cell_volumes[1], 1e-12);
        EXPECT_NEAR(rhs[2], a*geometry.cell_volumes[2], 1e-12);
        double rhs_sum = 0.0;
        double boundary_flux_sum = 0.0;
        for (const double value : rhs) rhs_sum += value;
        for (std::size_t face = 0; face < mesh.n_faces(); ++face)
            if (mesh.ownership().neighbour(face) < 0)
                boundary_flux_sum += face_flux[face];
        EXPECT_NEAR(rhs_sum, boundary_flux_sum, 1e-12);
    });

    run_case("variable_viscosity_transpose_stress_enters_momentum_rhs", [] {
        constexpr double a = 2.5;
        const Mesh mesh = make_four_cell_stress_channel();
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> U_y(mesh.n_cells(), "Uy", "m/s", 1);
        Field<double, Location::CELL> grad_p(mesh.n_cells(), "grad_p", "Pa/m", 3);
        Field<double, Location::CELL> body(mesh.n_cells(), "body", "N/m3", 1);
        Field<double, Location::FACE> mass_flux(mesh.n_faces(), "phi", "kg/s", 1);
        std::vector<double> mu(mesh.n_cells(), 0.0);
        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
            U.set(cell, a*geometry.cell_centres[cell].y, 0.0, 0.0);
            U_y(cell) = 0.0;
            grad_p.set(cell, 0.0, 0.0, 0.0);
            body(cell) = 0.0;
            mu[cell] = 1.0 + geometry.cell_centres[cell].x;
        }
        mass_flux.fill(0.0);
        const auto vector_bcs = affine_velocity_bcs(a);
        ScalarBoundaryConditions y_bcs;
        for (const auto& [name, bc] : vector_bcs)
            y_bcs[name] = {ScalarBoundaryType::FIXED_VALUE, bc.value.y, 0.0};

        const auto equation = assemble_momentum_component_with_transpose_diffusion(
            mesh, geometry, mass_flux, grad_p, body, mu, y_bcs, U,
            vector_bcs, 1, false, ConvectionScheme::UPWIND, &U_y);
        EXPECT_NEAR(equation.rhs(1), a*geometry.cell_volumes[1], 1e-12);
        EXPECT_NEAR(equation.rhs(2), a*geometry.cell_volumes[2], 1e-12);
    });


    run_case("variable_viscosity_transpose_stress_skewed_nonuniform_affine", [] {
        constexpr std::size_t n = 4;
        constexpr double shear = 0.65;
        constexpr double a = 2.5;
        const Mesh mesh = make_skewed_stress_cube(n, shear);
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        std::vector<double> mu(mesh.n_cells(), 0.0);
        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
            U.set(cell, a * geometry.cell_centres[cell].y, 0.0, 0.0);
            mu[cell] = 1.0 + geometry.cell_centres[cell].x;
        }
        VelocityBoundaryConditions bcs;
        bcs["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};

        const auto flux = momentum_transpose_diffusion_face_fluxes(
            mesh, geometry, U, bcs, mu, 1);
        for (std::size_t face = 0; face < mesh.n_faces(); ++face) {
            if (mesh.ownership().neighbour(face) < 0) continue;
            const double expected =
                (1.0 + geometry.face_centres[face].x) * a *
                geometry.face_area_vectors[face].x;
            EXPECT_NEAR(flux[face], expected, 2e-11);
        }

        const auto rhs = momentum_transpose_diffusion_rhs(
            mesh, geometry, U, bcs, mu, 1);
        const auto cell_id = [](std::size_t i, std::size_t j, std::size_t k) {
            return (k * n + j) * n + i;
        };
        for (std::size_t k = 1; k + 1 < n; ++k)
            for (std::size_t j = 1; j + 1 < n; ++j)
                for (std::size_t i = 1; i + 1 < n; ++i) {
                    const std::size_t cell = cell_id(i, j, k);
                    EXPECT_NEAR(rhs[cell], a * geometry.cell_volumes[cell], 2e-11);
                }
    });


    run_case("variable_viscosity_transpose_stress_smooth_mms_order", [] {
        const double pi = std::acos(-1.0);
        std::vector<double> errors;
        for (const std::size_t n : {4u, 8u, 16u}) {
            const Mesh mesh = make_skewed_stress_cube(n, 0.65);
            const auto geometry = build_fv_geometry(mesh);
            Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
            std::vector<double> mu(mesh.n_cells(), 0.0);
            for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
                const double y = geometry.cell_centres[cell].y;
                U.set(cell, std::sin(pi * y), 0.0, 0.0);
                mu[cell] = 1.0 + geometry.cell_centres[cell].x;
            }
            VelocityBoundaryConditions bcs;
            bcs["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
            const auto rhs = momentum_transpose_diffusion_rhs(
                mesh, geometry, U, bcs, mu, 1);
            double squared_error = 0.0;
            std::size_t samples = 0;
            const auto cell_id = [n](std::size_t i, std::size_t j, std::size_t k) {
                return (k * n + j) * n + i;
            };
            for (std::size_t k = 1; k + 1 < n; ++k)
                for (std::size_t j = 1; j + 1 < n; ++j)
                    for (std::size_t i = 1; i + 1 < n; ++i) {
                        const std::size_t cell = cell_id(i, j, k);
                        const double numerical = rhs[cell] / geometry.cell_volumes[cell];
                        const double exact = pi * std::cos(pi * geometry.cell_centres[cell].y);
                        const double error = numerical - exact;
                        squared_error += error * error;
                        ++samples;
                    }
            if (samples == 0)
                throw std::runtime_error("smooth transpose-stress MMS has no interior cells");
            errors.push_back(std::sqrt(squared_error / static_cast<double>(samples)));
        }
        EXPECT_TRUE(errors[2] < errors[0]);
        const double order_4_8 = std::log(errors[0] / errors[1]) / std::log(2.0);
        const double order_8_16 = std::log(errors[1] / errors[2]) / std::log(2.0);
        EXPECT_TRUE(std::isfinite(order_4_8));
        EXPECT_TRUE(std::isfinite(order_8_16));
        EXPECT_TRUE(order_4_8 > 1.0);
        EXPECT_TRUE(order_8_16 > 1.0);
    });

    return run_all();
}
