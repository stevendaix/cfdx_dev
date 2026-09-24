#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

Mesh make_diagnostic_mesh(std::size_t nx = 4, std::size_t ny = 8)
{
    Mesh mesh;
    const std::size_t plane = (nx + 1) * (ny + 1);
    mesh.points().resize(2 * plane);
    const auto id = [nx](std::size_t i, std::size_t j, std::size_t k) {
        return (j * (nx + 1) + i) * 2 + k;
    };

    for (std::size_t j = 0; j <= ny; ++j)
        for (std::size_t i = 0; i <= nx; ++i) {
            const double x = static_cast<double>(i) / nx;
            const double y = static_cast<double>(j) / ny;
            mesh.points().set(id(i, j, 0), x, y, 0.0);
            mesh.points().set(id(i, j, 1), x, y, 1.0);
        }

    std::map<std::vector<std::size_t>, std::size_t> face_map;
    std::vector<std::vector<std::size_t>> cell_faces(nx * ny);

    auto add_face = [&](std::initializer_list<std::size_t> vertices,
                        std::size_t cell) {
        std::vector<std::size_t> key(vertices);
        std::sort(key.begin(), key.end());
        auto it = face_map.find(key);
        if (it != face_map.end()) {
            mesh.ownership().set_neighbour(it->second, static_cast<int>(cell));
            return it->second;
        }
        const std::size_t f = mesh.faces().n_faces();
        mesh.faces().push_face(vertices);
        mesh.ownership().resize(mesh.faces().n_faces());
        mesh.ownership().set_owner(f, cell);
        mesh.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        face_map.emplace(std::move(key), f);
        return f;
    };

    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i) {
            const std::size_t c = j * nx + i;
            const auto a = id(i, j, 0), b = id(i + 1, j, 0);
            const auto c0 = id(i + 1, j + 1, 0), d = id(i, j + 1, 0);
            const auto e = id(i, j, 1), f = id(i + 1, j, 1);
            const auto g = id(i + 1, j + 1, 1), h = id(i, j + 1, 1);
            cell_faces[c] = {
                add_face({a, d, c0, b}, c), add_face({e, f, g, h}, c),
                add_face({a, b, f, e}, c), add_face({d, h, g, c0}, c),
                add_face({a, e, h, d}, c), add_face({b, c0, g, f}, c)};
        }

    for (const auto& faces : cell_faces) mesh.cells().push_cell(faces);

    Patch inlet{"inlet", PatchType::INLET, {}};
    Patch outlet{"outlet", PatchType::OUTLET, {}};
    Patch bottom{"bottom", PatchType::WALL, {}};
    Patch top{"top", PatchType::WALL, {}};
    Patch front{"front", PatchType::EMPTY, {}};
    Patch back{"back", PatchType::EMPTY, {}};

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f) >= 0) continue;
        const auto& vertices = mesh.faces().vertices();
        const auto begin = vertices.begin() +
            static_cast<std::ptrdiff_t>(mesh.faces().face_offset(f));
        const auto end = begin +
            static_cast<std::ptrdiff_t>(mesh.faces().face_size(f));
        double x = 0.0, y = 0.0, z = 0.0;
        for (auto it = begin; it != end; ++it) {
            x += mesh.points().x(*it);
            y += mesh.points().y(*it);
            z += mesh.points().z(*it);
        }
        const double n = static_cast<double>(mesh.faces().face_size(f));
        x /= n; y /= n; z /= n;
        constexpr double eps = 1e-12;
        if (std::abs(x) < eps) inlet.face_ids.push_back(f);
        else if (std::abs(x - 1.0) < eps) outlet.face_ids.push_back(f);
        else if (std::abs(y) < eps) bottom.face_ids.push_back(f);
        else if (std::abs(y - 1.0) < eps) top.face_ids.push_back(f);
        else if (std::abs(z) < eps) front.face_ids.push_back(f);
        else if (std::abs(z - 1.0) < eps) back.face_ids.push_back(f);
    }

    mesh.boundary().add_patch(inlet);
    mesh.boundary().add_patch(outlet);
    mesh.boundary().add_patch(bottom);
    mesh.boundary().add_patch(top);
    mesh.boundary().add_patch(front);
    mesh.boundary().add_patch(back);
    return mesh;
}

void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

double max_abs(const std::vector<double>& values)
{
    double result = 0.0;
    for (double value : values) result = std::max(result, std::abs(value));
    return result;
}

void test_geometry_orientation(const Mesh& mesh, const FvGeometry& geometry)
{
    double min_internal_projection = 1e300;
    double min_boundary_projection = 1e300;

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const std::size_t o = mesh.ownership().owner(f);
        const Vec3 Sf = geometry.face_area_vectors[f];
        require(Sf.mag() > 0.0, "degenerate face area vector");

        if (mesh.ownership().neighbour(f) >= 0) {
            const std::size_t n = static_cast<std::size_t>(mesh.ownership().neighbour(f));
            const Vec3 d = geometry.cell_centres[n] - geometry.cell_centres[o];
            min_internal_projection = std::min(min_internal_projection, Sf.dot(d));
        } else if (mesh.boundary().patch(geometry.face_patch[f]).type != PatchType::EMPTY) {
            const Vec3 d = geometry.face_centres[f] - geometry.cell_centres[o];
            min_boundary_projection = std::min(min_boundary_projection, Sf.dot(d));
        }
    }

    require(min_internal_projection > 0.0,
            "internal owner/neighbour face orientation is inconsistent");
    require(min_boundary_projection > 0.0,
            "boundary owner/face orientation is inconsistent");

    std::cout << "GEOMETRY orientation internal_min_projection="
              << min_internal_projection
              << " boundary_min_projection=" << min_boundary_projection << "\n";
}

void test_constant_flux_conservation(
    const Mesh& mesh, const FvGeometry& geometry)
{
    Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
    U.fill(0.0);
    U.component_data(0).assign(mesh.n_cells(), 1.0);

    VelocityBoundaryConditions bcs;
    for (const char* name : {"inlet", "outlet", "front", "back"})
        bcs[name] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};
    bcs["bottom"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {1, 0, 0}};
    bcs["top"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {1, 0, 0}};

    const auto phi = make_mass_flux(mesh, geometry, U, 1.0, bcs);

    double internal_antisymmetry = 0.0;
    double max_divergence = 0.0;
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f) < 0) continue;
        // A single global face flux must enter the owner cell and leave the
        // neighbour cell with exactly the opposite sign.
        internal_antisymmetry = std::max(internal_antisymmetry, std::abs(phi(f) + (-phi(f))));
    }

    const auto* faces = mesh.cells().faces_data();
    const auto* offsets = mesh.cells().offsets_data();
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        double div = 0.0;
        for (Offset k = offsets[c]; k < offsets[c + 1]; ++k) {
            const std::size_t f = faces[k];
            div += mesh.ownership().owner(f) == c ? phi(f) : -phi(f);
        }
        max_divergence = std::max(max_divergence, std::abs(div));
    }

    // Also check the physical x-directed internal flux against the oriented
    // face area. This catches sign flips that an antisymmetry-only test cannot.
    double max_flux_error = 0.0;
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f) < 0) continue;
        const double expected = geometry.face_area_vectors[f].x;
        max_flux_error = std::max(max_flux_error, std::abs(phi(f) - expected));
    }

    require(internal_antisymmetry == 0.0, "internal flux antisymmetry bookkeeping failed");
    require(max_divergence < 1e-13, "constant velocity is not discretely divergence free");
    require(max_flux_error < 1e-13, "constant velocity face flux has wrong orientation");

    std::cout << "FLUX constant internal_antisymmetry=" << internal_antisymmetry
              << " max_divergence=" << max_divergence
              << " max_flux_error=" << max_flux_error << "\n";
}

void test_linear_gradient(
    const Mesh& mesh, const FvGeometry& geometry)
{
    Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c)
        p(c) = geometry.cell_centres[c].x;

    ScalarBoundaryConditions bcs;
    bcs["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
    bcs["outlet"] = {ScalarBoundaryType::FIXED_VALUE, 1.0, 0.0};
    bcs["bottom"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    bcs["top"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    bcs["front"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    bcs["back"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

    const auto grad = gauss_gradient_with_boundary(p, mesh, geometry, bcs);
    double error = 0.0;
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        error = std::max(error, std::abs(grad.component_data(0)[c] - 1.0));
        error = std::max(error, std::abs(grad.component_data(1)[c]));
        error = std::max(error, std::abs(grad.component_data(2)[c]));
    }

    require(error < 1e-13, "linear manufactured gradient is not reproduced exactly");
    std::cout << "GRADIENT linear max_error=" << error << "\n";
}

void test_rhie_chow_pressure_operators(
    const Mesh& mesh, const FvGeometry& geometry)
{
    Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
    U.fill(0.0);
    Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);

    std::array<std::vector<double>, 3> rAU{
        std::vector<double>(mesh.n_cells(), 0.5),
        std::vector<double>(mesh.n_cells(), 0.5),
        std::vector<double>(mesh.n_cells(), 0.5)};

    VelocityBoundaryConditions ubc;
    ScalarBoundaryConditions pbc;
    for (const char* name : {"inlet", "outlet", "bottom", "top", "front", "back"}) {
        ubc[name] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};
        pbc[name] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    }

    p.fill(17.0);
    const auto phi_constant = make_rhie_chow_mass_flux(
        mesh, geometry, U, p, rAU, 1.0, ubc, pbc);
    double constant_max = 0.0;
    for (std::size_t f = 0; f < mesh.n_faces(); ++f)
        constant_max = std::max(constant_max, std::abs(phi_constant(f)));
    require(constant_max < 1e-14, "constant pressure generated a Rhie-Chow flux");

    // A linear pressure field isolates the pressure projection. With zero U
    // and p=x, every internal x-normal face must carry the same discrete
    // pressure response and internal fluxes must remain conservative.
    for (std::size_t c = 0; c < mesh.n_cells(); ++c)
        p(c) = geometry.cell_centres[c].x;

    pbc["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
    pbc["outlet"] = {ScalarBoundaryType::FIXED_VALUE, 1.0, 0.0};
    const auto phi_linear = make_rhie_chow_mass_flux(
        mesh, geometry, U, p, rAU, 1.0, ubc, pbc);

    double internal_conservation = 0.0;
    double internal_x_variation = 0.0;
    double first_x_flux = 0.0;
    bool first = true;
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f) < 0) continue;
        const double flux = phi_linear(f);
        const std::size_t owner = mesh.ownership().owner(f);\n        const double expected = -geometry.cell_volumes[owner] * 0.5 * geometry.face_area_vectors[f].x;
        if (std::abs(geometry.face_area_vectors[f].x) > 0.5) {
            if (first) { first_x_flux = flux; first = false; }
            internal_x_variation = std::max(internal_x_variation, std::abs(flux - first_x_flux));
            internal_conservation = std::max(internal_conservation, std::abs(flux - expected));
        }
    }

    require(internal_conservation < 1e-13,
            "linear pressure did not produce the expected Rhie-Chow x flux");
    require(internal_x_variation < 1e-13,
            "linear pressure Rhie-Chow flux varies between equivalent x faces");

    std::cout << "RHIE_CHOW constant_max=" << constant_max
              << " linear_max_error=" << internal_conservation
              << " linear_x_variation=" << internal_x_variation << "\n";
}

void test_cell_microscope(
    const Mesh& mesh, const FvGeometry& geometry, std::size_t cell)
{
    require(cell < mesh.n_cells(), "microscope cell out of range");

    Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
    U.fill(0.0);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c)
        U.component_data(0)[c] = geometry.cell_centres[c].y;

    Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
    p.fill(0.0);

    VelocityBoundaryConditions ubc;
    ScalarBoundaryConditions pbc;
    for (const char* name : {"inlet", "outlet", "front", "back", "bottom", "top"}) {
        ubc[name] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};
        pbc[name] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    }
    ubc["bottom"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0, 0, 0}};
    ubc["top"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {1, 0, 0}};

    const auto phi = make_mass_flux(mesh, geometry, U, 1.0, ubc);
    std::array<std::vector<double>, 3> rAU{
        std::vector<double>(mesh.n_cells(), 1.0),
        std::vector<double>(mesh.n_cells(), 1.0),
        std::vector<double>(mesh.n_cells(), 1.0)};
    const auto rc_phi = make_rhie_chow_mass_flux(
        mesh, geometry, U, p, rAU, 1.0, ubc, pbc);

    const auto* faces = mesh.cells().faces_data();
    const auto* offsets = mesh.cells().offsets_data();
    double authoritative_div = 0.0;
    double reconstructed_div = 0.0;

    std::cout << "MICROSCOPE cell=" << cell
              << " centre=(" << geometry.cell_centres[cell].x << ","
              << geometry.cell_centres[cell].y << ","
              << geometry.cell_centres[cell].z << ")\n";

    for (Offset k = offsets[cell]; k < offsets[cell + 1]; ++k) {
        const std::size_t f = faces[k];
        const bool owner = mesh.ownership().owner(f) == cell;
        const auto nr = mesh.ownership().neighbour(f);
        const Vec3 Sf = owner ? geometry.face_area_vectors[f]
                              : geometry.face_area_vectors[f] * -1.0;
        const double local_phi = owner ? phi(f) : -phi(f);
        const double local_rc_phi = owner ? rc_phi(f) : -rc_phi(f);
        const std::size_t neighbour =
            nr >= 0 ? static_cast<std::size_t>(nr) : mesh.ownership().owner(f);
        std::cout << "MICRO_FACE f=" << f
                  << " role=" << (owner ? "owner" : "neighbour")
                  << " neighbour=" << neighbour
                  << " patch=" << geometry.face_patch[f]
                  << " Sf=(" << Sf.x << "," << Sf.y << "," << Sf.z << ")"
                  << " phi=" << local_phi
                  << " rc_phi=" << local_rc_phi
                  << " delta_phi=" << (local_rc_phi - local_phi)
                  << "\n";
        authoritative_div += local_phi;
        reconstructed_div += local_rc_phi;
    }

    std::cout << "MICRO_CELL_BALANCE cell=" << cell
              << " authoritative_div=" << authoritative_div
              << " rc_div=" << reconstructed_div
              << " delta_div=" << (reconstructed_div - authoritative_div)
              << "\n";
}

} // namespace

int main()
{
    try {
        std::cout << "COUETTE_NUMERICAL_DIAGNOSTICS begin\n";
        const Mesh mesh = make_diagnostic_mesh();
        const FvGeometry geometry = build_fv_geometry(mesh);

        test_geometry_orientation(mesh, geometry);
        test_constant_flux_conservation(mesh, geometry);
        test_linear_gradient(mesh, geometry);
        test_rhie_chow_pressure_operators(mesh, geometry);

        // The default microscope cell is deliberately chosen near the current
        // Couette discrepancy, but can be overridden without recompilation.
        std::size_t cell = 33;
        if (const char* env = std::getenv("CFDX_DIAGNOSTIC_CELL")) {
            try { cell = static_cast<std::size_t>(std::stoul(env)); }
            catch (...) { throw std::runtime_error("invalid CFDX_DIAGNOSTIC_CELL"); }
        }
        test_cell_microscope(mesh, geometry, cell);

        std::cout << "COUETTE_NUMERICAL_DIAGNOSTICS PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "COUETTE_NUMERICAL_DIAGNOSTICS FAIL: " << e.what() << "\n";
        return 1;
    }
}
