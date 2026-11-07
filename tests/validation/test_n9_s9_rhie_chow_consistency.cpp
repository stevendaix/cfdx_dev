#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/physics/coupled_assembly_certificate.h"
#include "cfdx/physics/coupled_flux_equivalence_certificate.h"
#include "cfdx/physics/coupled_frozen_fv_certificate.h"
#include "cfdx/core/linalg/coupled_block_extraction.h"
#include "cfdx/core/linalg/small_dense_schur_oracle.h"
#include "cfdx/core/linalg/schur_inverse_quality_certificate.h"
#include "common/test_harness.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_two_cell_channel()
{
    Mesh m;
    m.points().resize(12);
    const double xyz[12][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i, xyz[i][0], xyz[i][1], xyz[i][2]);

    m.faces().push_face({0,3,7,4});       // x = 0
    m.faces().push_face({1,2,6,5});       // internal
    m.faces().push_face({8,9,11,10});     // x = 2
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({1,8,10,5});
    m.faces().push_face({2,6,11,9});
    m.faces().push_face({1,2,9,8});
    m.faces().push_face({5,10,11,6});

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 11; ++f) {
        const std::size_t owner = (f == 2 || f >= 7) ? 1u : 0u;
        m.ownership().set_owner(f, owner);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);

    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});

    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,2,3,4,5,6,7,8,9,10};
    m.boundary().add_patch(wall);
    return m;
}

static std::array<std::vector<double>, 3> unit_rAU(std::size_t n)
{
    return {
        std::vector<double>(n, 1.0),
        std::vector<double>(n, 1.0),
        std::vector<double>(n, 1.0)};
}

// Stretched two-cell channel: cell 0 spans x in [0, 0.5], cell 1 spans
// x in [0.5, 2] (stretch ratio 3). Same topology and face winding as
// make_two_cell_channel(); only the x-coordinates of the internal and
// right-boundary points change.
static Mesh make_two_cell_channel_stretched()
{
    Mesh m;
    m.points().resize(12);
    const double xyz[12][3] = {
        {0,0,0},{0.5,0,0},{0.5,1,0},{0,1,0},
        {0,0,1},{0.5,0,1},{0.5,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i, xyz[i][0], xyz[i][1], xyz[i][2]);

    m.faces().push_face({0,3,7,4});       // x = 0
    m.faces().push_face({1,2,6,5});       // internal, x = 0.5
    m.faces().push_face({8,9,11,10});     // x = 2
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({1,8,10,5});
    m.faces().push_face({2,6,11,9});
    m.faces().push_face({1,2,9,8});
    m.faces().push_face({5,10,11,6});

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 11; ++f) {
        const std::size_t owner = (f == 2 || f >= 7) ? 1u : 0u;
        m.ownership().set_owner(f, owner);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);

    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});

    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,2,3,4,5,6,7,8,9,10};
    m.boundary().add_patch(wall);
    return m;
}

static std::array<std::vector<double>, 3> scale_rAU(
    const std::array<std::vector<double>, 3>& rAU, double s)
{
    std::array<std::vector<double>, 3> out{};
    for (std::size_t k = 0; k < 3; ++k) {
        out[k].resize(rAU[k].size());
        for (std::size_t i = 0; i < rAU[k].size(); ++i)
            out[k][i] = s * rAU[k][i];
    }
    return out;
}

// The Rhie-Chow pressure correction of an internal face, extracted as
// phi_interpolated - phi_RC. Both fluxes come from the production paths
// (make_mass_flux / make_rhie_chow_mass_flux), so no formula is duplicated
// here: only the difference is taken.
static double rc_correction(std::size_t face,
    const Mesh& mesh,
    const FvGeometry& geometry,
    const Field<double, Location::CELL>& U,
    const Field<double, Location::CELL>& p,
    const std::array<std::vector<double>, 3>& rAU,
    double rho,
    const VelocityBoundaryConditions& ubc,
    const ScalarBoundaryConditions& pbc)
{
    const auto phi_interp = make_mass_flux(mesh, geometry, U, rho, ubc);
    const auto phi_rc = make_rhie_chow_mass_flux(
        mesh, geometry, U, p, rAU, rho, ubc, pbc);
    return phi_interp(face) - phi_rc(face);
}

static Mesh make_two_cell_channel_oblique()
{
    Mesh m;
    m.points().resize(12);
    constexpr double shear = 0.35;
    const double xyz[12][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i, xyz[i][0] + shear * xyz[i][1],
                       xyz[i][1], xyz[i][2]);

    m.faces().push_face({0,3,7,4});
    m.faces().push_face({1,2,6,5}); // Internal face is non-orthogonal to C_O C_N.
    m.faces().push_face({8,9,11,10});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({1,8,10,5});
    m.faces().push_face({2,6,11,9});
    m.faces().push_face({1,2,9,8});
    m.faces().push_face({5,10,11,6});

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 11; ++f) {
        const std::size_t owner = (f == 2 || f >= 7) ? 1u : 0u;
        m.ownership().set_owner(f, owner);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);
    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});

    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,2,3,4,5,6,7,8,9,10};
    m.boundary().add_patch(wall);
    return m;
}

// Four-cell channel so an alternating pressure mode is distinguishable
// from a genuinely affine pressure field. Interior cells have a centred
// zero gradient for p = (+1,-1,+1,-1), while adjacent-cell jumps remain nonzero.
static Mesh make_four_cell_channel(
    const std::array<double, 5>& x_positions = {0.0, 1.0, 2.0, 3.0, 4.0},
    bool split_boundary_patches = false)
{
    Mesh m;
    constexpr std::size_t nx = 4;
    m.points().resize((nx + 1) * 4);
    for (std::size_t i = 0; i <= nx; ++i) {
        const double x = x_positions[i];
        const double xyz[4][3] = {
            {x,0,0}, {x,1,0}, {x,1,1}, {x,0,1}};
        for (std::size_t j = 0; j < 4; ++j)
            m.points().set(4*i+j, xyz[j][0], xyz[j][1], xyz[j][2]);
    }

    std::vector<std::size_t> x_faces(nx + 1);
    std::vector<std::vector<std::size_t>> side_faces(nx);
    for (std::size_t i = 0; i <= nx; ++i) {
        const std::size_t b = 4*i;
        m.faces().push_face({b, b+1, b+2, b+3});
        x_faces[i] = i;
    }
    for (std::size_t i = 0; i < nx; ++i) {
        const std::size_t l = 4*i, r = 4*(i+1);
        const std::size_t first = nx + 1 + 4 * i;
        m.faces().push_face({l, r, r+3, l+3});       // y = 0
        m.faces().push_face({l+1, l+2, r+2, r+1});   // y = 1
        m.faces().push_face({l, l+1, r+1, r});       // z = 0
        m.faces().push_face({l+3, r+3, r+2, l+2});   // z = 1
        side_faces[i] = {first, first+1, first+2, first+3};
    }

    constexpr std::size_t n_faces = (nx + 1) + 4 * nx;
    m.ownership().resize(n_faces);
    for (std::size_t f = 0; f < n_faces; ++f)
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    for (std::size_t i = 0; i <= nx; ++i) {
        if (i == 0) m.ownership().set_owner(x_faces[i], 0);
        else m.ownership().set_owner(x_faces[i], i-1);
        if (i > 0 && i < nx)
            m.ownership().set_neighbour(x_faces[i], static_cast<std::int64_t>(i));
    }

    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    Patch inlet;
    inlet.name = "inlet";
    inlet.type = PatchType::WALL;
    Patch outlet;
    outlet.name = "outlet";
    outlet.type = PatchType::WALL;
    if (split_boundary_patches) {
        inlet.face_ids = {x_faces.front()};
        outlet.face_ids = {x_faces.back()};
    } else {
        wall.face_ids = {x_faces.front(), x_faces.back()};
    }
    for (std::size_t i = 0; i < nx; ++i) {
        const auto l = x_faces[i], r = x_faces[i+1];
        m.cells().push_cell({l, r, side_faces[i][0], side_faces[i][1],
                             side_faces[i][2], side_faces[i][3]});
        for (const auto f : side_faces[i]) {
            m.ownership().set_owner(f, i);
            wall.face_ids.push_back(f);
        }
    }
    m.boundary().add_patch(wall);
    if (split_boundary_patches) {
        m.boundary().add_patch(inlet);
        m.boundary().add_patch(outlet);
    }
    return m;
}

int main()
{
    run_case("n9_s9_rhie_chow_linear_pressure_consistency", [] {
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.set(0, 2.0, 0.0, 0.0);
        U.set(1, 2.0, 0.0, 0.0);
        p(0) = 1.0;
        p(1) = 3.0;

        // The internal pressure jump is 2 Pa across a unit centre distance.
        // The zero-gradient boundary condition intentionally defines the
        // discrete cell gradient used by the compensation term.

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const auto phi = make_rhie_chow_flux_from_physical_velocity(
            mesh, geometry, U, p, rAU, 1.0, ubc, pbc);
        const auto phi_linear = make_mass_flux(mesh, geometry, U, 1.0, ubc);
        const auto grad_p_discrete = gauss_gradient_with_boundary(
            p, mesh, geometry, pbc);

        // COUPLED supplies physical U, unlike the segregated H/A predictor.
        // Its flux therefore includes interp(dAU*grad(p)) before subtracting
        // the implicit pressure response. On this orthogonal internal face,
        // the response is dAU_f * (p_N-p_O) / d * Sf.x; this direct expression
        // is independent of the production Rhie-Chow flux routine.
        const double grad_face_x = 0.5 * (
            grad_p_discrete.component_data(0)[0] +
            grad_p_discrete.component_data(0)[1]);
        const double dAU_face_x = 0.5 * (
            geometry.cell_volumes[0] * rAU[0][0] +
            geometry.cell_volumes[1] * rAU[0][1]);
        const double pressure_jump_response =
            dAU_face_x * (p(1) - p(0));
        const double expected = phi_linear(1) +
            dAU_face_x * grad_face_x - pressure_jump_response;
        EXPECT_NEAR(phi(1), expected, 1e-12);
    });

    run_case("n9_s9_rhie_chow_pressure_offset_invariance", [] {
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p0(2, "p", "Pa", 1);
        Field<double, Location::CELL> p1(2, "p", "Pa", 1);
        U.fill(0.0);
        p0(0) = -2.0;
        p0(1) = 4.0;
        p1(0) = 1002.0;
        p1(1) = 1008.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        // Exercise the COUPLED physical-velocity contract, not the default
        // segregated H/A contract, for this pressure-gauge invariance check.
        const auto phi0 = make_rhie_chow_flux_from_physical_velocity(
            mesh, geometry, U, p0, rAU, 1.0, ubc, pbc);
        const auto phi1 = make_rhie_chow_flux_from_physical_velocity(
            mesh, geometry, U, p1, rAU, 1.0, ubc, pbc);

        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
            EXPECT_NEAR(phi0(f), phi1(f), 1e-12);
    });

    run_case("n9_s9_rhie_chow_affine_pressure_nonuniform_mesh", [] {
        const Mesh mesh = make_four_cell_channel(
            {0.0, 0.5, 2.0, 4.0, 7.0}, true);
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        U.fill(0.0);
        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
            U.set(cell, 2.0, 0.0, 0.0);
            p(cell) = 1.0 + 2.0 * geometry.cell_centres[cell].x;
        }

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ubc["inlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {2.0, 0.0, 0.0}};
        ubc["outlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {2.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        pbc["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 1.0, 0.0};
        pbc["outlet"] = {ScalarBoundaryType::FIXED_VALUE, 15.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const auto phi = make_rhie_chow_flux_from_physical_velocity(
            mesh, geometry, U, p, rAU, 1.0, ubc, pbc);
        const auto phi_hbya = make_rhie_chow_flux_from_hbya(
            mesh, geometry, U, p, rAU, 1.0, ubc, pbc);
        const auto phi_linear = make_mass_flux(mesh, geometry, U, 1.0, ubc);
        const auto grad_p = gauss_gradient_with_boundary(p, mesh, geometry, pbc);

        // Internal face IDs 1..3 are orthogonal but have unequal centre
        // distances and unequal cell volumes. The oracle is evaluated directly
        // from geometry, the pressure jump, and the independently inspected
        // cell gradients, rather than calling the production correction helper.
        for (std::size_t f = 1; f <= 3; ++f) {
            const auto owner = mesh.ownership().owner(f);
            const auto neighbour = static_cast<std::size_t>(mesh.ownership().neighbour(f));
            const double dAU_face_x = 0.5 * (
                geometry.cell_volumes[owner] * rAU[0][owner] +
                geometry.cell_volumes[neighbour] * rAU[0][neighbour]);
            const double dx = geometry.cell_centres[neighbour].x -
                              geometry.cell_centres[owner].x;
            const double grad_face_x = 0.5 * (
                grad_p.component_data(0)[owner] +
                grad_p.component_data(0)[neighbour]);
            const double sf_x = geometry.face_area_vectors[f].x;
            const double pressure_jump_response =
                dAU_face_x * (p(neighbour) - p(owner)) / dx * sf_x;
            const double compensation = dAU_face_x * grad_face_x * sf_x;
            const double expected_hbya = phi_linear(f) - pressure_jump_response;
            const double expected = expected_hbya + compensation;
            EXPECT_NEAR(grad_face_x, 2.0, 1e-12);
            EXPECT_NEAR(phi_hbya(f), expected_hbya, 1e-12);
            EXPECT_NEAR(phi(f), expected, 1e-12);
            EXPECT_NEAR(phi(f) - phi_hbya(f), compensation, 1e-12);
        }

        // Independent cellwise divergence of the authoritative reconstructed
        // flux. Constant U has zero divergence on this closed orthogonal mesh.
        std::vector<double> cell_divergence(mesh.n_cells(), 0.0);
        for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
            const auto owner = mesh.ownership().owner(f);
            cell_divergence[owner] += phi(f);
            const auto neighbour = mesh.ownership().neighbour(f);
            if (neighbour >= 0)
                cell_divergence[static_cast<std::size_t>(neighbour)] -= phi(f);
        }
        for (const double divergence : cell_divergence)
            EXPECT_NEAR(divergence, 0.0, 1e-12);
    });

    run_case("n9_s9_coupled_block_continuity_matches_flux_divergence", [] {
        const Mesh mesh = make_four_cell_channel(
            {0.0, 0.5, 2.0, 4.0, 7.0}, true);
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
            // Manufactured solution: U=(2,0,0), p=1+2x, body force=(2,0,0).
            // Start from the manufactured exact state: a one-step
            // frozen-deferred system cannot be required to recover a distinct
            // pressure field to exact precision.
            U.set(cell, 2.0, 0.0, 0.0);
            p(cell) = 1.0 + 2.0 * geometry.cell_centres[cell].x;
        }

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {2.0, 0.0, 0.0}};
        ubc["inlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {2.0, 0.0, 0.0}};
        ubc["outlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {2.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        pbc["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 1.0, 0.0};
        pbc["outlet"] = {ScalarBoundaryType::FIXED_VALUE, 15.0, 0.0};

        SparseMatrix captured_matrix;
        Vector captured_rhs;
        std::vector<double> captured_schur_diagonal;
        std::vector<double> captured_schur_diagonal_raw;
        bool captured_has_fixed_pressure = false;
        bool captured_structure_certificate = false;
        bool frozen_fv_certificate = false;
        bool mutated_frozen_certificate_rejected = false;
        bool captured = false;
        IncompressibleSolverControls controls;
        controls.algorithm = PressureVelocityAlgorithm::COUPLED;
        controls.pressure_gauge_policy = PressureGaugePolicy::NONE;
        controls.body_force = {2.0, 0.0, 0.0};
        controls.kinematic_viscosity = 1e-2;
        controls.coupling.alpha_u = 1.0;
        controls.coupling.alpha_p = 1.0;
        controls.coupling.n_outer_correctors = 1;
        controls.coupling.schur_model = CoupledSchurModel::BlockLocal;
        controls.coupling.coupled_linear_tolerance = 1e-12;
        controls.convergence.max_iterations = 1;
        controls.coupled_assembly_callback = [&](const CoupledAssemblySnapshot& snapshot) {
            if (snapshot.matrix == nullptr || snapshot.rhs == nullptr)
                throw std::runtime_error("coupled assembly evidence was incomplete");
            captured_matrix = *snapshot.matrix;
            captured_rhs = *snapshot.rhs;
            if (snapshot.pressure_schur_diagonal == nullptr ||
                snapshot.pressure_schur_diagonal_raw == nullptr)
                throw std::runtime_error("coupled Schur diagonal evidence was incomplete");
            captured_schur_diagonal = *snapshot.pressure_schur_diagonal;
            captured_schur_diagonal_raw = *snapshot.pressure_schur_diagonal_raw;
            captured_has_fixed_pressure = snapshot.has_fixed_pressure;
            if (snapshot.initial_guess == nullptr)
                throw std::runtime_error("coupled initial guess missing from evidence");
            const auto frozen = certify_coupled_frozen_fv_operator(
                snapshot, mesh, geometry, ubc, pbc, controls.density,
                controls.pressure_reference_value, *snapshot.initial_guess);
            frozen_fv_certificate = frozen.passed;
            Vector perturbed = *snapshot.initial_guess;
            for (std::size_t i = 0; i < perturbed.size(); ++i)
                perturbed(i) += 0.013 * static_cast<double>(i % 7 + 1);
            const auto perturbed_report = certify_coupled_frozen_fv_operator(
                snapshot, mesh, geometry, ubc, pbc, controls.density,
                controls.pressure_reference_value, perturbed);
            frozen_fv_certificate = frozen_fv_certificate && perturbed_report.passed;
            std::cerr << "N9_S9_FROZEN_PICARD_CERTIFICATE"
                      << " passed=" << frozen.passed
                      << " worst_row=" << frozen.worst_row
                      << " momentum_linf=" << frozen.momentum_linf
                      << " continuity_linf=" << frozen.continuity_linf << '\n';
            // A deliberately corrupted G entry MUST be detected by the
            // independent momentum/face-flux residual at the same state.
            SparseMatrix tampered = *snapshot.matrix;
            bool changed = false;
            for (std::size_t row = 0; row < 3*mesh.n_cells() && !changed; ++row)
                for (std::size_t k = tampered.row_offsets_data()[row];
                     k < tampered.row_offsets_data()[row+1]; ++k)
                    if (tampered.columns_data()[k] >= 3*mesh.n_cells()) {
                        tampered.values_data()[k] += 0.125;
                        changed = true;
                        break;
                    }
            CoupledAssemblySnapshot broken = snapshot;
            broken.matrix = &tampered;
            const auto negative = certify_coupled_frozen_fv_operator(
                broken, mesh, geometry, ubc, pbc, controls.density,
                controls.pressure_reference_value, *snapshot.initial_guess);
            mutated_frozen_certificate_rejected = changed && !negative.passed;
            CoupledAssemblyCertificateOptions audit_controls;
            audit_controls.pressure_schur_diagonal = snapshot.pressure_schur_diagonal;
            audit_controls.pressure_schur_diagonal_raw = snapshot.pressure_schur_diagonal_raw;
            // This case has fixed-pressure inlet/outlet: do NOT impose
            // a fake pinned gauge equation on genuine continuity rows.
            if (!snapshot.has_fixed_pressure)
                audit_controls.pressure_reference_cell = snapshot.reference_cell;
            const auto audit = certify_coupled_assembly_structure(
                *snapshot.matrix, *snapshot.rhs, *snapshot.initial_guess,
                3 * mesh.n_cells(), mesh.n_cells(), audit_controls);
            captured_structure_certificate = audit.passed &&
                audit.nnz_momentum > 0 && audit.nnz_gradient > 0 &&
                audit.nnz_divergence > 0 && audit.nnz_pressure > 0;
            std::cerr << "N9_S9_COUPLED_ASSEMBLY_CERTIFICATE"
                      << " structural_pass=" << audit.passed
                      << " gradient_nnz=" << audit.nnz_gradient
                      << " divergence_nnz=" << audit.nnz_divergence
                      << " linear_continuity_linf="
                      << audit.continuity_linear_residual_linf << '\n';
            captured = true;
        };

        const auto result = solve_steady_incompressible(mesh, U, p, ubc, pbc, controls);
        EXPECT_TRUE(captured);
        EXPECT_TRUE(captured_structure_certificate);
        EXPECT_TRUE(frozen_fv_certificate);
        EXPECT_TRUE(mutated_frozen_certificate_rejected);
        // Diagnostic for the manufactured fixed-pressure case: distinguish a
        // bad assembled operator/RHS from a Krylov solution error without
        // changing any acceptance threshold.
        if (captured) {
            Vector exact(4 * mesh.n_cells(), 0.0);
            for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
                exact(cell) = 2.0;
                exact(3 * mesh.n_cells() + cell) =
                    1.0 + 2.0 * geometry.cell_centres[cell].x;
            }
            double exact_residual_linf = 0.0;
            std::size_t exact_residual_row = 0;
            for (std::size_t row = 0; row < captured_matrix.n_rows(); ++row) {
                double ax = 0.0;
                for (std::size_t k = captured_matrix.row_offsets_data()[row];
                     k < captured_matrix.row_offsets_data()[row + 1]; ++k)
                    ax += captured_matrix.values_data()[k] *
                          exact(captured_matrix.columns_data()[k]);
                const double residual = std::abs(ax - captured_rhs(row));
                if (residual > exact_residual_linf) {
                    exact_residual_linf = residual;
                    exact_residual_row = row;
                }
            }
            std::cerr << "N9_S9_MANUFACTURED_FIXED_PRESSURE"
                      << " exact_matrix_residual_linf=" << exact_residual_linf
                      << " exact_matrix_residual_row=" << exact_residual_row;
            for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
                std::cerr << " cell" << cell
                          << "_Ux=" << U.component_data(0)[cell]
                          << "_p=" << p(cell)
                          << "_p_exact="
                          << (1.0 + 2.0 * geometry.cell_centres[cell].x);
            std::cerr << "\\n";
        }
        // V4: build the exact non-SPD-aware Schur oracle from the *captured
        // post-gauge production CSR* rather than an unrelated toy matrix.
        // The matrix includes a finite convective flux and Dirichlet pressure.
        const std::size_t n_capture = mesh.n_cells();
        const auto M = extract_coupled_block(captured_matrix,0,0,
            3*n_capture,3*n_capture);
        const auto G = extract_coupled_block(captured_matrix,0,1,
            3*n_capture,n_capture);
        const auto D = extract_coupled_block(captured_matrix,1,0,
            n_capture,3*n_capture);
        const auto C = extract_coupled_block(captured_matrix,1,1,
            n_capture,n_capture);
        const BlockOperator production_blocks(M,G,D,C);
        EXPECT_TRUE(production_blocks.is_valid());
        SmallDenseExactSchurInverse schur_oracle(64);
        EXPECT_TRUE(schur_oracle.setup(production_blocks));
        const auto probes = schur_inverse_quality_basis(n_capture);
        const auto schur_cert = certify_schur_inverse_action(
            production_blocks, schur_oracle,
            [&](const Vector& rhs, Vector& answer) {
                return schur_oracle.solve_velocity(rhs,answer);
            }, probes);
        std::cerr << "N9_S9_PRODUCTION_EXACT_SCHUR"
                  << " complete=" << schur_cert.complete
                  << " passed=" << schur_cert.passed
                  << " max_defect=" << schur_cert.worst_schur_defect
                  << " velocity_defect=" << schur_cert.worst_velocity_solve_defect
                  << '\n';
        EXPECT_TRUE(schur_cert.complete && schur_cert.passed);
        EXPECT_TRUE(schur_cert.probes_checked == n_capture + 1);
        EXPECT_TRUE(captured_has_fixed_pressure);
        EXPECT_TRUE(captured_schur_diagonal.size() == mesh.n_cells());
        EXPECT_TRUE(captured_schur_diagonal_raw.size() == mesh.n_cells());
        // A fixed-pressure boundary means no continuity row is replaced by
        // the pressure gauge. The Schur preconditioner must therefore retain
        // the actual assembled diagonal at every cell, including the nominal
        // reference cell; only the pure-Neumann gauge path may pin it to 1.
        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
            EXPECT_NEAR(captured_schur_diagonal[cell],
                        captured_schur_diagonal_raw[cell], 1e-14);
        EXPECT_TRUE(result.authoritative_mass_flux.size() == mesh.n_faces());
        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
            EXPECT_NEAR(U.component_data(0)[cell], 2.0, 1e-8);
            EXPECT_NEAR(U.component_data(1)[cell], 0.0, 1e-8);
            EXPECT_NEAR(U.component_data(2)[cell], 0.0, 1e-8);
            EXPECT_NEAR(p(cell), 1.0 + 2.0 * geometry.cell_centres[cell].x, 1e-8);
        }

        const std::size_t nc = mesh.n_cells();
        Vector x(4 * nc, 0.0);
        for (std::size_t cell = 0; cell < nc; ++cell) {
            x(cell) = U.component_data(0)[cell];
            x(nc + cell) = U.component_data(1)[cell];
            x(2 * nc + cell) = U.component_data(2)[cell];
            x(3 * nc + cell) = p(cell);
        }
        std::vector<double> cell_divergence(nc, 0.0);
        for (std::size_t face = 0; face < mesh.n_faces(); ++face) {
            const auto owner = mesh.ownership().owner(face);
            cell_divergence[owner] += result.authoritative_mass_flux(face);
            const auto neighbour = mesh.ownership().neighbour(face);
            if (neighbour >= 0)
                cell_divergence[static_cast<std::size_t>(neighbour)] -=
                    result.authoritative_mass_flux(face);
        }

        // Independent production flux-versus-matrix certificate: use the
        // AUTHORITATIVE FV flux and exact Krylov CSR capture, not a matrix-
        // derived proxy for the flux. The existing 1e-8 gate is unchanged.
        std::vector<std::size_t> face_owner(mesh.n_faces());
        std::vector<std::int64_t> face_neighbour(mesh.n_faces());
        std::vector<double> face_flux(mesh.n_faces());
        for (std::size_t face = 0; face < mesh.n_faces(); ++face) {
            face_owner[face] = mesh.ownership().owner(face);
            face_neighbour[face] = mesh.ownership().neighbour(face);
            face_flux[face] = result.authoritative_mass_flux(face);
        }
        CoupledFluxEquivalenceControls face_controls;
        face_controls.absolute_tolerance = 1e-8;
        face_controls.relative_tolerance = 1e-8;
        face_controls.continuity_linf_gate = 1e-8;
        const auto face_certificate = certify_coupled_face_flux_equivalence(
            captured_matrix, captured_rhs, x, nc, face_owner,
            face_neighbour, face_flux, face_controls);
        std::cerr << "N9_S9_COUPLED_FACE_FLUX_CERTIFICATE"
                  << " complete=" << face_certificate.complete
                  << " pass=" << face_certificate.passed
                  << " matrix_flux_linf=" << face_certificate.max_identity_error
                  << " continuity_linf=" <<
                     face_certificate.max_conservative_flux_divergence << '\n';
        EXPECT_TRUE(face_certificate.complete);
        EXPECT_TRUE(face_certificate.passed);

        // The post-solve continuity row is the exact algebraic residual of
        // the captured production block system. With the production sign
        // convention, A*x-b equals the outward divergence of the authoritative
        // face flux, independently reconstructed cell by cell.
        for (std::size_t cell = 0; cell < nc; ++cell) {
            const std::size_t row = 3 * nc + cell;
            double ax = 0.0;
            const auto rb = captured_matrix.row_offsets_data()[row];
            const auto re = captured_matrix.row_offsets_data()[row + 1];
            for (std::size_t k = rb; k < re; ++k)
                ax += captured_matrix.values_data()[k] *
                      x(captured_matrix.columns_data()[k]);
            const double block_residual = ax - captured_rhs(row);
            EXPECT_NEAR(block_residual, cell_divergence[cell], 1e-8);
        }
    });

    run_case("n9_s9_coupled_picard_frozen_nonzero_pressure_neumann", [] {
        const Mesh mesh = make_four_cell_channel(
            {0.0, 0.5, 2.0, 4.0, 7.0}, true);
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            U.set(c, 1.0, 0.0, 0.0);
            p(c) = 2.0 + 0.5 * geometry.cell_centres[c].x;
        }
        VelocityBoundaryConditions ubc;
        for (const char* name : {"wall", "inlet", "outlet"})
            ubc[name] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {1.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        pbc["inlet"] = {ScalarBoundaryType::FIXED_GRADIENT, 0.0, 0.5};
        pbc["outlet"] = {ScalarBoundaryType::FIXED_VALUE, 5.5, 0.0};
        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::COUPLED;
        c.pressure_gauge_policy = PressureGaugePolicy::NONE;
        c.kinematic_viscosity = 1e-2;
        c.coupling.schur_model = CoupledSchurModel::BlockLocal;
        c.coupling.coupled_linear_tolerance = 1e-10;
        c.convergence.max_iterations = 1;
        bool captured = false;
        c.coupled_assembly_callback = [&](const CoupledAssemblySnapshot& s) {
            captured = true;
            EXPECT_TRUE(s.initial_guess != nullptr && s.pressure_old != nullptr);
            if (!s.initial_guess || !s.pressure_old) return;
            EXPECT_TRUE(s.has_fixed_pressure);
            // Dirichlet pressure must NOT clobber this reference-cell
            // pressure in the Krylov initial guess.
            EXPECT_NEAR((*s.initial_guess)(3*mesh.n_cells() + s.reference_cell),
                        (*s.pressure_old)(s.reference_cell), 0.0);
            const auto report = certify_coupled_frozen_fv_operator(
                s, mesh, geometry, ubc, pbc, c.density,
                c.pressure_reference_value, *s.initial_guess);
            std::cerr << "N9_S9_NEUMANN_FROZEN"
                      << " passed=" << report.passed
                      << " worst_row=" << report.worst_row
                      << " error=" << report.max_absolute_error
                      << " ratio=" << report.max_tolerance_ratio << '\n';
            EXPECT_TRUE(report.passed);
            // Reproduce the old bug: delete the known boundary-gradient
            // source from the assembled momentum RHS.  The FV certificate
            // must fail even though the matrix G itself is unchanged.
            Vector broken_rhs = *s.rhs;
            const auto& inlet = mesh.boundary().patch(mesh.boundary().find("inlet"));
            for (std::size_t f : inlet.face_ids) {
                const std::size_t owner = mesh.ownership().owner(f);
                const auto Sf = geometry.face_area_vectors[f];
                const auto delta = geometry.face_centres[f] - geometry.cell_centres[owner];
                const double boundary_offset = pbc.at("inlet").gradient *
                    delta.dot(Sf * (1.0 / Sf.mag()));
                broken_rhs(owner) += Sf.x * boundary_offset;
                broken_rhs(mesh.n_cells()+owner) += Sf.y * boundary_offset;
                broken_rhs(2*mesh.n_cells()+owner) += Sf.z * boundary_offset;
            }
            CoupledAssemblySnapshot missing_gradient_source = s;
            missing_gradient_source.rhs = &broken_rhs;
            const auto negative = certify_coupled_frozen_fv_operator(
                missing_gradient_source, mesh, geometry, ubc, pbc,
                c.density, c.pressure_reference_value, *s.initial_guess);
            EXPECT_TRUE(!negative.passed);
            EXPECT_TRUE(negative.momentum_linf > 1e-3);
        };
        (void)solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
        EXPECT_TRUE(captured);
    });

    run_case("n9_s9_coupled_picard_xy_empty_spanwise_faces", [] {
        Mesh mesh = make_four_cell_channel(
            {0.0, 0.5, 2.0, 4.0, 7.0}, true);
        // Convert only z=constant faces to EMPTY, as on NACA XY.
        std::vector<FaceIndex> wall, front, back;
        const std::size_t wall_patch = mesh.boundary().find("wall");
        for (std::size_t f : mesh.boundary().patch(wall_patch).face_ids) {
            double z = 0.0;
            for (std::size_t k=mesh.faces().face_offset(f);
                 k<mesh.faces().face_offset(f)+mesh.faces().face_size(f); ++k)
                z += mesh.points().z(mesh.faces().vertices_data()[k]);
            z /= mesh.faces().face_size(f);
            if (std::abs(z) < 1e-12) front.push_back(f);
            else if (std::abs(z-1.0) < 1e-12) back.push_back(f);
            else wall.push_back(f);
        }
        mesh.boundary().patch(wall_patch).face_ids = std::move(wall);
        mesh.boundary().add_patch(Patch{"front", PatchType::EMPTY, front});
        mesh.boundary().add_patch(Patch{"back", PatchType::EMPTY, back});
        EXPECT_TRUE(front.size() == mesh.n_cells());
        EXPECT_TRUE(back.size() == mesh.n_cells());
        const auto geometry = build_fv_geometry(mesh);
        Field<double,Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double,Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        for(std::size_t cell=0; cell<mesh.n_cells(); ++cell) {
            U.set(cell, 1.0, 0.0, 0.0);
            p(cell) = 1.0 + geometry.cell_centres[cell].x;
        }
        VelocityBoundaryConditions ubc;
        for(const char* name:{"inlet","outlet","wall"})
            ubc[name]={VelocityBoundaryCondition::Type::FIXED_VALUE,{1.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        pbc["inlet"]={ScalarBoundaryType::FIXED_VALUE,1.0,0.0};
        pbc["outlet"]={ScalarBoundaryType::FIXED_VALUE,8.0,0.0};
        IncompressibleSolverControls c;
        c.algorithm=PressureVelocityAlgorithm::COUPLED;
        c.pressure_gauge_policy=PressureGaugePolicy::NONE;
        c.kinematic_viscosity=1e-2;
        c.coupling.schur_model=CoupledSchurModel::BlockLocal;
        c.coupling.coupled_linear_tolerance=1e-10;
        c.convergence.max_iterations=1;
        bool checked=false;
        c.coupled_assembly_callback=[&](const CoupledAssemblySnapshot& snap) {
            if (!snap.initial_guess) return;
            checked=true;
            Vector x=*snap.initial_guess;
            for (std::size_t i=0; i<x.size(); ++i)
                x(i)+=0.009*(i%5+1);
            const auto check=certify_coupled_frozen_fv_operator(
                snap,mesh,geometry,ubc,pbc,c.density,
                c.pressure_reference_value,x);
            std::cerr << "N9_S9_XY_EMPTY_FROZEN passed=" << check.passed
                      << " error=" << check.max_absolute_error << '\n';
            EXPECT_TRUE(check.passed);
        };
        (void)solve_steady_incompressible(mesh,U,p,ubc,pbc,c);
        EXPECT_TRUE(checked);
    });

    run_case("n9_s9_coupled_picard_frozen_pure_neumann_gauge", [] {
        const Mesh mesh = make_four_cell_channel();
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        for (std::size_t i = 0; i < mesh.n_cells(); ++i) {
            U.set(i, 0.5, 0.0, 0.0);
            p(i) = 1.0 + 0.2 * geometry.cell_centres[i].x;
        }
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,
                       {0.5,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::COUPLED;
        c.coupling.schur_model = CoupledSchurModel::BlockLocal;
        c.kinematic_viscosity = 1e-2;
        c.pressure_reference_cell = 1;
        c.pressure_reference_value = 0.0;
        c.coupling.coupled_linear_tolerance = 1e-10;
        c.convergence.max_iterations = 1;
        bool checked = false;
        c.coupled_assembly_callback = [&](const CoupledAssemblySnapshot& snap) {
            if (!snap.initial_guess) return;
            checked = true;
            EXPECT_TRUE(!snap.has_fixed_pressure);
            Vector perturbed = *snap.initial_guess;
            for (std::size_t i = 0; i < perturbed.size(); ++i)
                perturbed(i) += 0.008 * (i % 11 + 1);
            const auto report = certify_coupled_frozen_fv_operator(
                snap, mesh, geometry, ubc, pbc, c.density,
                c.pressure_reference_value, perturbed);
            std::cerr << "N9_S9_GAUGE_FROZEN passed=" << report.passed
                      << " worst_row=" << report.worst_row
                      << " error=" << report.max_absolute_error
                      << " ratio=" << report.max_tolerance_ratio << '\n';
            EXPECT_TRUE(report.passed);
            // The gauge row is p_ref=value, not a continuity equation.
            EXPECT_NEAR((*snap.initial_guess)(3*mesh.n_cells()+snap.reference_cell),
                        0.0, 0.0);
        };
        (void)solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
        EXPECT_TRUE(checked);
    });

    run_case("n9_s9_coupled_schwarz_optin_manufactured_four_cell", [] {
        // Exercise the actual production COUPLED dispatch, not just a
        // hand-made matrix. Two overlapping RAS subdomains of two core
        // cells each: correctness smoke, not a speed or scaling claim.
        const Mesh mesh = make_four_cell_channel({0.0, 0.5, 2.0, 4.0, 7.0}, true);
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            U.set(c, 2.0, 0.0, 0.0);
            p(c) = 1.0 + 2.0 * geometry.cell_centres[c].x;
        }
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{2.0,0.0,0.0}};
        ubc["inlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{2.0,0.0,0.0}};
        ubc["outlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{2.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        pbc["inlet"] = {ScalarBoundaryType::FIXED_VALUE,1.0,0.0};
        pbc["outlet"] = {ScalarBoundaryType::FIXED_VALUE,15.0,0.0};
        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::COUPLED;
        c.pressure_gauge_policy = PressureGaugePolicy::NONE;
        c.body_force = {2.0,0.0,0.0};
        c.kinematic_viscosity = 1e-2;
        c.coupling.alpha_u = 1.0;
        c.coupling.alpha_p = 1.0;
        c.coupling.schur_model = CoupledSchurModel::BlockLocal;
        c.coupling.coupled_schwarz_experimental = true;
        c.coupling.coupled_schwarz_cells_per_subdomain = 2;
        c.coupling.coupled_schwarz_overlap_hops = 1;
        c.coupling.coupled_schwarz_max_local_unknowns = 16;
        c.coupling.coupled_linear_tolerance = 1e-12;
        c.convergence.max_iterations = 1;
        std::size_t snapshot_count = 0;
        c.coupled_assembly_callback = [&](const CoupledAssemblySnapshot& snap) {
            EXPECT_TRUE(snap.matrix != nullptr && snap.rhs != nullptr);
            ++snapshot_count;
        };
        const auto result = solve_steady_incompressible(mesh,U,p,ubc,pbc,c);
        EXPECT_TRUE(snapshot_count > 0);
        EXPECT_TRUE(result.authoritative_mass_flux.size() == mesh.n_faces());
        for (std::size_t cell=0; cell<mesh.n_cells(); ++cell) {
            EXPECT_NEAR(U.component_data(0)[cell],2.0,1e-8);
            EXPECT_NEAR(U.component_data(1)[cell],0.0,1e-8);
            EXPECT_NEAR(U.component_data(2)[cell],0.0,1e-8);
            EXPECT_NEAR(p(cell),1.0 + 2.0*geometry.cell_centres[cell].x,1e-8);
        }
        double divergence_linf = 0.0;
        std::vector<double> divergence(mesh.n_cells(),0.0);
        for (std::size_t f=0; f<mesh.n_faces(); ++f) {
            divergence[mesh.ownership().owner(f)] += result.authoritative_mass_flux(f);
            const auto nb = mesh.ownership().neighbour(f);
            if (nb>=0) divergence[static_cast<std::size_t>(nb)] -= result.authoritative_mass_flux(f);
        }
        for (double div:divergence) divergence_linf = std::max(divergence_linf,std::abs(div));
        std::cout << "COUPLED_SCHWARZ_PRODUCTION continuity_linf=" << divergence_linf << '\n';
        EXPECT_TRUE(divergence_linf <= 1e-8);
    });

    run_case("n9_s9_rhie_chow_fixed_pressure_boundary_gauge_invariance", [] {
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p0(2, "p0", "Pa", 1);
        Field<double, Location::CELL> p1(2, "p1", "Pa", 1);
        U.fill(0.0);
        p0(0) = -2.0; p0(1) = 4.0;
        constexpr double gauge_shift = 1000.0;
        p1(0) = p0(0) + gauge_shift;
        p1(1) = p0(1) + gauge_shift;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc0, pbc1;
        pbc0["wall"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
        pbc1["wall"] = {ScalarBoundaryType::FIXED_VALUE, gauge_shift, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const auto phi0 = make_rhie_chow_flux_from_physical_velocity(
            mesh, geometry, U, p0, rAU, 1.0, ubc, pbc0);
        const auto phi1 = make_rhie_chow_flux_from_physical_velocity(
            mesh, geometry, U, p1, rAU, 1.0, ubc, pbc1);
        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
            EXPECT_NEAR(phi0(f), phi1(f), 1e-12);
    });

    run_case("n9_s9_rhie_chow_alternating_pressure_mode_is_coupled", [] {
        const Mesh mesh = make_four_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        U.fill(0.0);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            p(c) = (c % 2 == 0) ? 1.0 : -1.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const auto phi = make_rhie_chow_flux_from_physical_velocity(
            mesh, geometry, U, p, rAU, 1.0, ubc, pbc);
        const auto phi_hbya = make_rhie_chow_flux_from_hbya(
            mesh, geometry, U, p, rAU, 1.0, ubc, pbc);
        const auto phi_linear = make_mass_flux(mesh, geometry, U, 1.0, ubc);
        const auto grad_p = gauss_gradient_with_boundary(p, mesh, geometry, pbc);

        // Independent face-level oracle: on the central orthogonal face the
        // alternating mode has a nonzero pressure jump but a zero interpolated
        // cell-gradient contribution. The expected flux is the physical
        // velocity flux plus dAU_f*grad(p)_f.Sf minus the orthogonal pressure
        // jump response. Do not merely assert nonzero flux: that would also
        // pass for an incorrectly signed or incorrectly scaled correction.
        constexpr std::size_t f = 2;
        constexpr std::size_t owner = 1;
        constexpr std::size_t neighbour = 2;
        const double dAU_face_x = 0.5 * (
            geometry.cell_volumes[owner] * rAU[0][owner] +
            geometry.cell_volumes[neighbour] * rAU[0][neighbour]);
        const double grad_face_x = 0.5 * (
            grad_p.component_data(0)[owner] +
            grad_p.component_data(0)[neighbour]);
        const double dx = geometry.cell_centres[neighbour].x -
                          geometry.cell_centres[owner].x;
        const double sf_x = geometry.face_area_vectors[f].x;
        const double pressure_jump_response =
            dAU_face_x * (p(neighbour) - p(owner)) / dx * sf_x;
        const double compensation =
            dAU_face_x * grad_face_x * sf_x;
        const double expected = phi_linear(f) + compensation -
                                pressure_jump_response;
        EXPECT_NEAR(grad_face_x, 0.0, 1e-12);
        EXPECT_NEAR(phi_linear(f), 0.0, 1e-12);
        EXPECT_NEAR(phi_hbya(f), expected, 1e-12);
        EXPECT_NEAR(phi(f), expected, 1e-12);
        EXPECT_NEAR(phi(f) - phi_linear(f), -2.0, 1e-12);
        EXPECT_TRUE(std::abs(phi(f) - phi_linear(f)) > 1.0);

        // Grid-level pressure-mode response: compare the corrected authoritative
        // flux with the velocity-only interpolated flux for an alternating
        // pressure field. This is an operator-response metric, not a claim of
        // checkerboard damping in a converged pressure-velocity solution.
        std::vector<double> corrected_divergence(mesh.n_cells(), 0.0);
        std::vector<double> naive_divergence(mesh.n_cells(), 0.0);
        for (std::size_t face = 0; face < mesh.n_faces(); ++face) {
            const std::size_t face_owner = mesh.ownership().owner(face);
            corrected_divergence[face_owner] += phi(face);
            naive_divergence[face_owner] += phi_linear(face);
            const auto face_neighbour = mesh.ownership().neighbour(face);
            if (face_neighbour >= 0) {
                const auto neighbour_cell = static_cast<std::size_t>(face_neighbour);
                corrected_divergence[neighbour_cell] -= phi(face);
                naive_divergence[neighbour_cell] -= phi_linear(face);
            }
        }
        double corrected_l1 = 0.0;
        double naive_l1 = 0.0;
        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
            corrected_l1 += std::abs(corrected_divergence[cell]);
            naive_l1 += std::abs(naive_divergence[cell]);
        }
        EXPECT_TRUE(corrected_l1 > 1.0);
        EXPECT_NEAR(naive_l1, 0.0, 1e-12);
    });


    run_case("n9_s9_rhie_chow_correction_scales_linearly_with_rAU_majumdar", [] {
        // D^RC = V/a_P^r = alpha_u * V/a_P = alpha_u * D^pc: the production
        // RC term weights the face pressure flux with the RELAXED inverse
        // momentum diagonal. Scaling rAU by alpha_u must therefore scale the
        // RC correction by exactly alpha_u. This is the operator-level
        // statement of the Majumdar relation documented in
        // N9_S9_RHIE_CHOW_CONSISTENCY.md.
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p(0) = 1.0;
        p(1) = -1.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const double c1 = rc_correction(
            1, mesh, geometry, U, p, rAU, 1.0, ubc, pbc);
        EXPECT_TRUE(std::abs(c1) > 1e-12);
        for (const double alpha : {0.65, 0.8}) {
            const auto rAU_a = scale_rAU(rAU, alpha);
            const double ca = rc_correction(
                1, mesh, geometry, U, p, rAU_a, 1.0, ubc, pbc);
            EXPECT_NEAR(ca, alpha * c1, 1e-12);
        }
    });

    run_case("n9_s9_rhie_chow_correction_is_linear_in_pressure", [] {
        // The internal RC operator is linear in p (the Green-Gauss gradient
        // and the owner-to-neighbour difference are linear, and the
        // zero-gradient wall treatment preserves linearity):
        // corr(p_a + p_b) = corr(p_a) + corr(p_b).
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> pa(2, "pa", "Pa", 1);
        Field<double, Location::CELL> pb(2, "pb", "Pa", 1);
        Field<double, Location::CELL> psum(2, "psum", "Pa", 1);
        U.fill(0.0);
        pa(0) = 1.0;
        pa(1) = -1.0;      // alternating mode
        pb(0) = 1.0;
        pb(1) = 3.0;       // linear mode
        psum(0) = pa(0) + pb(0);
        psum(1) = pa(1) + pb(1);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const double ca = rc_correction(
            1, mesh, geometry, U, pa, rAU, 1.0, ubc, pbc);
        const double cb = rc_correction(
            1, mesh, geometry, U, pb, rAU, 1.0, ubc, pbc);
        const double cs = rc_correction(
            1, mesh, geometry, U, psum, rAU, 1.0, ubc, pbc);
        EXPECT_TRUE(std::abs(ca) > 1e-12);
        EXPECT_TRUE(std::abs(cb) > 1e-12);
        EXPECT_NEAR(cs, ca + cb, 1e-12);
    });

    run_case("n9_s9_rhie_chow_directional_rAU_face_normal_projection", [] {
        // rfn = cx*ex^2 + cy*ey^2 + cz*ez^2: only the face-normal projection
        // of the directional inverse diagonals enters the correction. The
        // internal face normal is +x, so anisotropic rAU in y/z must not
        // change the correction, while scaling the x component scales it
        // exactly.
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p(0) = 1.0;
        p(1) = -1.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU_iso = unit_rAU(mesh.n_cells());
        std::array<std::vector<double>, 3> rAU_aniso{};
        rAU_aniso[0] = std::vector<double>(mesh.n_cells(), 1.0);
        rAU_aniso[1] = std::vector<double>(mesh.n_cells(), 2.0);
        rAU_aniso[2] = std::vector<double>(mesh.n_cells(), 2.0);
        std::array<std::vector<double>, 3> rAU_x2{};
        rAU_x2[0] = std::vector<double>(mesh.n_cells(), 2.0);
        rAU_x2[1] = std::vector<double>(mesh.n_cells(), 1.0);
        rAU_x2[2] = std::vector<double>(mesh.n_cells(), 1.0);

        const double c_iso = rc_correction(
            1, mesh, geometry, U, p, rAU_iso, 1.0, ubc, pbc);
        const double c_aniso = rc_correction(
            1, mesh, geometry, U, p, rAU_aniso, 1.0, ubc, pbc);
        const double c_x2 = rc_correction(
            1, mesh, geometry, U, p, rAU_x2, 1.0, ubc, pbc);
        EXPECT_TRUE(std::abs(c_iso) > 1e-12);
        EXPECT_NEAR(c_aniso, c_iso, 1e-12);
        EXPECT_NEAR(c_x2, 2.0 * c_iso, 1e-12);
    });

    run_case("n9_s9_rhie_chow_stretched_mesh_volume_weighted_dAU", [] {
        // On a stretched mesh with unequal cell volumes and per-cell inverse
        // diagonals, the face coefficient must be the volume-weighted dAU
        // interpolation dAU_f = 0.5*(V_P*rAU_P + V_N*rAU_N) applied along the
        // owner-to-neighbour line:
        // corr = rho * dAU_f * (p_N - p_P)/d * (S_f . e).
        // The expected value is rebuilt from the same FvGeometry the
        // production path uses - no second implementation of the operator.
        const Mesh mesh = make_two_cell_channel_stretched();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p(0) = 1.0;
        p(1) = -1.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const double rAU0 = 0.5;
        const double rAU1 = 1.7;
        std::array<std::vector<double>, 3> rAU{};
        for (std::size_t k = 0; k < 3; ++k) {
            rAU[k] = std::vector<double>(mesh.n_cells(), 1.0);
            rAU[k][0] = rAU0;
            rAU[k][1] = rAU1;
        }

        const double corr = rc_correction(
            1, mesh, geometry, U, p, rAU, 1.0, ubc, pbc);

        const Vec3 c0 = geometry.cell_centres[0];
        const Vec3 c1 = geometry.cell_centres[1];
        const double dx = c1.x - c0.x;
        const double dy = c1.y - c0.y;
        const double dz = c1.z - c0.z;
        const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        const double ex = dx / d;
        const double ey = dy / d;
        const double ez = dz / d;
        const Vec3 Sf = geometry.face_area_vectors[1];
        const double Sf_dot_e = Sf.x * ex + Sf.y * ey + Sf.z * ez;
        // The mesh is orthogonal: S_non = 0, so the deferred non-orthogonal
        // part vanishes and the expected correction is exact.
        const double dAU_f = 0.5 * (geometry.cell_volumes[0] * rAU0 +
                                   geometry.cell_volumes[1] * rAU1);
        const double expected =
            1.0 * dAU_f * ((p(1) - p(0)) / d) * Sf_dot_e;

        EXPECT_TRUE(std::abs(geometry.cell_volumes[0] - 0.5) < 1e-12);
        EXPECT_TRUE(std::abs(geometry.cell_volumes[1] - 1.5) < 1e-12);
        EXPECT_NEAR(corr, expected, 1e-12);
    });

    run_case("n9_s9_rhie_chow_uniform_mesh_limit_exact_face_operator", [] {
        // On the uniform unit two-cell channel the production correction
        // must collapse to the analytic face operator
        // rho * V * rAU * (p_N - p_P)/d * (S_f . e) with V0 = V1, and the
        // geometry must be the unit construction the analytic constants
        // assume (V = d = S_f.e = 1).
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p(0) = 1.0;
        p(1) = 3.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const double corr = rc_correction(
            1, mesh, geometry, U, p, rAU, 1.0, ubc, pbc);

        const Vec3 c0 = geometry.cell_centres[0];
        const Vec3 c1 = geometry.cell_centres[1];
        const double dx = c1.x - c0.x;
        const double dy = c1.y - c0.y;
        const double dz = c1.z - c0.z;
        const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        const double ex = dx / d;
        const double ey = dy / d;
        const double ez = dz / d;
        const Vec3 Sf = geometry.face_area_vectors[1];
        const double Sf_dot_e = Sf.x * ex + Sf.y * ey + Sf.z * ez;

        EXPECT_NEAR(geometry.cell_volumes[0], 1.0, 1e-12);
        EXPECT_NEAR(geometry.cell_volumes[1], 1.0, 1e-12);
        EXPECT_NEAR(d, 1.0, 1e-12);
        EXPECT_NEAR(Sf_dot_e, 1.0, 1e-12);

        const double dAU_f = 0.5 * (geometry.cell_volumes[0] + geometry.cell_volumes[1]);
        const double expected =
            1.0 * dAU_f * ((p(1) - p(0)) / d) * Sf_dot_e;
        EXPECT_NEAR(corr, expected, 1e-12);
    });


    run_case("n9_s9_pressure_matrix_coefficient_matches_rhie_chow_on_skew_face", [] {
        // The implicit pressure matrix must use the same centre-line
        // projection as rhie_chow_pressure_flux_internal(). On this skew face
        // d points along x while Sf has a tangential y component. The
        // directional inverse diagonal is intentionally anisotropic so that
        // projecting it on Sf (rather than d) is also detected.
        const Vec3 dvec{1.0, 0.0, 0.0};
        const Vec3 Sf{1.0, 0.4, 0.0};
        const Vec3 dAU_diagonal{1.0, 4.0, 9.0};

        const double coefficient =
            pressure_face_orthogonal_coefficient(Sf, dvec, dAU_diagonal);
        EXPECT_NEAR(coefficient, 1.0, 1e-12);

        const double area = Sf.mag();
        const Vec3 face_normal{Sf.x / area, Sf.y / area, Sf.z / area};
        const double old_face_normal_rfn =
            dAU_diagonal.x * face_normal.x * face_normal.x +
            dAU_diagonal.y * face_normal.y * face_normal.y +
            dAU_diagonal.z * face_normal.z * face_normal.z;
        const double old_area_over_distance = old_face_normal_rfn * area / dvec.mag();

        // This is deliberately not a loose numerical tolerance: the old
        // area/d plus face-normal projection differs by more than 0.5 on this
        // analytic geometry, while the corrected coefficient is exactly 1.
        EXPECT_TRUE(std::abs(old_area_over_distance - coefficient) > 0.4);

        // Orthogonal limit: retain the previous coefficient exactly even for
        // anisotropic momentum response.
        const double orthogonal = pressure_face_orthogonal_coefficient(
            Vec3{1.0, 0.0, 0.0}, dvec, dAU_diagonal);
        EXPECT_NEAR(orthogonal, 1.0, 1e-12);
    });


    run_case("n9_s9_rhie_chow_skewed_face_uses_shared_implicit_coefficient", [] {
        const Mesh mesh = make_two_cell_channel();
        auto geometry = build_fv_geometry(mesh);
        geometry.face_area_vectors[1].y = 0.4;
        // Rotate the centre-line as well: anisotropic response then couples
        // the centre-line and tangential components in global coordinates.
        geometry.cell_centres[1].y += 0.3;

        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        Field<double, Location::CELL> grad_p(2, "grad_p", "Pa/m", 3);
        p(0) = 0.0;
        p(1) = 0.6;
        // A linear tangential pressure gradient with anisotropic dAU exposes
        // whether the deferred term preserves the full face flux.
        grad_p.set(0, 0.0, 2.0, 0.0);
        grad_p.set(1, 0.0, 2.0, 0.0);

        std::array<std::vector<double>, 3> rAU{
            std::vector<double>(mesh.n_cells(), 1.0),
            std::vector<double>(mesh.n_cells(), 4.0),
            std::vector<double>(mesh.n_cells(), 9.0)};
        const Vec3 dvec = geometry.cell_centres[1] - geometry.cell_centres[0];
        const Vec3 Sf = geometry.face_area_vectors[1];
        const Vec3 dAU_face{
            0.5 * (geometry.cell_volumes[0] * rAU[0][0] +
                   geometry.cell_volumes[1] * rAU[0][1]),
            0.5 * (geometry.cell_volumes[0] * rAU[1][0] +
                   geometry.cell_volumes[1] * rAU[1][1]),
            0.5 * (geometry.cell_volumes[0] * rAU[2][0] +
                   geometry.cell_volumes[1] * rAU[2][1])};
        const double d = dvec.mag();
        const Vec3 e{dvec.x / d, dvec.y / d, dvec.z / d};
        const double orthogonal_area = Sf.dot(e);
        const Vec3 gp{0.0, 2.0, 0.0};
        const double rfn = dAU_face.x * e.x * e.x +
                           dAU_face.y * e.y * e.y +
                           dAU_face.z * e.z * e.z;
        const double implicit = pressure_face_orthogonal_coefficient(
            Sf, dvec, dAU_face) * (p(1) - p(0));
        const double full_face_flux = dAU_face.x * gp.x * Sf.x +
                                      dAU_face.y * gp.y * Sf.y +
                                      dAU_face.z * gp.z * Sf.z;
        const double deferred = full_face_flux -
                               rfn * gp.dot(e) * orthogonal_area;
        const double expected = implicit + deferred;
        EXPECT_NEAR(expected, full_face_flux, 1e-12);
        const double actual = rhie_chow_pressure_flux_internal(
            mesh, geometry, 1, p, grad_p, rAU, 1.0);
        EXPECT_NEAR(actual, expected, 1e-12);
    });

    run_case("n9_s9_skew_boundary_deferred_flux_preserves_linear_pressure", [] {
        // A prescribed boundary pressure lies along an oblique centre-to-face
        // vector. The implicit jump plus deferred remainder must reproduce the
        // full anisotropic pressure-response flux for a linear pressure field.
        const Vec3 Sf{1.0, 0.4, 0.0};
        const Vec3 dvec{1.0, 0.3, 0.0};
        const Vec3 dAU{1.0, 4.0, 9.0};
        const Vec3 grad_p{0.0, 2.0, 0.0};
        const double pressure_jump = grad_p.dot(dvec);
        const double implicit = pressure_face_orthogonal_coefficient(
            Sf, dvec, dAU) * pressure_jump;
        const double deferred = pressure_face_nonorthogonal_flux(
            Sf, dvec, dAU, grad_p);
        const double full_flux = dAU.x * grad_p.x * Sf.x +
                                 dAU.y * grad_p.y * Sf.y +
                                 dAU.z * grad_p.z * Sf.z;
        EXPECT_NEAR(implicit + deferred, full_flux, 1e-12);
        EXPECT_TRUE(std::abs(deferred) > 1.0);
    });


    run_case("n9_s9_rhie_chow_oblique_3d_face_operator_and_relaxation_scaling", [] {
        // An affine shear creates a valid 3D hexahedral mesh whose internal
        // face normal is not parallel to the owner-neighbour centre vector.
        // For an exact linear pressure field, the production face operator
        // (#923 contract: implicit orthogonal coefficient plus deferred
        // non-orthogonal flux) must return the full anisotropic directional
        // flux rho * sum_d dAU_d * grad_d(p) * Sf_d, with
        // dAU_d = 0.5*(V_P*rAU_d[P] + V_N*rAU_d[N]). The complete pressure
        // flux must then scale linearly with the relaxed inverse momentum
        // diagonal; no extra alpha_u factor belongs in Rhie-Chow.
        const Mesh mesh = make_two_cell_channel_oblique();
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        Field<double, Location::CELL> grad_p(2, "grad_p", "Pa/m", 3);
        const Vec3 exact_gradient{2.0, -1.25, 0.75};
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            const Vec3 centre = geometry.cell_centres[c];
            p(c) = exact_gradient.dot(centre);
            grad_p.set(c, exact_gradient.x, exact_gradient.y, exact_gradient.z);
        }

        const Vec3 Sf = geometry.face_area_vectors[1];
        const Vec3 dvec = geometry.cell_centres[1] - geometry.cell_centres[0];
        const double d = dvec.mag();
        const Vec3 e{dvec.x / d, dvec.y / d, dvec.z / d};
        const double orthogonal_area = Sf.dot(e);
        const Vec3 Snon = Sf - e * orthogonal_area;
        EXPECT_TRUE(Snon.mag() > 1e-3);
        EXPECT_TRUE(std::abs(exact_gradient.dot(Snon)) > 1e-3);

        std::array<std::vector<double>, 3> rAU{
            std::vector<double>{0.7, 1.3},
            std::vector<double>{1.1, 0.6},
            std::vector<double>{0.9, 1.8}};
        constexpr double rho = 1.7;
        const double actual = rhie_chow_pressure_flux_internal(
            mesh, geometry, 1, p, grad_p, rAU, rho);

        const double cx = 0.5 * (geometry.cell_volumes[0] * rAU[0][0] +
                                 geometry.cell_volumes[1] * rAU[0][1]);
        const double cy = 0.5 * (geometry.cell_volumes[0] * rAU[1][0] +
                                 geometry.cell_volumes[1] * rAU[1][1]);
        const double cz = 0.5 * (geometry.cell_volumes[0] * rAU[2][0] +
                                 geometry.cell_volumes[1] * rAU[2][1]);
        // Oracle adapted to the #923 operator contract: the implicit
        // orthogonal coefficient plus the deferred non-orthogonal flux
        // equals the full directional tensor contraction for a linear
        // pressure field.
        const double expected = rho * (cx * exact_gradient.x * Sf.x +
                                       cy * exact_gradient.y * Sf.y +
                                       cz * exact_gradient.z * Sf.z);
        EXPECT_NEAR(actual, expected, 1e-12);

        constexpr double alpha_u = 0.4;
        const auto relaxed_rAU = scale_rAU(rAU, alpha_u);
        const double relaxed_actual = rhie_chow_pressure_flux_internal(
            mesh, geometry, 1, p, grad_p, relaxed_rAU, rho);
        EXPECT_NEAR(relaxed_actual, alpha_u * actual, 1e-12);
    });

    run_case("n9_s9_flux_velocity_mismatch_detects_known_inconsistency", [] {
        // Sensitivity proof for the restored independent audit channel:
        // phi_U (plain interpolation of the corrected velocity) and phi_RC
        // (authoritative Rhie-Chow flux) are two independent paths. For the
        // checkerboard pressure mode the channels MUST disagree by the known
        // Rhie-Chow correction (the channel is not blind); for the consistent
        // linear pressure field the internal-face mismatch MUST read zero
        // (the channel does not report an artifact when the identity holds).
        const Mesh mesh = make_four_cell_channel({0.0, 1.0, 2.0, 3.0, 4.0}, true);
        const auto geometry = build_fv_geometry(mesh);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ubc["inlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ubc["outlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        pbc["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
        pbc["outlet"] = {ScalarBoundaryType::FIXED_VALUE, 4.0, 0.0};
        const auto rAU = unit_rAU(mesh.n_cells());

        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        U.fill(0.0);
        Field<double, Location::CELL> p_check(mesh.n_cells(), "p", "Pa", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            p_check(c) = (c % 2 == 0) ? 1.0 : -1.0;
        // Compare against the physical cell velocity, not an H/A predictor.
        // The checkerboard pressure jump must remain visible in the physical
        // Rhie-Chow flux, while an affine pressure field must be exactly
        // compensated by the physical-velocity contract.
        const auto phi_rc_check = make_rhie_chow_mass_flux(
            mesh, geometry, U, p_check, rAU, 1.0, ubc, pbc,
            RhieChowVelocityContract::PHYSICAL_VELOCITY);
        const auto phi_u = make_mass_flux(mesh, geometry, U, 1.0, ubc);
        const double mismatch_check = std::abs(phi_rc_check(1) - phi_u(1));
        EXPECT_TRUE(mismatch_check > 1e-12);

        Field<double, Location::CELL> p_linear(mesh.n_cells(), "p", "Pa", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            p_linear(c) = geometry.cell_centres[c].x;
        const auto phi_rc_linear = make_rhie_chow_mass_flux(
            mesh, geometry, U, p_linear, rAU, 1.0, ubc, pbc,
            RhieChowVelocityContract::PHYSICAL_VELOCITY);
        for (std::size_t f = 1; f <= 3; ++f)
            EXPECT_NEAR(phi_rc_linear(f), phi_u(f), 1e-12);
    });

    run_case("n9_s9_flux_velocity_mismatch_boundary_channel_sees_expected_correction", [] {
        // The boundary channel must not be blind either: on a fixed-pressure
        // boundary face the two reconstruction paths differ by exactly the
        // Rhie-Chow boundary correction
        // rho*(implicit_coeff*(p_bc - p_P) + deferred), evaluated with the
        // owner-cell dAU = V*rAU and the owner-centre -> face-centre
        // distance, using the same production coefficient helpers.
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        U.fill(0.0);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        p(0) = 1.0;
        p(1) = 3.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::FIXED_VALUE, 0.5, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const auto phi_rc = make_rhie_chow_mass_flux(
            mesh, geometry, U, p, rAU, 1.0, ubc, pbc);
        const auto phi_u = make_mass_flux(mesh, geometry, U, 1.0, ubc);

        const std::size_t f = 0;  // x = 0 plane, owner = cell 0
        const std::size_t o = mesh.ownership().owner(f);
        const Vec3 Sf = geometry.face_area_vectors[f];
        const Vec3 dvec = geometry.face_centres[f] - geometry.cell_centres[o];
        const auto grad_p = gauss_gradient_with_boundary(p, mesh, geometry, pbc);
        const Vec3 dAU_boundary{
            geometry.cell_volumes[o] * rAU[0][o],
            geometry.cell_volumes[o] * rAU[1][o],
            geometry.cell_volumes[o] * rAU[2][o]};
        const double implicit_coefficient =
            pressure_face_orthogonal_coefficient(Sf, dvec, dAU_boundary);
        const Vec3 gp{
            grad_p.component_data(0)[o],
            grad_p.component_data(1)[o],
            grad_p.component_data(2)[o]};
        const double deferred =
            pressure_face_nonorthogonal_flux(Sf, dvec, dAU_boundary, gp);
        const double expected_correction =
            1.0 * (implicit_coefficient * (0.5 - p(o)) + deferred);

        EXPECT_TRUE(std::abs(expected_correction) > 1e-12);
        EXPECT_NEAR(phi_u(f), 0.0, 1e-12);
        EXPECT_NEAR(phi_rc(f), -expected_correction, 1e-12);
        EXPECT_NEAR(phi_rc(f) - phi_u(f), -expected_correction, 1e-12);
    });

    run_case("n9_s9_rhie_chow_identity_detects_stale_state", [] {
        // The identity channel re-derives the authoritative flux from the
        // current state: re-deriving from the same state reads exactly zero,
        // while a stored flux built from a different (stale) state MUST be
        // detected as a nonzero identity error.
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        U.fill(0.0);
        Field<double, Location::CELL> p_a(2, "p", "Pa", 1);
        Field<double, Location::CELL> p_b(2, "p", "Pa", 1);
        p_a(0) = 1.0; p_a(1) = -1.0;
        p_b(0) = 4.0; p_b(1) = -2.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        const auto rAU = unit_rAU(mesh.n_cells());

        const auto phi_stored = make_rhie_chow_mass_flux(
            mesh, geometry, U, p_a, rAU, 1.0, ubc, pbc);
        const auto phi_rederived_same = make_rhie_chow_mass_flux(
            mesh, geometry, U, p_a, rAU, 1.0, ubc, pbc);
        const auto phi_rederived_stale = make_rhie_chow_mass_flux(
            mesh, geometry, U, p_b, rAU, 1.0, ubc, pbc);

        double identity_same = 0.0;
        double identity_stale = 0.0;
        for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
            identity_same = std::max(
                identity_same, std::abs(phi_stored(f) - phi_rederived_same(f)));
            identity_stale = std::max(
                identity_stale, std::abs(phi_stored(f) - phi_rederived_stale(f)));
        }
        EXPECT_TRUE(identity_same == 0.0);
        EXPECT_TRUE(identity_stale > 1e-12);
    });

    run_case("n9_s9_coupled_symmetry_continuity_matches_authoritative_flux", [] {
        Mesh mesh = make_two_cell_channel();
        auto& walls = mesh.boundary().patch(0).face_ids;
        walls.erase(std::remove(walls.begin(), walls.end(), FaceIndex{8}), walls.end());
        Patch symmetry;
        symmetry.name = "sym";
        symmetry.type = PatchType::SYMMETRY;
        symmetry.face_ids = {8};
        mesh.boundary().add_patch(symmetry);
        const auto geometry = build_fv_geometry(mesh);
        const std::size_t nc = mesh.n_cells();
        Field<double, Location::CELL> U_old(nc, "U", "m/s", 3);
        Field<double, Location::CELL> p_old(nc, "p", "Pa", 1);
        U_old.fill(0.0); p_old.fill(0.0);
        U_old.component_data(1)[1] = 1.0; // Not yet tangent to symmetry: flux still MUST be zero.
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
        pbc["sym"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        auto diagonal_equation = [&]() {
            ScalarEquation eq;
            eq.matrix = SparseMatrix(nc,nc);
            eq.diagonal.assign(nc,2.0);
            eq.rhs = Vector(nc,0.0);
            for (std::size_t c=0; c<nc; ++c) eq.matrix.push_back(c,c,2.0);
            eq.matrix.finalize();
            return eq;
        };
        auto ex=diagonal_equation(), ey=diagonal_equation(), ez=diagonal_equation();
        auto U = U_old; auto p = p_old;
        IncompressibleSolverControls controls;
        controls.algorithm = PressureVelocityAlgorithm::COUPLED;
        controls.coupling.schur_model = CoupledSchurModel::BlockLocal;
        bool captured = false;
        controls.coupled_assembly_callback = [&](const CoupledAssemblySnapshot& snap) {
            captured = true;
            const auto& A = *snap.matrix;
            const auto& b = *snap.rhs;
            const auto& x = *snap.initial_guess;
            // The FV oracle must suppress impermeable SYMMETRY faces exactly
            // as the COUPLED assembly does, even for a non-tangent trial U.
            const auto frozen = certify_coupled_frozen_fv_operator(
                snap, mesh, geometry, ubc, pbc, controls.density,
                controls.pressure_reference_value, x);
            std::cerr << "SYMMETRY_FROZEN_FV passed=" << frozen.passed
                      << " continuity_error=" << frozen.continuity_linf << '\n';
            EXPECT_TRUE(frozen.passed);
            const std::size_t row = 3*nc + 1;
            double algebraic = -b(row);
            for (std::size_t k=A.row_offsets_data()[row];
                 k<A.row_offsets_data()[row+1]; ++k)
                algebraic += A.values_data()[k]*x(A.columns_data()[k]);
            std::array<std::vector<double>,3> rAU{
                std::vector<double>(nc,0.5),
                std::vector<double>(nc,0.5),
                std::vector<double>(nc,0.5)};
            const auto flux = make_rhie_chow_flux_from_physical_velocity(
                mesh,geometry,U_old,p_old,rAU,1.0,ubc,pbc);
            double independent = 0.0;
            for (std::size_t f=0; f<mesh.n_faces(); ++f) {
                if (mesh.ownership().owner(f)==1) independent += flux(f);
                if (mesh.ownership().neighbour(f)==1) independent -= flux(f);
            }
            std::cerr << "SYMMETRY_COUPLED_MISMATCH=" << (algebraic-independent) << '\n';
            EXPECT_NEAR(algebraic, independent, 1e-12);
        };
        (void)solve_coupled_momentum_continuity(
            mesh,geometry,ex,ey,ez,U_old,p_old,ubc,pbc,
            1.0,0.01,0,0.0,40,1e-10,controls.diagnostics,
            controls.coupling.schur_model,controls.coupled_linear_solver,controls,
            1.0,U,p);
        EXPECT_TRUE(captured);
    });

    run_case("n9_s9_coupled_fixed_gradient_boundary_momentum_matches_gauss", [] {
        Mesh mesh = make_two_cell_channel();
        auto& walls = mesh.boundary().patch(0).face_ids;
        walls.erase(std::remove(walls.begin(), walls.end(), FaceIndex{8}), walls.end());
        Patch gradient_patch;
        gradient_patch.name = "grad";
        gradient_patch.type = PatchType::WALL;
        gradient_patch.face_ids = {8};
        mesh.boundary().add_patch(gradient_patch);
        const auto geometry = build_fv_geometry(mesh);
        const std::size_t nc = mesh.n_cells();
        Field<double, Location::CELL> U_old(nc,"U","m/s",3);
        Field<double, Location::CELL> p_old(nc,"p","Pa",1);
        U_old.fill(0.0); p_old.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0,0.0,0.0}};
        ubc["grad"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        pbc["grad"] = {ScalarBoundaryType::FIXED_GRADIENT,0.0,3.0};
        const auto true_gradient = gauss_gradient_with_boundary(p_old,mesh,geometry,pbc);
        const double expected = true_gradient.component_data(1)[1]*geometry.cell_volumes[1];
        EXPECT_TRUE(std::abs(expected)>1e-8);
        auto diagonal_equation = [&](std::size_t component) {
            ScalarEquation eq;
            eq.matrix = SparseMatrix(nc,nc);
            eq.diagonal.assign(nc,2.0);
            eq.rhs = Vector(nc,0.0);
            for (std::size_t c=0; c<nc; ++c) {
                eq.matrix.push_back(c,c,2.0);
                eq.rhs(c) = -true_gradient.component_data(component)[c]*geometry.cell_volumes[c];
            }
            eq.matrix.finalize();
            return eq;
        };
        auto ex=diagonal_equation(0), ey=diagonal_equation(1), ez=diagonal_equation(2);
        auto U=U_old; auto p=p_old;
        IncompressibleSolverControls controls;
        controls.algorithm=PressureVelocityAlgorithm::COUPLED;
        controls.coupling.schur_model=CoupledSchurModel::BlockLocal;
        bool captured = false;
        controls.coupled_assembly_callback = [&](const CoupledAssemblySnapshot& snap) {
            captured = true;
            const auto& A=*snap.matrix; const auto& b=*snap.rhs;
            const auto& x=*snap.initial_guess;
            const std::size_t row = nc + 1;
            double algebraic = -b(row);
            for (std::size_t k=A.row_offsets_data()[row];
                 k<A.row_offsets_data()[row+1]; ++k)
                algebraic += A.values_data()[k]*x(A.columns_data()[k]);
            std::cerr << "FIXED_GRADIENT_GAUSS_MISMATCH=" << (algebraic-expected) << '\n';
            EXPECT_NEAR(algebraic,expected,1e-12);
        };
        (void)solve_coupled_momentum_continuity(
            mesh,geometry,ex,ey,ez,U_old,p_old,ubc,pbc,
            1.0,0.01,0,0.0,40,1e-10,controls.diagnostics,
            controls.coupling.schur_model,controls.coupled_linear_solver,controls,
            1.0,U,p);
        EXPECT_TRUE(captured);
    });

    run_case("n9_s9_coupled_empty_face_gradient_matches_gauss", [] {
        Mesh mesh = make_two_cell_channel();
        auto& walls = mesh.boundary().patch(0).face_ids;
        walls.erase(std::remove(walls.begin(),walls.end(),FaceIndex{8}),walls.end());
        Patch empty;
        empty.name="front_back"; empty.type=PatchType::EMPTY; empty.face_ids={8};
        mesh.boundary().add_patch(empty);
        const auto geometry=build_fv_geometry(mesh);
        const std::size_t nc=mesh.n_cells();
        Field<double,Location::CELL> U_old(nc,"U","m/s",3);
        Field<double,Location::CELL> p_old(nc,"p","Pa",1);
        U_old.fill(0.0); p_old.fill(0.0); p_old(1)=4.0;
        VelocityBoundaryConditions ubc;
        ubc["wall"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        const auto grad=gauss_gradient_with_boundary(p_old,mesh,geometry,pbc);
        auto mk_eq=[&](std::size_t component) {
            ScalarEquation eq;
            eq.matrix=SparseMatrix(nc,nc);
            eq.diagonal.assign(nc,2.0);
            eq.rhs=Vector(nc,0.0);
            for(std::size_t c=0;c<nc;++c) {
                eq.matrix.push_back(c,c,2.0);
                eq.rhs(c)=-geometry.cell_volumes[c]*grad.component_data(component)[c];
            }
            eq.matrix.finalize(); return eq;
        };
        auto ex=mk_eq(0),ey=mk_eq(1),ez=mk_eq(2);
        auto U=U_old,p=p_old;
        IncompressibleSolverControls controls;
        controls.algorithm=PressureVelocityAlgorithm::COUPLED;
        bool captured=false;
        controls.coupled_assembly_callback=[&](const CoupledAssemblySnapshot& snap) {
            captured=true;
            const std::size_t row=nc+1;
            const auto& A=*snap.matrix; const auto& b=*snap.rhs;
            const auto& x=*snap.initial_guess;
            double algebraic=-b(row);
            for(std::size_t k=A.row_offsets_data()[row]; k<A.row_offsets_data()[row+1];++k)
                algebraic+=A.values_data()[k]*x(A.columns_data()[k]);
            const double expected=geometry.cell_volumes[1]*grad.component_data(1)[1];
            std::cerr << "EMPTY_GAUSS_MISMATCH=" << algebraic-expected << '\n';
            EXPECT_NEAR(algebraic,expected,1e-12);
        };
        (void)solve_coupled_momentum_continuity(
            mesh,geometry,ex,ey,ez,U_old,p_old,ubc,pbc,1.0,0.01,
            0,0.0,40,1e-10,controls.diagnostics,controls.coupling.schur_model,
            controls.coupled_linear_solver,controls,1.0,U,p);
        EXPECT_TRUE(captured);
    });

    run_case("n9_s9_coupled_fixed_pressure_warm_start_is_not_gauge_pinned", [] {
        const Mesh mesh=make_two_cell_channel();
        const auto geometry=build_fv_geometry(mesh);
        const std::size_t nc=mesh.n_cells();
        Field<double,Location::CELL> U_old(nc,"U","m/s",3);
        Field<double,Location::CELL> p_old(nc,"p","Pa",1);
        U_old.fill(0.0); p_old(0)=3.0; p_old(1)=4.0;
        VelocityBoundaryConditions ubc;
        ubc["wall"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"]={ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
        const auto grad=gauss_gradient_with_boundary(p_old,mesh,geometry,pbc);
        auto mk_eq=[&](std::size_t component) {
            ScalarEquation eq;
            eq.matrix=SparseMatrix(nc,nc); eq.diagonal.assign(nc,2.0);
            eq.rhs=Vector(nc,0.0);
            for(std::size_t c=0;c<nc;++c) {
                eq.matrix.push_back(c,c,2.0);
                eq.rhs(c)=-geometry.cell_volumes[c]*grad.component_data(component)[c];
            }
            eq.matrix.finalize();return eq;
        };
        auto ex=mk_eq(0),ey=mk_eq(1),ez=mk_eq(2);
        auto U=U_old,p=p_old;
        IncompressibleSolverControls controls;
        controls.algorithm=PressureVelocityAlgorithm::COUPLED;
        bool captured=false;
        controls.coupled_assembly_callback=[&](const CoupledAssemblySnapshot& snap) {
            captured=true;
            EXPECT_TRUE(snap.has_fixed_pressure);
            EXPECT_NEAR((*snap.initial_guess)(3*nc),3.0,0.0);
            EXPECT_NEAR((*snap.initial_guess)(3*nc+1),4.0,0.0);
        };
        (void)solve_coupled_momentum_continuity(
            mesh,geometry,ex,ey,ez,U_old,p_old,ubc,pbc,1.0,0.01,
            0,0.0,40,1e-10,controls.diagnostics,controls.coupling.schur_model,
            controls.coupled_linear_solver,controls,1.0,U,p);
        EXPECT_TRUE(captured);
    });

    return run_all();
}
