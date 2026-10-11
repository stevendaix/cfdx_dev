#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

// Two-cell channel with a deliberately skewed internal face. The internal
// face corners are moved off x=1 and off the y=0/y=1 planes so that
//   * dc != dn across the internal face (the Gauss gradient then uses a
//     distance weight dn/(dc+dn) != 0.5, which is exactly the weighting the
//     coupled momentum rows must reproduce), and
//   * the internal face normal gains a y component (the y-momentum rows then
//     also see a non-trivial pressure operator, not only the x rows).
// Boundary pressure is Dirichlet on inlet/outlet and zero-gradient on the
// walls, so the probe exercises the implicit interior operator, the
// owner-extrapolation zero-gradient path AND the Dirichlet displacement
// recorded during assembly.
static Mesh make_skewed_two_cell_channel(double x_face_low, double x_face_high)
{
    Mesh m;
    m.points().resize(12);
    const double xyz[12][3] = {
        {0,0,0},{x_face_low,0,0},{x_face_high,1,0},{0,1,0},
        {0,0,1},{x_face_low,0,1},{x_face_high,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i, xyz[i][0], xyz[i][1], xyz[i][2]);

    m.faces().push_face({0,3,7,4});       // x = 0 (inlet)
    m.faces().push_face({1,2,6,5});       // internal, skewed
    m.faces().push_face({8,9,11,10});     // x = 2 (outlet)
    m.faces().push_face({0,1,5,4});       // y = 0, cell 0
    m.faces().push_face({3,7,6,2});       // y = 1, cell 0
    m.faces().push_face({0,3,2,1});       // z = 0, cell 0
    m.faces().push_face({4,5,6,7});       // z = 1, cell 0
    m.faces().push_face({1,8,10,5});       // y = 0, cell 1
    m.faces().push_face({2,6,11,9});       // y = 1, cell 1
    m.faces().push_face({1,2,9,8});       // z = 0, cell 1
    m.faces().push_face({5,10,11,6});     // z = 1, cell 1

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 11; ++f) {
        const auto owner = (f == 2 || f >= 7) ? 1 : 0;
        m.ownership().set_owner(f, owner);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);

    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});

    Patch inlet;
    inlet.name = "inlet";
    inlet.type = PatchType::WALL;
    inlet.face_ids = {0};
    m.boundary().add_patch(inlet);
    Patch outlet;
    outlet.name = "outlet";
    outlet.type = PatchType::WALL;
    outlet.face_ids = {2};
    m.boundary().add_patch(outlet);
    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {3,4,5,6,7,8,9,10};
    m.boundary().add_patch(wall);
    return m;
}

struct ProbeCase {
    double x_face_low;
    double x_face_high;
    double probe_alpha_u;
    double probe_alpha_p;
};

// Runs one coupled solve on the skewed channel with the operator probe
// enabled and returns the filled probe plus the solver status.
static CoupledPressureOperatorDiagnostics run_probe_case(
    const ProbeCase& probe_case,
    cfdx::core::SolverStatus& status_out,
    double& block_residual_inf_out)
{
    const Mesh mesh = make_skewed_two_cell_channel(
        probe_case.x_face_low, probe_case.x_face_high);
    const auto geometry = build_fv_geometry(mesh);

    const std::size_t nc = mesh.n_cells();
    Field<double, Location::CELL> U_old(nc, "U_old", "m/s", 3);
    Field<double, Location::CELL> p_old(nc, "p_old", "Pa", 1);
    for (std::size_t c = 0; c < nc; ++c) {
        U_old.set(c, 1.0, 0.0, 0.0);
        p_old(c) = 0.8 - 0.6 * static_cast<double>(c);
    }

    VelocityBoundaryConditions ubc;
    ubc["inlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {1.0, 0.0, 0.0}};
    ubc["outlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {1.0, 0.0, 0.0}};
    ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
    ScalarBoundaryConditions pbc;
    pbc["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 1.0, 0.0};
    pbc["outlet"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
    pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

    const double rho = 1.0;
    const double kinematic_viscosity = 0.1;
    const auto mass_flux = make_mass_flux(mesh, geometry, U_old, rho, ubc);
    const auto grad_p_old =
        gauss_gradient_with_boundary(p_old, mesh, geometry, pbc);

    Field<double, Location::CELL> body_x(nc, "body_x", "N/m3", 1);
    Field<double, Location::CELL> body_y(nc, "body_y", "N/m3", 1);
    Field<double, Location::CELL> body_z(nc, "body_z", "N/m3", 1);
    body_x.fill(0.0);
    body_y.fill(0.0);
    body_z.fill(0.0);

    Field<double, Location::CELL> u_x_field(nc, "Ux_reconstruction", "m/s", 1);
    Field<double, Location::CELL> u_y_field(nc, "Uy_reconstruction", "m/s", 1);
    Field<double, Location::CELL> u_z_field(nc, "Uz_reconstruction", "m/s", 1);
    for (std::size_t c = 0; c < nc; ++c) {
        u_x_field(c) = U_old.component_data(0)[c];
        u_y_field(c) = U_old.component_data(1)[c];
        u_z_field(c) = U_old.component_data(2)[c];
    }

    ScalarBoundaryConditions ubc_x, ubc_y, ubc_z;
    for (const auto& [name, bc] : ubc) {
        ubc_x[name] = {bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                           ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::ZERO_GRADIENT,
                       bc.value.x, 0.0};
        ubc_y[name] = {bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                           ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::ZERO_GRADIENT,
                       bc.value.y, 0.0};
        ubc_z[name] = {bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                           ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::ZERO_GRADIENT,
                       bc.value.z, 0.0};
    }

    const double mu_eff = rho * kinematic_viscosity;
    const auto ex = assemble_momentum_component(
        mesh, geometry, mass_flux, grad_p_old, body_x, mu_eff, ubc_x, 0,
        true, ConvectionScheme::UPWIND, &u_x_field);
    const auto ey = assemble_momentum_component(
        mesh, geometry, mass_flux, grad_p_old, body_y, mu_eff, ubc_y, 1,
        true, ConvectionScheme::UPWIND, &u_y_field);
    const auto ez = assemble_momentum_component(
        mesh, geometry, mass_flux, grad_p_old, body_z, mu_eff, ubc_z, 2,
        true, ConvectionScheme::UPWIND, &u_z_field);

    Field<double, Location::CELL> U(nc, "U", "m/s", 3);
    Field<double, Location::CELL> p(nc, "p", "Pa", 1);

    DiagnosticsControls diagnostics;
    const cfdx::core::LinearSolverRequest solver_request;
    const IncompressibleSolverControls controls;
    CoupledPressureOperatorDiagnostics probe;
    const auto result = solve_coupled_momentum_continuity(
        mesh, geometry, ex, ey, ez, U_old, p_old,
        ubc, pbc, rho, kinematic_viscosity,
        0, 0.0,
        500, 1e-12,
        diagnostics,
        CoupledSchurModel::BlockLocal,
        solver_request,
        controls,
        1.0,
        U, p, nullptr,
        &probe,
        probe_case.probe_alpha_u,
        probe_case.probe_alpha_p);
    status_out = result.status;
    block_residual_inf_out = probe.block_residual_inf;
    return probe;
}

// The probe tolerance is tied to double round-off and the assembled force
// scale, never calibrated against a case: a row sums O(10) coefficients of
// magnitude O(force_scale), so round-off noise is a few eps*scale; the bound
// below sits three orders of magnitude above that noise and orders of
// magnitude below any real operator discrepancy (an arithmetic 0.5 face
// weight where the Gauss gradient uses dn/(dc+dn) produces an O(0.1)
// relative error on the skewed face).
static double probe_absolute_tolerance(double force_scale)
{
    const double scale = std::max(1.0, std::abs(force_scale));
    return std::max(1e-12 * scale, 500.0 * 2.220446049250313e-16 * scale);
}

int main()
{
    run_case("coupled_operator_probe_matches_gauss_on_skewed_mesh", [] {
        cfdx::core::SolverStatus status = cfdx::core::SolverStatus::DIVERGED;
        double block_residual_inf = -1.0;
        const ProbeCase probe_case{0.65, 0.85, 0.7, 0.4};
        const auto probe =
            run_probe_case(probe_case, status, block_residual_inf);

        EXPECT_TRUE(status == cfdx::core::SolverStatus::CONVERGED);
        EXPECT_TRUE(probe.valid);

        // Non-vacuity: the solved pressure field must exercise the pressure
        // operator with forces far above round-off, otherwise dg=0 would be
        // trivially true and the equivalence untested.
        EXPECT_TRUE(probe.force_scale > 1e-6);

        // Operator equivalence (Steven P2): the assembled coupled
        // pressure-gradient operator must reproduce the independent Gauss
        // gradient on the same state, to round-off, not to a calibrated
        // tolerance.
        EXPECT_TRUE(probe.dg_max <= probe_absolute_tolerance(probe.force_scale));
        EXPECT_TRUE(probe.dg_rel_max <= 1e-10);

        // Post-relaxation state (Picard accepted fields): the same operator
        // identity must hold on the state the nonlinear iteration accepts.
        EXPECT_TRUE(probe.force_scale_post > 1e-6);
        EXPECT_TRUE(probe.dg_max_post <=
            probe_absolute_tolerance(probe.force_scale_post));
        EXPECT_TRUE(probe.dg_rel_max_post <= 1e-10);

        // The linear solution of the block system must be a real solution:
        // the probe's independent max|Ax-b| stays at the linear tolerance
        // scale before relaxation.
        EXPECT_TRUE(std::isfinite(block_residual_inf));
        EXPECT_TRUE(block_residual_inf <= 1e-9);
    });

    run_case("coupled_operator_probe_matches_gauss_on_heavier_skew", [] {
        cfdx::core::SolverStatus status = cfdx::core::SolverStatus::DIVERGED;
        double block_residual_inf = -1.0;
        // Stronger face skew (dc:dn = 0.4:1.2 along the y=0 corner) and unit
        // relaxation factors, so the post-relaxation state coincides with
        // the raw linear solution.
        const ProbeCase probe_case{0.4, 0.7, 1.0, 1.0};
        const auto probe =
            run_probe_case(probe_case, status, block_residual_inf);

        EXPECT_TRUE(status == cfdx::core::SolverStatus::CONVERGED);
        EXPECT_TRUE(probe.valid);
        EXPECT_TRUE(probe.force_scale > 1e-6);
        EXPECT_TRUE(probe.dg_max <= probe_absolute_tolerance(probe.force_scale));
        EXPECT_TRUE(probe.dg_rel_max <= 1e-10);
        EXPECT_TRUE(probe.dg_max_post <=
            probe_absolute_tolerance(probe.force_scale_post));
        EXPECT_TRUE(probe.dg_rel_max_post <= 1e-10);
    });

    return 0;
}
