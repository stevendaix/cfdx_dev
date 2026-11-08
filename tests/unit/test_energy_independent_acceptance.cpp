#include "audit/chain_fixture.h"
#include "cfdx/physics/energy_solver.h"
#include "common/test_harness.h"

#include <cmath>
#include <iostream>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {
struct EnergyFixture {
    Mesh mesh;
    FvGeometry geometry;
    Field<double, Location::FACE> flux;
    Field<double, Location::CELL> temperature;
    Field<double, Location::CELL> source;
    ScalarBoundaryConditions bcs;
    EnergySolverControls controls;

    explicit EnergyFixture(std::size_t cells)
        : mesh(make_1d_chain(cells)), geometry(build_fv_geometry(mesh)),
          flux(mesh.n_faces(), "mass_flux", "kg/s", 1),
          temperature(cells, "T", "K", 1),
          source(cells, "source", "W/m3", 1) {
        flux.fill(0.0);
        temperature.fill(250.0);
        source.fill(0.0);
        bcs["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 300.0, 0.0};
        bcs["outlet"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        bcs["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        controls.density=1.0;
        controls.cp=1.0;
        controls.conductivity=0.0;
        controls.relaxation=1.0;
        controls.max_iterations=3;
        controls.tolerance=1e-12;
    }

    void set_flux(double inlet, double internal, double outlet) {
        for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
            const auto n=mesh.ownership().neighbour(f);
            if (n >= 0) {
                flux(f)=internal;
            } else {
                const double x=geometry.face_centres[f].x;
                if (std::abs(x) < 1e-12) flux(f)=inlet;
                if (std::abs(x-1.0) < 1e-12) flux(f)=outlet;
            }
        }
    }
    EnergySolveResult solve() {
        return solve_energy(mesh, geometry, flux, temperature, source, controls, bcs);
    }
};
}

int main() {
    run_case("balanced_flow_is_accepted", [] {
        EnergyFixture f(1);
        f.set_flux(-1.0, 0.0, 1.0);
        const auto result=f.solve();
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(f.temperature(0), 300.0, 1e-10);
        EXPECT_NEAR(result.history.back().independent_conservation_normalized, 0.0, 1e-12);
    });

    run_case("nonconservative_flow_must_not_look_converged", [] {
        EnergyFixture f(1);
        f.set_flux(-1.0, 0.0, 0.5);
        const auto result=f.solve();
        // The bounded-convection correction can make Ax=b exact even when
        // the physical face-enthalpy flux has a nonzero divergence.
        EXPECT_TRUE(result.linear_solver_converged);
        EXPECT_TRUE(result.history.back().energy_imbalance < 1e-12);
        EXPECT_TRUE(result.history.back().independent_conservation_normalized > 0.1);
        EXPECT_FALSE(result.converged);
    });

    run_case("local_defects_cannot_cancel_globally", [] {
        EnergyFixture f(2);
        f.set_flux(-1.0, 2.0, 1.0);
        const auto result=f.solve();
        EXPECT_TRUE(result.linear_solver_converged);
        EXPECT_TRUE(result.history.back().energy_imbalance < 1e-12);
        EXPECT_NEAR(result.history.back().independent_conservation_normalized, 0.0, 1e-12);
        const auto local=reconstruct_energy_balance(
            f.mesh, f.geometry, f.flux, f.temperature, f.temperature,
            f.source, f.controls, f.bcs);
        EXPECT_TRUE(local.max_cell_residual_normalized > 0.1);
        EXPECT_FALSE(result.converged);
    });
    run_case("native_temperature_dirichlet_linear_profile_and_balance", [] {
        // Representative physical boundary check: a 1 m rod, constant lambda,
        // no flow, T(0)=300 K, T(1)=400 K. Exact T(x)=300+100*x.
        EnergyFixture f(16);
        f.controls.conductivity = 2.0;
        f.controls.max_iterations = 20;
        f.controls.tolerance = 1e-11;
        f.temperature.fill(280.0);
        f.bcs["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 300.0, 0.0};
        f.bcs["outlet"] = {ScalarBoundaryType::FIXED_VALUE, 400.0, 0.0};
        const auto result = f.solve();
        EXPECT_TRUE(result.converged);
        EXPECT_TRUE(result.history.back().independent_conservation_accepted);
        for (std::size_t i=0; i<16; ++i) {
            const double x = (static_cast<double>(i)+0.5)/16.0;
            EXPECT_NEAR(f.temperature(i), 300.0+100.0*x, 1e-8);
        }
        EXPECT_TRUE(result.history.back().independent_conservation_local_normalized <= 1e-10);
    });
    run_case("manufactured_conduction_requires_independent_closure", [] {
        // T = x^2, div(-grad T) = -2. A three-mesh thermal campaign
        // exposed matrix-relative convergence before the independent watts
        // balance closed; the solver may refine Krylov, but never the gate.
        EnergyFixture f(16);
        f.controls.conductivity = 1.0;
        f.controls.max_iterations = 16;
        f.controls.tolerance = 1e-11;
        f.temperature.fill(0.25);
        f.source.fill(-2.0);
        f.bcs["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
        f.bcs["outlet"] = {ScalarBoundaryType::FIXED_VALUE, 1.0, 0.0};
        const auto result = f.solve();
        EXPECT_TRUE(result.converged);
        EXPECT_TRUE(!result.history.empty());
        EXPECT_TRUE(result.history.back().independent_conservation_accepted);
        EXPECT_TRUE(result.history.back().independent_conservation_normalized <= 1e-11);
        EXPECT_TRUE(result.history.back().independent_conservation_local_normalized <= 1e-11);
        for (std::size_t i=0; i<16; ++i) {
            const double x = (static_cast<double>(i)+0.5)/16.0;
            EXPECT_NEAR(f.temperature(i), x*x, 0.00098);
        }
    });
    return run_all();
}
