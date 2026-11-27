#pragma once
// CoupledJfnkResidualV6N: monolithic 6N residual for the U+p+k+omega JFNK
// solver (Option D). Rebuilds strain rate, eddy viscosity, F1/F2 blending,
// momentum, Rhie-Chow flux and SST transport from the trial state alone.
// No accepted/lagged state is read. Residual layout:
//   [momentum_x(N), momentum_y(N), momentum_z(N), continuity(N), k(N), omega(N)]
#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/physics/physical_six_candidate_v51.h"
#include "cfdx/physics/state_only_rhie_chow_v49.h"
#include "cfdx/physics/state_only_six_residual_v50.h"
#include "cfdx/physics/sst_production_coupling.h"
#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/sst_strain_rate.h"
#include "cfdx/physics/mesh_wall_distance.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct CoupledJfnkResidualV6N {
    static PhysicalResidualV51 evaluate(
        const ::cfdx::core::Mesh& mesh,
        const FvGeometry& geometry,
        const ::cfdx::core::Field<double, ::cfdx::core::Location::CELL>& Ux,
        const ::cfdx::core::Field<double, ::cfdx::core::Location::CELL>& Uy,
        const ::cfdx::core::Field<double, ::cfdx::core::Location::CELL>& Uz,
        const ::cfdx::core::Field<double, ::cfdx::core::Location::CELL>& p,
        const ::cfdx::core::Field<double, ::cfdx::core::Location::CELL>& k,
        const ::cfdx::core::Field<double, ::cfdx::core::Location::CELL>& omega,
        const TurbulenceTransportControls& sst_controls,
        const ScalarBoundaryConditions& ubc_x,
        const ScalarBoundaryConditions& ubc_y,
        const ScalarBoundaryConditions& ubc_z,
        const ScalarBoundaryConditions& pbc,
        const ScalarBoundaryConditions& k_bc,
        const ScalarBoundaryConditions& omega_bc,
        const std::map<std::string, std::array<double, 3>>& fixed_velocity_values,
        const std::vector<std::string>& wall_patch_names,
        double density,
        double kinematic_viscosity)
    {
        using ::cfdx::core::Location;
        using ::cfdx::core::Field;
        const std::size_t nc = mesh.n_cells();
        const std::size_t nf = mesh.n_faces();
        if (nc == 0)
            throw std::invalid_argument("JFNK residual: empty mesh");
        const auto check_field = [&](const Field<double, Location::CELL>& f, const char* name) {
            if (f.size() != nc || f.dimension() != 1)
                throw std::invalid_argument(std::string("JFNK residual: ") + name + " field shape mismatch");
        };
        check_field(Ux, "Ux"); check_field(Uy, "Uy"); check_field(Uz, "Uz");
        check_field(p, "p"); check_field(k, "k"); check_field(omega, "omega");
        if (!(density > 0.0) || !std::isfinite(density))
            throw std::invalid_argument("JFNK residual: invalid density");
        if (!(kinematic_viscosity > 0.0) || !std::isfinite(kinematic_viscosity))
            throw std::invalid_argument("JFNK residual: invalid kinematic viscosity");

        // Build the 3-component candidate velocity from the trial components.
        Field<double, Location::CELL> cand_U(nc, "U_jfnk", "m/s", 3);
        for (std::size_t i = 0; i < nc; ++i) {
            cand_U.component_data(0)[i] = Ux(i);
            cand_U.component_data(1)[i] = Uy(i);
            cand_U.component_data(2)[i] = Uz(i);
            if (!std::isfinite(Ux(i)) || !std::isfinite(Uy(i)) || !std::isfinite(Uz(i)))
                throw std::domain_error("JFNK residual: nonfinite trial velocity");
        }
        Field<double, Location::CELL> cand_p(nc, "p_jfnk", "Pa", 1);
        for (std::size_t i = 0; i < nc; ++i) {
            cand_p(i) = p(i);
            if (!std::isfinite(p(i)))
                throw std::domain_error("JFNK residual: nonfinite trial pressure");
        }

        // Build candidate SST state from trial k/omega.
        SSTProductionState cand_sst(nc);
        for (std::size_t i = 0; i < nc; ++i) {
            cand_sst.k(i) = k(i);
            cand_sst.omega(i) = omega(i);
            if (!std::isfinite(k(i)) || k(i) < 0.0)
                throw std::domain_error("JFNK residual: invalid trial k");
            if (!std::isfinite(omega(i)) || omega(i) <= 0.0)
                throw std::domain_error("JFNK residual: invalid trial omega");
        }

        // Geometry-derived wall distance from the mesh.
        const auto wall_distance =
            compute_mesh_cell_wall_distance(mesh, geometry, wall_patch_names);
        for (std::size_t i = 0; i < nc; ++i)
            cand_sst.wall_distance(i) = wall_distance.distance[i];

        // Strain rate from the trial velocity (rebuilt every evaluation).
        cand_sst.strain_rate =
            compute_sst_strain_rate_with_boundary_values(cand_U, mesh, fixed_velocity_values);

        // F1/F2 blending from the trial SST state.
        Field<double, Location::CELL> F2(nc, "F2_jfnk", "1", 1);
        for (std::size_t i = 0; i < nc; ++i) {
            F2(i) = compute_sst_blending(
                cand_sst.k(i), cand_sst.omega(i), cand_sst.wall_distance(i),
                sst_controls.molecular_viscosity, sst_controls.beta_star, 0.0,
                sst_controls.density, sst_controls.sst_sigma_w2).second;
            if (!std::isfinite(F2(i)) || F2(i) < 0.0 || F2(i) > 1.0)
                throw std::domain_error("JFNK residual: invalid F2 blending");
        }

        // Eddy viscosity from the trial k/omega/strain_rate/F2.
        refresh_sst_production_eddy_viscosity(
            cand_sst.k, cand_sst.omega, cand_sst.strain_rate, F2,
            sst_controls, cand_sst.nu_t);
        validate_sst_production_state(cand_sst, nc);

        // Effective dynamic viscosity field: mu = rho * (nu + nu_t).
        std::vector<double> mu(nc);
        for (std::size_t i = 0; i < nc; ++i) {
            mu[i] = density * (kinematic_viscosity + cand_sst.nu_t(i));
            if (!(mu[i] > 0.0) || !std::isfinite(mu[i]))
                throw std::domain_error("JFNK residual: invalid effective viscosity");
        }

        // Merge component scalar BCs into velocity BCs.
        VelocityBoundaryConditions velocity_bcs;
        std::set<std::string> patch_names;
        for (const auto& [name, _] : ubc_x) patch_names.insert(name);
        for (const auto& [name, _] : ubc_y) patch_names.insert(name);
        for (const auto& [name, _] : ubc_z) patch_names.insert(name);
        for (const auto& name : patch_names) {
            VelocityBoundaryCondition vbc;
            vbc.type = VelocityBoundaryCondition::Type::ZERO_GRADIENT;
            vbc.value = ::cfdx::core::Vec3{0.0, 0.0, 0.0};
            const std::array<const ScalarBoundaryConditions*, 3> component_bcs{{
                &ubc_x, &ubc_y, &ubc_z}};
            std::array<bool, 3> has_fixed_value{};
            std::array<double, 3> prescribed_values{};
            for (std::size_t comp = 0; comp < component_bcs.size(); ++comp) {
                const auto it = component_bcs[comp]->find(name);
                if (it == component_bcs[comp]->end())
                    continue;
                has_fixed_value[comp] =
                    it->second.type == ScalarBoundaryType::FIXED_VALUE;
                prescribed_values[comp] = it->second.value;
            }
            // VelocityBoundaryCondition stores one type for the whole vector.
            // Do not collapse component-wise mixed Dirichlet/Neumann data into
            // a false all-components FIXED_VALUE condition used by mass-flux
            // reconstruction. The current vector BC API cannot represent it.
            if (has_fixed_value[0] != has_fixed_value[1] ||
                has_fixed_value[1] != has_fixed_value[2]) {
                throw std::invalid_argument(
                    "JFNK residual: mixed component velocity BCs on patch '" +
                    name + "' are not representable by VelocityBoundaryCondition");
            }
            if (has_fixed_value[0])
                vbc.type = VelocityBoundaryCondition::Type::FIXED_VALUE;
            vbc.value = ::cfdx::core::Vec3{
                prescribed_values[0], prescribed_values[1], prescribed_values[2]};
            velocity_bcs[name] = vbc;
        }

        // Pressure gradient from the trial pressure field.
        const auto grad_p = gauss_gradient_with_boundary(
            cand_p, mesh, geometry, pbc,
            /*skew_correct=*/false, /*projected=*/false);

        // Body force (zero for the residual callback).
        std::array<Field<double, Location::CELL>, 3> forces{{
            Field<double, Location::CELL>(nc, "body_x", "N/m3", 1),
            Field<double, Location::CELL>(nc, "body_y", "N/m3", 1),
            Field<double, Location::CELL>(nc, "body_z", "N/m3", 1)}};
        for (std::size_t i = 0; i < nc; ++i) {
            forces[0](i) = 0.0;
            forces[1](i) = 0.0;
            forces[2](i) = 0.0;
        }

        // Component fields for momentum assembly.
        std::array<Field<double, Location::CELL>, 3> comps{{
            Field<double, Location::CELL>(nc, "ux", "m/s", 1),
            Field<double, Location::CELL>(nc, "uy", "m/s", 1),
            Field<double, Location::CELL>(nc, "uz", "m/s", 1)}};
        for (std::size_t d = 0; d < 3; ++d)
            for (std::size_t i = 0; i < nc; ++i)
                comps[d](i) = cand_U.component_data(d)[i];

        const std::array<const ScalarBoundaryConditions*, 3> bcs{{&ubc_x, &ubc_y, &ubc_z}};

        // Deterministic seed flux from the trial velocity (not the accepted flux).
        const auto phi_seed = make_mass_flux(
            mesh, geometry, cand_U, density,
            velocity_bcs, pbc,
            /*projected_face_interpolation=*/false);
        std::vector<double> seed(nf);
        for (std::size_t f = 0; f < nf; ++f) seed[f] = phi_seed(f);

        // Momentum equation assembler: given a face flux and component,
        // assemble the momentum scalar equation.
        const auto momentum_equation = [&](
            const Field<double, Location::FACE>& flux, std::size_t d) -> ScalarEquation {
            return assemble_momentum_component(
                mesh, geometry, flux, grad_p, forces[d], mu,
                *bcs[d], d, /*bounded=*/false,
                ConvectionScheme::UPWIND, &comps[d], nullptr,
                density * kinematic_viscosity,
                /*high_order_fraction=*/1.0,
                /*face_conditions=*/nullptr,
                /*wall_face_mu_override=*/nullptr);
        };

        // Assemble rAU from a candidate flux: momentum diagonal -> inverse.
        const auto assemble = [&](const std::vector<double>& trial_phi)
            -> std::array<std::vector<double>, 3> {
            if (trial_phi.size() != nf)
                throw std::invalid_argument("JFNK residual: flux size mismatch");
            Field<double, Location::FACE> flux(nf, "phi_trial", "kg/s", 1);
            for (std::size_t f = 0; f < nf; ++f) flux(f) = trial_phi[f];
            std::array<std::vector<double>, 3> rAU;
            for (std::size_t d = 0; d < 3; ++d) {
                const auto eq = momentum_equation(flux, d);
                rAU[d].resize(nc);
                for (std::size_t i = 0; i < nc; ++i) {
                    if (!(eq.diagonal[i] > 0.0) || !std::isfinite(eq.diagonal[i]))
                        throw std::runtime_error("JFNK residual: invalid momentum diagonal");
                    rAU[d][i] = 1.0 / eq.diagonal[i];
                }
            }
            return rAU;
        };

        // Reconstruct Rhie-Chow flux from rAU.
        const auto reconstruct = [&](const std::array<std::vector<double>, 3>& rAU)
            -> std::vector<double> {
            const auto rc = make_rhie_chow_flux_from_physical_velocity(
                mesh, geometry, cand_U, cand_p, rAU,
                density, velocity_bcs, pbc,
                /*skew_correct=*/false);
            std::vector<double> values(nf);
            for (std::size_t f = 0; f < nf; ++f) values[f] = rc(f);
            return values;
        };

        // Solve the state-only Rhie-Chow fixed point for the trial state.
        const auto rc = evaluate_state_only_rhie_chow_v49(
            seed, nc, assemble, reconstruct, 64, 1e-10, 0.5);
        if (!rc.converged)
            throw std::runtime_error("JFNK residual: Rhie-Chow fixed point did not close");
        Field<double, Location::FACE> flux_field(nf, "phi_jfnk", "kg/s", 1);
        for (std::size_t f = 0; f < nf; ++f) flux_field(f) = rc.flux[f];

        // Raw momentum residuals: b - A*U for each component.
        std::array<std::vector<double>, 3> raw_momentum;
        std::array<double, 3> rhs_scale{};
        for (std::size_t d = 0; d < 3; ++d) {
            const auto eq = momentum_equation(flux_field, d);
            raw_momentum[d].resize(nc);
            const auto* offsets = eq.matrix.row_offsets_data();
            const auto* columns = eq.matrix.columns_data();
            const auto* values = eq.matrix.values_data();
            rhs_scale[d] = 1e-300;
            for (std::size_t i = 0; i < nc; ++i)
                rhs_scale[d] = std::max(rhs_scale[d], std::abs(eq.rhs(i)));
            for (std::size_t i = 0; i < nc; ++i) {
                long double sum = -static_cast<long double>(eq.rhs(i));
                for (auto j = offsets[i]; j < offsets[i + 1]; ++j) {
                    if (columns[j] >= nc)
                        throw std::logic_error("JFNK residual: bad momentum CSR");
                    sum += static_cast<long double>(values[j]) *
                        comps[d](columns[j]);
                }
                raw_momentum[d][i] = static_cast<double>(sum);
            }
        }

        // Continuity: divergence of the converged face flux.
        std::vector<double> div(nc, 0.0);
        double flux_l1 = 0.0;
        for (std::size_t f = 0; f < nf; ++f) {
            const double mf = rc.flux[f];
            flux_l1 += std::abs(mf);
            div[mesh.ownership().owner(f)] += mf;
            const auto adjacent = mesh.ownership().neighbour(f);
            if (adjacent >= 0)
                div[static_cast<std::size_t>(adjacent)] -= mf;
        }

        // SST k and omega residuals on the trial state.
        const auto sst = sst_coupled_residual_v22(
            mesh, geometry, flux_field, cand_sst.k, cand_sst.omega,
            cand_sst.strain_rate, cand_sst.wall_distance,
            sst_controls, k_bc, omega_bc,
            /*omega_face_values=*/nullptr,
            /*wall_production=*/nullptr);

        // Assemble the 6N raw residual.
        PhysicalResidualV51 residual;
        residual.signed_raw.resize(6 * nc);
        for (std::size_t d = 0; d < 3; ++d)
            std::copy(raw_momentum[d].begin(), raw_momentum[d].end(),
                      residual.signed_raw.begin() + static_cast<std::ptrdiff_t>(d * nc));
        std::copy(div.begin(), div.end(),
                  residual.signed_raw.begin() + static_cast<std::ptrdiff_t>(3 * nc));
        std::copy(sst.values.begin(), sst.values.end(),
                  residual.signed_raw.begin() + static_cast<std::ptrdiff_t>(4 * nc));

        // Normalized diagnostic (not differenced by J*v).
        residual.normalized.resize(6 * nc);
        double global_rhs = 0.0;
        for (double s : rhs_scale) global_rhs = std::max(global_rhs, s);
        for (std::size_t cell = 0; cell < nc; ++cell) {
            for (std::size_t d = 0; d < 6; ++d) {
                double raw = d < 3 ? raw_momentum[d][cell]
                    : d == 3 ? div[cell]
                    : sst.values[(d - 4) * nc + cell];
                double divisor = d < 3 ? std::max(rhs_scale[d], 1e-8 * global_rhs)
                    : d == 3 ? std::max(flux_l1, 1e-300)
                    : sst.row_norms[(d - 4) * nc + cell] * sst.local_scales[(d - 4) * nc + cell];
                residual.normalized[d * nc + cell] = raw / divisor;
            }
        }
        residual.face_flux = rc.flux;
        residual.inverse_momentum_diagonal.reserve(3 * nc);
        for (const auto& d : rc.inverse_diagonal)
            residual.inverse_momentum_diagonal.insert(
                residual.inverse_momentum_diagonal.end(), d.begin(), d.end());
        residual.eddy_viscosity.resize(nc);
        for (std::size_t i = 0; i < nc; ++i)
            residual.eddy_viscosity[i] = cand_sst.nu_t(i);
        residual.flux_fixed_point_defect = rc.relative_linf;
        residual.flux_closed = rc.converged;
        double merit = 0.0;
        for (double v : residual.signed_raw)
            merit = std::max(merit, std::abs(v));
        residual.merit = merit;
        return residual;
    }
};

} // namespace cfdx::physics
