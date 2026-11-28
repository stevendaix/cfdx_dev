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

// Convert the three scalar velocity BC maps into the vector BC contract used
// by mass-flux reconstruction. The current vector API supports only a single
// FIXED_VALUE/ZERO_GRADIENT type per patch; unsupported scalar types and
// component-wise mixed types must fail rather than be silently reinterpreted.
inline VelocityBoundaryConditions make_jfnk_velocity_boundary_conditions(
    const ScalarBoundaryConditions& ubc_x,
    const ScalarBoundaryConditions& ubc_y,
    const ScalarBoundaryConditions& ubc_z)
{
    VelocityBoundaryConditions velocity_bcs;
    std::set<std::string> patch_names;
    for (const auto& [name, _] : ubc_x) patch_names.insert(name);
    for (const auto& [name, _] : ubc_y) patch_names.insert(name);
    for (const auto& [name, _] : ubc_z) patch_names.insert(name);

    const std::array<const ScalarBoundaryConditions*, 3> component_bcs{{
        &ubc_x, &ubc_y, &ubc_z}};
    for (const auto& name : patch_names) {
        VelocityBoundaryCondition vbc;
        vbc.type = VelocityBoundaryCondition::Type::ZERO_GRADIENT;
        vbc.value = ::cfdx::core::Vec3{0.0, 0.0, 0.0};
        std::array<bool, 3> has_fixed_value{};
        std::array<double, 3> prescribed_values{};
        for (std::size_t comp = 0; comp < component_bcs.size(); ++comp) {
            const auto it = component_bcs[comp]->find(name);
            if (it == component_bcs[comp]->end())
                continue;
            switch (it->second.type) {
            case ScalarBoundaryType::FIXED_VALUE:
                has_fixed_value[comp] = true;
                prescribed_values[comp] = it->second.value;
                if (!std::isfinite(prescribed_values[comp]))
                    throw std::invalid_argument(
                        "JFNK residual: nonfinite fixed velocity on patch '" + name + "'");
                break;
            case ScalarBoundaryType::ZERO_GRADIENT:
                break;
            default:
                throw std::invalid_argument(
                    "JFNK residual: scalar velocity BC type is not representable "
                    "by VelocityBoundaryCondition on patch '" + name + "'");
            }
        }
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
    return velocity_bcs;
}

// The 6N residual contains continuity but no explicit pressure-gauge row.
// A pure-Neumann pressure system therefore retains its constant null mode.
inline bool jfnk_pressure_has_dirichlet_reference(
    const ScalarBoundaryConditions& pressure_bcs)
{
    for (const auto& [_, bc] : pressure_bcs)
        if (bc.type == ScalarBoundaryType::FIXED_VALUE)
            return true;
    return false;
}

// A face-valued pressure trace can depend on the trial SST state (for
// modified pressure p*=p+2 rho k/3). The current JFNK callback does not
// rebuild that trace for each trial vector, so using one here would freeze a
// state-dependent boundary term and produce an inconsistent matrix-free Jv.
inline bool jfnk_has_face_valued_pressure_snapshot(
    const ScalarBoundaryConditions& pressure_bcs)
{
    for (const auto& [_, bc] : pressure_bcs)
        if (bc.face_values != nullptr)
            return true;
    return false;
}

// Return a finite divisor without inventing a physical scale for a
// nonzero residual. A zero denominator is admissible only for an exactly
// zero residual; the unit sentinel then contributes exactly zero to the merit.
inline double jfnk_residual_scale_or_exact_zero(
    double raw_residual, double candidate_scale, const char* equation_name)
{
    if (!std::isfinite(raw_residual) || !std::isfinite(candidate_scale) ||
        candidate_scale < 0.0)
        throw std::runtime_error(std::string("JFNK residual: invalid ") + equation_name + " scale");
    if (candidate_scale == 0.0) {
        if (raw_residual != 0.0)
            throw std::runtime_error(std::string("JFNK residual: nonzero ") +
                equation_name + " residual has zero backward-error scale");
        return 1.0;
    }
    return candidate_scale;
}

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
        if (jfnk_has_face_valued_pressure_snapshot(pbc))
            throw std::invalid_argument(
                "JFNK residual: state-dependent face-valued pressure traces are not implemented");
        if (!jfnk_pressure_has_dirichlet_reference(pbc))
            throw std::invalid_argument(
                "JFNK residual: pure-Neumann pressure requires an explicit gauge, not implemented");
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

        // Convert component scalar BCs without losing their mathematical type.
        const auto velocity_bcs = make_jfnk_velocity_boundary_conditions(
            ubc_x, ubc_y, ubc_z);

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

        // Raw momentum residuals and componentwise backward-error scales.
        // For each row, |b_i| + sum_j |A_ij| |U_j| is invariant to cancellation
        // in b and remains meaningful when the assembled RHS is nearly zero.
        std::array<std::vector<double>, 3> raw_momentum;
        std::array<std::vector<double>, 3> momentum_row_scales;
        for (std::size_t d = 0; d < 3; ++d) {
            const auto eq = momentum_equation(flux_field, d);
            raw_momentum[d].resize(nc);
            momentum_row_scales[d].resize(nc);
            const auto* offsets = eq.matrix.row_offsets_data();
            const auto* columns = eq.matrix.columns_data();
            const auto* values = eq.matrix.values_data();
            for (std::size_t i = 0; i < nc; ++i) {
                long double sum = -static_cast<long double>(eq.rhs(i));
                long double row_scale = std::abs(static_cast<long double>(eq.rhs(i)));
                for (auto j = offsets[i]; j < offsets[i + 1]; ++j) {
                    if (columns[j] >= nc)
                        throw std::logic_error("JFNK residual: bad momentum CSR");
                    const long double coefficient = values[j];
                    const long double state_value = comps[d](columns[j]);
                    sum += coefficient * state_value;
                    row_scale += std::abs(coefficient) * std::abs(state_value);
                }
                raw_momentum[d][i] = static_cast<double>(sum);
                momentum_row_scales[d][i] = static_cast<double>(row_scale);
            }
        }

        // Continuity: divergence and local absolute-flux scale on the same
        // cell stencil. This is a local mass-imbalance ratio, not a global
        // flux L1 norm that can hide a poorly balanced low-flow cell.
        std::vector<double> div(nc, 0.0);
        std::vector<double> continuity_row_scales(nc, 0.0);
        for (std::size_t f = 0; f < nf; ++f) {
            const double mf = rc.flux[f];
            const auto owner = mesh.ownership().owner(f);
            div[owner] += mf;
            continuity_row_scales[owner] += std::abs(mf);
            const auto adjacent = mesh.ownership().neighbour(f);
            if (adjacent >= 0) {
                const auto neighbour = static_cast<std::size_t>(adjacent);
                div[neighbour] -= mf;
                continuity_row_scales[neighbour] += std::abs(mf);
            }
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

        // Per-equation physical residual scales. These are frozen by the
        // JFNK driver at the accepted base state and are never recomputed
        // inside a directional difference or candidate line search.
        residual.normalized.resize(6 * nc);
        residual.normalization_scales.resize(6 * nc);
        for (std::size_t cell = 0; cell < nc; ++cell) {
            for (std::size_t d = 0; d < 6; ++d) {
                const double raw = d < 3 ? raw_momentum[d][cell]
                    : d == 3 ? div[cell]
                    : sst.values[(d - 4) * nc + cell];
                const double candidate_scale = d < 3 ? momentum_row_scales[d][cell]
                    : d == 3 ? continuity_row_scales[cell]
                    : sst.row_norms[(d - 4) * nc + cell] * sst.local_scales[(d - 4) * nc + cell];
                const char* equation_name = d < 3 ? "momentum" : d == 3 ? "continuity" : "SST";
                const double divisor = jfnk_residual_scale_or_exact_zero(
                    raw, candidate_scale, equation_name);
                const auto index = d * nc + cell;
                residual.normalization_scales[index] = divisor;
                residual.normalized[index] = raw / divisor;
                if (!std::isfinite(residual.normalized[index]))
                    throw std::runtime_error(
                        "JFNK residual: normalized equation residual is nonfinite");
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
