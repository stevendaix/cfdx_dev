#pragma once

#include "cfdx/physics/momentum_face_budget.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

// Power contributions are formed from already integrated cell momentum
// contributions [N] paired with accepted velocity [m/s]. Values are [W].
// This is an aggregation layer: it does not reconstruct or re-evaluate fluxes.
struct MomentumCellPowerTerms {
    double convection = 0.0;
    double viscous = 0.0;
    double transpose_viscous = 0.0;
    double pressure = 0.0;
    double body_force = 0.0;
    double temporal = 0.0;
    double physical_residual = 0.0;
    // Algorithmic bounded-convection correction; deliberately not included
    // in physical_residual.
    double bounded_divergence = 0.0;
    double bounded_operator_residual = 0.0;
};

struct MomentumVectorPowerBudget {
    std::vector<MomentumCellPowerTerms> cells;
    MomentumCellPowerTerms total;
    // Preserve, rather than replace, the independent component-wise CSR gate.
    std::array<bool, 3> matrix_closure_checked{false, false, false};
    std::array<double, 3> matrix_closure_linf{};
    std::array<double, 3> matrix_closure_relative_linf{};
};

// Aggregate the three accepted-state scalar MomentumFaceBudget objects.
// Each budget entry is already an integrated signed cell contribution, so
// cell volume must NOT be applied again. The resulting velocity pairing is
// mechanical power. This function makes no stability or physical-sign claim.
inline MomentumVectorPowerBudget aggregate_momentum_vector_power(
    const std::array<const MomentumFaceBudget*, 3>& component_budgets,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity)
{
    if (velocity.dimension() != 3 || velocity.size() == 0)
        throw std::invalid_argument("momentum vector power: velocity must be a nonempty 3-component cell field");

    const std::size_t nc = velocity.size();
    for (const auto* budget : component_budgets) {
        if (budget == nullptr || budget->cells.size() != nc)
            throw std::invalid_argument("momentum vector power: missing component budget or cell-count mismatch");
    }

    MomentumVectorPowerBudget out;
    out.cells.resize(nc);
    for (std::size_t d = 0; d < 3; ++d) {
        const auto& budget = *component_budgets[d];
        out.matrix_closure_checked[d] = budget.matrix_closure_checked;
        out.matrix_closure_linf[d] = budget.matrix_closure_linf;
        out.matrix_closure_relative_linf[d] = budget.matrix_closure_relative_linf;
    }

    const auto add = [](double& target, double value) {
        target += value;
    };
    for (std::size_t cell = 0; cell < nc; ++cell) {
        auto& power = out.cells[cell];
        for (std::size_t d = 0; d < 3; ++d) {
            const double u = velocity(cell, d);
            const auto& residual = component_budgets[d]->cells[cell];
            const std::array<double, 11> inputs{
                u, residual.convection, residual.viscous,
                residual.transpose_viscous, residual.pressure,
                residual.body_force, residual.temporal,
                residual.physical_residual, residual.bounded_divergence,
                residual.bounded_operator_residual,
                component_budgets[d]->matrix_closure_checked
                    ? component_budgets[d]->matrix_closure_linf : 0.0
            };
            for (const double value : inputs) {
                if (!std::isfinite(value))
                    throw std::invalid_argument("momentum vector power: nonfinite velocity or residual contribution");
            }
            if (component_budgets[d]->matrix_closure_checked &&
                (!std::isfinite(component_budgets[d]->matrix_closure_relative_linf) ||
                 component_budgets[d]->matrix_closure_linf < 0.0 ||
                 component_budgets[d]->matrix_closure_relative_linf < 0.0))
                throw std::invalid_argument("momentum vector power: invalid checked CSR closure diagnostics");

            const double reconstructed_physical =
                residual.convection + residual.viscous +
                residual.transpose_viscous + residual.pressure +
                residual.body_force + residual.temporal;
            const double reconstructed_bounded =
                reconstructed_physical + residual.bounded_divergence;
            const double physical_scale = std::max(
                std::abs(residual.convection) + std::abs(residual.viscous) +
                std::abs(residual.transpose_viscous) + std::abs(residual.pressure) +
                std::abs(residual.body_force) + std::abs(residual.temporal),
                std::numeric_limits<double>::min());
            const double bounded_scale = std::max(
                physical_scale + std::abs(residual.bounded_divergence),
                std::numeric_limits<double>::min());
            constexpr double consistency_factor = 64.0;
            const double arithmetic_tolerance =
                consistency_factor * std::numeric_limits<double>::epsilon();
            if (std::abs(reconstructed_physical - residual.physical_residual) >
                    arithmetic_tolerance * physical_scale ||
                std::abs(reconstructed_bounded - residual.bounded_operator_residual) >
                    arithmetic_tolerance * bounded_scale)
                throw std::invalid_argument("momentum vector power: inconsistent physical/bounded residual decomposition");

            add(power.convection, u * residual.convection);
            add(power.viscous, u * residual.viscous);
            add(power.transpose_viscous, u * residual.transpose_viscous);
            add(power.pressure, u * residual.pressure);
            add(power.body_force, u * residual.body_force);
            add(power.temporal, u * residual.temporal);
            add(power.physical_residual, u * residual.physical_residual);
            add(power.bounded_divergence, u * residual.bounded_divergence);
            add(power.bounded_operator_residual,
                u * (residual.physical_residual + residual.bounded_divergence));
        }
        const std::array<double, 9> cell_values{
            power.convection, power.viscous, power.transpose_viscous,
            power.pressure, power.body_force, power.temporal,
            power.physical_residual, power.bounded_divergence,
            power.bounded_operator_residual
        };
        for (const double value : cell_values) {
            if (!std::isfinite(value))
                throw std::invalid_argument("momentum vector power: nonfinite cell power");
        }
        add(out.total.convection, power.convection);
        add(out.total.viscous, power.viscous);
        add(out.total.transpose_viscous, power.transpose_viscous);
        add(out.total.pressure, power.pressure);
        add(out.total.body_force, power.body_force);
        add(out.total.temporal, power.temporal);
        add(out.total.physical_residual, power.physical_residual);
        add(out.total.bounded_divergence, power.bounded_divergence);
        add(out.total.bounded_operator_residual, power.bounded_operator_residual);
    }

    const std::array<double, 9> totals{
        out.total.convection, out.total.viscous, out.total.transpose_viscous,
        out.total.pressure, out.total.body_force, out.total.temporal,
        out.total.physical_residual, out.total.bounded_divergence,
        out.total.bounded_operator_residual
    };
    for (const double value : totals) {
        if (!std::isfinite(value))
            throw std::invalid_argument("momentum vector power: nonfinite global power");
    }
    return out;
}


using MomentumVector3 = std::array<double, 3>;

// Independent raw-face oracle inputs. Face vectors are integrated forces [N]
// oriented owner -> neighbour; source vectors are already volume-integrated
// forces [N]. They must be taken from the accepted operators, not reconstructed
// from MomentumFaceBudget's cell-scattered output.
struct MomentumRawFacePowerInput {
    std::vector<MomentumVector3> convection_face_forces;
    std::vector<MomentumVector3> viscous_face_forces;
    std::vector<MomentumVector3> transpose_viscous_face_forces;
    std::vector<MomentumVector3> pressure_cell_forces;
    std::vector<MomentumVector3> body_force_cell_forces;
    std::vector<MomentumVector3> temporal_cell_forces;
};

struct MomentumRawFacePower {
    double convection = 0.0;
    double viscous = 0.0;
    double transpose_viscous = 0.0;
    double pressure = 0.0;
    double body_force = 0.0;
    double temporal = 0.0;
    double physical_residual = 0.0;
};

// Evaluate the finite-volume work identity directly from one shared force
// vector per face and integrated cell sources. This path does not consume
// MomentumFaceBudget and is intended as an independent audit oracle.
inline MomentumRawFacePower evaluate_momentum_raw_face_power(
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
    const MomentumRawFacePowerInput& input)
{
    const std::size_t nc = mesh.n_cells();
    const std::size_t nf = mesh.n_faces();
    if (nc == 0 || velocity.size() != nc || velocity.dimension() != 3 ||
        input.convection_face_forces.size() != nf ||
        input.viscous_face_forces.size() != nf ||
        input.transpose_viscous_face_forces.size() != nf ||
        input.pressure_cell_forces.size() != nc ||
        input.body_force_cell_forces.size() != nc ||
        input.temporal_cell_forces.size() != nc)
        throw std::invalid_argument("raw momentum face power: field shape mismatch");

    const auto finite_vector = [](const MomentumVector3& value) {
        return std::isfinite(value[0]) && std::isfinite(value[1]) &&
               std::isfinite(value[2]);
    };
    const auto dot = [](const MomentumVector3& a, const MomentumVector3& b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };
    const auto velocity_at = [&](std::size_t cell) {
        MomentumVector3 value{
            velocity(cell, 0), velocity(cell, 1), velocity(cell, 2)};
        if (!finite_vector(value))
            throw std::invalid_argument("raw momentum face power: nonfinite velocity");
        return value;
    };
    const auto face_work = [&](const std::vector<MomentumVector3>& forces) {
        double work = 0.0;
        const auto& ownership = mesh.ownership();
        for (std::size_t face = 0; face < nf; ++face) {
            const auto& force = forces[face];
            if (!finite_vector(force))
                throw std::invalid_argument("raw momentum face power: nonfinite face force");
            const std::size_t owner = ownership.owner(face);
            if (owner >= nc)
                throw std::invalid_argument("raw momentum face power: invalid owner");
            const auto owner_velocity = velocity_at(owner);
            const auto neighbour = ownership.neighbour(face);
            if (neighbour >= 0) {
                const auto neighbour_cell = static_cast<std::size_t>(neighbour);
                if (neighbour_cell >= nc || neighbour_cell == owner)
                    throw std::invalid_argument("raw momentum face power: invalid neighbour");
                const auto neighbour_velocity = velocity_at(neighbour_cell);
                MomentumVector3 jump{
                    owner_velocity[0] - neighbour_velocity[0],
                    owner_velocity[1] - neighbour_velocity[1],
                    owner_velocity[2] - neighbour_velocity[2]};
                work += dot(jump, force);
            } else {
                if (neighbour != cfdx::core::FaceOwnership::BOUNDARY)
                    throw std::invalid_argument("raw momentum face power: invalid boundary sentinel");
                work += dot(owner_velocity, force);
            }
            if (!std::isfinite(work))
                throw std::invalid_argument("raw momentum face power: nonfinite accumulated work");
        }
        return work;
    };
    const auto source_work = [&](const std::vector<MomentumVector3>& forces) {
        double work = 0.0;
        for (std::size_t cell = 0; cell < nc; ++cell) {
            if (!finite_vector(forces[cell]))
                throw std::invalid_argument("raw momentum face power: nonfinite cell source");
            work += dot(velocity_at(cell), forces[cell]);
            if (!std::isfinite(work))
                throw std::invalid_argument("raw momentum face power: nonfinite accumulated source work");
        }
        return work;
    };

    MomentumRawFacePower out;
    out.convection = face_work(input.convection_face_forces);
    out.viscous = face_work(input.viscous_face_forces);
    out.transpose_viscous = face_work(input.transpose_viscous_face_forces);
    out.pressure = source_work(input.pressure_cell_forces);
    out.body_force = source_work(input.body_force_cell_forces);
    out.temporal = source_work(input.temporal_cell_forces);
    out.physical_residual = out.convection + out.viscous +
        out.transpose_viscous + out.pressure + out.body_force + out.temporal;
    if (!std::isfinite(out.physical_residual))
        throw std::invalid_argument("raw momentum face power: nonfinite total");
    return out;
}

} // namespace cfdx::physics
