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
                1.0, std::abs(reconstructed_physical) +
                     std::abs(residual.physical_residual));
            const double bounded_scale = std::max(
                1.0, std::abs(reconstructed_bounded) +
                     std::abs(residual.bounded_operator_residual));
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

} // namespace cfdx::physics
