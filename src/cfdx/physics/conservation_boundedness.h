#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/numerics/conservation.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

struct TransportConservationReport {
    double boundary_flux = 0.0;
    double source_sum = 0.0;
    double accumulation = 0.0;
    double residual = 0.0;
    double normalized_residual = 0.0;
    double max_cell_residual = 0.0;
    std::size_t worst_cell = 0;
    std::size_t nonfinite_cells = 0;
    std::size_t nonfinite_sources = 0;

    bool finite() const noexcept {
        return nonfinite_cells == 0 && nonfinite_sources == 0 &&
               std::isfinite(boundary_flux) && std::isfinite(source_sum) &&
               std::isfinite(accumulation) && std::isfinite(residual);
    }
};

inline TransportConservationReport audit_transport_balance(
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& face_flux,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& source,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>* accumulation = nullptr)
{
    if (source.size() != mesh.n_cells() || source.dimension() != 1)
        throw std::invalid_argument(
            "audit_transport_balance: source field does not match mesh");
    const auto r = cfdx::core::audit_cell_balance(
        mesh, face_flux, source, accumulation);

    TransportConservationReport out;
    out.boundary_flux =
        cfdx::core::audit_face_flux_conservation(mesh, face_flux).global_boundary_flux;
    out.source_sum = r.source_integral;
    out.accumulation = r.accumulation;
    out.residual = r.residual;
    out.normalized_residual = r.normalized_residual;
    out.max_cell_residual = r.max_cell_residual;
    out.worst_cell = r.worst_cell;
    out.nonfinite_cells =
        r.nonfinite_balance + r.nonfinite_accumulation;
    out.nonfinite_sources = r.nonfinite_source;
    return out;
}

struct ScalarBoundReport {
    double minimum=std::numeric_limits<double>::infinity();
    double maximum=-std::numeric_limits<double>::infinity();
    double lower=0.0;
    double upper=std::numeric_limits<double>::infinity();
    std::size_t nonfinite=0;
    std::size_t below=0;
    std::size_t above=0;
    std::size_t worst_cell=0;
    double worst_violation=0.0;

    bool bounded() const noexcept {
        return nonfinite==0 && below==0 && above==0;
    }
};

inline ScalarBoundReport audit_scalar_bounds(
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& field,
    double lower, double upper=std::numeric_limits<double>::infinity(),
    double tolerance=0.0)
{
    const auto d = cfdx::core::audit_boundedness(
        field, lower, upper, tolerance);
    ScalarBoundReport out;
    out.minimum = d.minimum;
    out.maximum = d.maximum;
    out.lower = lower;
    out.upper = upper;
    out.nonfinite = d.nonfinite;
    out.below = d.below_lower;
    out.above = d.above_upper;
    out.worst_cell = d.worst_cell;
    out.worst_violation = d.worst_violation;
    return out;
}


inline ScalarBoundReport audit_positive_scalar(
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& field,
    double floor=0.0)
{
    if (!std::isfinite(floor) || floor < 0.0)
        throw std::invalid_argument("audit_positive_scalar: invalid floor");
    return audit_scalar_bounds(field, floor);
}


} // namespace cfdx::physics
