#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::core {

struct ConservationDiagnostics {
    double global_boundary_flux = 0.0;
    double global_abs_boundary_flux = 0.0;
    double max_local_imbalance = 0.0;
    double l1_local_imbalance = 0.0;
    double l2_local_imbalance = 0.0;
    double global_cell_balance = 0.0;
    std::size_t worst_cell = 0;
    double worst_cell_imbalance = 0.0;
    std::size_t internal_faces = 0;
    std::size_t boundary_faces = 0;
    std::size_t nonfinite_faces = 0;
};

struct BoundednessDiagnostics {
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    std::size_t below_lower = 0;
    std::size_t above_upper = 0;
    std::size_t nonfinite = 0;
    double minimum_violation = 0.0;
    double maximum_violation = 0.0;
    std::size_t worst_cell = 0;
    double worst_violation = 0.0;

    bool finite() const noexcept { return nonfinite == 0; }
    bool bounded() const noexcept {
        return finite() && below_lower == 0 && above_upper == 0;
    }
};


struct FaceAssemblyDiagnostics {
    std::size_t internal_faces = 0;
    std::size_t boundary_faces = 0;
    std::size_t mismatched_faces = 0;
    std::size_t nonfinite_faces = 0;
    std::size_t worst_face = 0;
    double max_abs_mismatch = 0.0;
    double l1_mismatch = 0.0;
    double tolerance = 0.0;

    bool finite() const noexcept { return nonfinite_faces == 0; }
    bool antisymmetric() const noexcept {
        return finite() && mismatched_faces == 0;
    }
};

// Audit independently assembled owner/neighbour face contributions. Unlike
// audit_face_flux_conservation(), this does not reconstruct the neighbour
// contribution from a single face flux, so it can expose an assembly bug.
inline FaceAssemblyDiagnostics audit_face_assembly(
    const Mesh& mesh,
    const std::vector<double>& owner_contribution,
    const std::vector<double>& neighbour_contribution,
    double tolerance = 0.0)
{
    const std::size_t nf = mesh.n_faces();
    if (owner_contribution.size() != nf ||
        neighbour_contribution.size() != nf)
        throw std::invalid_argument(
            "audit_face_assembly: contribution size does not match mesh");
    if (!std::isfinite(tolerance) || tolerance < 0.0)
        throw std::invalid_argument("audit_face_assembly: invalid tolerance");

    FaceAssemblyDiagnostics out;
    out.tolerance = tolerance;
    const auto& own = mesh.ownership();

    for (std::size_t f = 0; f < nf; ++f) {
        const double owner = owner_contribution[f];
        const double neighbour = neighbour_contribution[f];
        if (!std::isfinite(owner) || !std::isfinite(neighbour)) {
            ++out.nonfinite_faces;
            continue;
        }

        const bool internal = own.neighbour(f) >= 0;
        if (internal) {
            ++out.internal_faces;
        } else {
            ++out.boundary_faces;
        }

        // Internal face contributions must cancel. A boundary face has no
        // neighbour equation contribution and therefore expects zero there.
        const double mismatch = internal ? owner + neighbour : neighbour;
        const double abs_mismatch = std::abs(mismatch);
        out.l1_mismatch += abs_mismatch;
        if (abs_mismatch > out.max_abs_mismatch) {
            out.max_abs_mismatch = abs_mismatch;
            out.worst_face = f;
        }
        if (abs_mismatch > tolerance)
            ++out.mismatched_faces;
    }
    return out;
};

inline std::vector<double> reconstruct_cell_balance(
    const Mesh& mesh,
    const Field<double, Location::FACE>& face_flux)
{
    const std::size_t nf = mesh.n_faces();
    const std::size_t nc = mesh.n_cells();
    if (face_flux.size() != nf || face_flux.dimension() != 1)
        throw std::invalid_argument(
            "reconstruct_cell_balance: face flux size/dimension mismatch");

    std::vector<double> balance(nc, 0.0);
    const auto& own = mesh.ownership();
    for (std::size_t f = 0; f < nf; ++f) {
        const double flux = face_flux(f);
        if (!std::isfinite(flux))
            continue;
        const std::size_t owner = own.owner(f);
        if (owner >= nc)
            throw std::runtime_error("reconstruct_cell_balance: invalid owner");
        balance[owner] += flux;
        const auto neighbour = own.neighbour(f);
        if (neighbour >= 0) {
            const std::size_t n = static_cast<std::size_t>(neighbour);
            if (n >= nc)
                throw std::runtime_error(
                    "reconstruct_cell_balance: invalid neighbour");
            balance[n] -= flux;
        }
    }
    return balance;
};

inline ConservationDiagnostics audit_face_flux_conservation(
    const Mesh& mesh,
    const Field<double, Location::FACE>& face_flux)
{
    const std::size_t nf = mesh.n_faces();
    const std::size_t nc = mesh.n_cells();
    if (face_flux.size() != nf || face_flux.dimension() != 1)
        throw std::invalid_argument("audit_face_flux_conservation: face flux size/dimension mismatch");

    ConservationDiagnostics out;
    std::vector<double> cell_balance(nc, 0.0);
    const auto& own = mesh.ownership();

    for (std::size_t f = 0; f < nf; ++f) {
        const double flux = face_flux(f);
        if (!std::isfinite(flux)) {
            ++out.nonfinite_faces;
            continue;
        }

        const std::size_t owner = own.owner(f);
        if (owner >= nc)
            throw std::runtime_error("audit_face_flux_conservation: invalid owner");

        const auto neighbour = own.neighbour(f);
        cell_balance[owner] += flux;

        if (neighbour >= 0) {
            const std::size_t n = static_cast<std::size_t>(neighbour);
            if (n >= nc)
                throw std::runtime_error("audit_face_flux_conservation: invalid neighbour");
            cell_balance[n] -= flux;
            ++out.internal_faces;
        } else {
            out.global_boundary_flux += flux;
            out.global_abs_boundary_flux += std::abs(flux);
            ++out.boundary_faces;
        }
    }

    double sum_sq = 0.0;
    for (std::size_t c = 0; c < cell_balance.size(); ++c) {
        const double b = cell_balance[c];
        const double a = std::abs(b);
        out.global_cell_balance += b;
        if (a > out.worst_cell_imbalance) {
            out.worst_cell_imbalance = a;
            out.worst_cell = c;
        }
        out.max_local_imbalance = std::max(out.max_local_imbalance, a);
        out.l1_local_imbalance += a;
        sum_sq += b * b;
    }
    out.l2_local_imbalance = std::sqrt(sum_sq);

    // An internal face is stored once; owner/neighbour antisymmetry is therefore
    // enforced by construction. Pairwise mismatch must be audited at assembly
    // time if both owner and neighbour contributions are stored independently.

    return out;
}

inline BoundednessDiagnostics audit_boundedness(
    const std::vector<double>& values,
    double lower,
    double upper,
    double tolerance = 0.0)
{
    if (values.empty())
        throw std::invalid_argument("audit_boundedness: empty field");
    if (!std::isfinite(lower) || (std::isfinite(upper) && lower > upper) ||
        (std::isinf(upper) && upper < 0.0) ||
        !std::isfinite(tolerance) || tolerance < 0.0)
        throw std::invalid_argument("audit_boundedness: invalid bounds/tolerance");

    BoundednessDiagnostics out;
    const double lo = lower - tolerance;
    const double hi = upper + tolerance;
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double v = values[i];
        if (!std::isfinite(v)) {
            ++out.nonfinite;
            continue;
        }
        out.minimum = std::min(out.minimum, v);
        out.maximum = std::max(out.maximum, v);
        if (v < lo) {
            ++out.below_lower;
            const double violation = lo - v;
            out.minimum_violation = std::max(out.minimum_violation, violation);
            if (violation > out.worst_violation) { out.worst_violation = violation; out.worst_cell = i; }
        }
        if (std::isfinite(hi) && v > hi) {
            ++out.above_upper;
            const double violation = v - hi;
            out.maximum_violation = std::max(out.maximum_violation, violation);
            if (violation > out.worst_violation) { out.worst_violation = violation; out.worst_cell = i; }
        }
    }
    return out;
}

inline BoundednessDiagnostics audit_positive_field(
    const std::vector<double>& values,
    double tolerance = 0.0)
{
    return audit_boundedness(values, 0.0,
                             std::numeric_limits<double>::infinity(), tolerance);
}

inline BoundednessDiagnostics audit_boundedness(
    const Field<double, Location::CELL>& field,
    double lower,
    double upper,
    double tolerance = 0.0)
{
    if (field.dimension() != 1)
        throw std::invalid_argument("audit_boundedness: cell field must be scalar");
    if (field.size() == 0)
        throw std::invalid_argument("audit_boundedness: empty field");
    BoundednessDiagnostics out;
    const double lo = lower;
    const double hi = upper;
    if (!std::isfinite(lo) || (std::isfinite(hi) && lo > hi) ||
        (std::isinf(hi) && hi < 0.0) || !std::isfinite(tolerance) || tolerance < 0.0)
        throw std::invalid_argument("audit_boundedness: invalid bounds/tolerance");
    const double effective_lower = lo - tolerance;
    const double effective_upper = hi + tolerance;
    for (std::size_t i = 0; i < field.size(); ++i) {
        const double v = field(i);
        if (!std::isfinite(v)) { ++out.nonfinite; continue; }
        out.minimum = std::min(out.minimum, v);
        out.maximum = std::max(out.maximum, v);
        if (v < effective_lower) {
            ++out.below_lower;
            const double violation = effective_lower - v;
            out.minimum_violation = std::max(out.minimum_violation, violation);
            if (violation > out.worst_violation) { out.worst_violation = violation; out.worst_cell = i; }
        }
        if (std::isfinite(effective_upper) && v > effective_upper) {
            ++out.above_upper;
            const double violation = v - effective_upper;
            out.maximum_violation = std::max(out.maximum_violation, violation);
            if (violation > out.worst_violation) { out.worst_violation = violation; out.worst_cell = i; }
        }
    }
    return out;
}

inline BoundednessDiagnostics audit_positive_field(
    const Field<double, Location::CELL>& field,
    double tolerance = 0.0)
{
    return audit_boundedness(field, 0.0,
                             std::numeric_limits<double>::infinity(), tolerance);
}



/**
 * Generic cell-wise balance audit.
 *
 * The face-flux reconstruction is deliberately separated from the equation
 * balance: callers provide the independently reconstructed net face
 * contribution and integrated source/accumulation terms. This makes the same
 * contract usable for mass, momentum components, scalar transport, energy,
 * species and turbulence equations without embedding physics-specific logic.
 *
 * Per cell:
 *   residual[c] = face_balance[c] + source_integral[c] - accumulation[c]
 */
struct CellBalanceDiagnostics {
    double source_integral = 0.0;
    double accumulation = 0.0;
    double residual = 0.0;
    double normalized_residual = 0.0;
    double max_cell_residual = 0.0;
    double l1_cell_residual = 0.0;
    double l2_cell_residual = 0.0;
    std::size_t worst_cell = 0;
    std::size_t nonfinite_balance = 0;
    std::size_t nonfinite_source = 0;
    std::size_t nonfinite_accumulation = 0;

    bool finite() const noexcept {
        return nonfinite_balance == 0 && nonfinite_source == 0 &&
               nonfinite_accumulation == 0 &&
               std::isfinite(source_integral) &&
               std::isfinite(accumulation) &&
               std::isfinite(residual);
    }

    bool closed(double tolerance = 0.0) const noexcept {
        return finite() && std::abs(residual) <= tolerance;
    }
};

inline CellBalanceDiagnostics audit_cell_balance(
    const std::vector<double>& face_balance,
    const std::vector<double>& source_integral,
    const std::vector<double>* accumulation = nullptr)
{
    if (face_balance.empty() || source_integral.size() != face_balance.size())
        throw std::invalid_argument(
            "audit_cell_balance: balance/source sizes must match and be non-empty");
    if (accumulation && accumulation->size() != face_balance.size())
        throw std::invalid_argument(
            "audit_cell_balance: accumulation size does not match balance");

    CellBalanceDiagnostics out;
    double sum_sq = 0.0;
    double scale = 1.0;

    for (std::size_t c = 0; c < face_balance.size(); ++c) {
        const double balance = face_balance[c];
        const double source = source_integral[c];
        if (!std::isfinite(balance)) {
            ++out.nonfinite_balance;
            continue;
        }
        if (!std::isfinite(source)) {
            ++out.nonfinite_source;
            continue;
        }

        double acc = 0.0;
        if (accumulation) {
            acc = (*accumulation)[c];
            if (!std::isfinite(acc)) {
                ++out.nonfinite_accumulation;
                continue;
            }
        }

        const double residual = balance + source - acc;
        if (!std::isfinite(residual)) {
            ++out.nonfinite_balance;
            continue;
        }

        out.source_integral += source;
        out.accumulation += acc;
        out.residual += residual;
        out.l1_cell_residual += std::abs(residual);
        sum_sq += residual * residual;
        scale = std::max({scale, std::abs(balance),
                          std::abs(source), std::abs(acc)});
        if (std::abs(residual) > out.max_cell_residual) {
            out.max_cell_residual = std::abs(residual);
            out.worst_cell = c;
        }
    }

    out.l2_cell_residual = std::sqrt(sum_sq);
    out.normalized_residual = std::abs(out.residual) / scale;
    return out;
}

inline CellBalanceDiagnostics audit_cell_balance(
    const Mesh& mesh,
    const Field<double, Location::FACE>& face_flux,
    const Field<double, Location::CELL>& source_integral,
    const Field<double, Location::CELL>* accumulation = nullptr)
{
    if (source_integral.size() != mesh.n_cells() ||
        source_integral.dimension() != 1)
        throw std::invalid_argument(
            "audit_cell_balance: source field does not match mesh");

    const auto face_balance = reconstruct_cell_balance(mesh, face_flux);
    std::vector<double> source(source_integral.size());
    for (std::size_t c = 0; c < source.size(); ++c)
        source[c] = source_integral(c);

    std::vector<double> acc;
    if (accumulation) {
        if (accumulation->size() != mesh.n_cells() ||
            accumulation->dimension() != 1)
            throw std::invalid_argument(
                "audit_cell_balance: accumulation field does not match mesh");
        acc.resize(accumulation->size());
        for (std::size_t c = 0; c < acc.size(); ++c)
            acc[c] = (*accumulation)(c);
    }

    return audit_cell_balance(face_balance, source,
                              accumulation ? &acc : nullptr);
}

struct IntegratedBalanceDiagnostics {
    double boundary_flux = 0.0;
    double volume_source = 0.0;
    double accumulation = 0.0;
    double residual = 0.0;
    double normalized_residual = 0.0;
    double max_cell_residual = 0.0;
    double l1_cell_residual = 0.0;
    double l2_cell_residual = 0.0;
    std::size_t worst_cell = 0;
    std::size_t nonfinite_faces = 0;
    std::size_t nonfinite_source = 0;
    std::size_t nonfinite_accumulation = 0;

    bool finite() const noexcept {
        return nonfinite_faces == 0 && nonfinite_source == 0 &&
               nonfinite_accumulation == 0 &&
               std::isfinite(boundary_flux) &&
               std::isfinite(volume_source) &&
               std::isfinite(accumulation) &&
               std::isfinite(residual);
    }
    bool closed(double tolerance = 0.0) const noexcept {
        return finite() && std::abs(residual) <= tolerance;
    }
};

// Independent integral balance reconstruction:
//   sum(boundary face flux) + sum(volume source) - sum(accumulation).
// Internal faces are accounted for only once and cancel pairwise. This API
// deliberately consumes final fields/fluxes, not a matrix residual.
inline IntegratedBalanceDiagnostics audit_integrated_balance(
    const Mesh& mesh,
    const Field<double, Location::FACE>& face_flux,
    const Field<double, Location::CELL>& volume_source,
    const std::vector<double>& cell_volumes,
    const Field<double, Location::CELL>* accumulation = nullptr)
{
    const std::size_t nf = mesh.n_faces();
    const std::size_t nc = mesh.n_cells();
    if (face_flux.size() != nf || face_flux.dimension() != 1 ||
        volume_source.size() != nc || volume_source.dimension() != 1 ||
        cell_volumes.size() != nc)
        throw std::invalid_argument("audit_integrated_balance: dimensions do not match mesh");
    if (accumulation && (accumulation->size() != nc || accumulation->dimension() != 1))
        throw std::invalid_argument("audit_integrated_balance: invalid accumulation field");

    IntegratedBalanceDiagnostics out;
    const auto& own = mesh.ownership();
    std::vector<double> cell(nc, 0.0);
    for (std::size_t f = 0; f < nf; ++f) {
        const double flux = face_flux(f);
        if (!std::isfinite(flux)) { ++out.nonfinite_faces; continue; }
        const std::size_t o = own.owner(f);
        if (o >= nc) throw std::runtime_error("audit_integrated_balance: invalid owner");
        cell[o] += flux;
        const auto n = own.neighbour(f);
        if (n >= 0) {
            const std::size_t ni = static_cast<std::size_t>(n);
            if (ni >= nc) throw std::runtime_error("audit_integrated_balance: invalid neighbour");
            cell[ni] -= flux;
        } else {
            out.boundary_flux += flux;
        }
    }

    double sum_sq = 0.0;
    for (std::size_t c = 0; c < nc; ++c) {
        const double source = volume_source(c);
        if (!std::isfinite(source) || !std::isfinite(cell_volumes[c])) {
            ++out.nonfinite_source;
            continue;
        }
        const double source_integral = source * cell_volumes[c];
        out.volume_source += source_integral;
        double acc = 0.0;
        if (accumulation) {
            acc = (*accumulation)(c);
            if (!std::isfinite(acc)) { ++out.nonfinite_accumulation; continue; }
        }
        out.accumulation += acc;
        const double rc = cell[c] + source_integral - acc;
        out.residual += rc;
        const double a = std::abs(rc);
        out.l1_cell_residual += a;
        sum_sq += rc * rc;
        if (a > out.max_cell_residual) {
            out.max_cell_residual = a;
            out.worst_cell = c;
        }
    }
    out.l2_cell_residual = std::sqrt(sum_sq);
    const double scale = std::max({1.0, std::abs(out.boundary_flux),
                                   std::abs(out.volume_source),
                                   std::abs(out.accumulation)});
    out.normalized_residual = std::abs(out.residual) / scale;
    return out;
}

} // namespace cfdx::core
