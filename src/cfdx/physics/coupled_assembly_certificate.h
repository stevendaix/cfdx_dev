#pragma once

#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

// Integrity certificate of the EXACT monolithic CSR system handed to Krylov.
// A valid certificate asserts finite coefficients, a well-formed pressure
// gauge (if configured), and consistent pressure Schur evidence arrays.
// The residual magnitudes below are Ax-b at the INITIAL GUESS, and are
// diagnostics only. They are NOT the conservative FV continuity residual.
struct CoupledAssemblyCertificate {
    bool passed = false;
    std::size_t unknowns = 0;
    std::size_t nnz_momentum = 0;
    std::size_t nnz_gradient = 0;
    std::size_t nnz_divergence = 0;
    std::size_t nnz_pressure = 0;
    std::size_t zero_rows = 0;
    double momentum_linear_residual_linf = std::numeric_limits<double>::infinity();
    double continuity_linear_residual_linf = std::numeric_limits<double>::infinity();
    double gauge_linear_residual = std::numeric_limits<double>::infinity();
};

struct CoupledAssemblyCertificateOptions {
    // For a pure-Neumann case, the local cell index of the imposed reference
    // pressure. For fixed-pressure boundaries, leave nullopt (no gauge row).
    std::optional<std::size_t> pressure_reference_cell = std::nullopt;
    double gauge_coefficient_tolerance = 1e-12;
    const std::vector<double>* pressure_schur_diagonal = nullptr;
    const std::vector<double>* pressure_schur_diagonal_raw = nullptr;
};

inline CoupledAssemblyCertificate certify_coupled_assembly_structure(
    const cfdx::core::SparseMatrix& A,
    const cfdx::core::Vector& rhs,
    const cfdx::core::Vector& initial_guess,
    std::size_t velocity_unknowns,
    std::size_t pressure_unknowns,
    const CoupledAssemblyCertificateOptions& controls = {})
{
    CoupledAssemblyCertificate result;
    if (velocity_unknowns == 0 || pressure_unknowns == 0 ||
        velocity_unknowns > std::numeric_limits<std::size_t>::max() - pressure_unknowns)
        return result;
    const std::size_t n = velocity_unknowns + pressure_unknowns;
    result.unknowns = n;
    if (A.n_rows() != n || A.n_cols() != n || rhs.size() != n ||
        initial_guess.size() != n || !rhs.is_valid() ||
        !initial_guess.is_valid() ||
        !std::isfinite(controls.gauge_coefficient_tolerance) ||
        controls.gauge_coefficient_tolerance < 0.0 ||
        (controls.pressure_reference_cell &&
         *controls.pressure_reference_cell >= pressure_unknowns))
        return result;
    const auto* diag = controls.pressure_schur_diagonal;
    const auto* raw = controls.pressure_schur_diagonal_raw;
    if ((diag && diag->size() != pressure_unknowns) ||
        (raw && raw->size() != pressure_unknowns)) return result;
    // Gauge replacement in the actual monolithic CSR sets its diagonal to
    // one. The adjusted Schur evidence must represent that same gauge row;
    // merely exempting it from the raw-vs-adjusted comparison was fail-open.
    if (diag && controls.pressure_reference_cell &&
        std::abs((*diag)[*controls.pressure_reference_cell] - 1.0) >
            controls.gauge_coefficient_tolerance)
        return result;
    for (std::size_t c = 0; c < pressure_unknowns; ++c) {
        if ((diag && !std::isfinite((*diag)[c])) ||
            (raw && !std::isfinite((*raw)[c]))) return result;
    const auto* raw = controls.pressure_schur_diagonal_raw;
    if ((diag && diag->size() != pressure_unknowns) ||
        (raw && raw->size() != pressure_unknowns)) return result;
    for (std::size_t c = 0; c < pressure_unknowns; ++c) {
        if ((diag && !std::isfinite((*diag)[c])) ||
            (raw && !std::isfinite((*raw)[c]))) return result;
        // The reference cell is the ONLY row whose diagnostic Schur
        // diagonal may differ after gauge replacement.
        if (diag && raw && (!controls.pressure_reference_cell ||
            c != *controls.pressure_reference_cell) &&
            std::abs((*diag)[c] - (*raw)[c]) >
                controls.gauge_coefficient_tolerance)
            return result;
    }
    const std::size_t gauge_row = controls.pressure_reference_cell
        ? velocity_unknowns + *controls.pressure_reference_cell : n;

    double momentum = 0.0, continuity = 0.0, gauge = 0.0;
    for (std::size_t row = 0; row < n; ++row) {
        double ax = 0.0;
        const std::size_t first = A.row_offsets_data()[row];
        const std::size_t last = A.row_offsets_data()[row+1];
        if (first == last) ++result.zero_rows;
        for (std::size_t k = first; k < last; ++k) {
            const std::size_t col = A.columns_data()[k];
            const double val = A.values_data()[k];
            if (col >= n || !std::isfinite(val)) return result;
            ax += val * initial_guess(col);
            if (row < velocity_unknowns) {
                if (col < velocity_unknowns) ++result.nnz_momentum;
                else ++result.nnz_gradient;
            } else {
                if (col < velocity_unknowns) ++result.nnz_divergence;
                else ++result.nnz_pressure;
            }
            if (row == gauge_row) {
                if (col == gauge_row) {
                    if (std::abs(val - 1.0) > controls.gauge_coefficient_tolerance)
                        return result;
                } else if (std::abs(val) > controls.gauge_coefficient_tolerance)
                    return result;
            }
            // The gauge pressure column is eliminated everywhere except
            // the diagonal of its replacement row.
            if (col == gauge_row && row != gauge_row &&
                std::abs(val) > controls.gauge_coefficient_tolerance)
                return result;
        }
        if (!std::isfinite(ax)) return result;
        const double defect = std::abs(ax - rhs(row));
        if (!std::isfinite(defect)) return result;
        if (row < velocity_unknowns) momentum = std::max(momentum, defect);
        else if (row == gauge_row) gauge = defect;
        else continuity = std::max(continuity, defect);
    }
    if (result.zero_rows != 0 ||
        (gauge_row < n && std::abs(A(gauge_row, gauge_row) - 1.0) >
                                controls.gauge_coefficient_tolerance))
        return result;
    result.momentum_linear_residual_linf = momentum;
    result.continuity_linear_residual_linf = continuity;
    result.gauge_linear_residual = gauge_row < n ? gauge : 0.0;
    result.passed = true;
    return result;
}

} // namespace cfdx::physics
