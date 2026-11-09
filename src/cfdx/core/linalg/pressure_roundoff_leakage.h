#pragma once

// Diagnostic for the *rounded* pressure matrix, never a convergence waiver.
// On a finite-volume Laplacian with exact Neumann faces, row sums should be
// zero. Rounding its diagonal in FP64 leaves a small artificial row sum; on
// poorly conditioned systems, multiplying by the large pressure correction
// creates a false boundary leakage comparable with the strict residual gate.
// This function only reports the magnitude on the ORIGINAL CSR and solution;
// it does not modify the matrix, selected Krylov solver, fluxes or tolerances.
// All CSR rows are independent and can be processed on CPU threads/GPU.

#include "sparse_matrix.h"
#include "vector.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace cfdx::core {

struct PressureRoundoffLeakageAudit {
    bool valid = false;
    std::size_t rows = 0;
    std::size_t pressure_anchor_rows = 0;
    std::size_t negligible_roundoff_rows = 0;
    std::size_t worst_row = 0;
    double spurious_leakage_relative_l2 =
        std::numeric_limits<double>::infinity();
    double spurious_leakage_max_abs =
        std::numeric_limits<double>::infinity();
    // True means only an algebraic artifact was measured; the FV physical
    // face flux (Rhie-Chow predictor/nonorthogonal) has NOT been certified.
    bool physical_face_flux_certified = false;
};

inline PressureRoundoffLeakageAudit diagnose_pressure_roundoff_leakage(
    const SparseMatrix& matrix, const Vector& rhs, const Vector& pressure,
    double relative_anchor_cutoff = 1e-11)
{
    PressureRoundoffLeakageAudit report;
    const auto n = matrix.n_rows();
    if (n == 0 || matrix.n_cols() != n || rhs.size() != n ||
        pressure.size() != n ||
        !(relative_anchor_cutoff > 0.0 && relative_anchor_cutoff < 1.0) ||
        !std::isfinite(relative_anchor_cutoff)) return report;
    const auto* offsets = matrix.row_offsets_data();
    const auto* columns = matrix.columns_data();
    const auto* values = matrix.values_data();
    long double bnorm2 = 0.0L;
    long double noise2 = 0.0L;
    long double maximum = 0.0L;
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(rhs(i)) || !std::isfinite(pressure(i)))
            return report;
        bnorm2 += static_cast<long double>(rhs(i)) * rhs(i);
        long double row_sum = 0.0L, row_abs = 0.0L;
        for (auto k = offsets[i]; k < offsets[i + 1]; ++k) {
            if (columns[k] >= n || !std::isfinite(values[k])) return report;
            row_sum += static_cast<long double>(values[k]);
            row_abs += std::abs(static_cast<long double>(values[k]));
        }
        if (!(row_abs > 0.0L)) return report;
        if (std::abs(row_sum) > relative_anchor_cutoff * row_abs) {
            ++report.pressure_anchor_rows;
            continue;  // Real pressure Dirichlet leakage is not roundoff.
        }
        ++report.negligible_roundoff_rows;
        const long double noise = row_sum * static_cast<long double>(pressure(i));
        const long double magnitude = std::abs(noise);
        noise2 += noise * noise;
        if (magnitude > maximum) {
            maximum = magnitude;
            report.worst_row = i;
        }
    }
    if (!(bnorm2 > 0.0L)) return report;
    report.rows = n;
    report.spurious_leakage_relative_l2 =
        static_cast<double>(std::sqrt(noise2 / bnorm2));
    report.spurious_leakage_max_abs = static_cast<double>(maximum);
    report.valid = std::isfinite(report.spurious_leakage_relative_l2) &&
                   std::isfinite(report.spurious_leakage_max_abs);
    return report;
}

}  // namespace cfdx::core
