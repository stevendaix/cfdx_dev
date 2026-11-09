#pragma once

// Read-only, original-operator pressure balance in extended precision.
//
// This kernel is deliberately separate from the CFD acceptance path. A small
// algebraic balance does not certify a double-precision pressure field, its
// Rhie-Chow face fluxes or a converged FV calculation. The high/low interface
// provides the data contract needed by a future compensated FV face flux path.
// Each CSR row can be evaluated independently (CPU threading/GPU reduction).

#include "sparse_matrix.h"
#include "vector.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::core {

struct CompensatedPressureBalance {
    bool valid = false;
    bool symmetric_m_matrix = false;
    double original_relative_l2 = std::numeric_limits<double>::infinity();
    double algebraic_face_relative_l2 = std::numeric_limits<double>::infinity();
    double max_row_face_identity = std::numeric_limits<double>::infinity();
    double global_internal_pair_defect = std::numeric_limits<double>::infinity();
    // Never set by the algebraic replay: this requires FV face reconstruction.
    bool physical_face_flux_certified = false;
};

// Error-free double + double accumulation (Knuth's TwoSum). The returned pair
// represents the real sum of the INPUT doubles up to under/overflow limits.
inline bool compensated_add(double& high, double& low, double increment) {
    if (!std::isfinite(high) || !std::isfinite(low) ||
        !std::isfinite(increment)) return false;
    const double s = high + increment;
    const double bb = s - high;
    const double error = (high - (s - bb)) + (increment - bb);
    const double t = low + error;
    const double updated = s + t;
    const double correction = t - (updated - s);
    if (!std::isfinite(updated) || !std::isfinite(correction)) return false;
    high = updated;
    low = correction;
    return true;
}

// The row-wise pressure operator is evaluated before any FP64 rounding of
// (high+low). Pure Neumann is supported as an algebraic system, but requires a
// separate gauge/compatibility contract for any physics promotion.
inline CompensatedPressureBalance audit_compensated_pressure_balance(
    const SparseMatrix& matrix, const Vector& rhs,
    const Vector& high, const Vector& low,
    double symmetry_relative_tolerance = 1e-12)
{
    CompensatedPressureBalance report;
    const std::size_t n = matrix.n_rows();
    if (!n || matrix.n_cols() != n || rhs.size() != n ||
        high.size() != n || low.size() != n ||
        !(symmetry_relative_tolerance >= 0.0) ||
        !std::isfinite(symmetry_relative_tolerance)) return report;
    const auto* offsets = matrix.row_offsets_data();
    const auto* columns = matrix.columns_data();
    const auto* values = matrix.values_data();
    std::vector<long double> direct(n), graph(n, 0.0L);
    long double rhs_norm2 = 0.0L, norm2 = 0.0L;
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(high(i)) || !std::isfinite(low(i)) ||
            !std::isfinite(rhs(i))) return report;
        rhs_norm2 += static_cast<long double>(rhs(i)) * rhs(i);
    }
    if (rhs_norm2 <= 0.0L) return report;  // no scale to certify
    bool symmetric = true;
    long double boundary_sum = 0.0L;
    for (std::size_t i = 0; i < n; ++i) {
        long double off_high = 0.0L, off_low = 0.0L, diagonal = 0.0L;
        long double off_diagonal_sum = 0.0L;
        std::size_t diagonal_count = 0;
        for (std::uint32_t k = offsets[i]; k < offsets[i + 1]; ++k) {
            const std::size_t j = columns[k];
            const double a = values[k];
            if (j >= n || !std::isfinite(a)) return report;
            if (j == i) {
                ++diagonal_count;
                diagonal += a;
            } else {
                off_diagonal_sum += a;
                off_high += static_cast<long double>(a) *
                    (static_cast<long double>(high(j)) - high(i));
                off_low += static_cast<long double>(a) *
                    (static_cast<long double>(low(j)) - low(i));
                if (!(a <= 0.0)) symmetric = false;
                // Sparse row lookup is only used for the diagnostic of
                // reciprocity; do not assume a sorted-CSR binary search.
                long double transpose = 0.0L;
                bool found = false;
                for (std::uint32_t t = offsets[j]; t < offsets[j + 1]; ++t) {
                    if (columns[t] == i) {
                        transpose += values[t];
                        found = true;
                    }
                }
                if (!found || std::abs(static_cast<long double>(a) - transpose) >
                    symmetry_relative_tolerance *
                    std::max({1.0e-300L, std::abs(static_cast<long double>(a)),
                              std::abs(transpose)})) symmetric = false;
                if (i < j) {
                    const long double face = -static_cast<long double>(a) * (
                        (static_cast<long double>(high(i)) - high(j)) +
                        (static_cast<long double>(low(i)) - low(j)));
                    graph[i] += face;
                    graph[j] -= face;
                }
            }
        }
        if (diagonal_count != 1 || !(diagonal > 0.0L)) symmetric = false;
        // Boundary leakage corresponds to row sums, not to a made-up face.
        const long double row_sum = diagonal + off_diagonal_sum;
        const long double leakage = row_sum * high(i) + row_sum * low(i);
        boundary_sum += leakage;
        graph[i] += leakage;
        direct[i] = (off_high + off_low + leakage) -
                    static_cast<long double>(rhs(i));
        norm2 += direct[i] * direct[i];
    }
    report.valid = true;
    report.original_relative_l2 =
        static_cast<double>(std::sqrt(norm2 / rhs_norm2));
    report.symmetric_m_matrix = symmetric;
    if (symmetric) {
        long double face_norm2 = 0.0L, worst = 0.0L;
        for (std::size_t i = 0; i < n; ++i) {
            const long double defect = graph[i] - static_cast<long double>(rhs(i));
            face_norm2 += defect * defect;
            worst = std::max(worst, std::abs(graph[i] -
                (direct[i] + static_cast<long double>(rhs(i)))));
        }
        report.algebraic_face_relative_l2 =
            static_cast<double>(std::sqrt(face_norm2 / rhs_norm2));
        report.max_row_face_identity = static_cast<double>(worst);
        long double total = -boundary_sum;
        for (const long double value : graph) total += value;
        report.global_internal_pair_defect = static_cast<double>(std::abs(total));
    }
    return report;
}

} // namespace cfdx::core
