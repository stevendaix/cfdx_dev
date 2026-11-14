#pragma once
// V42 evidence-only precision decomposition of the EXACT returned SST FV row.
// The existing convergence gates remain in sst_solver.h; this report may
// neither set converged nor substitute a compensated residual for a gate.
#include "cfdx/physics/finite_volume_transport.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {

struct SSTStallRowV42 {
    std::size_t cell = 0;
    double rhs = 0.0;
    double diagonal_term = 0.0;
    double off_diagonal_terms = 0.0;
    double native_signed = 0.0;
    double compensated_signed = 0.0;
    double row_l1 = 0.0;
    double local_field_scale = 0.0;
    double native_row_scaled = 0.0;
    double compensated_row_scaled = 0.0;
    double native_local_relative = 0.0;
    double compensated_local_relative = 0.0;
    double cancellation_ratio = 0.0;
};

inline SSTStallRowV42 sst_stall_row_v42(
    const ScalarEquation& equation, const cfdx::core::Vector& x,
    std::size_t row, double floor, bool verified_csr = false)
{
    // Overview verifies CSR once: revalidating every row is O(n * nnz).
    // Direct callers still validate by default.
    if (x.empty() || x.size() != equation.rhs.size() ||
        (!verified_csr && !equation.matrix.is_consistent()) ||
        equation.matrix.n_rows() != x.size() ||
        equation.matrix.n_cols() != x.size() || row >= x.size() ||
        !(floor > 0.0) || !std::isfinite(floor))
        throw std::invalid_argument("V42 SST row: invalid assembled FV dimensions");
    SSTStallRowV42 d;
    d.cell = row;
    d.rhs = equation.rhs(row);
    if (!std::isfinite(d.rhs) || !std::isfinite(x(row)))
        throw std::invalid_argument("V42 SST row: nonfinite right-hand side/state");
    long double sum = -static_cast<long double>(d.rhs);
    long double correction = 0.0L;
    long double row_l1 = 0.0L;
    long double scale = std::max(static_cast<long double>(floor),
                                 std::abs(static_cast<long double>(x(row))));
    long double largest_term = std::abs(static_cast<long double>(d.rhs));
    double native = -d.rhs;
    double native_row_l1 = 0.0;
    double native_scale = std::max(floor, std::abs(x(row)));
    long double diag = 0.0L, others = 0.0L;
    const auto* offsets = equation.matrix.row_offsets_data();
    const auto* columns = equation.matrix.columns_data();
    const auto* values = equation.matrix.values_data();
    for (std::size_t j=offsets[row]; j<offsets[row+1]; ++j) {
        const std::size_t col = columns[j];
        if (col >= x.size() || !std::isfinite(values[j]) ||
            !std::isfinite(x(col)))
            throw std::invalid_argument("V42 SST row: nonfinite CSR coefficient/state");
        const long double term = static_cast<long double>(values[j]) *
                                 static_cast<long double>(x(col));
        const long double next = sum + term;
        correction += std::abs(sum) >= std::abs(term)
            ? (sum-next)+term : (term-next)+sum;
        sum = next;
        if (col==row) diag+=term;
        else others+=term;
        row_l1+=std::abs(static_cast<long double>(values[j]));
        scale=std::max(scale,std::abs(static_cast<long double>(x(col))));
        largest_term=std::max(largest_term,std::abs(term));
        native+=values[j]*x(col);
        native_row_l1+=std::abs(values[j]);
        native_scale=std::max(native_scale,std::abs(x(col)));
    }
    if (!(row_l1>0.0L) || !(scale>0.0L) ||
        !(native_row_l1>0.0) || !(native_scale>0.0) ||
        !std::isfinite(sum+correction) || !std::isfinite(native))
        throw std::domain_error("V42 SST row: invalid FV row normalization");
    d.diagonal_term=static_cast<double>(diag);
    d.off_diagonal_terms=static_cast<double>(others);
    d.native_signed=native;
    d.compensated_signed=static_cast<double>(sum+correction);
    d.row_l1=static_cast<double>(row_l1);
    d.local_field_scale=static_cast<double>(scale);
    d.native_row_scaled=std::abs(native)/native_row_l1;
    d.compensated_row_scaled=static_cast<double>(std::abs(sum+correction)/row_l1);
    d.native_local_relative=d.native_row_scaled/native_scale;
    d.compensated_local_relative=static_cast<double>(
        std::abs(sum+correction)/(row_l1*scale));
    // Indicator of cancellation severity; not an error estimate.
    const long double residual_magnitude=std::abs(sum+correction);
    // An exactly balanced row has infinite mathematical cancellation;
    // represent it by the largest finite diagnostic, not a rejected state.
    d.cancellation_ratio=residual_magnitude==0.0L
        ? (largest_term>0.0L ? std::numeric_limits<double>::max() : 1.0)
        : static_cast<double>(std::min(
            static_cast<long double>(std::numeric_limits<double>::max()),
            largest_term/residual_magnitude));
    if (!std::isfinite(d.native_local_relative) ||
        !std::isfinite(d.compensated_local_relative) ||
        !std::isfinite(d.native_row_scaled) ||
        !std::isfinite(d.compensated_row_scaled) ||
        !std::isfinite(d.cancellation_ratio))
        throw std::domain_error("V42 SST row: diagnostic overflow");
    return d;
}

struct SSTStallOverviewV42 {
    SSTStallRowV42 native_row_scaled_worst;
    SSTStallRowV42 compensated_row_scaled_worst;
    SSTStallRowV42 native_local_worst;
    SSTStallRowV42 compensated_local_worst;
};

inline SSTStallOverviewV42 sst_stall_overview_v42(
    const ScalarEquation& equation, const cfdx::core::Vector& x, double floor)
{
    if (x.empty()) throw std::invalid_argument("V42 SST overview: empty field");
    if (!equation.matrix.is_consistent())
        throw std::invalid_argument("V42 SST overview: malformed CSR");
    SSTStallOverviewV42 result;
    for (std::size_t row=0; row<x.size(); ++row) {
        const auto d=sst_stall_row_v42(equation,x,row,floor,true);
        if (row==0 || d.native_row_scaled > result.native_row_scaled_worst.native_row_scaled)
            result.native_row_scaled_worst=d;
        if (row==0 || d.compensated_row_scaled > result.compensated_row_scaled_worst.compensated_row_scaled)
            result.compensated_row_scaled_worst=d;
        if (row==0 || d.native_local_relative > result.native_local_worst.native_local_relative)
            result.native_local_worst=d;
        if (row==0 || d.compensated_local_relative > result.compensated_local_worst.compensated_local_relative)
            result.compensated_local_worst=d;
    }
    return result;
}

} // namespace cfdx::physics
