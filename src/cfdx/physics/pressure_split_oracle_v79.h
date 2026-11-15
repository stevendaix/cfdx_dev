#pragma once
// V79: independent pressure split oracle. READ-ONLY; no solver state mutation,
// no production certificate, no pressure/flux acceptance. A small CSR residual
// alone never proves conservation of the full Rhie-Chow FV flux.
#include "cfdx/core/linalg/pressure_compensated_refinement.h"
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace cfdx::physics {

struct PressureSplitOracleV79 {
    cfdx::core::Vector high;
    cfdx::core::Vector low;
    explicit PressureSplitOracleV79(std::size_t size)
        : high(size,0.0),low(size,0.0) {}
};

inline PressureSplitOracleV79 read_pressure_split_oracle_v79(
    const std::string& path, std::size_t nrows, std::size_t nnz)
{
    if (path.empty() || nrows == 0 || nnz == 0)
        throw std::invalid_argument("V79: invalid split read parameters");
    std::ifstream input(path);
    if (!input) throw std::runtime_error("V79: cannot open pressure oracle");
    std::string magic;
    if (!std::getline(input,magic) || magic!="CFDX_PRESSURE_SPLIT_V79")
        throw std::runtime_error("V79: invalid split format");
    std::size_t file_rows=0,file_nnz=0;
    if (!(input >> file_rows >> file_nnz) || file_rows != nrows || file_nnz != nnz)
        throw std::runtime_error("V79: matrix/split dimensions or nnz mismatch");
    PressureSplitOracleV79 out(nrows);
    for (std::size_t i=0;i<nrows;++i) {
        if (!(input >> out.high(i) >> out.low(i)) ||
            !std::isfinite(out.high(i)) || !std::isfinite(out.low(i)))
            throw std::runtime_error("V79: truncated or nonfinite split data");
    }
    std::string extra;
    if (input >> extra)
        throw std::runtime_error("V79: extra split data");
    return out;
}

struct PressureSplitEvidenceV79 {
    double csr_relative=std::numeric_limits<double>::infinity();
    double double_only_relative=std::numeric_limits<double>::infinity();
    bool csr_closed=false;
    bool production_commit_authorized=false;
};

inline PressureSplitEvidenceV79 audit_pressure_split_oracle_v79(
    const cfdx::core::SparseMatrix& a, const cfdx::core::Vector& b,
    const PressureSplitOracleV79& trial, double original_tolerance)
{
    if (!(original_tolerance > 0.) || !std::isfinite(original_tolerance))
        throw std::invalid_argument("V79: invalid unmodified pressure tolerance");
    PressureSplitEvidenceV79 result;
    cfdx::core::Vector zero(trial.low.size(),0.0);
    result.csr_relative = cfdx::core::compensated_pressure_residual(
        a,b,trial.high,trial.low).second;
    result.double_only_relative = cfdx::core::compensated_pressure_residual(
        a,b,trial.high,zero).second;
    result.csr_closed = result.csr_relative <= original_tolerance;
    // The face proof is a SEPARATE physics calculation in the caller.
    result.production_commit_authorized = false;
    return result;
}

} // namespace cfdx::physics