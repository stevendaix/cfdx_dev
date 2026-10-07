#pragma once

#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::core {

/// Matrix-quality evidence used by the numerical robustness layer.
///
/// These diagnostics are intentionally descriptive: they do not modify the
/// matrix, relax tolerances, or select a different solver.
struct MatrixDiagnostics {
    std::size_t rows = 0;
    std::size_t columns = 0;
    std::size_t nnz = 0;

    bool finite = true;
    bool structurally_consistent = true;

    std::size_t empty_rows = 0;
    std::size_t empty_columns = 0;
    std::size_t missing_diagonal = 0;
    std::size_t near_zero_diagonal = 0;
    std::size_t isolated_dofs = 0;
    std::size_t duplicate_entries = 0;

    double max_abs_value = 0.0;
    double min_nonzero_abs_value = std::numeric_limits<double>::infinity();
    double max_row_sum = 0.0;
    double min_nonzero_row_sum = std::numeric_limits<double>::infinity();
    double max_column_sum = 0.0;
    double min_nonzero_column_sum = std::numeric_limits<double>::infinity();

    double max_abs_diagonal = 0.0;
    double min_nonzero_abs_diagonal = std::numeric_limits<double>::infinity();

    // These are condition indicators, not claims of an exact condition number.
    double diagonal_dynamic_range = 1.0;
    double row_sum_dynamic_range = 1.0;
    double column_sum_dynamic_range = 1.0;

    // For square matrices, this is the minimum Gershgorin lower margin
    // |a_ii| - sum_{j != i}|a_ij|. It is useful evidence for diagonal
    // dominance but is not an eigenvalue estimate for general matrices.
    double minimum_gershgorin_margin = std::numeric_limits<double>::infinity();

    std::size_t connected_components = 0;
};

enum class MatrixScaling {
    None,
    Row,
    Column,
    RowColumn,
    SymmetricDiagonal
};

struct MatrixScalingResult {
    SparseMatrix matrix;
    std::vector<double> row_scale;
    std::vector<double> column_scale;
    std::size_t zero_rows = 0;
    std::size_t zero_columns = 0;
    bool applied = false;
};

enum class SolverFailureClass {
    None,
    NonFinite,
    MatrixPathology,
    IncompatibleRhs,
    Divergence,
    Stagnation,
    MaxIterations
};

inline const char* to_string(SolverFailureClass failure) {
    switch (failure) {
        case SolverFailureClass::None: return "none";
        case SolverFailureClass::NonFinite: return "non_finite";
        case SolverFailureClass::MatrixPathology: return "matrix_pathology";
        case SolverFailureClass::IncompatibleRhs: return "incompatible_rhs";
        case SolverFailureClass::Divergence: return "divergence";
        case SolverFailureClass::Stagnation: return "stagnation";
        case SolverFailureClass::MaxIterations: return "max_iterations";
    }
    return "unknown";
}

enum class MatrixPathology {
    None,
    NonFiniteCoefficient,
    EmptyRow,
    EmptyColumn,
    MissingDiagonal,
    NearZeroDiagonal,
    IsolatedDof,
    DisconnectedMatrix
};

inline const char* to_string(MatrixPathology pathology) {
    switch (pathology) {
        case MatrixPathology::None: return "none";
        case MatrixPathology::NonFiniteCoefficient: return "non_finite_coefficient";
        case MatrixPathology::EmptyRow: return "empty_row";
        case MatrixPathology::EmptyColumn: return "empty_column";
        case MatrixPathology::MissingDiagonal: return "missing_diagonal";
        case MatrixPathology::NearZeroDiagonal: return "near_zero_diagonal";
        case MatrixPathology::IsolatedDof: return "isolated_dof";
        case MatrixPathology::DisconnectedMatrix: return "disconnected_matrix";
    }
    return "unknown";
}

inline MatrixDiagnostics diagnose_matrix(
    const SparseMatrix& matrix,
    double diagonal_relative_tolerance = 1e-14) {
    if (!(diagonal_relative_tolerance >= 0.0) ||
        !std::isfinite(diagonal_relative_tolerance))
        throw std::invalid_argument("diagonal relative tolerance must be finite and non-negative");

    MatrixDiagnostics d;
    d.rows = matrix.n_rows();
    d.columns = matrix.n_cols();
    d.nnz = matrix.nnz();
    d.structurally_consistent = matrix.is_consistent();

    if (!d.structurally_consistent) {
        d.finite = false;
        return d;
    }

    std::vector<double> row_sum(d.rows, 0.0);
    std::vector<double> column_sum(d.columns, 0.0);
    std::vector<bool> row_offdiag(d.rows, false);
    std::vector<bool> column_offdiag(d.columns, false);

    const auto* row = matrix.row_offsets_data();
    const auto* col = matrix.columns_data();
    const auto* val = matrix.values_data();

    for (std::size_t i = 0; i < d.rows; ++i) {
        const auto begin = row[i];
        const auto end = row[i + 1];
        if (begin == end) ++d.empty_rows;

        bool diagonal_present = false;
        double diagonal = 0.0;
        std::size_t previous_column = std::numeric_limits<std::size_t>::max();

        for (std::size_t k = begin; k < end; ++k) {
            const std::size_t j = col[k];
            const double a = val[k];
            if (!std::isfinite(a)) {
                d.finite = false;
                continue;
            }

            const double magnitude = std::abs(a);
            d.max_abs_value = std::max(d.max_abs_value, magnitude);
            if (magnitude > 0.0)
                d.min_nonzero_abs_value = std::min(d.min_nonzero_abs_value, magnitude);
            row_sum[i] += magnitude;
            column_sum[j] += magnitude;

            if (j == i) {
                diagonal_present = true;
                diagonal += a;
            } else if (magnitude > 0.0) {
                row_offdiag[i] = true;
                column_offdiag[j] = true;
            }

            if (j == previous_column) ++d.duplicate_entries;
            previous_column = j;
        }

        if (!diagonal_present && i < d.columns) ++d.missing_diagonal;
        const double diagonal_abs = std::abs(diagonal);
        if (diagonal_abs > 0.0) {
            d.max_abs_diagonal = std::max(d.max_abs_diagonal, diagonal_abs);
            d.min_nonzero_abs_diagonal =
                std::min(d.min_nonzero_abs_diagonal, diagonal_abs);
        }

        if (i < d.columns && diagonal_abs <=
                diagonal_relative_tolerance * std::max(1.0, d.max_abs_value))
            ++d.near_zero_diagonal;

        if (i < d.columns) {
            double offdiag_sum = 0.0;
            for (std::size_t k = begin; k < end; ++k)
                if (col[k] != i) offdiag_sum += std::abs(val[k]);
            d.minimum_gershgorin_margin =
                std::min(d.minimum_gershgorin_margin, diagonal_abs - offdiag_sum);
        }
    }

    for (std::size_t j = 0; j < d.columns; ++j)
        if (column_sum[j] == 0.0) ++d.empty_columns;

    for (std::size_t i = 0; i < std::min(d.rows, d.columns); ++i)
        if (!row_offdiag[i] && !column_offdiag[i]) ++d.isolated_dofs;

    for (const double value : row_sum) {
        if (value > 0.0) {
            d.max_row_sum = std::max(d.max_row_sum, value);
            d.min_nonzero_row_sum = std::min(d.min_nonzero_row_sum, value);
        }
    }
    for (const double value : column_sum) {
        if (value > 0.0) {
            d.max_column_sum = std::max(d.max_column_sum, value);
            d.min_nonzero_column_sum = std::min(d.min_nonzero_column_sum, value);
        }
    }

    if (d.min_nonzero_abs_diagonal < std::numeric_limits<double>::infinity())
        d.diagonal_dynamic_range =
            d.max_abs_diagonal / d.min_nonzero_abs_diagonal;
    if (d.min_nonzero_row_sum < std::numeric_limits<double>::infinity())
        d.row_sum_dynamic_range =
            d.max_row_sum / d.min_nonzero_row_sum;
    if (d.min_nonzero_column_sum < std::numeric_limits<double>::infinity())
        d.column_sum_dynamic_range =
            d.max_column_sum / d.min_nonzero_column_sum;

    if (d.rows != d.columns)
        d.minimum_gershgorin_margin = std::numeric_limits<double>::quiet_NaN();

    // Connected components of the undirected sparsity graph. For a rectangular
    // matrix, rows and columns are mapped into one bipartite graph.
    const std::size_t nodes = d.rows + d.columns;
    std::vector<std::vector<std::size_t>> graph(nodes);
    for (std::size_t i = 0; i < d.rows; ++i) {
        for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
            const std::size_t j = col[k];
            if (std::abs(val[k]) == 0.0) continue;
            graph[i].push_back(d.rows + j);
            graph[d.rows + j].push_back(i);
        }
    }
    std::vector<bool> visited(nodes, false);
    for (std::size_t root = 0; root < nodes; ++root) {
        if (visited[root]) continue;
        ++d.connected_components;
        std::vector<std::size_t> stack{root};
        visited[root] = true;
        while (!stack.empty()) {
            const auto node = stack.back();
            stack.pop_back();
            for (const auto neighbour : graph[node]) {
                if (!visited[neighbour]) {
                    visited[neighbour] = true;
                    stack.push_back(neighbour);
                }
            }
        }
    }

    return d;
}


/// Convert measured CSR diagnostics into the dispatcher-facing matrix profile.
///
/// This is descriptive classification only. It does not modify the matrix and
/// none of the thresholds below are solver convergence tolerances.
inline MatrixCharacteristics measure_matrix_characteristics(
    const SparseMatrix& matrix,
    bool saddle_point = false,
    double symmetry_relative_tolerance = 1e-12) {
    if (!(symmetry_relative_tolerance >= 0.0) ||
        !std::isfinite(symmetry_relative_tolerance))
        throw std::invalid_argument(
            "matrix symmetry relative tolerance must be finite and non-negative");

    const auto d = diagnose_matrix(matrix);
    MatrixCharacteristics c;
    c.equations = d.rows;
    c.average_nnz_per_row =
        d.rows == 0 ? 0.0 : static_cast<double>(d.nnz) / static_cast<double>(d.rows);
    c.coefficient_range =
        d.min_nonzero_abs_value == std::numeric_limits<double>::infinity()
            ? 1.0
            : d.max_abs_value / d.min_nonzero_abs_value;
    c.square = d.rows == d.columns;
    c.diagonally_dominant =
        c.square && std::isfinite(d.minimum_gershgorin_margin) &&
        d.minimum_gershgorin_margin >=
            -1e-14 * std::max(1.0, d.max_abs_diagonal);
    c.strongly_scaled = std::isfinite(c.coefficient_range) &&
                        c.coefficient_range >= 1.0e8;
    c.anisotropic = std::isfinite(d.diagonal_dynamic_range) &&
                    d.diagonal_dynamic_range >= 1.0e4;
    c.saddle_point = saddle_point;

    if (c.square && matrix.n_rows() > 0) {
        const auto* row = matrix.row_offsets_data();
        const auto* col = matrix.columns_data();
        const auto* val = matrix.values_data();
        const double scale = std::max(1.0, d.max_abs_value);
        c.numerically_symmetric = true;
        for (std::size_t i = 0; i < matrix.n_rows() && c.numerically_symmetric; ++i) {
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                const std::size_t j = col[k];
                const double aji = matrix(j, i);
                if (std::abs(val[k] - aji) > symmetry_relative_tolerance * scale) {
                    c.numerically_symmetric = false;
                    break;
                }
            }
        }
    }

    return c;
}

inline std::vector<MatrixPathology> classify_matrix_pathologies(
    const MatrixDiagnostics& d) {
    std::vector<MatrixPathology> result;
    if (!d.finite || !d.structurally_consistent)
        result.push_back(MatrixPathology::NonFiniteCoefficient);
    if (d.empty_rows) result.push_back(MatrixPathology::EmptyRow);
    if (d.empty_columns) result.push_back(MatrixPathology::EmptyColumn);
    if (d.missing_diagonal) result.push_back(MatrixPathology::MissingDiagonal);
    if (d.near_zero_diagonal) result.push_back(MatrixPathology::NearZeroDiagonal);
    if (d.isolated_dofs) result.push_back(MatrixPathology::IsolatedDof);
    if (d.connected_components > 1)
        result.push_back(MatrixPathology::DisconnectedMatrix);
    if (result.empty()) result.push_back(MatrixPathology::None);
    return result;
}

inline SolverFailureClass classify_solver_failure(
    const SolverResult& result,
    const MatrixDiagnostics& diagnostics,
    bool rhs_compatible = true,
    double stagnation_ratio = 0.95) {
    if (!(stagnation_ratio >= 0.0 && stagnation_ratio <= 1.0) ||
        !std::isfinite(stagnation_ratio))
        throw std::invalid_argument("stagnation ratio must be finite and in [0,1]");

    if (result.status == SolverStatus::CONVERGED)
        return SolverFailureClass::None;
    if (!std::isfinite(result.residual) ||
        !std::isfinite(result.residual_relative))
        return SolverFailureClass::NonFinite;

    const auto pathologies = classify_matrix_pathologies(diagnostics);
    if (!(pathologies.size() == 1 && pathologies.front() == MatrixPathology::None))
        return SolverFailureClass::MatrixPathology;
    if (!rhs_compatible)
        return SolverFailureClass::IncompatibleRhs;
    if (result.status == SolverStatus::DIVERGED)
        return SolverFailureClass::Divergence;

    if (result.status == SolverStatus::MAX_ITER_REACHED &&
        std::isfinite(result.min_true_residual) &&
        std::isfinite(result.max_true_residual) &&
        result.max_true_residual > 0.0 &&
        result.min_true_residual / result.max_true_residual >= stagnation_ratio)
        return SolverFailureClass::Stagnation;

    return SolverFailureClass::MaxIterations;
}

inline MatrixScalingResult scale_matrix(
    const SparseMatrix& matrix,
    MatrixScaling mode,
    double zero_row_tolerance = 0.0) {
    if (!(zero_row_tolerance >= 0.0) || !std::isfinite(zero_row_tolerance))
        throw std::invalid_argument("zero-row tolerance must be finite and non-negative");

    MatrixScalingResult result;
    result.matrix = SparseMatrix(matrix.n_rows(), matrix.n_cols());
    result.row_scale.assign(matrix.n_rows(), 1.0);
    result.column_scale.assign(matrix.n_cols(), 1.0);

    if (mode == MatrixScaling::None) {
        for (std::size_t i = 0; i < matrix.n_rows(); ++i)
            for (std::size_t k = matrix.row_offsets_data()[i];
                 k < matrix.row_offsets_data()[i + 1]; ++k)
                result.matrix.push_back(i, matrix.columns_data()[k],
                                       matrix.values_data()[k]);
        result.matrix.finalize();
        return result;
    }

    std::vector<double> row_max(matrix.n_rows(), 0.0);
    std::vector<double> col_max(matrix.n_cols(), 0.0);
    for (std::size_t i = 0; i < matrix.n_rows(); ++i) {
        for (std::size_t k = matrix.row_offsets_data()[i];
             k < matrix.row_offsets_data()[i + 1]; ++k) {
            const std::size_t j = matrix.columns_data()[k];
            const double magnitude = std::abs(matrix.values_data()[k]);
            row_max[i] = std::max(row_max[i], magnitude);
            col_max[j] = std::max(col_max[j], magnitude);
        }
    }

    if (mode == MatrixScaling::Row || mode == MatrixScaling::RowColumn) {
        for (std::size_t i = 0; i < row_max.size(); ++i) {
            if (row_max[i] <= zero_row_tolerance) {
                ++result.zero_rows;
            } else {
                result.row_scale[i] = 1.0 / row_max[i];
            }
        }
    }

    if (mode == MatrixScaling::Column || mode == MatrixScaling::RowColumn) {
        for (std::size_t j = 0; j < col_max.size(); ++j) {
            if (col_max[j] <= zero_row_tolerance) {
                ++result.zero_columns;
            } else {
                result.column_scale[j] = 1.0 / col_max[j];
            }
        }
    }

    if (mode == MatrixScaling::SymmetricDiagonal) {
        if (matrix.n_rows() != matrix.n_cols())
            throw std::invalid_argument("symmetric diagonal scaling requires a square matrix");
        for (std::size_t i = 0; i < matrix.n_rows(); ++i) {
            const double diagonal = std::abs(matrix(i, i));
            if (diagonal <= zero_row_tolerance) {
                ++result.zero_rows;
                ++result.zero_columns;
            } else {
                const double scale = 1.0 / std::sqrt(diagonal);
                result.row_scale[i] = scale;
                result.column_scale[i] = scale;
            }
        }
    }

    for (std::size_t i = 0; i < matrix.n_rows(); ++i) {
        for (std::size_t k = matrix.row_offsets_data()[i];
             k < matrix.row_offsets_data()[i + 1]; ++k) {
            const std::size_t j = matrix.columns_data()[k];
            result.matrix.push_back(
                i, j,
                matrix.values_data()[k] * result.row_scale[i] *
                    result.column_scale[j]);
        }
    }
    result.matrix.finalize();
    result.applied = mode != MatrixScaling::None;
    return result;
}

inline double rhs_compatibility_norm(
    const NullSpaceProjector& null_space,
    const Vector& rhs) {
    return null_space.component_norm(rhs);
}

} // namespace cfdx::core
