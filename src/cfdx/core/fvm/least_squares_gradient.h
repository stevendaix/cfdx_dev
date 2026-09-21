#pragma once

#include "cfdx/core/field/field.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::core {

/**
 * Compute a cell-centred least-squares gradient from neighbouring cell values.
 *
 * For cell P, the gradient g minimises
 *
 *   sum_N w_N [phi_N - phi_P - g . (x_N - x_P)]^2
 *
 * with inverse-distance-squared weights. The implementation solves the
 * symmetric 3x3 normal equations with partial pivoting and detects rank
 * deficiency. A rank-deficient component is left at zero; this is required
 * for valid 2-D/1-D meshes embedded in 3-D.
 */
inline Vec3 least_squares_gradient(
    const Vec3& centre,
    double value,
    const std::vector<Vec3>& neighbour_centres,
    const std::vector<double>& neighbour_values)
{
    if (neighbour_centres.size() != neighbour_values.size())
        throw std::invalid_argument(
            "least_squares_gradient: neighbour dimensions do not match");

    double a[3][3] = {};
    double b[3] = {};

    for (std::size_t k = 0; k < neighbour_centres.size(); ++k) {
        const Vec3 d = neighbour_centres[k] - centre;
        const double r2 = d.mag2();
        if (!(r2 > 0.0) || !std::isfinite(r2))
            continue;

        const double w = 1.0 / r2;
        const double dv = neighbour_values[k] - value;

        a[0][0] += w * d.x * d.x;
        a[0][1] += w * d.x * d.y;
        a[0][2] += w * d.x * d.z;
        a[1][0] += w * d.y * d.x;
        a[1][1] += w * d.y * d.y;
        a[1][2] += w * d.y * d.z;
        a[2][0] += w * d.z * d.x;
        a[2][1] += w * d.z * d.y;
        a[2][2] += w * d.z * d.z;

        b[0] += w * d.x * dv;
        b[1] += w * d.y * dv;
        b[2] += w * d.z * dv;
    }

    double scale = 0.0;
    for (int i = 0; i < 3; ++i)
        scale = std::max(scale, std::abs(a[i][i]));

    if (!(scale > 0.0) || !std::isfinite(scale))
        return {};

    const double rank_tol =
        128.0 * std::numeric_limits<double>::epsilon() * scale;

    // Gaussian elimination with partial pivoting. Work on a copy so the
    // assembled normal matrix remains conceptually symmetric.
    double m[3][4] = {};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            m[i][j] = a[i][j];
        m[i][3] = b[i];
    }

    bool pivoted[3] = {false, false, false};
    int rank = 0;
    for (int col = 0; col < 3; ++col) {
        int pivot = col;
        double pivot_abs = std::abs(m[col][col]);
        for (int row = col + 1; row < 3; ++row) {
            const double candidate = std::abs(m[row][col]);
            if (candidate > pivot_abs) {
                pivot_abs = candidate;
                pivot = row;
            }
        }

        if (!(pivot_abs > rank_tol))
            continue;

        if (pivot != col) {
            for (int j = col; j < 4; ++j)
                std::swap(m[col][j], m[pivot][j]);
        }

        pivoted[col] = true;
        ++rank;

        for (int row = col + 1; row < 3; ++row) {
            const double factor = m[row][col] / m[col][col];
            if (factor == 0.0)
                continue;
            for (int j = col; j < 4; ++j)
                m[row][j] -= factor * m[col][j];
        }
    }

    if (rank == 0)
        return {};

    double x[3] = {};
    for (int i = 2; i >= 0; --i) {
        if (!pivoted[i])
            continue;
        double rhs = m[i][3];
        for (int j = i + 1; j < 3; ++j)
            rhs -= m[i][j] * x[j];
        if (std::abs(m[i][i]) > rank_tol)
            x[i] = rhs / m[i][i];
    }

    if (!std::isfinite(x[0]) || !std::isfinite(x[1]) || !std::isfinite(x[2]))
        throw std::runtime_error(
            "least_squares_gradient: non-finite gradient");

    return {x[0], x[1], x[2]};
}

} // namespace cfdx::core
