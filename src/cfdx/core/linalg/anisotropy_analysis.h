#pragma once

#include "cfdx/core/linalg/sparse_matrix.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::core {

// Measurements on the *assembled equation*, not an inferred mesh aspect ratio.
// Optional geometric inputs are kept separate: one mesh can induce very
// different pressure, momentum, and thermal operators.
struct AnisotropyGeometry {
    std::vector<std::array<double, 3>> cell_centers;
    std::vector<double> cell_aspect_ratios;
    std::vector<double> face_nonorthogonality_degrees;
    std::vector<double> face_skewness;
};

struct AnisotropyReport {
    std::size_t rows = 0;
    std::size_t nonzeros = 0;
    std::size_t anisotropic_rows = 0;
    std::size_t strong_directed_edges = 0;
    std::size_t branching_rows = 0;
    std::size_t asymmetric_strong_edges = 0;
    std::size_t nonpositive_diagonals = 0;
    std::size_t geometrically_directional_rows = 0;
    std::size_t geometric_dimension = 0;
    double worst_offdiagonal_contrast = 1.0;
    double mean_strong_coupling_fraction = 0.0;
    double max_cell_aspect_ratio = 0.0;
    double max_face_nonorthogonality_degrees = 0.0;
    double max_face_skewness = 0.0;
    double mean_neighbor_direction_alignment = 0.0;
    bool has_geometry = false;
    bool symmetric_strong_graph = true;
    std::vector<std::size_t> worst_rows;

    // JSON has no NaN or infinity. All reported numbers are finite.
    std::string json(const std::string& equation = "unspecified") const {
        std::ostringstream out;
        out << std::setprecision(17);
        out << "{\"equation\":\"";
        for (char ch : equation) {
            if (ch == '"' || ch == '\\') out << '\\';
            if (static_cast<unsigned char>(ch) < 32) out << '_';
            else out << ch;
        }
        out << "\",\"rows\":" << rows
            << ",\"nonzeros\":" << nonzeros
            << ",\"anisotropic_rows\":" << anisotropic_rows
            << ",\"strong_directed_edges\":" << strong_directed_edges
            << ",\"branching_rows\":" << branching_rows
            << ",\"asymmetric_strong_edges\":" << asymmetric_strong_edges
            << ",\"nonpositive_diagonals\":" << nonpositive_diagonals
            << ",\"worst_offdiagonal_contrast\":" << worst_offdiagonal_contrast
            << ",\"mean_strong_coupling_fraction\":" << mean_strong_coupling_fraction
            << ",\"has_geometry\":" << (has_geometry ? "true" : "false")
            << ",\"geometric_dimension\":" << geometric_dimension
            << ",\"geometrically_directional_rows\":" << geometrically_directional_rows
            << ",\"max_cell_aspect_ratio\":" << max_cell_aspect_ratio
            << ",\"max_face_nonorthogonality_degrees\":" << max_face_nonorthogonality_degrees
            << ",\"max_face_skewness\":" << max_face_skewness
            << ",\"mean_neighbor_direction_alignment\":" << mean_neighbor_direction_alignment
            << ",\"symmetric_strong_graph\":"
            << (symmetric_strong_graph ? "true" : "false")
            << ",\"worst_rows\":[";
        for (std::size_t i = 0; i < worst_rows.size(); ++i) {
            if (i) out << ',';
            out << worst_rows[i];
        }
        out << "]}";
        return out.str();
    }
};

namespace anisotropy_detail {
inline std::array<double, 3> principal_direction(
    const std::array<double, 9>& covariance, std::size_t dimension,
    double& eigen_ratio) {
    // Symmetric Jacobi rotations (2x2 or 3x3), deterministic and stable for
    // weak couplings. Eigenvectors are sorted by decreasing eigenvalue.
    double a[3][3] = {}, v[3][3] = {};
    for (std::size_t i = 0; i < dimension; ++i) {
        v[i][i] = 1.0;
        for (std::size_t j = 0; j < dimension; ++j)
            a[i][j] = covariance[3 * i + j];
    }
    for (int sweep = 0; sweep < 24; ++sweep) {
        std::size_t p = 0, q = 1;
        double maximum = 0.0;
        for (std::size_t i = 0; i < dimension; ++i)
            for (std::size_t j = i + 1; j < dimension; ++j)
                if (std::abs(a[i][j]) > maximum) {
                    p = i; q = j; maximum = std::abs(a[i][j]);
                }
        if (maximum <= 1e-14 * std::max(1.0, std::abs(a[0][0]))) break;
        const double tau = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
        const double t = (tau >= 0.0 ? 1.0 : -1.0) /
            (std::abs(tau) + std::sqrt(1.0 + tau * tau));
        const double c = 1.0 / std::sqrt(1.0 + t * t), s = t * c;
        const double app = a[p][p], aqq = a[q][q], apq = a[p][q];
        a[p][p] = app - t * apq;
        a[q][q] = aqq + t * apq;
        a[p][q] = a[q][p] = 0.0;
        for (std::size_t k = 0; k < dimension; ++k) {
            if (k != p && k != q) {
                const double akp = a[k][p], akq = a[k][q];
                a[k][p] = a[p][k] = c * akp - s * akq;
                a[k][q] = a[q][k] = s * akp + c * akq;
            }
            const double vkp = v[k][p], vkq = v[k][q];
            v[k][p] = c * vkp - s * vkq;
            v[k][q] = s * vkp + c * vkq;
        }
    }
    std::size_t best = 0;
    double lowest = a[0][0];
    for (std::size_t k = 1; k < dimension; ++k) {
        if (a[k][k] > a[best][best]) best = k;
        lowest = std::min(lowest, a[k][k]);
    }
    eigen_ratio = a[best][best] / std::max(1e-30, lowest);
    return {v[0][best], v[1][best], v[2][best]};
}
} // namespace anisotropy_detail

inline AnisotropyReport analyze_anisotropy(
    const SparseMatrix& A, const AnisotropyGeometry* geometry = nullptr,
    double strength_fraction = 0.25, double directional_contrast = 32.0) {
    if (!A.is_consistent() || A.n_rows() == 0 || A.n_rows() != A.n_cols() ||
        !(strength_fraction > 0.0 && strength_fraction < 1.0) ||
        !std::isfinite(strength_fraction) ||
        !(directional_contrast > 1.0) || !std::isfinite(directional_contrast))
        throw std::invalid_argument("analyze_anisotropy: invalid CSR or thresholds");
    const std::size_t n = A.n_rows();
    if (geometry && !geometry->cell_centers.empty() &&
        geometry->cell_centers.size() != n)
        throw std::invalid_argument("analyze_anisotropy: cell-center size mismatch");
    if (geometry && !geometry->cell_aspect_ratios.empty() &&
        geometry->cell_aspect_ratios.size() != n)
        throw std::invalid_argument("analyze_anisotropy: aspect-ratio size mismatch");
    AnisotropyReport report;
    report.rows = n;
    report.nonzeros = A.nnz();
    const auto* ro = A.row_offsets_data();
    const auto* ci = A.columns_data();
    const auto* av = A.values_data();
    const bool has_centers = geometry && !geometry->cell_centers.empty();
    report.has_geometry = geometry && (
        has_centers || !geometry->cell_aspect_ratios.empty() ||
        !geometry->face_nonorthogonality_degrees.empty() ||
        !geometry->face_skewness.empty());
    report.geometric_dimension = has_centers ? 2 : 0;
    if (has_centers) {
        const double z0 = geometry->cell_centers[0][2];
        for (const auto& x : geometry->cell_centers) {
            for (const double v : x)
                if (!std::isfinite(v)) throw std::invalid_argument("nonfinite cell center");
            if (x[2] != z0) report.geometric_dimension = 3;
        }
    }
    const auto read_quality = [](const std::vector<double>& values,
                                 double& maximum) {
        for (const double v : values) {
            if (!std::isfinite(v) || v < 0.0)
                throw std::invalid_argument("nonfinite or negative mesh quality");
            maximum = std::max(maximum, v);
        }
    };
    if (geometry) {
        read_quality(geometry->cell_aspect_ratios, report.max_cell_aspect_ratio);
        read_quality(geometry->face_nonorthogonality_degrees,
                     report.max_face_nonorthogonality_degrees);
        read_quality(geometry->face_skewness, report.max_face_skewness);
    }
    std::vector<double> maxima(n, 0.0), diagonals(n, 0.0);
    std::vector<std::size_t> strong_count(n, 0);
    std::vector<double> contrast(n, 1.0);
    for (std::size_t i = 0; i < n; ++i) {
        double minimum = std::numeric_limits<double>::infinity();
        for (std::size_t k = ro[i]; k < ro[i+1]; ++k) {
            if (!std::isfinite(av[k]))
                throw std::invalid_argument("analyze_anisotropy: nonfinite matrix");
            if (ci[k] == i) diagonals[i] += av[k];
            else if (av[k] < 0.0) maxima[i] = std::max(maxima[i], -av[k]);
        }
        if (!(diagonals[i] > 0.0)) ++report.nonpositive_diagonals;
        if (!(maxima[i] > 0.0)) continue;
        for (std::size_t k = ro[i]; k < ro[i+1]; ++k) {
            if (ci[k] == i || av[k] >= 0.0 ||
                -av[k] <= maxima[i] * 1e-12) continue;
            minimum = std::min(minimum, -av[k]);
        }
        if (minimum < std::numeric_limits<double>::infinity())
            contrast[i] = maxima[i] / minimum;
        report.worst_offdiagonal_contrast =
            std::max(report.worst_offdiagonal_contrast, contrast[i]);
        if (contrast[i] >= directional_contrast) ++report.anisotropic_rows;
    }
    double strong_fraction_sum = 0.0, alignment_sum = 0.0;
    std::size_t valid_alignments = 0;
    std::vector<std::array<double, 3>> directions(n);
    std::vector<bool> directional(n, false);
    for (std::size_t i = 0; i < n; ++i) {
        if (!(maxima[i] > 0.0)) continue;
        std::array<double, 9> covariance{};
        double strong_weight = 0.0, total_weight = 0.0;
        for (std::size_t k = ro[i]; k < ro[i+1]; ++k) {
            const auto j = ci[k];
            if (j == i || av[k] >= 0.0) continue;
            const double coupling = -av[k];
            total_weight += coupling;
            if (coupling < strength_fraction * maxima[i]) continue;
            ++strong_count[i];
            ++report.strong_directed_edges;
            strong_weight += coupling;
            const double scale = std::max(maxima[j], 1e-300);
            if (coupling < strength_fraction * scale)
                ++report.asymmetric_strong_edges;
            else {
                // Confirm the reverse entry exists and is negative.
                const auto begin = ci + ro[j], end = ci + ro[j+1];
                const auto it = std::lower_bound(begin, end,
                    static_cast<SparseMatrix::Index>(i));
                if (it == end || *it != i ||
                    av[it - ci] >= 0.0)
                    ++report.asymmetric_strong_edges;
            }
            if (has_centers) {
                std::array<double, 3> delta{};
                double length2 = 0.0;
                for (std::size_t a = 0; a < report.geometric_dimension; ++a) {
                    delta[a] = geometry->cell_centers[j][a] - geometry->cell_centers[i][a];
                    length2 += delta[a] * delta[a];
                }
                if (!(length2 > 0.0) || !std::isfinite(length2))
                    throw std::invalid_argument("zero/invalid neighbor distance");
                for (std::size_t a = 0; a < report.geometric_dimension; ++a)
                    for (std::size_t b = 0; b < report.geometric_dimension; ++b)
                        covariance[3*a+b] += coupling * delta[a] * delta[b] / length2;
            }
        }
        if (total_weight > 0.0) strong_fraction_sum += strong_weight / total_weight;
        if (strong_count[i] > 2) ++report.branching_rows;
        if (has_centers && strong_count[i] >= 1) {
            double ratio = 1.0;
            directions[i] = anisotropy_detail::principal_direction(
                covariance, report.geometric_dimension, ratio);
            directional[i] = std::isfinite(ratio) && ratio >= 3.0;
            if (directional[i]) ++report.geometrically_directional_rows;
        }
    }
    report.mean_strong_coupling_fraction = strong_fraction_sum / static_cast<double>(n);
    report.symmetric_strong_graph = (report.asymmetric_strong_edges == 0);
    if (has_centers) {
        for (std::size_t i = 0; i < n; ++i) {
            if (!directional[i]) continue;
            for (std::size_t k = ro[i]; k < ro[i+1]; ++k) {
                const std::size_t j = ci[k];
                if (j <= i || av[k] >= 0.0 || !directional[j] ||
                    -av[k] < strength_fraction * maxima[i]) continue;
                double cosine = 0.0;
                for (std::size_t a = 0; a < report.geometric_dimension; ++a)
                    cosine += directions[i][a] * directions[j][a];
                alignment_sum += std::abs(cosine);
                ++valid_alignments;
            }
        }
    }
    report.mean_neighbor_direction_alignment = valid_alignments
        ? alignment_sum / static_cast<double>(valid_alignments) : 0.0;
    std::vector<std::size_t> rank(n);
    std::iota(rank.begin(), rank.end(), 0);
    std::stable_sort(rank.begin(), rank.end(), [&](std::size_t i, std::size_t j) {
        return contrast[i] > contrast[j];
    });
    rank.resize(std::min<std::size_t>(8, rank.size()));
    report.worst_rows = std::move(rank);
    return report;
}

} // namespace cfdx::core