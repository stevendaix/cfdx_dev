// M0.3-T01 — Mesh validator (topology + geometry + quality + conservation)
//
// Spécification CFDX v0.7 §19 :
//   Le validator doit être indépendant d'OpenFOAM.
//   Il vérifie :
//     Topologie : indices valides, pas de faces orphelines,
//                 cohérence owner/neighbour, cohérence cellules/faces,
//                 doublons, patches valides.
//     Géométrie : volumes positifs, surfaces non nulles,
//                 NaN, Inf, centres valides, normales cohérentes.
//     Qualité   : skewness, non-orthogonalité, aspect ratio,
//                 volume ratio, cellules dégénérées.
//     Conservation : Σ Sf ≈ 0, Σ flux = 0.
//
// N10 contract:
//   Internal-face non-orthogonality/skewness must be measured from the
//   owner-neighbour centre line, not from owner->face.  The latter is only
//   appropriate for boundary faces.

#pragma once

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/geometry/mesh_quality.h"
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <set>

namespace cfdx {
namespace core {

struct MeshQualityReport {
    bool ok = true;
    std::vector<std::string> errors;
    std::vector<std::string> warnings;

    // Global metrics.
    double max_skewness = 0.0;
    double max_non_orthogonality_rad = 0.0;
    double max_non_orthogonality_deg = 0.0;
    double max_aspect_ratio = 0.0;
    double min_volume_ratio = 1.0;
    double min_cell_volume = 0.0;
    double max_cell_volume = 0.0;
    double surface_closure_error = 0.0;  // |Σ Sf|

    void add_error(const std::string& msg) {
        ok = false;
        errors.push_back(msg);
    }
    void add_warning(const std::string& msg) {
        warnings.push_back(msg);
    }
};

inline double mesh_cell_aspect_ratio(const Mesh& m, std::size_t cell)
{
    const auto offset = m.cells().cell_offset(cell);
    const auto size = m.cells().cell_size(cell);

    std::set<VertexIndex> vertices;
    for (std::size_t k = 0; k < size; ++k) {
        const FaceIndex f = m.cells().faces_data()[offset + k];
        const auto fo = m.faces().face_offset(f);
        const auto fn = m.faces().face_size(f);
        for (std::size_t q = 0; q < fn; ++q)
            vertices.insert(m.faces().vertices_data()[fo + q]);
    }

    if (vertices.empty()) return std::numeric_limits<double>::infinity();

    Vec3 lo{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    Vec3 hi{
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()};

    for (const VertexIndex v : vertices) {
        lo.x = std::min(lo.x, m.points().x(v));
        lo.y = std::min(lo.y, m.points().y(v));
        lo.z = std::min(lo.z, m.points().z(v));
        hi.x = std::max(hi.x, m.points().x(v));
        hi.y = std::max(hi.y, m.points().y(v));
        hi.z = std::max(hi.z, m.points().z(v));
    }

    const double ex = hi.x - lo.x;
    const double ey = hi.y - lo.y;
    const double ez = hi.z - lo.z;
    const double emin = std::min({ex, ey, ez});
    const double emax = std::max({ex, ey, ez});

    if (!(emin > 0.0) || !std::isfinite(emin) || !std::isfinite(emax))
        return std::numeric_limits<double>::infinity();

    return emax / emin;
}

// Valide un maillage complet (§19).
inline MeshQualityReport validate_mesh(const Mesh& m) {
    MeshQualityReport report;

    // --- 1. Topologie ---
    auto topo = m.topo_validate();
    if (!topo.ok) {
        for (const auto& e : topo.errors) report.add_error(e);
    }

    const std::size_t n_points = m.n_points();
    const std::size_t n_faces = m.n_faces();
    const std::size_t n_cells = m.n_cells();

    if (n_points == 0) report.add_error("mesh has no points");
    if (n_faces == 0) report.add_error("mesh has no faces");
    if (n_cells == 0) report.add_error("mesh has no cells");

    // --- 2. Face geometry ---
    std::vector<Vec3> face_centres(n_faces);
    std::vector<Vec3> face_Sf(n_faces);
    std::vector<double> face_areas(n_faces);

    for (std::size_t f = 0; f < n_faces; ++f) {
        const auto offset = m.faces().face_offset(f);
        const auto size = m.faces().face_size(f);
        try {
            auto g = compute_face_geometry(
                m.points().x_data(), m.points().y_data(), m.points().z_data(),
                m.faces().vertices_data(), offset, size);
            face_centres[f] = g.centre;
            face_Sf[f] = g.Sf;
            face_areas[f] = g.area;

            if (g.area <= 0.0) {
                report.add_error("face " + std::to_string(f) + " has non-positive area");
            }
            if (!std::isfinite(g.centre.x) || !std::isfinite(g.centre.y) ||
                !std::isfinite(g.centre.z)) {
                report.add_error("face " + std::to_string(f) + " has non-finite centre");
            }
        } catch (const std::exception& e) {
            report.add_error(std::string("face ") + std::to_string(f) + ": " + e.what());
        }
    }

    // --- 3. Cell geometry ---
    std::vector<Vec3> cell_centres(n_cells);
    std::vector<double> cell_volumes(n_cells);

    for (std::size_t c = 0; c < n_cells; ++c) {
        const auto offset = m.cells().cell_offset(c);
        const auto size = m.cells().cell_size(c);
        try {
            auto geom = compute_cell_geometry_oriented(
                face_centres.data(), face_Sf.data(),
                m.cells().faces_data() + offset, size,
                static_cast<CellIndex>(c), m.ownership());
            cell_centres[c] = geom.centre;
            cell_volumes[c] = geom.volume;

            if (!(geom.signed_volume > 0.0) || !std::isfinite(geom.signed_volume)) {
                report.add_error("cell " + std::to_string(c) +
                                 " has non-positive or inverted signed volume");
            }
            if (!(geom.volume > 0.0) || !std::isfinite(geom.volume)) {
                report.add_error("cell " + std::to_string(c) + " has invalid volume");
            }

            const double aspect = mesh_cell_aspect_ratio(m, c);
            if (!std::isfinite(aspect) || aspect <= 0.0) {
                report.add_error("cell " + std::to_string(c) +
                                 " has invalid aspect ratio");
            } else {
                report.max_aspect_ratio = std::max(report.max_aspect_ratio, aspect);
            }

            if (c == 0 || geom.volume < report.min_cell_volume)
                report.min_cell_volume = geom.volume;
            if (c == 0 || geom.volume > report.max_cell_volume)
                report.max_cell_volume = geom.volume;
        } catch (const std::exception& e) {
            report.add_error(std::string("cell ") + std::to_string(c) + ": " + e.what());
        }
    }

    // --- 4. Face quality ---
    for (std::size_t f = 0; f < n_faces; ++f) {
        const auto owner = m.ownership().owner(f);
        if (owner >= n_cells) continue;

        const auto neighbour = m.ownership().neighbour(f);
        FaceQuality q;
        if (neighbour >= 0 && static_cast<std::size_t>(neighbour) < n_cells) {
            // Internal-face quality must use the owner-neighbour centre line.
            q = compute_face_quality(
                face_centres[f], cell_centres[owner],
                cell_centres[static_cast<std::size_t>(neighbour)], face_Sf[f]);

            const double vo = cell_volumes[owner];
            const double vn = cell_volumes[static_cast<std::size_t>(neighbour)];
            if (vo > 0.0 && vn > 0.0) {
                report.min_volume_ratio = std::min(
                    report.min_volume_ratio,
                    std::min(vo / vn, vn / vo));
            }
        } else {
            q = compute_face_quality(
                face_centres[f], cell_centres[owner], face_Sf[f]);
        }

        if (q.skewness > report.max_skewness) report.max_skewness = q.skewness;
        if (q.non_orthogonality > report.max_non_orthogonality_rad) {
            report.max_non_orthogonality_rad = q.non_orthogonality;
            report.max_non_orthogonality_deg = q.non_orthogonality_deg;
        }

        if (q.skewness > 0.5) {
            report.add_warning("face " + std::to_string(f) +
                               " has high skewness (" + std::to_string(q.skewness) + ")");
        }
        if (q.non_orthogonality_deg > 70.0) {
            report.add_warning("face " + std::to_string(f) +
                               " has high non-orthogonality (" +
                               std::to_string(q.non_orthogonality_deg) + " deg)");
        }
    }

    if (report.max_aspect_ratio > 100.0) {
        report.add_warning("mesh has very high cell aspect ratio (" +
                           std::to_string(report.max_aspect_ratio) + ")");
    }
    if (report.min_volume_ratio < 1e-2) {
        report.add_warning("mesh has strong cell-volume disparity (minimum ratio " +
                           std::to_string(report.min_volume_ratio) + ")");
    }

    // --- 5. Conservation: Σ Sf ≈ 0 for each cell ---
    for (std::size_t c = 0; c < n_cells; ++c) {
        const auto offset = m.cells().cell_offset(c);
        const auto size = m.cells().cell_size(c);
        Vec3 sum_Sf;
        for (std::size_t k = 0; k < size; ++k) {
            const FaceIndex f = m.cells().faces_data()[offset + k];
            const bool owner = m.ownership().owner(f) == c;
            sum_Sf = sum_Sf + (owner ? face_Sf[f] : face_Sf[f] * (-1.0));
        }
        const double closure_error = sum_Sf.mag();
        if (closure_error > report.surface_closure_error)
            report.surface_closure_error = closure_error;

        if (cell_volumes[c] > 0.0) {
            const double rel = closure_error / (3.0 * cell_volumes[c]);
            if (rel > 1e-6) {
                report.add_warning("cell " + std::to_string(c) +
                                   " surface closure error " +
                                   std::to_string(rel) + " (relative)");
            }
        }
    }

    return report;
}

}  // namespace core
}  // namespace cfdx
