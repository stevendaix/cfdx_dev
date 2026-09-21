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

#pragma once

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/geometry/mesh_quality.h"
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

namespace cfdx {
namespace core {

struct MeshQualityReport {
    bool ok = true;
    std::vector<std::string> errors;
    std::vector<std::string> warnings;

    // Métriques globales.
    double max_skewness = 0.0;
    double max_non_orthogonality_rad = 0.0;
    double max_non_orthogonality_deg = 0.0;
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

// Valide un maillage complet (§19).
inline MeshQualityReport validate_mesh(const Mesh& m) {
    MeshQualityReport report;

    // --- 1. Topologie (déjà dans Mesh::topo_validate) ---
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

    // --- 2. Géométrie des faces ---
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

    // --- 3. Géométrie des cellules ---
    std::vector<Vec3> cell_centres(n_cells);
    std::vector<double> cell_volumes(n_cells);
    compute_area_weighted_cell_centres(m, face_centres.data(), face_Sf.data(), cell_centres.data());
    orient_mesh_face_vectors(m, face_centres, cell_centres, face_Sf);

    for (std::size_t c = 0; c < n_cells; ++c) {
        const auto offset = m.cells().cell_offset(c);
        const auto size = m.cells().cell_size(c);
        try {
            auto geom = compute_cell_geometry(m, face_centres.data(), face_Sf.data(), m.cells().faces_data() + offset, c, size);
            cell_centres[c] = geom.centre;
            cell_volumes[c] = geom.volume;

            if (!(geom.signed_volume > 0.0) || !std::isfinite(geom.signed_volume)) {
                report.add_error("cell " + std::to_string(c) +
                                 " has non-positive or inverted signed volume");
            }
            if (!(geom.volume > 0.0) || !std::isfinite(geom.volume)) {
                report.add_error("cell " + std::to_string(c) + " has invalid volume");
            }

            // Stats.
            if (c == 0 || geom.volume < report.min_cell_volume)
                report.min_cell_volume = geom.volume;
            if (c == 0 || geom.volume > report.max_cell_volume)
                report.max_cell_volume = geom.volume;
        } catch (const std::exception& e) {
            report.add_error(std::string("cell ") + std::to_string(c) + ": " + e.what());
        }
    }

    // --- 4. Qualité (skewness, non-orthogonalité) ---
    for (std::size_t f = 0; f < n_faces; ++f) {
        const auto owner = m.ownership().owner(f);
        if (owner >= n_cells) continue;  // shouldn't happen if topology is valid
        const auto q = compute_face_quality(face_centres[f], cell_centres[owner], face_Sf[f]);

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

    // --- 5. Conservation : Σ Sf ≈ 0 pour chaque cellule ---
    for (std::size_t c = 0; c < n_cells; ++c) {
        const auto offset = m.cells().cell_offset(c);
        const auto size = m.cells().cell_size(c);
        Vec3 sum_Sf;
        for (std::size_t k = 0; k < size; ++k) {
            const FaceIndex f = m.cells().faces_data()[offset + k];
            sum_Sf = sum_Sf + face_Sf[f];
        }
        const double closure_error = sum_Sf.mag();
        if (closure_error > report.surface_closure_error)
            report.surface_closure_error = closure_error;

        // Tolérance relative au volume de la cellule.
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
