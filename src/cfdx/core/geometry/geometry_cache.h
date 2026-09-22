// M0.11 — Geometry Cache
//
// Spécification CFDX v0.7 §15, §17 :
//   La géométrie est dérivée de la topologie et mise en cache.
//   /mesh/geometry/
//     face_centres, face_area_vectors, face_areas,
//     delta_coeffs, non_orthogonal_correction,
//     skewness, non_orthogonality
//     cell_centres, cell_volumes
//
//   La topologie reste la source de vérité (§15).
//   GeometryCache est calculé une seule fois par maillage (§17).

#pragma once

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/geometry/mesh_quality.h"
#include <vector>
#include <stdexcept>
#include <cmath>

namespace cfdx {
namespace core {

struct GeometryCache {
    // Face geometry. face_Sf is always oriented owner -> exterior.
    std::vector<Vec3> face_centres;
    std::vector<Vec3> face_Sf;
    std::vector<double> face_areas;
    std::vector<Vec3> face_normals;

    // Face quality metrics.
    std::vector<double> face_skewness;
    std::vector<double> face_non_orthogonality;
    std::vector<double> face_non_orthogonality_deg;

    // Cell geometry.
    std::vector<Vec3> cell_centres;
    std::vector<double> cell_volumes;

    // Distance from owner cell centre to face centre. Kept under the existing
    // public name for compatibility; future code should prefer explicit
    // owner/neighbour distances when an internal-face coefficient is needed.
    std::vector<double> delta_coeffs;

    // Surface closure error per cell: |sum_f Sf_cell|.
    std::vector<double> surface_closure_error;

    // Quality stats.
    double max_skewness = 0.0;
    double max_non_orthogonality_deg = 0.0;
    double min_cell_volume = 0.0;
    double max_cell_volume = 0.0;
    double global_surface_closure_error = 0.0;

    bool valid = false;
};

// Calcule et remplit le cache géométrique à partir du maillage.
//
// Important: face orientation is resolved only after cell centres are known.
// This avoids assuming that every importer uses the same polygon winding while
// preserving the topology as the source of truth.
inline void compute_geometry_cache(const Mesh& m, GeometryCache& cache) {
    const std::size_t n_faces = m.n_faces();
    const std::size_t n_cells = m.n_cells();

    cache = GeometryCache{};
    cache.face_centres.resize(n_faces);
    cache.face_Sf.resize(n_faces);
    cache.face_areas.resize(n_faces);
    cache.face_normals.resize(n_faces);
    cache.face_skewness.resize(n_faces);
    cache.face_non_orthogonality.resize(n_faces);
    cache.face_non_orthogonality_deg.resize(n_faces);
    cache.cell_centres.resize(n_cells);
    cache.cell_volumes.resize(n_cells);
    cache.delta_coeffs.resize(n_faces, 0.0);
    cache.surface_closure_error.resize(n_cells, 0.0);

    const PointCloud& pts = m.points();
    const double* px = pts.x_data();
    const double* py = pts.y_data();
    const double* pz = pts.z_data();
    const auto* verts = m.faces().vertices_data();
    const auto* offsets = m.faces().offsets_data();

    // --- 1. Raw face geometry ---
    for (std::size_t f = 0; f < n_faces; ++f) {
        const VertexIndex off = offsets[f];
        const VertexIndex n = offsets[f + 1] - off;
        const FaceGeometry fg = compute_face_geometry(px, py, pz, verts, off, n);
        if (!(fg.area > 0.0) || !std::isfinite(fg.area) ||
            !std::isfinite(fg.centre.x) || !std::isfinite(fg.centre.y) ||
            !std::isfinite(fg.centre.z)) {
            throw std::runtime_error("GeometryCache: invalid face geometry at face " +
                                     std::to_string(f));
        }
        cache.face_centres[f] = fg.centre;
        cache.face_Sf[f] = fg.Sf;
        cache.face_areas[f] = fg.area;
        cache.face_normals[f] = fg.normal;
    }

    // --- 2. Cell geometry using owner-relative face orientation ---
    const CellConnectivity& cells = m.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();
    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;
        const CellGeometry cg = compute_cell_geometry_oriented(
            cache.face_centres.data(), cache.face_Sf.data(),
            cell_faces + off, n, static_cast<CellIndex>(c), m.ownership());
        if (!(cg.volume > 0.0) || !std::isfinite(cg.volume) ||
            !std::isfinite(cg.signed_volume)) {
            throw std::runtime_error("GeometryCache: non-positive or non-finite cell volume at cell " +
                                     std::to_string(c));
        }
        cache.cell_centres[c] = cg.centre;
        cache.cell_volumes[c] = cg.volume;
    }

    // --- 3. Canonical owner-oriented face vectors ---
    const FaceOwnership& own = m.ownership();
    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::size_t owner = own.owner(f);
        if (owner >= n_cells) {
            throw std::runtime_error("GeometryCache: owner index out of range at face " +
                                     std::to_string(f));
        }

        const std::int64_t neighbour = own.neighbour(f);
        const Vec3* neighbour_centre = nullptr;
        if (neighbour >= 0) {
            const std::size_t nb = static_cast<std::size_t>(neighbour);
            if (nb >= n_cells) {
                throw std::runtime_error("GeometryCache: neighbour index out of range at face " +
                                         std::to_string(f));
            }
            neighbour_centre = &cache.cell_centres[nb];
        }

        ensure_face_orientation(
            cache.face_Sf[f], cache.face_centres[f], cache.cell_centres[owner],
            neighbour_centre);
        cache.face_normals[f] = cache.face_Sf[f].normalized();

        const FaceQuality q = neighbour_centre
            ? compute_face_quality(cache.face_centres[f], cache.cell_centres[owner],
                                   *neighbour_centre, cache.face_Sf[f])
            : compute_face_quality(cache.face_centres[f], cache.cell_centres[owner],
                                   cache.face_Sf[f]);
        cache.face_skewness[f] = q.skewness;
        cache.face_non_orthogonality[f] = q.non_orthogonality;
        cache.face_non_orthogonality_deg[f] = q.non_orthogonality_deg;
        cache.max_skewness = std::max(cache.max_skewness, q.skewness);
        cache.max_non_orthogonality_deg =
            std::max(cache.max_non_orthogonality_deg, q.non_orthogonality_deg);

        cache.delta_coeffs[f] =
            face_cell_distance(cache.face_centres[f], cache.cell_centres[owner]);
    }

    // --- 4. Surface closure with the correct orientation for each cell ---
    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;
        Vec3 sum_Sf;
        for (Offset k = 0; k < n; ++k) {
            const FaceIndex f = cell_faces[off + k];
            const Vec3 Sf_cell = (own.owner(f) == c)
                ? cache.face_Sf[f]
                : cache.face_Sf[f] * (-1.0);
            sum_Sf = sum_Sf + Sf_cell;
        }
        const double closure = sum_Sf.mag();
        cache.surface_closure_error[c] = closure;
        cache.global_surface_closure_error =
            std::max(cache.global_surface_closure_error, closure);
    }

    for (std::size_t c = 0; c < n_cells; ++c) {
        const double vol = cache.cell_volumes[c];
        if (c == 0 || vol < cache.min_cell_volume) cache.min_cell_volume = vol;
        if (c == 0 || vol > cache.max_cell_volume) cache.max_cell_volume = vol;
    }

    cache.valid = true;
}

inline GeometryCache make_geometry_cache(const Mesh& m) {
    GeometryCache cache;
    compute_geometry_cache(m, cache);
    return cache;
}

inline bool is_valid(const GeometryCache& cache, const Mesh& m) {
    if (!cache.valid) return false;
    if (cache.face_centres.size() != m.n_faces()) return false;
    if (cache.face_Sf.size() != m.n_faces()) return false;
    if (cache.cell_centres.size() != m.n_cells()) return false;
    if (cache.cell_volumes.size() != m.n_cells()) return false;
    for (double v : cache.cell_volumes) {
        if (!(v > 0.0) || !std::isfinite(v)) return false;
    }
    return true;
}

}  // namespace core
}  // namespace cfdx
