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

namespace cfdx {
namespace core {

struct GeometryCache {
    // Face geometry
    std::vector<Vec3> face_centres;
    std::vector<Vec3> face_Sf;
    std::vector<double> face_areas;
    std::vector<Vec3> face_normals;
    
    // Face quality metrics
    std::vector<double> face_skewness;
    std::vector<double> face_non_orthogonality;
    std::vector<double> face_non_orthogonality_deg;
    
    // Cell geometry
    std::vector<Vec3> cell_centres;
    std::vector<double> cell_volumes;
    
    // Delta coefficients (for non-orthogonal correction)
    std::vector<double> delta_coeffs;
    
    // Surface closure error per cell
    std::vector<double> surface_closure_error;
    
    // Quality stats
    double max_skewness = 0.0;
    double max_non_orthogonality_deg = 0.0;
    double min_cell_volume = 0.0;
    double max_cell_volume = 0.0;
    double global_surface_closure_error = 0.0;
    
    bool valid = false;
};

// Calcule et remplit le cache géométrique à partir du maillage.
inline void compute_geometry_cache(const Mesh& m, GeometryCache& cache) {
    const std::size_t n_points = m.n_points();
    const std::size_t n_faces = m.n_faces();
    const std::size_t n_cells = m.n_cells();
    
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
    cache.surface_closure_error.resize(n_cells);
    
    // --- 1. Face geometry ---
    const PointCloud& pts = m.points();
    const double* px = pts.x_data();
    const double* py = pts.y_data();
    const double* pz = pts.z_data();
    const auto* verts = m.faces().vertices_data();
    const auto* offsets = m.faces().offsets_data();
    
    for (std::size_t f = 0; f < n_faces; ++f) {
        const VertexIndex off = offsets[f];
        const VertexIndex n = offsets[f + 1] - off;
        const FaceGeometry fg = compute_face_geometry(px, py, pz, verts, off, n);
        cache.face_centres[f] = fg.centre;
        cache.face_Sf[f] = fg.Sf;
        cache.face_areas[f] = fg.area;
        cache.face_normals[f] = fg.normal;
    }
    
    // --- 2. Establish the authoritative global face orientation ---
    // The provisional centres are orientation-independent because they use only
    // scalar face areas. They are then used to enforce owner -> neighbour (and
    // owner -> boundary) orientation before any signed cell geometry is built.
    compute_area_weighted_cell_centres(
        m, cache.face_centres.data(), cache.face_Sf.data(), cache.cell_centres.data());
    orient_mesh_face_vectors(m, cache.face_centres, cache.cell_centres, cache.face_Sf);
    for (std::size_t f = 0; f < n_faces; ++f)
        cache.face_normals[f] = cache.face_Sf[f].normalized();
    
    // --- 3. Cell geometry ---
    const CellConnectivity& cells = m.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();
    
    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;
        const CellGeometry cg = compute_cell_geometry(m, cache.face_centres.data(), cache.face_Sf.data(), cell_faces + off, c, n);
        cache.cell_centres[c] = cg.centre;
        cache.cell_volumes[c] = cg.volume;
        if (!(cg.signed_volume > 0.0))
            throw std::runtime_error("compute_geometry_cache: inverted cell orientation");
    }
    
    // --- 4. Quality metrics (skewness, non-orthogonality) ---
    const FaceOwnership& own = m.ownership();
    cache.max_skewness = 0.0;
    cache.max_non_orthogonality_deg = 0.0;
    cache.min_cell_volume = 0.0;
    cache.max_cell_volume = 0.0;
    cache.global_surface_closure_error = 0.0;
    
    for (std::size_t c = 0; c < n_cells; ++c) {
        const double vol = cache.cell_volumes[c];
        if (c == 0 || vol < cache.min_cell_volume) cache.min_cell_volume = vol;
        if (c == 0 || vol > cache.max_cell_volume) cache.max_cell_volume = vol;
    }
    
    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::size_t owner = own.owner(f);
        if (owner >= n_cells) continue;
        
        const Vec3* neigh_centre = nullptr;
        const std::int64_t neighbour = own.neighbour(f);
        if (neighbour >= 0) {
            neigh_centre = &cache.cell_centres[static_cast<std::size_t>(neighbour)];
        }
        
        const FaceQuality q = compute_face_quality(
            cache.face_centres[f], cache.cell_centres[owner], cache.face_Sf[f]);
        cache.face_skewness[f] = q.skewness;
        cache.face_non_orthogonality[f] = q.non_orthogonality;
        cache.face_non_orthogonality_deg[f] = q.non_orthogonality_deg;
        
        if (q.skewness > cache.max_skewness) cache.max_skewness = q.skewness;
        if (q.non_orthogonality_deg > cache.max_non_orthogonality_deg) {
            cache.max_non_orthogonality_deg = q.non_orthogonality_deg;
        }
    }
    
    // --- 5. Surface closure (conservation) ---
    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;
        Vec3 sum_Sf;
        for (Offset k = 0; k < n; ++k) {
            const FaceIndex f = cell_faces[off + k];
            const std::size_t owner = own.owner(f);
            const Vec3 local_Sf = owner == c ? cache.face_Sf[f] : cache.face_Sf[f] * -1.0;
            sum_Sf = sum_Sf + local_Sf;
        }
        const double closure = sum_Sf.mag();
        cache.surface_closure_error[c] = closure;
        if (closure > cache.global_surface_closure_error) {
            cache.global_surface_closure_error = closure;
        }
    }
    
    // --- 6. Delta coefficients (distance from cell centre to face centre) ---
    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::size_t owner = own.owner(f);
        if (owner < n_cells) {
            cache.delta_coeffs[f] = face_cell_distance(cache.face_centres[f], cache.cell_centres[owner]);
        }
    }
    
    cache.valid = true;
}

// Interface pratique : retourne un cache prêt à l'emploi.
inline GeometryCache make_geometry_cache(const Mesh& m) {
    GeometryCache cache;
    compute_geometry_cache(m, cache);
    return cache;
}

// Vérifie la validité du cache.
inline bool is_valid(const GeometryCache& cache, const Mesh& m) {
    if (!cache.valid) return false;
    if (cache.face_centres.size() != m.n_faces()) return false;
    if (cache.cell_centres.size() != m.n_cells()) return false;
    return true;
}

}  // namespace core
}  // namespace cfdx
