// M0.1-T06 — Mesh topology + validator
//
// Spécification CFDX v0.7 §9, §19 :
//   Le Mesh regroupe les données fondamentales :
//     points, face_vertices, face_vertices_offsets,
//     face_owner, face_neighbour, cell_faces, cell_faces_offsets,
//     boundary patches
//
//   Le Mesh Validator vérifie :
//     Topologie : indices valides, pas de faces orphelines,
//                 cohérence owner/neighbour, cohérence cellules/faces,
//                 doublons, patches valides.

#pragma once

#include "point.h"
#include "face.h"
#include "ownership.h"
#include "cell.h"
#include "boundary.h"
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>

namespace cfdx {
namespace core {

struct MeshStats {
    std::size_t n_points = 0;
    std::size_t n_faces = 0;
    std::size_t n_cells = 0;
    std::size_t n_boundary_faces = 0;
    std::size_t n_internal_faces = 0;
    std::size_t n_patches = 0;
};

class Mesh {
public:
    Mesh() = default;

    // --- Accès aux composants ---

    PointCloud& points() { return points_; }
    const PointCloud& points() const { return points_; }

    FaceConnectivity& faces() { return faces_; }
    const FaceConnectivity& faces() const { return faces_; }

    FaceOwnership& ownership() { return ownership_; }
    const FaceOwnership& ownership() const { return ownership_; }

    CellConnectivity& cells() { return cells_; }
    const CellConnectivity& cells() const { return cells_; }

    BoundaryPatches& boundary() { return boundary_; }
    const BoundaryPatches& boundary() const { return boundary_; }
    void set_boundary(const BoundaryPatches& bp) { boundary_ = bp; }

    // --- Dimensions ---

    std::size_t n_points() const noexcept { return points_.size(); }
    std::size_t n_faces() const noexcept { return faces_.n_faces(); }
    std::size_t n_cells() const noexcept { return cells_.n_cells(); }

    MeshStats stats() const {
        MeshStats s;
        s.n_points = n_points();
        s.n_faces = n_faces();
        s.n_cells = n_cells();
        s.n_boundary_faces = ownership_.n_boundary_faces();
        s.n_internal_faces = ownership_.n_internal_faces();
        s.n_patches = boundary_.n_patches();
        return s;
    }

    // --- Validation complète (§19) ---

    struct ValidationResult {
        bool ok = true;
        std::vector<std::string> errors;
        std::vector<std::string> warnings;

        void add_error(const std::string& msg) {
            ok = false;
            errors.push_back(msg);
        }
        void add_warning(const std::string& msg) {
            warnings.push_back(msg);
        }
    };

    ValidationResult validate() const {
        ValidationResult result;

        // 1. Dimensions cohérentes
        if (ownership_.size() != faces_.n_faces()) {
            result.add_error("owner/neighbour size mismatch: " +
                             std::to_string(ownership_.size()) + " vs " +
                             std::to_string(faces_.n_faces()) + " faces");
        }

        // 2. Validité des points (NaN/Inf)
        if (!points_.is_valid()) {
            result.add_error("points contain NaN or Inf values");
        }

        // 3. Cohérence des offsets CSR (faces)
        if (!faces_.is_consistent()) {
            result.add_error("face CSR offsets inconsistent");
        }

        // 4. Validité des indices de sommets
        if (!faces_.indices_valid(n_points())) {
            result.add_error("face vertex indices out of range");
        }

        // 5. Cohérence owner/neighbour
        if (!ownership_.is_consistent(n_cells())) {
            result.add_error("owner/neighbour inconsistent with cell count");
        }

        // 6. Cohérence des offsets CSR (cells)
        if (!cells_.is_consistent()) {
            result.add_error("cell CSR offsets inconsistent");
        }

        // 7. Validité des face_ids dans cellules
        if (!cells_.face_ids_valid(n_faces())) {
            result.add_error("cell face ids out of range");
        }

        // 8. Validité des patches
        if (!boundary_.is_consistent(n_faces())) {
            result.add_error("boundary patches inconsistent (overlap or out-of-range)");
        }

        // 9. Chaque face interne doit être référencée par 2 cellules,
        //    chaque face de frontière par 1. On vérifie le total des références
        //    de cellules contre la somme attendue.
        const std::size_t expected_refs =
            ownership_.n_internal_faces() * 2 + ownership_.n_boundary_faces() * 1;
        if (cells_.n_face_refs() != expected_refs) {
            result.add_error("cell face reference count mismatch: got " +
                             std::to_string(cells_.n_face_refs()) +
                             ", expected " + std::to_string(expected_refs));
        }

        return result;
    }

    // --- Nettoyage ---

    void clear() {
        points_.clear();
        faces_.clear();
        ownership_.clear();
        cells_.clear();
        boundary_.clear();
    }

private:
    PointCloud points_;
    FaceConnectivity faces_;
    FaceOwnership ownership_;
    CellConnectivity cells_;
    BoundaryPatches boundary_;
};

}  // namespace core
}  // namespace cfdx