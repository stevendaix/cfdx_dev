// M0.1-T06 — Mesh topology + topological validator
//
// Spécification CFDX v0.7 §9, §19 :
//   Le Mesh regroupe les données fondamentales :
//     points, face_vertices, face_vertices_offsets,
//     face_owner, face_neighbour, cell_faces, cell_faces_offsets,
//     boundary patches
//
//   Validation complète (géométrie + qualité + conservation) :
//     cfdx::core::validate_mesh()  — voir mesh_validator.h

#pragma once

#include "point.h"
#include "face.h"
#include "ownership.h"
#include "cell.h"
#include "boundary.h"
#include "mesh_topology.h"
#include "index_types.h"
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

struct TopoValidation {
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

class Mesh {
public:
    Mesh() = default;

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

    std::size_t n_points() const noexcept { return points_.size(); }
    std::size_t n_faces() const noexcept { return faces_.n_faces(); }
    std::size_t n_cells() const noexcept { return cells_.n_cells(); }

    // Topology-only view for algorithms that should not depend on geometry.
    MeshTopology topology() const noexcept {
        return MeshTopology(faces_, ownership_, cells_, boundary_);
    }

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

    // --- Validation topologique (§19) ---
    TopoValidation topo_validate() const {
        TopoValidation result;

        if (ownership_.size() != faces_.n_faces()) {
            result.add_error("owner/neighbour size mismatch: " +
                             std::to_string(ownership_.size()) + " vs " +
                             std::to_string(faces_.n_faces()) + " faces");
        }

        if (!points_.is_valid()) {
            result.add_error("points contain NaN or Inf values");
        }

        if (!faces_.is_consistent()) {
            result.add_error("face CSR offsets inconsistent");
        }

        if (!faces_.indices_valid(n_points())) {
            result.add_error("face vertex indices out of range");
        }

        if (!ownership_.is_consistent(n_cells())) {
            result.add_error("owner/neighbour inconsistent with cell count");
        }

        if (!cells_.is_consistent()) {
            result.add_error("cell CSR offsets inconsistent");
        }

        if (!cells_.face_ids_valid(n_faces())) {
            result.add_error("cell face ids out of range");
        }

        if (!boundary_.is_consistent(n_faces())) {
            result.add_error("boundary patches inconsistent (overlap or out-of-range)");
        }
        
        // Vérifier la cohérence complète patches <-> topologie (§19)
        if (!boundary_.is_consistent_with_mesh(*this)) {
            result.add_error("boundary patches do not match mesh topology (orphan or internal faces in patches)");
        }

        // Chaque face interne doit être référencée par 2 cellules,
        // chaque face de frontière par 1.
        const std::size_t expected_refs =
            ownership_.n_internal_faces() * 2 + ownership_.n_boundary_faces() * 1;
        if (cells_.n_face_refs() != expected_refs) {
            result.add_error("cell face reference count mismatch: got " +
                             std::to_string(cells_.n_face_refs()) +
                             ", expected " + std::to_string(expected_refs));
        }

        // Vérification de cohérence owner/neighbour <-> cell faces
        // (plus lent mais plus complet)
        const auto* cell_faces = cells_.faces_data();
        const auto* cell_offsets = cells_.offsets_data();
        for (std::size_t f = 0; f < n_faces(); ++f) {
            const CellIndex owner = ownership_.owner(f);
            const std::int64_t neighbour = ownership_.neighbour(f);
            if (owner >= n_cells()) continue;
            
            bool found_owner = false;
            const Offset off = cell_offsets[owner];
            const Offset n = cell_offsets[owner + 1] - off;
            for (Offset k = 0; k < n; ++k) {
                if (cell_faces[off + k] == static_cast<FaceIndex>(f)) {
                    found_owner = true;
                    break;
                }
            }
            if (!found_owner) {
                result.add_error("face " + std::to_string(f) +
                                 " not found in owner cell " + std::to_string(owner));
            }
            if (neighbour >= 0 && static_cast<std::size_t>(neighbour) < n_cells()) {
                const CellIndex nb = static_cast<CellIndex>(neighbour);
                bool found_neighbour = false;
                const Offset off_nb = cell_offsets[nb];
                const Offset n_nb = cell_offsets[nb + 1] - off_nb;
                for (Offset k = 0; k < n_nb; ++k) {
                    if (cell_faces[off_nb + k] == static_cast<FaceIndex>(f)) {
                        found_neighbour = true;
                        break;
                    }
                }
                if (!found_neighbour) {
                    result.add_error("face " + std::to_string(f) +
                                     " not found in neighbour cell " + std::to_string(nb));
                }
            }
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
