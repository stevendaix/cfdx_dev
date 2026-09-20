// M0.1-T05 — Boundary patches implementation

#include "cfdx/core/mesh/boundary.h"
#include "cfdx/core/mesh/mesh.h"

namespace cfdx {
namespace core {

bool BoundaryPatches::is_consistent_with_mesh(const Mesh& m) const {
    if (!is_consistent(m.n_faces())) return false;

    const FaceOwnership& own = m.ownership();
    std::vector<bool> is_boundary_face(m.n_faces(), false);

    // Marquer les faces de frontière selon owner/neighbour
    for (std::size_t f = 0; f < m.n_faces(); ++f) {
        is_boundary_face[f] = (own.neighbour(f) == FaceOwnership::BOUNDARY);
    }

    // Vérifier que chaque face dans un patch est bien une face de frontière
    for (const auto& p : patches_) {
        for (const auto fid : p.face_ids) {
            if (!is_boundary_face[fid]) {
                return false;  // Face interne dans un patch de frontière
            }
        }
    }

    // Vérifier que toutes les faces de frontière sont dans un patch
    for (std::size_t f = 0; f < m.n_faces(); ++f) {
        if (is_boundary_face[f]) {
            bool found = false;
            for (const auto& p : patches_) {
                for (const auto fid : p.face_ids) {
                    if (fid == static_cast<FaceIndex>(f)) {
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            if (!found) return false;  // Face de frontière orpheline
        }
    }

    return true;
}

}  // namespace core
}  // namespace cfdx
