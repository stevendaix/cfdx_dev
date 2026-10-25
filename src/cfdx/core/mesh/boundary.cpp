// M0.1-T05 — Boundary patches implementation

#include "cfdx/core/mesh/boundary.h"
#include "cfdx/core/mesh/mesh.h"

namespace cfdx {
namespace core {

bool BoundaryPatches::is_consistent_with_mesh(const Mesh& m) const {
    if (!is_consistent(m.n_faces())) return false;

    const FaceOwnership& own = m.ownership();
    std::vector<bool> covered(m.n_faces(), false);

    // is_consistent() has already rejected duplicate/out-of-range ids.
    // Mark each listed boundary face once, retaining the internal-face gate.
    for (const auto& p : patches_) {
        for (const auto fid : p.face_ids) {
            if (own.neighbour(fid) != FaceOwnership::BOUNDARY) {
                return false;
            }
            covered[fid] = true;
        }
    }

    // Membership must be tested in O(1), not by rescanning every patch for
    // every boundary face. In a one-cell-thick 2-D extrusion, the two EMPTY
    // patches contain O(n_cells) faces and the former search was quadratic.
    for (std::size_t f = 0; f < m.n_faces(); ++f) {
        if (own.neighbour(f) == FaceOwnership::BOUNDARY && !covered[f]) {
            return false;
        }
    }

    return true;
}

}  // namespace core
}  // namespace cfdx
