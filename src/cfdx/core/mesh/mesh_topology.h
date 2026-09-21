#pragma once

#include "cfdx/core/mesh/boundary.h"
#include "cfdx/core/mesh/cell.h"
#include "cfdx/core/mesh/face.h"
#include "cfdx/core/mesh/ownership.h"

namespace cfdx::core {

// Non-geometric mesh topology view. Geometry, fields and solver state remain
// outside this object so topology-only algorithms can depend on a narrow API.
class MeshTopology {
public:
    MeshTopology(const FaceConnectivity& faces,
                 const FaceOwnership& ownership,
                 const CellConnectivity& cells,
                 const BoundaryPatches& boundary)
        : faces_(&faces), ownership_(&ownership), cells_(&cells), boundary_(&boundary) {}

    const FaceConnectivity& faces() const noexcept { return *faces_; }
    const FaceOwnership& ownership() const noexcept { return *ownership_; }
    const CellConnectivity& cells() const noexcept { return *cells_; }
    const BoundaryPatches& boundary() const noexcept { return *boundary_; }

    std::size_t n_faces() const noexcept { return faces_->n_faces(); }
    std::size_t n_cells() const noexcept { return cells_->n_cells(); }

private:
    const FaceConnectivity* faces_;
    const FaceOwnership* ownership_;
    const CellConnectivity* cells_;
    const BoundaryPatches* boundary_;
};

} // namespace cfdx::core
