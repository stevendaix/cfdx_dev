#pragma once
// Experimental production adapter: make a SINGLE complete triplet for one
// U/k/omega/nut candidate. OpenFOAM's omegaWallFunction additionally uses
// cell constraints/corner handling which are not represented exactly here.
// Therefore this is not an equivalence certificate nor a validation gate.
#include "cfdx/physics/sst_openfoam_wall_triplet_v15.h"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {

struct SSTWallOperatorV17 {
    ScalarBoundaryFaceValues omega_faces;
    SSTWallProductionV17 production;
    std::vector<double> momentum_face_mu;
    std::size_t faces=0;
    std::size_t owners=0;
    bool openfoam_equivalent=false;
};

inline SSTWallOperatorV17 make_sst_wall_operator_v17(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const std::vector<std::string>& wall_patches,
    const SstOpenFoamWallTripletV15& triplet,
    const TurbulenceTransportControls& controls)
{
    const auto nf=mesh.n_faces(), nc=mesh.n_cells();
    if (nf==0 || nc==0 || wall_patches.empty() ||
        triplet.omega_by_face.size()!=nf || triplet.nut_by_face.size()!=nf ||
        triplet.has_wall_face.size()!=nf ||
        triplet.production_by_cell.size()!=nc || triplet.has_wall_cell.size()!=nc ||
        triplet.omega_by_cell.size()!=nc ||
        geometry.face_centres.size()!=nf ||
        geometry.face_area_vectors.size()!=nf ||
        geometry.cell_centres.size()!=nc ||
        !(controls.molecular_viscosity>0.) || !(controls.density>0.) ||
        !std::isfinite(controls.molecular_viscosity) || !std::isfinite(controls.density))
        throw std::invalid_argument("V17 wall operator: invalid shape or physical coefficients");
    SSTWallOperatorV17 result;
    result.momentum_face_mu.assign(nf,-1.);
    result.production.g_by_cell=triplet.production_by_cell;
    result.production.has_wall_cell=triplet.has_wall_cell;
    result.production.validate(nc);
    for (const auto& name : wall_patches) {
        if (result.omega_faces.values.count(name))
            throw std::invalid_argument("V17 wall operator: duplicated wall patch");
        const auto& patch=mesh.boundary().patch(mesh.boundary().find(name));
        if (patch.type!=cfdx::core::PatchType::WALL || patch.face_ids.empty())
            throw std::invalid_argument("V17 wall operator: missing WALL patch");
        auto& face_values=result.omega_faces.values[name];
        face_values.assign(nf,std::numeric_limits<double>::quiet_NaN());
        for (const auto face:patch.face_ids) {
            if (face>=nf || !triplet.has_wall_face[face] ||
                mesh.ownership().owner(face)>=nc ||
                !triplet.has_wall_cell[mesh.ownership().owner(face)] ||
                mesh.ownership().neighbour(face)>=0 ||
                !(triplet.omega_by_face[face]>0.) ||
                !std::isfinite(triplet.omega_by_face[face]) ||
                !(triplet.nut_by_face[face]>=0.) ||
                !std::isfinite(triplet.nut_by_face[face]) ||
                result.momentum_face_mu[face]>=0.)
                throw std::invalid_argument("V17 wall operator: incomplete or repeated face");
            face_values[face]=triplet.omega_by_face[face];
            result.momentum_face_mu[face]=controls.density *
                (controls.molecular_viscosity+triplet.nut_by_face[face]);
            if (!std::isfinite(result.momentum_face_mu[face]) ||
                !(result.momentum_face_mu[face]>0.))
                throw std::invalid_argument("V17 wall operator: invalid momentum viscosity");
            ++result.faces;
        }
    }
    if (result.faces != triplet.wall_faces)
        throw std::invalid_argument("V17 wall operator: patch coverage is incomplete");
    for (std::size_t i=0;i<nc;++i) {
        if (triplet.has_wall_cell[i]) {
            if (!(triplet.omega_by_cell[i]>0.) || !std::isfinite(triplet.omega_by_cell[i]))
                throw std::invalid_argument("V17 wall operator: invalid corner omega");
            ++result.owners;
        }
    }
    if (result.owners!=triplet.wall_cells)
        throw std::invalid_argument("V17 wall operator: owner coverage mismatch");
    return result;
}
} // namespace cfdx::physics
