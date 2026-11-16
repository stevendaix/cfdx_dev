#pragma once

// Pure evidence path: compare CFDX's prescribed omega FACE value against
// OpenFOAM Foundation v13's local wall function at the SAME U/k/mesh state.
// Does not claim equivalence of the owner-cell omega/G/nut patch equations.
#include "cfdx/physics/sst_openfoam_wall_reference_v7.h"
#include "cfdx/physics/sst_production_coupling.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct OpenFoamOmegaWallAudit {
    std::size_t face_count=0;
    std::size_t viscous_faces=0;
    std::size_t buffer_faces=0;
    std::size_t logarithmic_faces=0;
    double minimum_k_based_y_plus=std::numeric_limits<double>::infinity();
    double maximum_k_based_y_plus=0.0;
    double worst_relative_omega_gap=0.0;
    std::size_t worst_face=0;
    double worst_cfdx_face_value=0.0;
    double worst_openfoam_local_value=0.0;
};

inline OpenFoamOmegaWallAudit audit_openfoam_omega_wall_reference(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const std::vector<std::string>& wall_patches,
    const ScalarBoundaryFaceValues& actual_omega_face_values,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    const std::map<std::string,std::array<double,3>>& fixed_wall_velocity,
    const TurbulenceTransportControls& controls,
    ReferenceOmegaWallBlend blend=ReferenceOmegaWallBlend::Stepwise)
{
    if (U.size()!=mesh.n_cells() || U.dimension()!=3 ||
        k.size()!=mesh.n_cells() || k.dimension()!=1 ||
        geometry.face_patch.size()!=mesh.n_faces())
        throw std::invalid_argument("omega wall comparison: invalid fields/geometry");
    OpenFoamOmegaWallAudit result;
    for (const auto& name:wall_patches) {
        const auto it=actual_omega_face_values.values.find(name);
        if (it==actual_omega_face_values.values.end() ||
            it->second.size()!=mesh.n_faces())
            throw std::invalid_argument("omega wall comparison: missing CFDX omega faces");
        const auto bc=fixed_wall_velocity.find(name);
        if (bc==fixed_wall_velocity.end())
            throw std::invalid_argument(
                "omega wall comparison: no fixed wall velocity (incomparable)");
        const std::size_t patch=mesh.boundary().find(name);
        const auto& wall=mesh.boundary().patch(patch);
        if (wall.type!=cfdx::core::PatchType::WALL)
            throw std::invalid_argument("omega wall comparison: non-wall patch");
        for (const auto face:wall.face_ids) {
            if (mesh.ownership().neighbour(face)>=0)
                throw std::invalid_argument("omega wall comparison: internal face");
            const std::size_t owner=mesh.ownership().owner(face);
            if (owner>=mesh.n_cells() || geometry.face_patch[face]!=patch)
                throw std::invalid_argument("omega wall comparison: invalid patch incidence");
            const auto Sf=geometry.face_area_vectors.at(face);
            const double area=Sf.mag();
            if (!(area>0.0))
                throw std::invalid_argument("omega wall comparison: zero wall area");
            const auto normal=Sf*(1.0/area);
            const double y=(geometry.face_centres.at(face)-
                geometry.cell_centres.at(owner)).dot(normal);
            if (!(y>0.0) || !std::isfinite(y))
                throw std::invalid_argument("omega wall comparison: invalid normal distance");
            const cfdx::core::Vec3 delta_u{
                U.component_data(0)[owner]-bc->second[0],
                U.component_data(1)[owner]-bc->second[1],
                U.component_data(2)[owner]-bc->second[2]};
            const double grad=delta_u.mag()/y;
            const auto reference=openfoam_omega_wall_reference(
                k(owner),controls.molecular_viscosity,y,grad,blend,
                0.09,controls.beta1);
            const double actual=it->second.at(face);
            if (!std::isfinite(actual) || !(actual>0.0))
                throw std::domain_error("omega wall comparison: invalid CFDX omega face");
            ++result.face_count;
            result.minimum_k_based_y_plus=std::min(
                result.minimum_k_based_y_plus,reference.y_plus_from_k);
            result.maximum_k_based_y_plus=std::max(
                result.maximum_k_based_y_plus,reference.y_plus_from_k);
            if (reference.y_plus_from_k<=5.0) ++result.viscous_faces;
            else if (reference.y_plus_from_k<30.0) ++result.buffer_faces;
            else ++result.logarithmic_faces;
            const double gap=std::abs(actual-reference.reference_omega)/
                std::max(actual,reference.reference_omega);
            if (gap>result.worst_relative_omega_gap) {
                result.worst_relative_omega_gap=gap;
                result.worst_face=face;
                result.worst_cfdx_face_value=actual;
                result.worst_openfoam_local_value=reference.reference_omega;
            }
        }
    }
    if (result.face_count==0)
        throw std::invalid_argument("omega wall comparison: no measured faces");
    return result;
}

} // namespace cfdx::physics
