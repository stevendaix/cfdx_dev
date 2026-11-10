#pragma once

// Atomic, comparable OpenFOAM-style wall REFERENCE staged on ONE returned
// U,k,omega,nu_t state. Not yet an OpenFOAM-equivalent production boundary:
// the solver must additionally consume the omega and G constraints in BOTH
// turbulent equations, and nut at the momentum wall face, before claiming it.
// The three arrays are NEVER published separately. No partial opt-in changes
// the physical finite-volume operator.
#include "cfdx/physics/sst_openfoam_wall_coupling_v8.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct SstOpenFoamWallTripletV15 {
    std::vector<double> omega_by_face;
    std::vector<double> nut_by_face;
    std::vector<unsigned char> has_wall_face;
    std::vector<double> production_by_cell;
    std::vector<unsigned char> has_wall_cell;
    std::size_t wall_faces = 0;
    std::size_t wall_cells = 0;
    std::size_t corner_cells = 0;
    // Deliberately false until the production operators actually consume
    // all three quantities (nu_t, G and omega) on the SAME outer state.
    bool applied_to_momentum = false;
    bool applied_to_sst = false;
    bool openfoam_equivalent = false;
};

inline SstOpenFoamWallTripletV15 prepare_openfoam_wall_triplet_v15(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const std::vector<std::string>& wall_patches,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& velocity,
    const SSTProductionState& state,
    const std::map<std::string,std::array<double,3>>& fixed_wall_velocity,
    const TurbulenceTransportControls& controls)
{
    if (velocity.size()!=mesh.n_cells() || velocity.dimension()!=3 ||
        wall_patches.empty())
        throw std::invalid_argument("V15 wall triplet: missing wall or velocity field");
    validate_sst_production_state(state,mesh.n_cells());
    SstOpenFoamWallTripletV15 result;
    result.omega_by_face.assign(mesh.n_faces(),0.);
    result.nut_by_face.assign(mesh.n_faces(),0.);
    result.has_wall_face.assign(mesh.n_faces(),0);
    result.production_by_cell.assign(mesh.n_cells(),0.);
    result.has_wall_cell.assign(mesh.n_cells(),0);
    std::vector<long double> weighted_g(mesh.n_cells(),0.L);
    std::vector<long double> total_area(mesh.n_cells(),0.L);
    std::vector<std::size_t> faces_per_cell(mesh.n_cells(),0);
    for (const auto& name:wall_patches) {
        const auto it=fixed_wall_velocity.find(name);
        if (it==fixed_wall_velocity.end())
            throw std::invalid_argument("V15 wall triplet: no fixed wall velocity");
        const auto patch_index=mesh.boundary().find(name);
        const auto& patch=mesh.boundary().patch(patch_index);
        if (patch.type!=cfdx::core::PatchType::WALL)
            throw std::invalid_argument("V15 wall triplet: patch is not WALL");
        for (const auto face:patch.face_ids) {
            const auto owner=mesh.ownership().owner(face);
            if (mesh.ownership().neighbour(face)>=0 ||
                face>=mesh.n_faces() || owner>=mesh.n_cells() ||
                result.has_wall_face[face])
                throw std::invalid_argument("V15 wall triplet: repeated/non-wall face");
            const auto Sf=geometry.face_area_vectors.at(face);
            const double area=Sf.mag();
            if (!(area>0.) || !std::isfinite(area))
                throw std::invalid_argument("V15 wall triplet: invalid face area");
            const auto normal=Sf*(1./area);
            const auto normal_distance=(geometry.face_centres.at(face)-
                geometry.cell_centres.at(owner)).dot(normal);
            if (!(normal_distance>0.) || !std::isfinite(normal_distance))
                throw std::invalid_argument("V15 wall triplet: invalid wall distance");
            const cfdx::core::Vec3 jump{
                velocity.component_data(0)[owner]-it->second[0],
                velocity.component_data(1)[owner]-it->second[1],
                velocity.component_data(2)[owner]-it->second[2]};
            const double tangential=(jump-normal*jump.dot(normal)).mag();
            const double k=state.k(owner), omega=state.omega(owner);
            const double production_proxy=std::min(
                controls.density*state.nu_t(owner)*state.strain_rate(owner)*state.strain_rate(owner),
                controls.sst_production_limiter*controls.beta_star*controls.density*k*omega);
            const auto point=openfoam_wall_coupling_point_v8(
                k,controls.molecular_viscosity,normal_distance,tangential,
                production_proxy);
            if (!(point.omega>0.) || !std::isfinite(point.production) ||
                !(point.production>=0.) || !std::isfinite(point.nut_spalding) ||
                !(point.nut_spalding>=0.))
                throw std::runtime_error("V15 wall triplet: nonphysical reference");
            result.omega_by_face[face]=point.omega;
            result.nut_by_face[face]=point.nut_spalding;
            result.has_wall_face[face]=1;
            weighted_g[owner]+=static_cast<long double>(point.production)*area;
            total_area[owner]+=area;
            ++faces_per_cell[owner];
            ++result.wall_faces;
        }
    }
    if (!result.wall_faces)
        throw std::invalid_argument("V15 wall triplet: no wall faces");
    for (std::size_t owner=0;owner<mesh.n_cells();++owner) {
        if (!faces_per_cell[owner]) continue;
        if (!(total_area[owner]>0.))
            throw std::invalid_argument("V15 wall triplet: zero weighted owner area");
        result.production_by_cell[owner]=static_cast<double>(weighted_g[owner]/total_area[owner]);
        if (!std::isfinite(result.production_by_cell[owner]))
            throw std::runtime_error("V15 wall triplet: overflow in weighted G");
        result.has_wall_cell[owner]=1;
        ++result.wall_cells;
        if (faces_per_cell[owner]>1) ++result.corner_cells;
    }
    return result;
}

} // namespace cfdx::physics
