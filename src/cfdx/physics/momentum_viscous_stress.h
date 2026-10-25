#pragma once

#include "cfdx/physics/scalar_transport_gradient.h"
#include <array>

namespace cfdx::physics {

// Owner-oriented flux of mu_eff * transpose(grad(U)). The scalar diffusion
// operator already supplies mu_eff * grad(U). This complementary conservative
// flux is required when mu_eff varies, even for divergence-free velocity.
// Internal mu uses the same normal-distance harmonic interpolation as scalar
// diffusion; boundary mu and face gradients use their owner-cell values.
// Boundary U values enter the Green--Gauss reconstruction explicitly. No
// non-orthogonal/skew correction or higher-order wall extrapolation is implied.
inline cfdx::core::Field<double,cfdx::core::Location::FACE>
transpose_viscous_stress_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& velocity,
    const std::vector<double>& cell_mu,
    const std::array<ScalarBoundaryConditions,3>& velocity_bcs)
{
    using namespace cfdx::core;
    const auto nc=mesh.n_cells(), nf=mesh.n_faces();
    if (velocity.size()!=nc || velocity.dimension()!=3 || cell_mu.size()!=nc)
        throw std::invalid_argument("transpose viscous stress: invalid field shape");
    for (const double mu : cell_mu)
        if (!std::isfinite(mu) || mu<0.0)
            throw std::invalid_argument("transpose viscous stress: invalid dynamic viscosity");
    std::array<Field<double,Location::CELL>,3> gradients;
    for (std::size_t component=0; component<3; ++component) {
        Field<double,Location::CELL> scalar(nc,"U_stress_component","m/s",1);
        for (std::size_t cell=0; cell<nc; ++cell) scalar(cell)=velocity(cell,component);
        gradients[component]=scalar_gradient_with_transport_boundaries(
            scalar,mesh,geometry,velocity_bcs[component]);
    }
    Field<double,Location::FACE> flux(nf,"transpose_viscous_stress_flux","N",3);
    flux.fill(0.0);
    for (std::size_t face=0; face<nf; ++face) {
        const auto owner=mesh.ownership().owner(face);
        const auto neighbour=mesh.ownership().neighbour(face);
        const auto patch=geometry.face_patch[face];
        if (neighbour<0 && patch<mesh.boundary().n_patches() &&
            mesh.boundary().patch(patch).type==PatchType::EMPTY)
            continue;
        double mu=cell_mu[owner], weight=0.0;
        const auto& Sf=geometry.face_area_vectors[face];
        if (neighbour>=0) {
            const auto nb=static_cast<std::size_t>(neighbour);
            const double area=Sf.mag();
            const double dp=std::abs((geometry.face_centres[face]-geometry.cell_centres[owner]).dot(Sf))/area;
            const double dn=std::abs((geometry.cell_centres[nb]-geometry.face_centres[face]).dot(Sf))/area;
            if (!std::isfinite(dp) || !std::isfinite(dn) || !(dp>0.0) || !(dn>0.0))
                throw std::invalid_argument("transpose viscous stress: invalid normal half-distance");
            mu=mu>0.0 && cell_mu[nb]>0.0 ? (dp+dn)/(dp/mu+dn/cell_mu[nb]) : 0.0;
            weight=geometric_interpolation_weight(
                geometry.face_centres[face],geometry.cell_centres[owner],&geometry.cell_centres[nb]);
        }
        const std::array<double,3> areas{Sf.x,Sf.y,Sf.z};
        for (std::size_t component=0; component<3; ++component) {
            double contraction=0.0;
            for (std::size_t direction=0; direction<3; ++direction) {
                double gradient=gradients[direction](owner,component);
                if (neighbour>=0)
                    gradient=(1.0-weight)*gradient+
                        weight*gradients[direction](static_cast<std::size_t>(neighbour),component);
                contraction+=gradient*areas[direction];
            }
            flux(face,component)=mu*contraction;
            if (!std::isfinite(flux(face,component)))
                throw std::runtime_error("transpose viscous stress: non-finite face flux");
        }
    }
    return flux;
}

inline cfdx::core::Field<double,cfdx::core::Location::CELL>
viscous_face_flux_divergence(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& flux)
{
    using namespace cfdx::core;
    if (flux.size()!=mesh.n_faces() || flux.dimension()!=3 ||
        geometry.cell_volumes.size()!=mesh.n_cells())
        throw std::invalid_argument("viscous face divergence: invalid field shape");
    Field<double,Location::CELL> result(mesh.n_cells(),"transpose_viscous_force","N/m3",3);
    result.fill(0.0);
    for (std::size_t face=0; face<mesh.n_faces(); ++face) {
        const auto owner=mesh.ownership().owner(face);
        const auto neighbour=mesh.ownership().neighbour(face);
        if (owner>=mesh.n_cells() ||
            (neighbour>=0 && static_cast<std::size_t>(neighbour)>=mesh.n_cells()))
            throw std::invalid_argument("viscous face divergence: invalid face ownership");
        for (std::size_t component=0; component<3; ++component) {
            if (!std::isfinite(flux(face,component)))
                throw std::invalid_argument("viscous face divergence: non-finite flux");
            result(owner,component)+=flux(face,component);
            if (neighbour>=0) result(static_cast<std::size_t>(neighbour),component)-=flux(face,component);
        }
    }
    for (std::size_t cell=0; cell<mesh.n_cells(); ++cell) {
        if (!(geometry.cell_volumes[cell]>0.0) || !std::isfinite(geometry.cell_volumes[cell]))
            throw std::invalid_argument("viscous face divergence: invalid cell volume");
        for (std::size_t component=0; component<3; ++component)
            result(cell,component)/=geometry.cell_volumes[cell];
    }
    return result;
}

} // namespace cfdx::physics
