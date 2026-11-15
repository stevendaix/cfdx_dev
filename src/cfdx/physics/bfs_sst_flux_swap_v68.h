@@ -0,0 +1,162 @@
#pragma once
// V68: read-only SST carrier/physical mass-flux substitution at EXACTLY one
// U,p,k,omega state. The SST source, diffusion, wall and strain operators are
// evaluated by the caller ONCE per state; only face flux may differ.
// The native V22 FV residual is the oracle, not the face formula below.
#include "cfdx/physics/sst_solver.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct BfsSstFluxSwapV68 {
    std::size_t cells=0, faces=0, focus_face=0;
    std::size_t omega_worst_physical=0, omega_worst_carrier=0;
    std::size_t omega_worst_delta=0, k_worst_delta=0;
    double omega_physical_linf=0.,omega_carrier_linf=0.;
    double k_physical_linf=0.,k_carrier_linf=0.;
    double omega_delta_linf=0.,k_delta_linf=0.;
    double omega_swap_relative_effect=0.; // positive = carrier phi reduces omega raw norm
    double omega_physical_worst_signed=0.,omega_carrier_worst_signed=0.;
    double omega_physical_cell0=0.,omega_carrier_cell0=0.;
    double omega_physical_cell12=0.,omega_carrier_cell12=0.;
    double face_physical=0.,face_carrier=0.;
    std::array<double,2> face_omega_delta{}; // owner, neighbour (if internal)
    double max_face_identity_gap=0.;
    double boundary_only_flux_delta=0.; // signed net boundary flux difference
    bool exact_single_state=true, convective_identity=true;
};

// Independent per-face FV advection identity, including the same bounded
// upwind selection and the -div(phi)*psi owner-row term as the SST assembler.
// At flux sign changes the upwind state is selected independently for BOTH
// fluxes. Never approximate delta(phi*psi) by delta(phi)*psi(old sign).
inline BfsSstFluxSwapV68 analyze_bfs_sst_flux_swap_v68(
    const cfdx::core::Mesh& mesh,const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& carrier,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& physical,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const ScalarBoundaryConditions& k_bc,
    const ScalarBoundaryConditions& omega_bc,
    const ScalarBoundaryFaceValues* omega_faces,
    const SSTWallProductionV17* wall_production,
    const SSTCoupledResidualV22& carrier_residual,
    const SSTCoupledResidualV22& physical_residual,
    std::size_t focus_face)
{
    const auto n=mesh.n_cells(),nf=mesh.n_faces();
    if (!n || !nf || k.size()!=n || omega.size()!=n ||
        carrier.size()!=nf || physical.size()!=nf || focus_face>=nf ||
        carrier_residual.values.size()!=2*n ||
        physical_residual.values.size()!=2*n)
        throw std::invalid_argument("V68 SST swap: incomplete identical-state inputs");
    if(wall_production)wall_production->validate(n);
    using ScalarField=cfdx::core::Field<double,cfdx::core::Location::CELL>;
    std::array<const ScalarField*,2> fields{{&k,&omega}};
    std::array<const ScalarBoundaryConditions*,2> bcs{{&k_bc,&omega_bc}};
    std::array<const ScalarBoundaryFaceValues*,2> face_bcs{{nullptr,omega_faces}};
    std::vector<long double> predicted(2*n,0.L);
    BfsSstFluxSwapV68 out;
    out.cells=n;out.faces=nf;out.focus_face=focus_face;
    for(std::size_t f=0;f<nf;++f) {
        const auto owner=mesh.ownership().owner(f);
        const auto neighbour=mesh.ownership().neighbour(f);
        if(owner>=n || neighbour < -1 ||
           (neighbour>=0 && static_cast<std::size_t>(neighbour)>=n))
            throw std::invalid_argument("V68 invalid face ownership");
        const double fc=carrier(f),fp=physical(f);
        if(!std::isfinite(fc) || !std::isfinite(fp))
            throw std::invalid_argument("V68 nonfinite face flux");
        if (f==focus_face){out.face_carrier=fc;out.face_physical=fp;}
        if(neighbour<0) {
            const auto patch=geometry.face_patch[f];
            if(patch>=mesh.boundary().n_patches())
                throw std::invalid_argument("V68 unassigned boundary face");
            if(mesh.boundary().patch(patch).type==cfdx::core::PatchType::EMPTY)
                continue;
            out.boundary_only_flux_delta+=fp-fc;
        }
        for(std::size_t b=0;b<2;++b) {
            const auto& psi=*fields[b];
            const double po=psi(owner);
            if(!std::isfinite(po))throw std::invalid_argument("V68 nonfinite cell state");
            double vc=po,vp=po;
            if(neighbour>=0){
                const auto pn=psi(static_cast<std::size_t>(neighbour));
                if(!std::isfinite(pn))throw std::invalid_argument("V68 nonfinite neighbour");
                vc=fc>=0. ?po:pn;
                vp=fp>=0. ?po:pn;
            } else {
                const auto bc=resolve_scalar_boundary_condition(
                    mesh,geometry,f,*bcs[b],face_bcs[b]);
                if(bc.type==ScalarBoundaryType::FIXED_VALUE){
                    if(fc<0.)vc=bc.value;
                    if(fp<0.)vp=bc.value;
                } else if(bc.type==ScalarBoundaryType::FIXED_GRADIENT){
                    const auto dist=boundary_normal_distance(geometry,f,owner);
                    vc=vp=po+bc.gradient*dist;
                }
            }
            const auto owner_delta=static_cast<long double>(fp)*(vp-po)-
                                   static_cast<long double>(fc)*(vc-po);
            predicted[b*n+owner]+=owner_delta;
            if(f==focus_face && b==1)out.face_omega_delta[0]=static_cast<double>(owner_delta);
            if(neighbour>=0){
                const auto ni=static_cast<std::size_t>(neighbour);
                const double pn=psi(ni);
                const auto neighbour_delta=static_cast<long double>(fc)*(vc-pn)-
                                           static_cast<long double>(fp)*(vp-pn);
                predicted[b*n+ni]+=neighbour_delta;
                if(f==focus_face && b==1)out.face_omega_delta[1]=static_cast<double>(neighbour_delta);
            }
        }
    }
    for(std::size_t b=0;b<2;++b)for(std::size_t i=0;i<n;++i) {
        const auto j=b*n+i;
        if(b==1 && wall_production && wall_production->has_owner_omega_constraints() &&
           wall_production->has_omega_owner_cell[i])
            predicted[j]=0.;
        const double full=physical_residual.values[j];
        const double base=carrier_residual.values[j];
        const double delta=full-base;
        if(!std::isfinite(full) || !std::isfinite(base))
            throw std::invalid_argument("V68 nonfinite SST residual");
        const auto gap=std::abs(static_cast<double>(predicted[j])-delta);
        const auto scale=std::max({1.,std::abs(full),std::abs(base),
                                   std::abs(static_cast<double>(predicted[j]))});
        out.max_face_identity_gap=std::max(out.max_face_identity_gap,gap/scale);
        if(gap>1e-9*scale)
            throw std::runtime_error("V68 SST swap violates exact convection-only FV identity");
        if(b==0) {
            if(std::abs(full)>out.k_physical_linf)out.k_physical_linf=std::abs(full);
            if(std::abs(base)>out.k_carrier_linf)out.k_carrier_linf=std::abs(base);
            if(std::abs(delta)>out.k_delta_linf){out.k_delta_linf=std::abs(delta);out.k_worst_delta=i;}
        } else {
            if(std::abs(full)>out.omega_physical_linf) {
                out.omega_physical_linf=std::abs(full);out.omega_worst_physical=i;
            }
            if(std::abs(base)>out.omega_carrier_linf){
                out.omega_carrier_linf=std::abs(base);out.omega_worst_carrier=i;
            }
            if(std::abs(delta)>out.omega_delta_linf){
                out.omega_delta_linf=std::abs(delta);out.omega_worst_delta=i;
            }
            if(i==0){out.omega_physical_cell0=full;out.omega_carrier_cell0=base;}
            if(i==12){out.omega_physical_cell12=full;out.omega_carrier_cell12=base;}
        }
    }
    out.omega_physical_worst_signed=
        physical_residual.values[n+out.omega_worst_physical];
    out.omega_carrier_worst_signed=
        carrier_residual.values[n+out.omega_worst_physical];
    out.omega_swap_relative_effect=(out.omega_physical_linf-out.omega_carrier_linf)/
        std::max(out.omega_physical_linf,1e-300);
    return out;
}
}
