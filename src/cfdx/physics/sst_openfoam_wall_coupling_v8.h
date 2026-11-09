#pragma once

// Diagnostic-only independent OpenFOAM Foundation v13 near-wall reference.
// This implements the published Spalding relation and omegaWallFunction's
// local G/omega formula, not their shared wall-cell corner weighting,
// constraint application, or turbulence model solve. It MUST NOT silently
// replace CFDX's low-Re SST model or certify OpenFOAM equivalence.
#include "cfdx/physics/sst_openfoam_wall_reference_v7.h"
#include "cfdx/physics/sst_production_coupling.h"
#include "cfdx/physics/sst_wall_triplet_cell_audit_v13.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct OpenFoamWallCouplingPointV8 {
    double omega = 0.;
    double production = 0.;
    double nut_spalding = 0.;
    double yplus_spalding = 0.;
    double yplus_from_k = 0.;
};

inline double spalding_wall_nut_v8(double speed_tangential, double wall_distance,
                                    double nu, double kappa = 0.41, double E = 9.8)
{
    if (!(speed_tangential >= 0.) || !std::isfinite(speed_tangential) ||
        !(wall_distance > 0.) || !std::isfinite(wall_distance) ||
        !(nu > 0.) || !std::isfinite(nu) ||
        !(kappa > 0.) || !std::isfinite(kappa) ||
        !(E > 0.) || !std::isfinite(E))
        throw std::invalid_argument("Spalding wall: invalid physical inputs");
    if (speed_tangential == 0.) return 0.;
    const auto law = [&](double uplus) {
        const double q = kappa*uplus;
        if (q > 700.) return std::numeric_limits<double>::infinity();
        // Subtracting the first three exponential series terms at q << 1
        // loses all significance. Use the positive fourth-order tail there.
        const double tail = q < .02
            ? (q*q*q*q/24.)*(1. + q/5. + q*q/30. + q*q*q/210.)
            : std::expm1(q) - q - q*q/2. - q*q*q/6.;
        return uplus + std::max(0., tail) / E;
    };
    const auto defect = [&](double ustar) {
        const double uplus = speed_tangential/ustar;
        return wall_distance*ustar/nu - law(uplus);
    };
    // Viscous u_tau is a strict lower bound (Spalding y+ >= u+).
    const double viscous = std::sqrt(nu*speed_tangential/wall_distance);
    double low=viscous, high=2.*viscous;
    if (!(high>0.) || !std::isfinite(high))
        throw std::domain_error("Spalding wall: unrepresentable friction velocity");
    for (int i=0; i<60 && defect(high)<0.; ++i) high*=2.;
    if (!std::isfinite(high) || defect(high)<0.)
        throw std::domain_error("Spalding wall: could not bracket root");
    for (int i=0; i<90; ++i) {
        const double middle=low + 0.5*(high-low);
        if (defect(middle) < 0.) low=middle;
        else high=middle;
    }
    const double friction=low + 0.5*(high-low);
    const double wall_nut=friction*friction*wall_distance/speed_tangential - nu;
    if (!std::isfinite(wall_nut))
        throw std::domain_error("Spalding wall: nonfinite viscosity");
    return std::max(0.,wall_nut);
}

inline OpenFoamWallCouplingPointV8 openfoam_wall_coupling_point_v8(
    double k, double viscosity, double y, double tangential_u,
    double production_viscous,
    ReferenceOmegaWallBlend blend = ReferenceOmegaWallBlend::Stepwise,
    double beta1 = 0.075, double kappa = 0.41, double E = 9.8)
{
    if (!std::isfinite(production_viscous) || production_viscous < 0.)
        throw std::invalid_argument("OpenFOAM wall G: invalid production");
    const double grad=tangential_u/y;
    const auto omega=openfoam_omega_wall_reference(
        k,viscosity,y,grad,blend,.09,beta1,kappa);
    const double cmu5=std::sqrt(.09);
    const double ustar = blend==ReferenceOmegaWallBlend::Stepwise
        ? std::sqrt(cmu5*k)
        : std::sqrt(omega.laminar_fraction*viscosity*grad +
                     (1.-omega.laminar_fraction)*cmu5*k);
    double gturb=0.;
    if (omega.y_plus_from_k > 0.) {
        const double uplus = std::log(E*omega.y_plus_from_k)/kappa;
        if (uplus>0. && std::isfinite(uplus)) {
            const double temp=ustar*grad*y/uplus;
            gturb=temp*temp/(viscosity*kappa*omega.y_plus_from_k);
        }
    }
    const double g=omega.laminar_fraction*production_viscous+
        (1.-omega.laminar_fraction)*gturb;
    OpenFoamWallCouplingPointV8 result{
        omega.reference_omega,g,spalding_wall_nut_v8(tangential_u,y,viscosity,kappa,E),
        0.,omega.y_plus_from_k};
    if (!std::isfinite(g))
        throw std::domain_error("OpenFOAM wall G: nonfinite reference");
    // A Spalding viscosity defines a velocity-based friction velocity.
    if (tangential_u>0.) {
        const double friction=std::sqrt(
            (result.nut_spalding+viscosity)*tangential_u/y);
        result.yplus_spalding=y*friction/viscosity;
    }
    return result;
}

struct OpenFoamWallCouplingAuditV8 {
    std::size_t faces=0;
    double max_omega_owner_gap=0.;
    double max_nut_wall_gap=0.;
    double max_production_gap=0.;
    double max_spalding_yplus=0.;
    double max_k_yplus=0.;
    std::size_t worst_face=0;
    // Area-weighted MULTI-FACE wall-cell reference is diagnostic only.
    // OpenFOAM corner weights and patch fields are not yet implemented.
    std::size_t wall_owner_cells=0;
    std::size_t multi_face_wall_owner_cells=0;
    double weighted_omega_owner_gap=0.;
    double weighted_nut_owner_proxy_gap=0.;
    double weighted_G_owner_proxy_gap=0.;
    // This MUST remain false until the shared wall-cell constraints,
    // nut patch field, and G replacement are integrated and validated.
    bool openfoam_equivalent=false;
};

inline OpenFoamWallCouplingAuditV8 audit_openfoam_wall_coupling_v8(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const std::vector<std::string>& wall_patches,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& U,
    const SSTProductionState& state,
    const std::map<std::string,std::array<double,3>>& fixed_wall_velocity,
    const TurbulenceTransportControls& c)
{
    if (U.size()!=mesh.n_cells() || U.dimension()!=3)
        throw std::invalid_argument("V8 wall audit: U shape mismatch");
    validate_sst_production_state(state,mesh.n_cells());
    OpenFoamWallCouplingAuditV8 report;
    std::vector<SstWallFaceTripletV13> wall_samples;
    for (const auto& name:wall_patches) {
        const auto it=fixed_wall_velocity.find(name);
        if (it==fixed_wall_velocity.end())
            throw std::invalid_argument("V8 wall audit: wall velocity missing");
        const auto& patch=mesh.boundary().patch(mesh.boundary().find(name));
        if (patch.type!=cfdx::core::PatchType::WALL)
            throw std::invalid_argument("V8 wall audit: non-wall patch");
        for (auto face:patch.face_ids) {
            if (mesh.ownership().neighbour(face)>=0)
                throw std::invalid_argument("V8 wall audit: nonboundary face");
            const auto owner=mesh.ownership().owner(face);
            const auto Sf=geometry.face_area_vectors.at(face);
            const double area=Sf.mag();
            if (!(area>0.)) throw std::invalid_argument("V8 wall audit: zero area");
            const auto normal=Sf*(1./area);
            const double y=(geometry.face_centres.at(face)-
                geometry.cell_centres.at(owner)).dot(normal);
            if (!(y>0.) || !std::isfinite(y))
                throw std::invalid_argument("V8 wall audit: invalid wall spacing");
            const cfdx::core::Vec3 jump{
                U.component_data(0)[owner]-it->second[0],
                U.component_data(1)[owner]-it->second[1],
                U.component_data(2)[owner]-it->second[2]};
            const double tangential=(jump-normal*jump.dot(normal)).mag();
            // The Gvf input is a CFDX cell-production *proxy*; OpenFOAM's
            // internal G field and corner-weighted wall constraints are not
            // identical to this unweighted face-local reference.
            const double k=state.k(owner), omega=state.omega(owner);
            const double g_local=std::min(c.density*state.nu_t(owner)*
                    state.strain_rate(owner)*state.strain_rate(owner),
                c.sst_production_limiter*c.beta_star*c.density*k*omega);
            const auto ref=openfoam_wall_coupling_point_v8(
                k,c.molecular_viscosity,y,tangential,g_local);
            wall_samples.push_back({owner,face,area,ref.omega,ref.nut_spalding,
                ref.production,omega,state.nu_t(owner),g_local});
            const double omega_gap=std::abs(omega-ref.omega)/
                std::max(omega,ref.omega);
            // Current SST momentum uses molecular-only WALL diffusion.
            const double nut_gap=ref.nut_spalding/
                std::max(c.molecular_viscosity,ref.nut_spalding);
            // Both numerator and scale are production terms. The molecular
            // kinematic viscosity has different physical units and must
            // NEVER be used as an arbitrary denominator for the G deficit.
            const double g_gap=std::abs(g_local-ref.production)/
                std::max({std::numeric_limits<double>::min(),
                          std::abs(g_local),std::abs(ref.production)});
            if (std::max({omega_gap,nut_gap,g_gap}) >
                std::max({report.max_omega_owner_gap,report.max_nut_wall_gap,
                          report.max_production_gap})) report.worst_face=face;
            report.max_omega_owner_gap=std::max(report.max_omega_owner_gap,omega_gap);
            report.max_nut_wall_gap=std::max(report.max_nut_wall_gap,nut_gap);
            report.max_production_gap=std::max(report.max_production_gap,g_gap);
            report.max_spalding_yplus=std::max(report.max_spalding_yplus,ref.yplus_spalding);
            report.max_k_yplus=std::max(report.max_k_yplus,ref.yplus_from_k);
            ++report.faces;
        }
    }
    if (report.faces==0)
        throw std::invalid_argument("V8 wall audit: no wall faces");
    const auto cells = audit_sst_wall_triplet_cells_v13(wall_samples);
    report.wall_owner_cells=cells.size();
    for (const auto& sample:cells) {
        report.multi_face_wall_owner_cells += (sample.face_count>1);
        report.weighted_omega_owner_gap = std::max(
            report.weighted_omega_owner_gap,sample.omega_owner_gap);
        report.weighted_nut_owner_proxy_gap = std::max(
            report.weighted_nut_owner_proxy_gap,sample.nut_owner_proxy_gap);
        report.weighted_G_owner_proxy_gap = std::max(
            report.weighted_G_owner_proxy_gap,sample.G_owner_proxy_gap);
    }
    return report;
}

} // namespace cfdx::physics
