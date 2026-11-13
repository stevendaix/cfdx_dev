#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/radiation.h"
#include "cfdx/physics/energy_solver.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <limits>
#include <optional>
#include <map>

namespace cfdx::physics {

    double absorption = 0.0;
    double scattering = 0.0;
    double intensity_relaxation = 1.0;
    std::size_t max_iterations = 100;
    double tolerance = 1e-8;
    std::size_t linear_max_iterations = 2000;
    double linear_tolerance = 1e-10;
};

struct RadiationIteration {
    std::size_t iteration = 0;
    double intensity_residual = 0.0;
    double energy_source_residual = 0.0;
};

struct RadiationSolveResult {
    bool converged = false;
    std::size_t iterations = 0;
    std::vector<RadiationIteration> history;
    // V41 opt-in independent wall-flux diagnostic, computed from the FINAL
    // accepted directional intensities. Positive = radiant power leaving gas.
    std::map<std::string,double> wall_outward_W;
    // V42 opt-in: defect of the FINAL nonlinear DOM FV operator, not a
    // difference between successive Picard iterates. Empty in legacy mode.
    bool accepted_state_verified_v42 = false;
    double max_cell_transport_residual_v42 =
        std::numeric_limits<double>::quiet_NaN();
    double relative_energy_closure_v42 =
        std::numeric_limits<double>::quiet_NaN();
};

inline void validate_radiation_transport_controls(const RadiationTransportControls& c)
    if(c.absorption<0.0 || c.scattering<0.0 ||
       c.max_iterations==0 || c.tolerance<=0.0 ||
       c.linear_max_iterations==0 || c.linear_tolerance<=0.0 ||
       !(c.intensity_relaxation>0.0 && c.intensity_relaxation<=1.0))
        throw std::invalid_argument("invalid radiation transport controls");
}

struct P1RadiationSolveResult {
    bool converged = false;
    std::size_t iterations = 0;
    double residual = 0.0;
};

// Solve the gray isotropic P1 equation
// div(D grad G) - kappa_a G + 4 kappa_a sigma T^4 = 0,
// with D = 1/[3(kappa_a+kappa_s)]. The supplied scalar BCs are
// boundary conditions on irradiation G [W/m2].
inline P1RadiationSolveResult solve_p1_radiation(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& irradiation,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& radiation_source,
    double absorption,
    double scattering,
    const ScalarBoundaryConditions& irradiation_bcs = {},
    std::size_t max_iterations = 2000,
    double tolerance = 1e-10)
{
    if(!std::isfinite(absorption) || !std::isfinite(scattering) ||
       absorption < 0.0 || scattering < 0.0 ||
       absorption + scattering <= 0.0 ||
       max_iterations == 0 || !std::isfinite(tolerance) || tolerance <= 0.0 ||
       temperature.size()!=mesh.n_cells() ||
       irradiation.size()!=mesh.n_cells() ||
       radiation_source.size()!=mesh.n_cells())
        throw std::invalid_argument("invalid P1 radiation solve inputs");

    const double D = 1.0/(3.0*(absorption+scattering));
    cfdx::core::Field<double,cfdx::core::Location::FACE> zero_flux(
        mesh.n_faces(),"p1_flux","m2/s",1);
    zero_flux.fill(0.0);
    cfdx::core::Field<double,cfdx::core::Location::CELL> su(
        mesh.n_cells(),"p1_source","W/m3",1);
    cfdx::core::Field<double,cfdx::core::Location::CELL> sp(
        mesh.n_cells(),"p1_sp","1/m",1);
    for(std::size_t cell=0; cell<mesh.n_cells(); ++cell) {
        const double T=temperature(cell);
        if(!std::isfinite(T) || T<0.0)
            throw std::invalid_argument("P1 temperature must be finite and non-negative");
        su(cell)=4.0*absorption*blackbody_emissive_power(T);
        sp(cell)=-absorption;
    }

    auto eq=assemble_scalar_equation(
        mesh,geometry,zero_flux,D,su,sp,irradiation_bcs,true,
        nullptr,nullptr,nullptr,nullptr,ConvectionScheme::UPWIND,&irradiation);
    cfdx::core::Vector G(mesh.n_cells(),0.0);
    for(std::size_t cell=0; cell<mesh.n_cells(); ++cell)
        G(cell)=irradiation(cell);
    ScalarSolveControls controls;
    controls.max_iterations=max_iterations;
    controls.tolerance=tolerance;
    controls.relaxation=1.0;
    const auto linear=solve_scalar_equation(eq,G,controls);

    P1RadiationSolveResult result;
    result.converged=(linear.status==cfdx::core::SolverStatus::CONVERGED);
    result.iterations=linear.iterations;
    result.residual=linear.residual;
    if(!result.converged)
        return result;

    for(std::size_t cell=0; cell<mesh.n_cells(); ++cell) {
        if(!std::isfinite(G(cell)) || G(cell)<-1e-10)
            throw std::runtime_error("P1 produced non-physical irradiation");
        irradiation(cell)=std::max(0.0,G(cell));
        radiation_source(cell)=absorption*
            (4.0*M_PI*blackbody_intensity(temperature(cell))-irradiation(cell));
    }
    return result;
}

struct DomWallBoundaryCondition {
    double emissivity = 1.0;
    double temperature = 300.0;
};

using DomWallBoundaryConditions = std::map<std::string, DomWallBoundaryCondition>;

// Direction-aware diffuse-gray DOM wall treatment. For each boundary face,
// incoming ordinates are set from the emitted blackbody intensity plus the
// diffusely reflected outgoing irradiation:
// I_in = epsilon I_b(T_w) + (1-epsilon) G_out/H_in,
// G_out = sum_{m: s_m.n>0} w_m I_m (s_m.n),
// H_in = sum_{m: s_m.n<0} w_m |s_m.n|.
// Important: for a finite angular quadrature H_in is generally not pi.
// Dividing by the CONTINUOUS hemispherical integral pi violates the
// isothermal diffuse-gray equilibrium and discrete wall-power conservation.
// Outgoing ordinates use zero-gradient/extrapolation from the owner cell.
inline void validate_dom_wall_conditions(
    const cfdx::core::Mesh& mesh,
{
    for (std::size_t p=0; p<mesh.boundary().n_patches(); ++p) {
        const auto& patch=mesh.boundary().patch(p);
        const auto it=walls.find(patch.name);
        if (it==walls.end())
            throw std::invalid_argument(
                "DOM diffuse-gray walls: missing condition for patch " + patch.name);
        const auto& bc=it->second;
        if (!std::isfinite(bc.emissivity) || bc.emissivity<0.0 ||
            bc.emissivity>1.0 || !std::isfinite(bc.temperature) ||
            bc.temperature<0.0)
            throw std::invalid_argument(
                "DOM diffuse-gray walls: invalid emissivity/temperature on patch " +
                patch.name);
    }
}

// V38: cache the same outgoing diffuse reflection used by the V37 wall
// operator. Rebuilt on each outer iteration and updated after every accepted
// angular solve so the historical Gauss-Seidel order is unchanged.
// H_in/H_out remain the DISCRETE angular measures from V37; never replace
// either by pi without a separate quadrature qualification.
struct DomDiffuseWallCache {
    std::vector<long double> outgoing_irradiation;
    std::vector<long double> incoming_measure;
    std::vector<long double> outgoing_measure;
    std::size_t direction_count=0;
};

inline DomDiffuseWallCache make_dom_diffuse_wall_cache(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const std::vector<DiscreteDirection>& directions,
    const std::vector<cfdx::core::Field<double,cfdx::core::Location::CELL>>& intensities)
{
    if (directions.empty() || intensities.size()!=directions.size() ||
        geometry.face_area_vectors.size()!=mesh.n_faces())
        throw std::invalid_argument("DOM cache: dimensions inconsistent");
    DomDiffuseWallCache cache;
    cache.direction_count=directions.size();
    cache.outgoing_irradiation.assign(mesh.n_faces(),0.0L);
    cache.incoming_measure.assign(mesh.n_faces(),0.0L);
    cache.outgoing_measure.assign(mesh.n_faces(),0.0L);
    for (const auto& intensity:intensities)
        if (intensity.size()!=mesh.n_cells())
            throw std::invalid_argument("DOM cache: intensity size mismatch");
    for (std::size_t p=0;p<mesh.boundary().n_patches();++p)
        for (const auto f:mesh.boundary().patch(p).face_ids) {
            const auto Sf=geometry.face_area_vectors[f];
            const double area=Sf.mag();
            if (!(area>0.0) || !std::isfinite(area))
                throw std::runtime_error("DOM cache: invalid wall area");
            const std::size_t owner=mesh.ownership().owner(f);
            if (owner>=mesh.n_cells())
                throw std::invalid_argument("DOM cache: invalid owner");
            for (std::size_t m=0;m<directions.size();++m) {
                const auto& d=directions[m];
                const double mu=(d.dx*Sf.x+d.dy*Sf.y+d.dz*Sf.z)/area;
                if (!std::isfinite(mu) || !std::isfinite(d.weight) || d.weight<0.0)
                    throw std::invalid_argument("DOM cache: invalid quadrature");
                if (mu>0.0) {
                    const double I=intensities[m](owner);
                    if (!(I>=0.0) || !std::isfinite(I))
                        throw std::invalid_argument("DOM cache: nonphysical outgoing intensity");
                    const long double coefficient=static_cast<long double>(d.weight)*mu;
                    cache.outgoing_measure[f]+=coefficient;
                    cache.outgoing_irradiation[f]+=coefficient*I;
                } else if (mu<0.0) {
                    cache.incoming_measure[f]-=static_cast<long double>(d.weight)*mu;
                }
            }
            const long double Hminus=cache.incoming_measure[f];
            const long double Hplus=cache.outgoing_measure[f];
            if (!(Hminus>0.0L) || !(Hplus>0.0L) ||
                !std::isfinite(Hminus) || !std::isfinite(Hplus) ||
                !std::isfinite(cache.outgoing_irradiation[f]) ||
                std::abs(Hminus-Hplus)>1e-10L*std::max(Hminus,Hplus))
                throw std::invalid_argument("DOM cache: unbalanced wall hemispheres");
        }
    return cache;
}

// Delta is the ACCEPTED change in one angle's cell intensities. The caller
// must not update the cache for a trial that is later rolled back.
inline void update_dom_diffuse_wall_cache(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const DiscreteDirection& direction,
    const std::vector<double>& accepted_intensity_delta,
    DomDiffuseWallCache& cache)
{
    if (accepted_intensity_delta.size()!=mesh.n_cells() ||
        cache.outgoing_irradiation.size()!=mesh.n_faces() ||
        cache.incoming_measure.size()!=mesh.n_faces() ||
        cache.outgoing_measure.size()!=mesh.n_faces())
        throw std::invalid_argument("DOM cache update: dimensions inconsistent");
    for (std::size_t p=0;p<mesh.boundary().n_patches();++p)
        for (const auto f:mesh.boundary().patch(p).face_ids) {
            const auto Sf=geometry.face_area_vectors[f];
            const double area=Sf.mag();
            if (!(area>0.0) || !std::isfinite(area))
                throw std::runtime_error("DOM cache update: invalid wall area");
            const double mu=(direction.dx*Sf.x+direction.dy*Sf.y+
                             direction.dz*Sf.z)/area;
            if (mu>0.0) {
                const double delta=accepted_intensity_delta[mesh.ownership().owner(f)];
                if (!std::isfinite(delta))
                    throw std::invalid_argument("DOM cache update: invalid intensity delta");
                cache.outgoing_irradiation[f]+=
                    static_cast<long double>(direction.weight)*mu*delta;
                if (!std::isfinite(cache.outgoing_irradiation[f]))
                    throw std::runtime_error("DOM cache update: nonfinite irradiance");
            }
        }
}

// Independently report whether a quadrature integrates a hemispherical flux
// close to pi: global 0th/1st/2nd sphere moments do NOT imply this property.
struct DomHemisphereMoments {
    double outgoing=0.0, incoming=0.0, relative_pi_error=0.0;
};
inline DomHemisphereMoments dom_hemisphere_moments(
    const std::vector<DiscreteDirection>& directions,
    const std::array<double,3>& outward_normal)
{
    validate_discrete_directions(directions);
    const double norm=std::hypot(outward_normal[0],outward_normal[1],outward_normal[2]);
    if (!(norm>0.0) || !std::isfinite(norm))
        throw std::invalid_argument("DOM hemisphere: invalid normal");
    long double Hp=0.0L,Hm=0.0L;
    for (const auto& d:directions) {
        const double mu=(d.dx*outward_normal[0]+d.dy*outward_normal[1]+
                         d.dz*outward_normal[2])/norm;
        if (mu>0) Hp+=static_cast<long double>(d.weight)*mu;
        else if (mu<0) Hm-=static_cast<long double>(d.weight)*mu;
    }
    if (!(Hp>0.0L) || !(Hm>0.0L) ||
        std::abs(Hp-Hm)>1e-10L*std::max(Hp,Hm))
        throw std::invalid_argument("DOM hemisphere: unbalanced angular measure");
    return {static_cast<double>(Hp),static_cast<double>(Hm),
        std::abs(static_cast<double>(Hp)/M_PI-1.0)};
}

// The discrete-isothermal balance (H+ == H-) is not the same as the
// physical hemispherical condition (H+ == pi). A separate preflight reports
// the maximum error over actual 3D wall-face normals. It never silently
// rescales physical emission, reflection or the angular quadrature.
struct DomWallQuadratureCertificate {
    std::size_t faces_checked=0, worst_face=0;
    double maximum_relative_pi_error=0.0;
};
inline DomWallQuadratureCertificate audit_dom_wall_quadrature(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const std::vector<DiscreteDirection>& directions)
{
    if (geometry.face_area_vectors.size()!=mesh.n_faces())
        throw std::invalid_argument("DOM quadrature audit: geometry size mismatch");
    DomWallQuadratureCertificate result;
    for (std::size_t p=0;p<mesh.boundary().n_patches();++p)
        for (const auto f:mesh.boundary().patch(p).face_ids) {
            const auto Sf=geometry.face_area_vectors[f];
            const auto m=dom_hemisphere_moments(directions,{Sf.x,Sf.y,Sf.z});
            ++result.faces_checked;
            if (m.relative_pi_error>result.maximum_relative_pi_error) {
                result.maximum_relative_pi_error=m.relative_pi_error;
                result.worst_face=f;
            }
        }
    if (result.faces_checked==0)
        throw std::invalid_argument("DOM quadrature audit: no boundary faces");
    return result;
}
inline DomWallQuadratureCertificate require_dom_wall_quadrature_accuracy(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const std::vector<DiscreteDirection>& directions,
    double maximum_allowed_relative_pi_error)
{
    if (!(maximum_allowed_relative_pi_error>=0.0) ||
        !std::isfinite(maximum_allowed_relative_pi_error))
        throw std::invalid_argument("DOM quadrature audit: invalid tolerance");
    auto cert=audit_dom_wall_quadrature(mesh,geometry,directions);
    if (cert.maximum_relative_pi_error>maximum_allowed_relative_pi_error)
        throw std::runtime_error("DOM quadrature fails physical hemisphere moment gate");
    return cert;
}

inline ScalarBoundaryFaceValues build_dom_diffuse_gray_face_values(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const std::vector<DiscreteDirection>& directions,
    const std::vector<cfdx::core::Field<double,cfdx::core::Location::CELL>>& intensities,
    const DomWallBoundaryConditions& walls,
    std::size_t direction_index,
    const DomDiffuseWallCache* cache = nullptr)
{
    if (cache && (cache->direction_count!=directions.size() ||
        cache->outgoing_irradiation.size()!=mesh.n_faces() ||
        cache->incoming_measure.size()!=mesh.n_faces() ||
        cache->outgoing_measure.size()!=mesh.n_faces()))
        throw std::invalid_argument("DOM wall operator: cache size mismatch");
    if (intensities.size()!=directions.size() || direction_index>=directions.size())
        throw std::invalid_argument("DOM wall operator: direction/intensity size mismatch");
    const std::size_t nf=mesh.n_faces();
    const auto& target_direction=directions[direction_index];
    for (std::size_t p=0; p<mesh.boundary().n_patches(); ++p) {
        const auto& patch=mesh.boundary().patch(p);
        const auto& wall=walls.at(patch.name);
        auto& pv=values.values[patch.name];
        pv.assign(nf,std::numeric_limits<double>::quiet_NaN());
        const double Ib=blackbody_intensity(wall.temperature);
        for (const auto f:patch.face_ids) {
            const auto Sf=geometry.face_area_vectors[f];
            const double area=Sf.mag();
            if (!(area>0.0) || !std::isfinite(area))
                throw std::runtime_error("DOM wall operator: invalid boundary face area");
            const double nx=Sf.x/area, ny=Sf.y/area, nz=Sf.z/area;
            const double target_mu=target_direction.dx*nx+
                                  target_direction.dy*ny+
                                  target_direction.dz*nz;
            if (target_mu>=0.0)
                continue; // outgoing ordinate: zero-gradient/extrapolated

            long double Gout=cache ? cache->outgoing_irradiation[f] : 0.0L;
            long double incoming_measure=cache ? cache->incoming_measure[f] : 0.0L;
            long double outgoing_measure=cache ? cache->outgoing_measure[f] : 0.0L;
            if (!cache) for (std::size_t m=0; m<directions.size(); ++m) {
                const auto& d=directions[m];
                const double mu=d.dx*nx+d.dy*ny+d.dz*nz;
                if (mu>0.0) {
                    Gout += static_cast<long double>(d.weight)*
                        intensities[m](mesh.ownership().owner(f))*mu;
                    outgoing_measure+=static_cast<long double>(d.weight)*mu;
                } else if (mu<0.0) {
                    incoming_measure-=static_cast<long double>(d.weight)*mu;
                }
            }
            if (!(incoming_measure>0.0L) || !(outgoing_measure>0.0L) ||
                !std::isfinite(incoming_measure) || !std::isfinite(outgoing_measure) ||
                !std::isfinite(Gout) || Gout < -1e-12L ||
                std::abs(incoming_measure-outgoing_measure) >
                    1e-10L*std::max(incoming_measure,outgoing_measure))
                throw std::invalid_argument(
                    "DOM wall quadrature has no balanced incoming/outgoing hemisphere");
            pv[f]=wall.emissivity*Ib +
                  (1.0-wall.emissivity)*static_cast<double>(Gout/incoming_measure);
        }
    }
    return values;
}

// V42 independent accepted-state audit of the isotropic-scattering DOM
// equation: div(s I) - I div(s) = ka Ib + ks G/(4*pi) - (ka+ks) I.
// Fluxes are reconstructed face by face (UPWIND) from the ACCEPTED field;
// diffuse wall reflection is rebuilt from all accepted ordinates. The
// bounded-convection term is included to match the production FV contract.
// This checks the DISCRETE transport balance; a separate quadrature oracle is
// still mandatory before claiming physical hemisphere accuracy.
struct DomAcceptedStateDefectV42 {
    double max_cell_relative=0.0;
    double max_cell_absolute=0.0;
    double relative_global_energy=0.0;
    double wall_outward_W=0.0;
    double gas_emission_minus_absorption_W=0.0;
    std::size_t worst_direction=0, worst_cell=0;
};
inline DomAcceptedStateDefectV42 audit_dom_accepted_state_v42(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    const std::vector<cfdx::core::Field<double,cfdx::core::Location::CELL>>& intensities,
    const std::vector<DiscreteDirection>& directions,
    const RadiationTransportControls& controls,
    const DomWallBoundaryConditions& walls,
    std::vector<double>* directional_cell_defects=nullptr)
{
    validate_radiation_transport_controls(controls);
    validate_discrete_directions(directions);
    validate_dom_wall_conditions(mesh,walls);
    const std::size_t nc=mesh.n_cells(), nf=mesh.n_faces();
    if (temperature.size()!=nc || intensities.size()!=directions.size() ||
        geometry.cell_volumes.size()!=nc ||
        geometry.face_area_vectors.size()!=nf ||
        geometry.face_patch.size()!=nf)
        throw std::invalid_argument("DOM V42 audit: incompatible field/geometry sizes");
    std::vector<long double> G(nc,0.0L);
    for(std::size_t m=0;m<directions.size();++m) {
        if (intensities[m].size()!=nc)
            throw std::invalid_argument("DOM V42 audit: invalid angle field size");
        for(std::size_t c=0;c<nc;++c) {
            const double I=intensities[m](c);
            if (!(I>=0.0) || !std::isfinite(I))
                throw std::invalid_argument("DOM V42 audit: negative or nonfinite intensity");
            G[c]+=static_cast<long double>(directions[m].weight)*I;
        }
    }
    DomAcceptedStateDefectV42 audit;
    std::vector<double> local_defects;
    if (directional_cell_defects)
        local_defects.assign(nc*directions.size(),0.0);
    long double total_wall=0.0L;
    long double total_gas=0.0L, gas_l1=0.0L, emission_absorption=0.0L;
    for(std::size_t c=0;c<nc;++c) {
        if (!(geometry.cell_volumes[c]>0.0) ||
            !std::isfinite(geometry.cell_volumes[c]))
            throw std::invalid_argument("DOM V42 audit: invalid cell volume");
        const double Ib=blackbody_intensity(temperature(c));
        const long double power=static_cast<long double>(geometry.cell_volumes[c])*
            controls.absorption*(4.0L*M_PI*Ib-G[c]);
        total_gas+=power;
        gas_l1+=std::abs(power);
        emission_absorption+=static_cast<long double>(geometry.cell_volumes[c])*
            controls.absorption*(4.0L*M_PI*Ib+G[c]);
    }
    const auto cache=make_dom_diffuse_wall_cache(mesh,geometry,directions,intensities);
    const auto& own=mesh.ownership();
    for(std::size_t m=0;m<directions.size();++m) {
        const auto& d=directions[m];
        const auto wall_values=build_dom_diffuse_gray_face_values(
            mesh,geometry,directions,intensities,walls,m,&cache);
        std::vector<long double> divergence(nc,0.0L), div_s(nc,0.0L);
        std::vector<long double> face_l1(nc,0.0L);
        for(std::size_t f=0;f<nf;++f) {
            const auto patch=geometry.face_patch[f];
            const auto neighbour=own.neighbour(f);
            if (neighbour<0 && patch<mesh.boundary().n_patches() &&
                mesh.boundary().patch(patch).type==cfdx::core::PatchType::EMPTY)
                continue;
            const auto Sf=geometry.face_area_vectors[f];
            const long double F=static_cast<long double>(d.dx)*Sf.x+
                static_cast<long double>(d.dy)*Sf.y+
                static_cast<long double>(d.dz)*Sf.z;
            if (!std::isfinite(F))
                throw std::invalid_argument("DOM V42 audit: nonfinite direction flux");
            const auto o=own.owner(f);
            if (o>=nc) throw std::invalid_argument("DOM V42 audit: invalid owner");
            double I_face=intensities[m](o);
            if (neighbour>=0) {
                const auto n=static_cast<std::size_t>(neighbour);
                if (n>=nc) throw std::invalid_argument("DOM V42 audit: invalid neighbour");
                if (F<0.0L) I_face=intensities[m](n);
            } else if (F<0.0L) {
                if (patch>=mesh.boundary().n_patches())
                    throw std::invalid_argument("DOM V42 audit: unassigned boundary face");
                const auto& name=mesh.boundary().patch(patch).name;
                I_face=wall_values.values.at(name).at(f);
                if (!(I_face>=0.0) || !std::isfinite(I_face))
                    throw std::invalid_argument("DOM V42 audit: invalid incoming wall intensity");
            }
            const long double power=F*I_face;
            divergence[o]+=power;div_s[o]+=F;face_l1[o]+=std::abs(power);
            if(neighbour>=0) {
                const auto n=static_cast<std::size_t>(neighbour);
                divergence[n]-=power;div_s[n]-=F;face_l1[n]+=std::abs(power);
            } else {
                const long double weighted=static_cast<long double>(d.weight)*power;
                total_wall+=weighted;
            }
        }
        for(std::size_t c=0;c<nc;++c) {
            const long double I=intensities[m](c);
            const long double Ib=blackbody_intensity(temperature(c));
            const long double V=geometry.cell_volumes[c];
            const long double emission=V*controls.absorption*Ib;
            const long double scattering=V*controls.scattering*G[c]/(4.0L*M_PI);
            const long double extinction=V*(controls.absorption+controls.scattering)*I;
            const long double defect=divergence[c]-I*div_s[c]-
                emission-scattering+extinction;
            const long double scale=std::max(1e-25L,face_l1[c]+
                std::abs(I*div_s[c])+std::abs(emission)+
                std::abs(scattering)+std::abs(extinction));
            if (directional_cell_defects)
                local_defects[m*nc+c]=static_cast<double>(defect);
            const double relative=static_cast<double>(std::abs(defect)/scale);
            if(!std::isfinite(relative))
                throw std::runtime_error("DOM V42 audit: nonfinite cell defect");
            if (relative>audit.max_cell_relative) {
                audit.max_cell_relative=relative;
                audit.worst_cell=c;audit.worst_direction=m;
            }
            audit.max_cell_absolute=std::max(audit.max_cell_absolute,
                static_cast<double>(std::abs(defect)));
        }
    }
    audit.wall_outward_W=static_cast<double>(total_wall);
    audit.gas_emission_minus_absorption_W=static_cast<double>(total_gas);
    audit.relative_global_energy=static_cast<double>(
        std::abs(total_wall-total_gas)/
        // Normalize by actual emitting/absorbing powers, NOT by the much
        // larger streaming gross wall intensity that cancels at equilibrium.
        std::max(1e-20L,gas_l1+std::abs(total_wall)+emission_absorption));
    if (!std::isfinite(audit.relative_global_energy))
        throw std::runtime_error("DOM V42 audit: nonfinite energy closure");
    if (directional_cell_defects)
        directional_cell_defects->swap(local_defects);
    return audit;
}

// DOM solve with physically coupled diffuse-gray wall boundary conditions.
// The wall reflection is lagged one transport iteration (Picard), so the
// directional linear systems remain independent and strictly diagonally
inline RadiationSolveResult solve_participating_radiation_diffuse_gray_walls(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& irradiation,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& radiation_source,
    const std::vector<DiscreteDirection>& directions,
    const RadiationTransportControls& controls,
    const DomWallBoundaryConditions& walls,
    bool use_wall_cache = true,
    bool compute_wall_power_v41 = false,
    bool verify_accepted_state_v42 = false)
{
    validate_radiation_transport_controls(controls);
    validate_discrete_directions(directions);
    if (temperature.size()!=mesh.n_cells() ||
        irradiation.size()!=mesh.n_cells() ||
        radiation_source.size()!=mesh.n_cells())
        throw std::invalid_argument("radiation field size mismatch");

    const std::size_t nc=mesh.n_cells();
    // Keep the caller state untouched until ALL strict criteria pass. The
    // legacy path performs no field copies and retains its prior behaviour.
    using CellField=cfdx::core::Field<double,cfdx::core::Location::CELL>;
    std::optional<CellField> staged_G, staged_source;
    if (verify_accepted_state_v42) {
        staged_G.emplace(irradiation);
        staged_source.emplace(radiation_source);
        for(std::size_t c=0;c<nc;++c)
            if (!std::isfinite(irradiation(c)) ||
                !std::isfinite(radiation_source(c)))
                throw std::invalid_argument("DOM V42: nonfinite initial field");
    }
    auto& active_G=verify_accepted_state_v42 ? *staged_G : irradiation;
    auto& active_source=verify_accepted_state_v42 ? *staged_source : radiation_source;
    std::vector<cfdx::core::Field<double,cfdx::core::Location::CELL>> intensities;
    intensities.reserve(directions.size());
    for (std::size_t m=0;m<directions.size();++m) {
        intensities.back().fill(0.0);
    }

    ScalarBoundaryConditions extrapolation_bcs;
    for (std::size_t p=0;p<mesh.boundary().n_patches();++p)
        extrapolation_bcs[mesh.boundary().patch(p).name]=
            {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

    RadiationSolveResult result;
    std::vector<double> accepted_delta(nc,0.0);
    for (std::size_t iter=1;iter<=controls.max_iterations;++iter) {
        auto old_source=active_source;
        DomDiffuseWallCache wall_cache;
        if (use_wall_cache)
            wall_cache=make_dom_diffuse_wall_cache(
                mesh,geometry,directions,intensities);
        double max_delta=0.0;

        std::vector<double> J(nc,0.0);
            for (std::size_t c=0;c<nc;++c)
                J[c]+=directions[m].weight*intensities[m](c)/(4.0*M_PI);

        for (std::size_t m=0;m<directions.size();++m) {
            const auto wall_values=build_dom_diffuse_gray_face_values(
                mesh,geometry,directions,intensities,walls,m,
                use_wall_cache ? &wall_cache : nullptr);
            cfdx::core::Field<double,cfdx::core::Location::FACE> directional_flux(
                mesh.n_faces(),"sI","W/m2",1);
            const auto& d=directions[m];
                directional_flux(f)=d.dx*geometry.face_area_vectors[f].x+
                                    d.dy*geometry.face_area_vectors[f].y+
                                    d.dz*geometry.face_area_vectors[f].z;

            cfdx::core::Field<double,cfdx::core::Location::CELL> source(
                nc,"radiation_source","W/m3/sr",1);
            cfdx::core::Field<double,cfdx::core::Location::CELL> sp(
                nc,"radiation_sp","1/m",1);
            for (std::size_t c=0;c<nc;++c) {
                source(c)=controls.absorption*blackbody_intensity(temperature(c))+
                          controls.scattering*J[c];
                sp(c)=-(controls.absorption+controls.scattering);
            }

            auto eq=assemble_scalar_equation(
                mesh,geometry,directional_flux,0.0,source,sp,
                extrapolation_bcs,true,&wall_values,nullptr,nullptr,nullptr,
                ConvectionScheme::UPWIND,&intensities[m]);

            cfdx::core::Vector intensity(nc,0.0);
            for (std::size_t c=0;c<nc;++c) intensity(c)=intensities[m](c);
            ScalarSolveControls sc;
            sc.max_iterations=controls.linear_max_iterations;
            sc.tolerance=controls.linear_tolerance;
            sc.relaxation=controls.intensity_relaxation;
            const auto lr=solve_scalar_equation(eq,intensity,sc);
            if (lr.status!=cfdx::core::SolverStatus::CONVERGED)
                throw std::runtime_error(
                    "DOM diffuse-gray wall intensity solve did not converge: residual="+
                    std::to_string(lr.residual));
            for (std::size_t c=0;c<nc;++c) {
                if (!std::isfinite(intensity(c)) || intensity(c)<-1e-12)
                    throw std::runtime_error(
                        "DOM diffuse-gray wall solve produced non-physical intensity");
                const double bounded=std::max(0.0,intensity(c));
                accepted_delta[c]=bounded-intensities[m](c);
                max_delta=std::max(max_delta,std::abs(accepted_delta[c]));
                intensities[m](c)=bounded;
            }
            if (use_wall_cache)
                update_dom_diffuse_wall_cache(
                    mesh,geometry,d,accepted_delta,wall_cache);
        }

        for (std::size_t c=0;c<nc;++c) {
            double G=0.0;
            for (std::size_t m=0;m<directions.size();++m)
                G+=directions[m].weight*intensities[m](c);
            active_G(c)=G;
            active_source(c)=controls.absorption*
                (4.0*M_PI*blackbody_intensity(temperature(c))-G);
        }

        double source_delta=0.0;
        for (std::size_t c=0;c<nc;++c)
            source_delta=std::max(source_delta,
                std::abs(active_source(c)-old_source(c)));

        double scale=1.0, source_scale=1.0;
        for (std::size_t c=0;c<nc;++c) {
            scale=std::max(scale,std::abs(active_G(c)));
            source_scale=std::max(source_scale,std::abs(active_source(c)));
        }
        const double rel_i=max_delta/scale;
        const double rel_s=source_delta/source_scale;
        result.history.push_back({iter,max_delta,source_delta});
        result.iterations=iter;
        if (rel_i<=controls.tolerance && rel_s<=controls.tolerance) {
            if (verify_accepted_state_v42) {
                const auto defect=audit_dom_accepted_state_v42(
                    mesh,geometry,temperature,intensities,directions,controls,walls);
                result.max_cell_transport_residual_v42=defect.max_cell_relative;
                result.relative_energy_closure_v42=defect.relative_global_energy;
                if (defect.max_cell_relative>controls.tolerance ||
                    defect.relative_global_energy>controls.tolerance)
                    continue; // Never accept a small Picard delta alone.
                result.accepted_state_verified_v42=true;
            }
            result.converged=true;
            break;
        }
    }
    if(result.converged && compute_wall_power_v41) {
        // One flux evaluation on the accepted angular field. The discrete
        // outgoing measure and incoming diffuse-gray radiation are evaluated
        // using identical quadrature; NO hidden rescale to physical pi.
        const auto cache=make_dom_diffuse_wall_cache(
            mesh,geometry,directions,intensities);
        for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
            const auto& patch=mesh.boundary().patch(p);
            if(patch.type==cfdx::core::PatchType::EMPTY) continue;
            const auto& wall=walls.at(patch.name);
            const long double Ib=blackbody_intensity(wall.temperature);
            long double watts=0.0L;
            for(const auto f:patch.face_ids) {
                const long double Hminus=cache.incoming_measure[f];
                const long double incoming_intensity=
                    wall.emissivity*Ib + (1.0-wall.emissivity)*
                    cache.outgoing_irradiation[f]/Hminus;
                watts+=static_cast<long double>(geometry.face_area_vectors[f].mag())*
                    (cache.outgoing_irradiation[f]-Hminus*incoming_intensity);
            }
            if(!std::isfinite(watts))
                throw std::runtime_error("DOM V41 nonfinite final wall power");
            result.wall_outward_W[patch.name]=static_cast<double>(watts);
        }
    }
    // Strict mode atomically publishes both physical fields only after the
    // accepted nonlinear state satisfies local and global FV balances.
    if(result.converged && verify_accepted_state_v42) {
        irradiation=std::move(*staged_G);
        radiation_source=std::move(*staged_source);
    }
    return result;
}

    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& irradiation,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& radiation_source,
    const std::vector<DiscreteDirection>& directions,
    const RadiationTransportControls& controls = {},
    const ScalarBoundaryConditions& wall_intensity_bcs = {})
{
    validate_radiation_transport_controls(controls);
    validate_discrete_directions(directions);
    if(temperature.size()!=mesh.n_cells() ||
       irradiation.size()!=mesh.n_cells() ||
       radiation_source.size()!=mesh.n_cells())
        throw std::invalid_argument("radiation field size mismatch");

    const std::size_t nc=mesh.n_cells();
    std::vector<cfdx::core::Field<double,cfdx::core::Location::CELL>> intensities;
    intensities.reserve(directions.size());
    for(std::size_t m=0;m<directions.size();++m) {
        intensities.emplace_back(nc,"I_"+std::to_string(m),"W/m2/sr",1);
        intensities.back().fill(0.0);
    }

    RadiationSolveResult result;
    for(std::size_t iter=1;iter<=controls.max_iterations;++iter) {
        auto old_radiation_source=radiation_source;
        std::vector<double> J(nc,0.0);
        for(std::size_t m=0;m<directions.size();++m)
            for(std::size_t c=0;c<nc;++c)
                J[c]+=directions[m].weight*intensities[m](c)/(4.0*M_PI);

        double max_delta=0.0;
        for(std::size_t m=0;m<directions.size();++m) {
            cfdx::core::Field<double,cfdx::core::Location::FACE> directional_flux(
                mesh.n_faces(),"sI","W/m2",1);
            const auto& d=directions[m];
            for(std::size_t f=0;f<mesh.n_faces();++f)
                directional_flux(f)=d.dx*geometry.face_area_vectors[f].x+
                                    d.dy*geometry.face_area_vectors[f].y+
                                    d.dz*geometry.face_area_vectors[f].z;

            cfdx::core::Field<double,cfdx::core::Location::CELL> source(
                nc,"radiation_source","W/m3/sr",1);
            cfdx::core::Field<double,cfdx::core::Location::CELL> sp(
                nc,"radiation_sp","1/m",1);

            for(std::size_t c=0;c<nc;++c) {
                const double Ib=blackbody_intensity(temperature(c));
                source(c)=controls.absorption*Ib +
                          controls.scattering*J[c];
                sp(c)=-(controls.absorption+controls.scattering);
            }

            auto eq=assemble_scalar_equation(
                mesh,geometry,directional_flux,0.0,source,sp,
                wall_intensity_bcs,true,nullptr,nullptr,nullptr,nullptr,
                ConvectionScheme::UPWIND,&intensities[m]);

            auto old=intensities[m];
            cfdx::core::Vector intensity(nc,0.0);
            for(std::size_t c=0;c<nc;++c) intensity(c)=intensities[m](c);
            ScalarSolveControls sc;
            sc.max_iterations=controls.linear_max_iterations;
            sc.tolerance=controls.linear_tolerance;
            sc.relaxation=controls.intensity_relaxation;
            const auto lr=solve_scalar_equation(eq,intensity,sc);
            for(std::size_t c=0;c<nc;++c) intensities[m](c)=intensity(c);
            if(lr.status!=cfdx::core::SolverStatus::CONVERGED) {
                throw std::runtime_error(
                    "radiation intensity linear solve did not converge: residual=" +
                    std::to_string(lr.residual) +
                    " relative=" + std::to_string(lr.residual_relative) +
                    " iterations=" + std::to_string(lr.iterations));
            }

            for(std::size_t c=0;c<nc;++c)
                max_delta=std::max(max_delta,
                    std::abs(intensities[m](c)-old(c)));
        }

        for(std::size_t c=0;c<nc;++c) {
            double G=0.0;
            for(std::size_t m=0;m<directions.size();++m)
                G+=directions[m].weight*intensities[m](c);
            irradiation(c)=G;
            radiation_source(c)=controls.absorption*
                (4.0*M_PI*blackbody_intensity(temperature(c))-G);
        }

        double max_source_delta = 0.0;
        for(std::size_t c=0;c<nc;++c)
            max_source_delta=std::max(max_source_delta,
                std::abs(radiation_source(c)-old_radiation_source(c)));
        result.history.push_back({iter,max_delta,max_source_delta});
        result.iterations=iter;
        double max_relative_delta=0.0;
        double max_relative_source_delta=0.0;
        for(std::size_t c=0;c<nc;++c) {
            double old_scale=1.0, source_scale=1.0;
            for(std::size_t m=0;m<directions.size();++m) {
                old_scale=std::max(old_scale,std::abs(intensities[m](c)));
            }
            old_scale=std::max(old_scale,std::abs(irradiation(c)));
            source_scale=std::max(source_scale,std::abs(radiation_source(c)));
            max_relative_delta=std::max(max_relative_delta,
                max_delta/old_scale);
            max_relative_source_delta=std::max(max_relative_source_delta,
                max_source_delta/source_scale);
        }
        if (mesh.n_cells() <= 64) {
            std::cerr << "RADIATION_RESIDUAL: iteration=" << iter
                      << " intensity_delta=" << max_delta
                      << " intensity_relative=" << max_relative_delta
                      << " source_delta=" << max_source_delta
                      << " source_relative=" << max_relative_source_delta
                      << '\n';
        }
        if(max_relative_delta<=controls.tolerance &&
           max_relative_source_delta<=controls.tolerance) {
            result.converged=true;
            break;
        }
    }
    return result;
}

struct RadiationEnergyCouplingControls {
    RadiationTransportControls radiation;
    EnergySolverControls energy;
    std::size_t max_outer_iterations = 100;
    // Relax the nonlinear radiation/energy fixed-point update separately from
    // the linear energy solver relaxation.
    double outer_relaxation = 0.5;
    // Legacy common tolerance; dedicated criteria override it when > 0.
    double tolerance = 1e-8;
    double temperature_tolerance = -1.0;
    double source_tolerance = -1.0;
    double energy_balance_tolerance = -1.0;
};

struct RadiationEnergyCouplingResult {
    bool converged = false;
    std::size_t iterations = 0;
    std::vector<double> source_residuals;
    std::vector<double> energy_balance_residuals;
};

inline RadiationEnergyCouplingResult solve_radiation_energy_coupled(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& non_radiative_source,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& irradiation,
    const std::vector<DiscreteDirection>& directions,
    const RadiationEnergyCouplingControls& controls = {},
    const ScalarBoundaryConditions& radiation_bcs = {},
    const ScalarBoundaryConditions& thermal_bcs = {})
{
    if(controls.max_outer_iterations==0 || controls.tolerance<=0.0 ||
       (controls.temperature_tolerance<=0.0 && controls.temperature_tolerance!=-1.0) ||
       (controls.source_tolerance<=0.0 && controls.source_tolerance!=-1.0) ||
       (controls.energy_balance_tolerance<=0.0 && controls.energy_balance_tolerance!=-1.0) ||
       !(controls.outer_relaxation>0.0 && controls.outer_relaxation<=1.0))
        throw std::invalid_argument("invalid radiation-energy coupling controls");

    const std::size_t nc=mesh.n_cells();
    cfdx::core::Field<double,cfdx::core::Location::CELL> qrad(
        nc,"qrad","W/m3",1);
    cfdx::core::Field<double,cfdx::core::Location::CELL> qtot(
        nc,"qtot","W/m3",1);

    RadiationEnergyCouplingResult result;
    // Keep the physical time level fixed across the nonlinear radiation/energy
    // outer iterations. The outer iterate is allowed to move independently.
    const auto T_previous_time=temperature;
    for(std::size_t iter=1;iter<=controls.max_outer_iterations;++iter) {
        cfdx::core::Field<double,cfdx::core::Location::CELL> oldT=temperature;
        cfdx::core::Field<double,cfdx::core::Location::CELL> old_qrad=qrad;
        cfdx::core::Field<double,cfdx::core::Location::CELL> source(
            nc,"radiation_source","W/m3",1);

        auto rr=solve_participating_radiation(
            mesh,geometry,temperature,irradiation,qrad,directions,
            controls.radiation,radiation_bcs);
        if(!rr.converged)
            throw std::runtime_error("radiation inner solve did not converge after " + std::to_string(rr.iterations) + " iterations");

        cfdx::core::Field<double,cfdx::core::Location::CELL> source_implicit(
            nc,"radiation_source_sp","W/m3/K",1);
        // Linearize the radiation source about the current outer iterate:
        // S(T) = non_radiative_source - q_rad(T),
        // Sp = dS/dT = -16 * absorption * sigma * T^3 <= 0,
        // Su = S(T*) - Sp*T*. This gives the energy solver the local Jacobian
        // while preserving the exact source value at the current iterate.
        for(std::size_t c=0;c<nc;++c) {
            const double Tstar=temperature(c);
            const double sigma=5.670374419e-8;
            const double sp=-16.0*controls.radiation.absorption*sigma*Tstar*Tstar*Tstar;
            source_implicit(c)=sp;
            source(c)=non_radiative_source(c)-qrad(c)-sp*Tstar;
        }

        // Solve the energy equation into a predictor, then relax the
        // nonlinear radiation/energy coupling update. The inner energy
        // relaxation controls linear/nonlinear convergence of that solve;
        // it must not be confused with the outer Picard relaxation.
        auto predictedT=temperature;
        auto er=solve_energy(
            mesh,geometry,mass_flux,predictedT,source,
            controls.energy,thermal_bcs,nullptr,&source_implicit,&T_previous_time);
        if(!er.converged)
            throw std::runtime_error("energy inner solve did not converge after " + std::to_string(er.iterations) + " iterations");

        double max_delta=0.0;
        for(std::size_t c=0;c<nc;++c) {
            const double candidate=predictedT(c);
            if(!std::isfinite(candidate) || candidate<0.0)
                throw std::runtime_error(
                    "radiation-energy coupling produced non-physical predicted temperature");
            const double updated =
                oldT(c) + controls.outer_relaxation*(candidate-oldT(c));
            if(!std::isfinite(updated) || updated<0.0)
                throw std::runtime_error(
                    "radiation-energy coupling produced non-physical temperature");
            temperature(c)=updated;
            max_delta=std::max(max_delta,std::abs(updated-oldT(c)));
        }

        double qrad_delta=0.0;
        double qrad_scale=1.0;
        for(std::size_t c=0;c<nc;++c) {
            qrad_delta=std::max(qrad_delta,std::abs(qrad(c)-old_qrad(c)));
            qrad_scale=std::max(qrad_scale,std::abs(qrad(c)));
            qrad_scale=std::max(qrad_scale,std::abs(old_qrad(c)));
        }
        const double source_relative=qrad_delta/qrad_scale;

        // Reassemble the thermal equation at the actually accepted outer
        // state. This avoids reporting the predictor balance after relaxation.
        for(std::size_t c=0;c<nc;++c)
            source(c)=non_radiative_source(c)-qrad(c);
        cfdx::core::Field<double,cfdx::core::Location::CELL> accepted_source(
            nc,"accepted_radiation_source","W/m3",1);
        for(std::size_t c=0;c<nc;++c)
            accepted_source(c)=non_radiative_source(c)-qrad(c);
        auto accepted_eq=assemble_energy_equation(
            mesh,geometry,mass_flux,accepted_source,T_previous_time,controls.energy,
            thermal_bcs,nullptr,nullptr,nullptr,nullptr,&temperature);
        cfdx::core::Vector accepted_vec(nc,0.0);
        for(std::size_t c=0;c<nc;++c) accepted_vec(c)=temperature(c);
        const double energy_balance_residual =
            energy_balance_relative_from_equation(accepted_eq,accepted_vec);

        result.source_residuals.push_back(qrad_delta);
        result.energy_balance_residuals.push_back(energy_balance_residual);
        result.iterations=iter;
        if (mesh.n_cells() <= 64) {
            std::cerr << "THERMAL_RADIATION_RESIDUAL: iteration=" << iter
                      << " dT=" << max_delta
                      << " qrad_delta=" << qrad_delta
                      << " qrad_relative=" << source_relative
                      << " energy_balance=" << energy_balance_residual
                      << '\n';
        }
        const double source_tolerance =
            controls.source_tolerance > 0.0 ? controls.source_tolerance : controls.tolerance;
        const double temperature_tolerance =
            controls.temperature_tolerance > 0.0 ? controls.temperature_tolerance : controls.tolerance;
        const double energy_balance_tolerance =
            controls.energy_balance_tolerance > 0.0 ? controls.energy_balance_tolerance : controls.tolerance;
        const double temperature_scale = [&]() {
            double scale=1.0;
            for(std::size_t c=0;c<nc;++c) scale=std::max(scale,std::abs(temperature(c)));
            return scale;
        }();
        if(max_delta/temperature_scale<=temperature_tolerance &&
           source_relative<=source_tolerance &&
           energy_balance_residual<=energy_balance_tolerance) {
            result.converged=true;
            break;
        }
    }
    return result;
}

} // namespace cfdx::physics
