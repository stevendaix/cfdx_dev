#pragma once
// V43: research-only diffusion synthetic acceleration candidate for the
// EXACT UPWIND finite-volume scalar flux fixed point. This deliberately does
// not change the production DOM solver and supports BLACK WALLS ONLY.
//
// The low-order diffusion matrix is SPD and conservative, but its boundary
// extrapolation and non-orthogonal metric are not a proven discretely
// consistent Schur complement of CFDX's high-order angular transport.
// Thus acceptance is based on the ACTUAL reevaluated transport fixed-point
// residual. No acceleration is claimed in situations where the gate rejects.
#include "cfdx/physics/radiation_solver.h"
#include <array>
#include <chrono>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics {

struct DomDsaStudyV43Result {
    bool converged=false;
    std::size_t iterations=0;
    std::size_t transport_maps=0;
    std::size_t accepted_accelerations=0;
    std::size_t rejected_accelerations=0;
    double fixed_point_relative=std::numeric_limits<double>::infinity();
    double local_transport_relative=std::numeric_limits<double>::infinity();
    double global_energy_relative=std::numeric_limits<double>::infinity();
    double elapsed_s=0.;
    std::vector<double> scalar_mean_intensity;
};

inline cfdx::core::SparseMatrix build_dom_dsa_diffusion_matrix_v43(
    const cfdx::core::Mesh& mesh, const FvGeometry& geo,
    double absorption, double scattering)
{
    const std::size_t nc=mesh.n_cells();
    if (geo.cell_volumes.size()!=nc ||
        geo.cell_centres.size()!=nc ||
        geo.face_area_vectors.size()!=mesh.n_faces() ||
        geo.face_centres.size()!=mesh.n_faces() ||
        !(absorption>0.) || !(scattering>=0.) ||
        !std::isfinite(absorption) || !std::isfinite(scattering))
        throw std::invalid_argument("DOM DSA V43: invalid geometry or optical coefficients");
    const double D=1.0/(3.0*(absorption+scattering));
    cfdx::core::SparseMatrix B(nc,nc);
    std::vector<double> diagonal(nc,0.0);
    const auto& own=mesh.ownership();
    for (std::size_t c=0;c<nc;++c) {
        if (!(geo.cell_volumes[c]>0.) || !std::isfinite(geo.cell_volumes[c]))
            throw std::invalid_argument("DOM DSA V43: non-positive cell volume");
        diagonal[c]+=absorption*geo.cell_volumes[c];
    }
    for (std::size_t f=0;f<mesh.n_faces();++f) {
        const auto sf=geo.face_area_vectors[f];
        const double area=sf.mag();
        if (!(area>0.) || !std::isfinite(area))
            throw std::invalid_argument("DOM DSA V43: non-positive face area");
        const std::size_t o=own.owner(f);
        if (o>=nc) throw std::invalid_argument("DOM DSA V43: invalid owner");
        const auto nb=own.neighbour(f);
        if (nb>=0) {
            const auto n=static_cast<std::size_t>(nb);
            if(n>=nc || n==o) throw std::invalid_argument("DOM DSA V43: invalid neighbour");
            const auto d=geo.cell_centres[n]-geo.cell_centres[o];
            const double dnormal=d.dot(sf)/area;
            if (!(dnormal>0.) || !std::isfinite(dnormal))
                throw std::invalid_argument("DOM DSA V43: inverted internal face");
            const double conductance=D*area/dnormal;
            diagonal[o]+=conductance; diagonal[n]+=conductance;
            B.push_back(o,n,-conductance);
            B.push_back(n,o,-conductance);
        } else {
            // Homogeneous black-wall transport correction: Marshak-like
            // extrapolation length 2D. Geometrically conservative face flux.
            const auto d=geo.face_centres[f]-geo.cell_centres[o];
            const double dnormal=d.dot(sf)/area;
            if (!(dnormal>0.) || !std::isfinite(dnormal))
                throw std::invalid_argument("DOM DSA V43: inverted boundary face");
            diagonal[o]+=area/(dnormal/D+2.0);
        }
    }
    for (std::size_t c=0;c<nc;++c)
        B.push_back(c,c,diagonal[c]);
    B.finalize();
    return B;
}

// Solve a sequence of REAL angular finite-volume systems using the same
// upwind, bounded convection, scalar solver and face metrics as production.
// The fixed-point residual is F(J)-J; a candidate DSA step is accepted only
// if that complete transport residual is below the unaccelerated candidate's.
inline DomDsaStudyV43Result study_dom_dsa_v43(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geo,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    const std::vector<DiscreteDirection>& directions,
    const DomWallBoundaryConditions& walls,
    const RadiationTransportControls& controls,
    bool use_dsa,
    std::size_t max_correction_trials=5)
{
    validate_radiation_transport_controls(controls);
    validate_discrete_directions(directions);
    validate_dom_wall_conditions(mesh,walls);
    const std::size_t nc=mesh.n_cells();
    if (temperature.size()!=nc || nc==0 || !(controls.absorption>0.) ||
        !std::isfinite(controls.absorption) ||
        !(controls.scattering>=0.) || !std::isfinite(controls.scattering) ||
        max_correction_trials==0 || max_correction_trials>12)
        throw std::invalid_argument("DOM DSA V43: invalid study inputs");
    for (std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& patch=mesh.boundary().patch(p);
        if (patch.type==cfdx::core::PatchType::EMPTY ||
            walls.at(patch.name).emissivity!=1.0)
            throw std::invalid_argument("DOM DSA V43: prototype supports black walls only");
    }
    const auto start=std::chrono::steady_clock::now();
    const double four_pi=4.0*M_PI;
    using CellField=cfdx::core::Field<double,cfdx::core::Location::CELL>;
    using FaceField=cfdx::core::Field<double,cfdx::core::Location::FACE>;
    std::vector<ScalarEquation> transport;
    transport.reserve(directions.size());
    CellField zero(nc,"zero","W/m3/sr",1),sp(nc,"sp","1/m",1);
    CellField old_intensity(nc,"I","W/m2/sr",1);
    zero.fill(0.0);sp.fill(-(controls.absorption+controls.scattering));
    old_intensity.fill(0.0);
    ScalarBoundaryConditions bcs;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p)
        bcs[mesh.boundary().patch(p).name]={ScalarBoundaryType::ZERO_GRADIENT,0.,0.};
    for(const auto& d:directions) {
        FaceField flux(mesh.n_faces(),"s dot Sf","m2",1);
        ScalarBoundaryFaceValues values;
        for (std::size_t f=0;f<mesh.n_faces();++f) {
            const auto sf=geo.face_area_vectors[f];
            flux(f)=d.dx*sf.x+d.dy*sf.y+d.dz*sf.z;
        }
        for (std::size_t p=0;p<mesh.boundary().n_patches();++p) {
            const auto& patch=mesh.boundary().patch(p);
            auto& v=values.values[patch.name];
            v.assign(mesh.n_faces(),std::numeric_limits<double>::quiet_NaN());
            for (auto f:patch.face_ids)
                v[f]=blackbody_intensity(walls.at(patch.name).temperature);
        }
        transport.push_back(assemble_scalar_equation(
            mesh,geo,flux,0.0,zero,sp,bcs,true,&values,
            nullptr,nullptr,nullptr,ConvectionScheme::UPWIND,&old_intensity));
    }
    std::vector<cfdx::core::Vector> rhs_base;
    rhs_base.reserve(transport.size());
    for(const auto& eq:transport) rhs_base.push_back(eq.rhs);
    std::vector<double> Ib(nc,0.);
    for(std::size_t c=0;c<nc;++c) {
        const double T=temperature(c);
        if (!(T>=0.) || !std::isfinite(T))
            throw std::invalid_argument("DOM DSA V43: invalid gas temperature");
        Ib[c]=blackbody_intensity(T);
    }
    cfdx::core::SparseMatrix diffusion;
    if(use_dsa) diffusion=build_dom_dsa_diffusion_matrix_v43(
        mesh,geo,controls.absorption,controls.scattering);
    ScalarEquation low_order;
    if(use_dsa) {
        low_order.matrix=std::move(diffusion);
        low_order.rhs=cfdx::core::Vector(nc,0.);
    }
    DomDsaStudyV43Result result;
    std::vector<double> J=Ib;
    std::vector<CellField> accepted_intensities;
    const auto map=[&](const std::vector<double>& input,
                       std::vector<CellField>* capture=nullptr) {
        if (input.size()!=nc) throw std::invalid_argument("DOM DSA V43: J dimension mismatch");
        std::vector<double> out(nc,0.);
        if(capture) capture->clear();
        for (std::size_t m=0;m<directions.size();++m) {
            auto& eq=transport[m];
            for(std::size_t c=0;c<nc;++c)
                eq.rhs(c)=rhs_base[m](c)+geo.cell_volumes[c]*(
                    controls.absorption*Ib[c]+controls.scattering*input[c]);
            cfdx::core::Vector intensity(nc,0.);
            ScalarSolveControls solver;
            solver.max_iterations=controls.linear_max_iterations;
            solver.tolerance=controls.linear_tolerance;
            solver.relaxation=1.;
            const auto lr=solve_scalar_equation(eq,intensity,solver);
            if(lr.status!=cfdx::core::SolverStatus::CONVERGED)
                throw std::runtime_error("DOM DSA V43: angular linear solve failed");
            if(capture) {
                capture->emplace_back(nc,"I","W/m2/sr",1);
            }
            for(std::size_t c=0;c<nc;++c) {
                if (!std::isfinite(intensity(c)) || intensity(c)<-1e-9)
                    throw std::runtime_error("DOM DSA V43: nonphysical angular result");
                out[c]+=directions[m].weight*intensity(c)/four_pi;
                if(capture) capture->back()(c)=std::max(0.,intensity(c));
            }
        }
        ++result.transport_maps;
        return out;
    };
    const auto norm=[&](const std::vector<double>& input,
                        const std::vector<double>& output) {
        double defect=0.,scale=1.;
        for(std::size_t c=0;c<nc;++c) {
            defect=std::max(defect,std::abs(output[c]-input[c]));
            scale=std::max(scale,std::abs(output[c]));
        }
        return defect/scale;
    };
    for(std::size_t k=0;k<controls.max_iterations;++k) {
        const auto FJ=map(J);
        result.iterations=k+1;
        const double r=norm(J,FJ);
        result.fixed_point_relative=r;
        if(r<=controls.tolerance) {
            result.iterations=k+1;
            // Recompute directional intensity from the FINAL J for independent
            // exact FV physical audit. The map's J discrepancy is also gated.
            (void)map(J,&accepted_intensities);
            const auto physical=audit_dom_accepted_state_v42(
                mesh,geo,temperature,accepted_intensities,directions,controls,walls);
            result.local_transport_relative=physical.max_cell_relative;
            result.global_energy_relative=physical.relative_global_energy;
            if(result.local_transport_relative<=controls.tolerance &&
               result.global_energy_relative<=controls.tolerance) {
                result.converged=true;
                break;
            }
        }
        if(!use_dsa) { J=FJ;continue; }
        cfdx::core::Vector delta(nc,0.);
        for (std::size_t c=0;c<nc;++c)
            low_order.rhs(c)=controls.scattering*geo.cell_volumes[c]*(FJ[c]-J[c]);
        ScalarSolveControls solver;
        solver.tolerance=std::min(controls.linear_tolerance,1e-11);
        solver.max_iterations=controls.linear_max_iterations;
        solver.relaxation=1.;
        const auto low=solve_scalar_equation(low_order,delta,solver);
        if(low.status!=cfdx::core::SolverStatus::CONVERGED)
            throw std::runtime_error("DOM DSA V43: low order diffusion solve failed");
        // Baseline true fixed-point residual after the usual SI update.
        const auto FFJ=map(FJ);
        const double base_r=norm(FJ,FFJ);
        bool accepted=false;
        for(std::size_t trial=0;trial<max_correction_trials;++trial) {
            const double alpha=std::ldexp(1.0,-static_cast<int>(trial));
            std::vector<double> candidate(nc,0.);
            bool physical=true;
            for(std::size_t c=0;c<nc;++c) {
                candidate[c]=FJ[c]+alpha*delta(c);
                if(!std::isfinite(candidate[c]) || candidate[c]<0.) physical=false;
            }
            if(!physical) continue;
            const auto next=map(candidate);
            const double next_r=norm(candidate,next);
            if(next_r<base_r*(1.-1e-8)) {
                J=std::move(candidate);++result.accepted_accelerations;
                accepted=true;break;
            }
        }
        if(!accepted) {J=FJ;++result.rejected_accelerations;}
    }
    result.scalar_mean_intensity=std::move(J);
    result.elapsed_s=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-start).count();
    return result;
}
} // namespace cfdx::physics
