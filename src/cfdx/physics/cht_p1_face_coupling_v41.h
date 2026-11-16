#pragma once
// CFDX V41 opt-in: face-resolved coupled CHT/P1 interface. The gas receives
// only its physical volume emission/absorption source; the wall radiation
// enters the solid through the SAME geometric interface overlap on the face
// RHS, not a surrogate volume source. Conduction is shared, radiative power
// is nonlocal in the gas. Transactional nonlinear convergence gates preserved.
#include "cfdx/physics/cht_nonconformal_v37.h"
#include "cfdx/physics/p1_marshak_v39.h"
#include <algorithm>
#include <cmath>
#include <limits>

#include <vector>

namespace cfdx::physics::v41 {
using namespace cfdx::core;
using namespace cfdx::physics::v37;
struct ChtP1FaceControls {
    double emissivity=0.85;
    double absorption=0.7;
    double scattering=1.3;
    double coupling_tolerance=1.e-10;
    double outer_relaxation=0.1;
    std::size_t max_outer_iterations=80;
    P1MarshakControls p1;
};
struct ChtP1FaceReport {
    bool converged=false;
    std::size_t outer_iterations=0;
    double radiation_power_to_solid_W=0;
    double radiation_power_from_fluid_W=0;
    double radiative_balance=std::numeric_limits<double>::infinity();
    double max_temperature_relative=std::numeric_limits<double>::infinity();
    double max_radiative_source_relative=std::numeric_limits<double>::infinity();
    double cht_max_cell_residual=std::numeric_limits<double>::infinity();
    double cht_global_balance=std::numeric_limits<double>::infinity();
    std::size_t interface_overlaps=0;
    double wall_face_power_mismatch_W=std::numeric_limits<double>::infinity();
    double interface_conduction_into_solid_W=0.0;
    double interface_total_into_solid_W=0.0;
    double accepted_global_energy_residual=std::numeric_limits<double>::infinity();
};
// Independently re-evaluate the FINAL accepted nonlinear (T,G) state from
// conservative FV face fluxes, not the CHT linear predictor's frozen RHS.
// Radiation is represented exactly once on solid faces and exactly once as
// gas volume emission/absorption. This is not a fictitious paired cell flux.
struct ChtP1FacePhysicalResidual {double maximum_cell=0,global=0;};
inline ChtP1FacePhysicalResidual evaluate_cht_p1_face_residual_v41(
    const Mesh& m1,const FvGeometry& g1,const Mesh& m2,const FvGeometry& g2,
    const InterfaceMap& map,
    const Field<double,Location::FACE>& phi1,
    const Field<double,Location::FACE>& phi2,
    const Field<double,Location::CELL>& T1,
    const Field<double,Location::CELL>& T2,
    const Field<double,Location::CELL>& old1,
    const Field<double,Location::CELL>& old2,
    const Field<double,Location::CELL>& source1,
    const Field<double,Location::CELL>& source2,
    const EnergySolverControls& c1,const EnergySolverControls& c2,
    const ChtInterfaceControls& interface,
    const ScalarBoundaryConditions& raw1,const ScalarBoundaryConditions& raw2,
    const std::vector<double>& incoming_solid_overlap_watts)
{
    if(incoming_solid_overlap_watts.size()!=map.overlaps.size())
        throw std::invalid_argument("CHT-P1 V41 independent audit overlap size mismatch");
    std::vector<double> conduction1(m1.n_faces(),0.),conduction2(m2.n_faces(),0.),
                        radiation2(m2.n_faces(),0.);
    for(std::size_t j=0;j<map.overlaps.size();++j) {
        const auto& edge=map.overlaps[j];
        const double cond=edge.conductance*(T1(edge.cell1)-T2(edge.cell2));
        conduction1[edge.face1]+=cond;
        conduction2[edge.face2]-=cond;
        radiation2[edge.face2]+=incoming_solid_overlap_watts[j];
    }
    ScalarBoundaryFaceConditions fb1(m1.n_faces()),fb2(m2.n_faces());
    for(const auto f:m1.boundary().patch(patch_id(m1,interface.region1_patch)).face_ids) {
        fb1.specified[f]=1;
        fb1.conditions[f]={ScalarBoundaryFaceCondition::Type::FIXED_GRADIENT,
            -conduction1[f]/(c1.conductivity*g1.face_area_vectors[f].mag())};
    }
    for(const auto f:m2.boundary().patch(patch_id(m2,interface.region2_patch)).face_ids) {
        fb2.specified[f]=1;
        fb2.conditions[f]={ScalarBoundaryFaceCondition::Type::FIXED_GRADIENT,
            -(conduction2[f]-radiation2[f])/(c2.conductivity*g2.face_area_vectors[f].mag())};
    }
    auto bc1=raw1,bc2=raw2;
    bc1[interface.region1_patch]={ScalarBoundaryType::ZERO_GRADIENT,0.,0.};
    bc2[interface.region2_patch]={ScalarBoundaryType::ZERO_GRADIENT,0.,0.};
    const auto q1=reconstruct_energy_face_flux(m1,g1,phi1,T1,c1,bc1,
                                                nullptr,nullptr,&fb1);
    const auto q2=reconstruct_energy_face_flux(m2,g2,phi2,T2,c2,bc2,
                                                nullptr,nullptr,&fb2);
    std::vector<long double> residual1(m1.n_cells(),0.),gross1(m1.n_cells(),0.),
                             residual2(m2.n_cells(),0.),gross2(m2.n_cells(),0.);
    const auto add_external=[&](const Mesh& m,const FvGeometry& g,
          const Field<double,Location::FACE>& q,
          const Field<double,Location::CELL>& T,
          const Field<double,Location::CELL>& prev,
          const Field<double,Location::CELL>& source,
          const EnergySolverControls& c,const std::string& skip,
          std::vector<long double>& residual,std::vector<long double>& gross) {
        const auto pid=patch_id(m,skip);
        for(std::size_t f=0;f<m.n_faces();++f) {
            if(m.ownership().neighbour(f)<0 && g.face_patch[f]==pid)continue;
            const auto owner=m.ownership().owner(f);
            if(!std::isfinite(q(f)))throw std::runtime_error("CHT-P1 V41 nonfinite physical face flux");
            residual[owner]+=q(f);gross[owner]+=std::abs(q(f));
            if(m.ownership().neighbour(f)>=0) {
                const auto n=static_cast<std::size_t>(m.ownership().neighbour(f));
                residual[n]-=q(f);gross[n]+=std::abs(q(f));
            }
        }
        for(std::size_t i=0;i<m.n_cells();++i) {
            const long double unsteady=c.dt>0.0 ?
                static_cast<long double>(c.density)*c.cp*g.cell_volumes[i]*
                (static_cast<long double>(T(i))-prev(i))/c.dt : 0.0L;
            const long double src=static_cast<long double>(source(i))*g.cell_volumes[i];
            residual[i]+=unsteady-src;
            gross[i]+=std::abs(unsteady)+std::abs(src);
        }
    };
    add_external(m1,g1,q1,T1,old1,source1,c1,interface.region1_patch,residual1,gross1);
    add_external(m2,g2,q2,T2,old2,source2,c2,interface.region2_patch,residual2,gross2);
    for(std::size_t j=0;j<map.overlaps.size();++j){
        const auto& edge=map.overlaps[j];
        const long double cond=static_cast<long double>(edge.conductance)*
            (static_cast<long double>(T1(edge.cell1))-T2(edge.cell2));
        const long double rad=incoming_solid_overlap_watts[j];
        residual1[edge.cell1]+=cond;
        residual2[edge.cell2]-=cond+rad;
        gross1[edge.cell1]+=std::abs(cond);
        gross2[edge.cell2]+=std::abs(cond)+std::abs(rad);
    }
    long double net=0.,scale=0.;
    ChtP1FacePhysicalResidual result;
    for(std::size_t i=0;i<residual1.size();++i) {
        net+=residual1[i];scale+=gross1[i];
        result.maximum_cell=std::max(result.maximum_cell,
            static_cast<double>(std::abs(residual1[i])/std::max(1.0L,gross1[i])));
    }
    for(std::size_t i=0;i<residual2.size();++i) {
        net+=residual2[i];scale+=gross2[i];
        result.maximum_cell=std::max(result.maximum_cell,
            static_cast<double>(std::abs(residual2[i])/std::max(1.0L,gross2[i])));
    }
    result.global=static_cast<double>(std::abs(net)/std::max(1.0L,scale));
    if(!std::isfinite(result.maximum_cell)||!std::isfinite(result.global))
        throw std::runtime_error("CHT-P1 V41 nonfinite accepted physical residual");
    return result;
}

inline ChtP1FaceReport solve_nonconformal_cht_p1_face_v41(
    const Mesh& mfluid,const FvGeometry& gfluid,
    const Mesh& msolid,const FvGeometry& gsolid,
    const Field<double,Location::FACE>& phi_fluid,
    const Field<double,Location::FACE>& phi_solid,
    Field<double,Location::CELL>& Tfluid,
    Field<double,Location::CELL>& Tsolid,
    const Field<double,Location::CELL>& base_source_fluid,
    const Field<double,Location::CELL>& base_source_solid,
    Field<double,Location::CELL>& irradiation,
    Field<double,Location::CELL>& emission_minus_absorption,
    const EnergySolverControls& energyfluid,
    const EnergySolverControls& energysolid,
    const ChtInterfaceControls& interface,
    const ChtP1FaceControls& controls={},
    const ScalarBoundaryConditions& thermal_fluid={},
    const ScalarBoundaryConditions& thermal_solid={})
{
    if(controls.max_outer_iterations==0 || !(controls.emissivity>0&&controls.emissivity<=1) ||
       !(controls.absorption>0) || !(controls.scattering>=0) ||
       !(controls.coupling_tolerance>0) ||
       !(controls.outer_relaxation>0&&controls.outer_relaxation<=1) ||
       !std::isfinite(controls.emissivity)||!std::isfinite(controls.absorption)||
       !std::isfinite(controls.scattering)||!std::isfinite(controls.coupling_tolerance))
        throw std::invalid_argument("CHT-P1 V41 invalid controls");
    if(Tfluid.size()!=mfluid.n_cells() || Tsolid.size()!=msolid.n_cells() ||
       irradiation.size()!=mfluid.n_cells() ||
       emission_minus_absorption.size()!=mfluid.n_cells() ||
       base_source_fluid.size()!=mfluid.n_cells() ||
       base_source_solid.size()!=msolid.n_cells())
        throw std::invalid_argument("CHT-P1 V41 incompatible fields");
    const auto map=build_nonconformal_interface(mfluid,gfluid,msolid,gsolid,interface);
    const auto iface=patch_id(mfluid,interface.region1_patch);
    const auto& interface_patch=mfluid.boundary().patch(iface);
    const auto original_fluid=Tfluid,original_solid=Tsolid;
    auto current_fluid=Tfluid,current_solid=Tsolid;
    auto current_G=irradiation,current_Q=emission_minus_absorption;
    std::vector<double> aa(mfluid.n_cells(),controls.absorption),
                        ss(mfluid.n_cells(),controls.scattering);
    P1MarshakWalls walls;
    for(std::size_t p=0;p<mfluid.boundary().n_patches();++p) {
        const auto& patch=mfluid.boundary().patch(p);
        walls[patch.name]={patch.name==interface.region1_patch ? controls.emissivity:0.0,300.};
    }
    P1MarshakControls pc=controls.p1;
    pc.skew_marshak_robin_v40=true;
    ChtP1FaceReport report;report.interface_overlaps=map.overlaps.size();
    // The region pair is solved once for every staged radiation load, but the
    // original physical time level must remain frozen throughout this loop.
    auto evaluate_radiation=[&](const Field<double,Location::CELL>& f,
                                const Field<double,Location::CELL>& s,
                                Field<double,Location::CELL>& G,
                                Field<double,Location::CELL>& Q,
                                Field<double,Location::CELL>& srcf,
                                Field<double,Location::CELL>& srcs,
                                std::vector<double>& overlap_watts,
                                double& outward,
                                double& radiative_balance){
        P1MarshakFaceWalls face_walls;
        std::vector<long double> wt(mfluid.n_faces(),0.0L),
                                 area(mfluid.n_faces(),0.0L);
        for(const auto& e:map.overlaps) {
            const double d1=boundary_normal_distance(gfluid,e.face1,e.cell1);
            const double d2=boundary_normal_distance(gsolid,e.face2,e.cell2);
            const double gf=energyfluid.conductivity/d1,gs=energysolid.conductivity/d2;
            const double T=(gf*f(e.cell1)+gs*s(e.cell2))/(gf+gs);
            wt[e.face1]+=static_cast<long double>(e.area)*T;
            area[e.face1]+=e.area;
        }
        for(const auto face:interface_patch.face_ids) {
            if(!(area[face]>0)) throw std::runtime_error("CHT-P1 uncovered radiative interface");
            const double T=static_cast<double>(wt[face]/area[face]);
            if(!(T>0.) || !std::isfinite(T))
                throw std::runtime_error("CHT-P1 nonphysical interface wall temperature");
            face_walls.emplace(face,P1MarshakWall{controls.emissivity,T});
        }
        const auto rad=solve_p1_radiation_marshak(
            mfluid,gfluid,f,G,Q,aa,ss,walls,pc,&face_walls);
        if(!rad.converged)return false;
        srcf=base_source_fluid;srcs=base_source_solid;
        overlap_watts.assign(map.overlaps.size(),0.0);
        long double total_radiation_cell=0.0L, gross_radiation_cell=0.0L;
        for(std::size_t i=0;i<mfluid.n_cells();++i){
            srcf(i)-=Q(i); // physical matter energy source is -emission
            const long double power=static_cast<long double>(Q(i))*gfluid.cell_volumes[i];
            total_radiation_cell+=power;
            gross_radiation_cell+=std::abs(power);
        }
        long double watt_solid=0.0L, gross_wall=0.0L;
        for(std::size_t j=0;j<map.overlaps.size();++j){
            const auto& e=map.overlaps[j];
            const double Af=gfluid.face_area_vectors[e.face1].mag();
            const double watt=rad.face_outward_W.at(e.face1)*(e.area/Af);
            if(!std::isfinite(watt))
                throw std::runtime_error("CHT-P1 V41 nonfinite overlap wall power");
            // One shared watt rate per geometric overlap; the solid FV face
            // will receive this through its Neumann boundary term in V37.
            overlap_watts[j]=watt;
            watt_solid+=watt;
            gross_wall+=std::abs(static_cast<long double>(watt));
        }
        outward=static_cast<double>(watt_solid);
        // Same gross-power normalization as the independent P1 physical gate.
        // Replacing gross by net powers spuriously doubles this measure for
        // equal/opposite exchanges; the numerical tolerance itself is unchanged.
        radiative_balance=static_cast<double>(
            std::abs(watt_solid-total_radiation_cell)/
            std::max(1.0L,gross_wall+gross_radiation_cell));
        return radiative_balance<=pc.balance_tolerance &&
            rad.relative_global_balance<=pc.balance_tolerance &&
            rad.maximum_normalized_cell_residual<=pc.physical_tolerance;
    };
    for(std::size_t iter=1;iter<=controls.max_outer_iterations;++iter) {
        auto srcf=base_source_fluid,srcs=base_source_solid;
        auto Gpred=current_G,Qpred=current_Q;
        std::vector<double> overlap_old;
        double watt_old=0,balance_old=0;
        if(!evaluate_radiation(current_fluid,current_solid,Gpred,Qpred,
                               srcf,srcs,overlap_old,watt_old,balance_old))return report;
        auto trial_fluid=current_fluid,trial_solid=current_solid;
        const auto cht=solve_nonconformal_cht(
            mfluid,gfluid,msolid,gsolid,phi_fluid,phi_solid,
            trial_fluid,trial_solid,srcf,srcs,energyfluid,energysolid,
            interface,thermal_fluid,thermal_solid,nullptr,
            &original_fluid,&original_solid,&overlap_old);
        if(!cht.converged)return report;
        // Damped nonlinear step: the monolithic regional temperature predictor
        // can be far outside the basin of the stiff fourth-power radiation
        // source. The *accepted* stage is damped, then P1 is reassembled there.
        for(std::size_t i=0;i<trial_fluid.size();++i) {
            trial_fluid(i)=current_fluid(i)+controls.outer_relaxation*
                (trial_fluid(i)-current_fluid(i));
            if(!(trial_fluid(i)>0.)||!std::isfinite(trial_fluid(i))) return report;
        }
        for(std::size_t i=0;i<trial_solid.size();++i) {
            trial_solid(i)=current_solid(i)+controls.outer_relaxation*
                (trial_solid(i)-current_solid(i));
            if(!(trial_solid(i)>0.)||!std::isfinite(trial_solid(i))) return report;
        }
        // Physically RE-EVALUATE the radiation at the damped candidate
        // temperature, not the previous/frozen CHT temperature.
        auto Gnew=Gpred,Qnew=Qpred;
        auto srcf_new=base_source_fluid,srcs_new=base_source_solid;
        std::vector<double> overlap_new;
        double watt_new=0,balance_new=0;
        if(!evaluate_radiation(trial_fluid,trial_solid,Gnew,Qnew,
                               srcf_new,srcs_new,overlap_new,watt_new,balance_new))return report;
        double change_t=0,scale_t=1,change_source_relative=0;
        for(std::size_t i=0;i<current_fluid.size();++i) {
            change_t=std::max(change_t,std::abs(trial_fluid(i)-current_fluid(i)));
            scale_t=std::max({scale_t,std::abs(trial_fluid(i)),std::abs(current_fluid(i))});
            const double v1=static_cast<double>(gfluid.cell_volumes[i])*srcf(i);
            const double v2=static_cast<double>(gfluid.cell_volumes[i])*srcf_new(i);
            change_source_relative=std::max(change_source_relative,
                std::abs(v1-v2)/std::max({1.0,std::abs(v1),std::abs(v2)}));
        }
        for(std::size_t i=0;i<current_solid.size();++i){
            change_t=std::max(change_t,std::abs(trial_solid(i)-current_solid(i)));
            scale_t=std::max({scale_t,std::abs(trial_solid(i)),std::abs(current_solid(i))});
            const double v1=static_cast<double>(gsolid.cell_volumes[i])*srcs(i);
            const double v2=static_cast<double>(gsolid.cell_volumes[i])*srcs_new(i);
            change_source_relative=std::max(change_source_relative,
                std::abs(v1-v2)/std::max({1.0,std::abs(v1),std::abs(v2)}));
        }
        // The candidate P1 wall flux must converge face-by-face as well
        // as the cell volume emission. Global cancellation is insufficient.
        double max_face_power_mismatch=0.0;
        for(std::size_t j=0;j<overlap_old.size();++j) {
            const double a=overlap_old[j],b=overlap_new[j];
            const double relative=std::abs(b-a)/std::max({1.0,std::abs(a),std::abs(b)});
            change_source_relative=std::max(change_source_relative,relative);
            max_face_power_mismatch=std::max(max_face_power_mismatch,std::abs(b-a));
        }
        report.wall_face_power_mismatch_W=max_face_power_mismatch;
        report.interface_conduction_into_solid_W=cht.interface_conduction_into_region2_W;
        report.interface_total_into_solid_W=cht.interface_total_into_region2_W;
        // This is the explicitly recomputed (not predictor-only) global
        // conservation measure: P1 volume radiation equals the accepted
        // wall face power; the independent CHT residual is checked below.
        const auto accepted_audit=evaluate_cht_p1_face_residual_v41(
            mfluid,gfluid,msolid,gsolid,map,phi_fluid,phi_solid,
            trial_fluid,trial_solid,original_fluid,original_solid,
            srcf_new,srcs_new,energyfluid,energysolid,interface,
            thermal_fluid,thermal_solid,overlap_new);
        report.accepted_global_energy_residual=std::max(balance_new,accepted_audit.global);
        report.outer_iterations=iter;
        report.max_temperature_relative=change_t/scale_t;
        report.max_radiative_source_relative=change_source_relative;
        report.radiation_power_to_solid_W=watt_new;
        report.radiation_power_from_fluid_W=static_cast<double>([&](){
            long double v=0;for(std::size_t i=0;i<Qnew.size();++i)
                v+=static_cast<long double>(Qnew(i))*gfluid.cell_volumes[i];
            return v;
        }());
        report.radiative_balance=balance_new;
        report.cht_max_cell_residual=accepted_audit.maximum_cell;
        // Guard against bypassing the new face assembly by inserting the
        // radiation a second time as a region-2 volume source.
        if(std::abs(cht.interface_radiation_into_region2_W-watt_old)>
           controls.p1.balance_tolerance*std::max(1.0,std::abs(watt_old)))
            return report;
        report.cht_global_balance=accepted_audit.global;
        current_fluid=std::move(trial_fluid);
        current_solid=std::move(trial_solid);
        current_G=std::move(Gnew);
        current_Q=std::move(Qnew);
        // A residual on the frozen radiation RHS is not enough: source drift
        // from re-evaluating P1 at the accepted temperature must ALSO close.
        if(report.max_temperature_relative<=controls.coupling_tolerance &&
           report.max_radiative_source_relative<=controls.coupling_tolerance &&
           report.radiative_balance<=controls.p1.balance_tolerance &&
           report.cht_max_cell_residual<=interface.tolerance &&
           report.cht_global_balance<=interface.tolerance &&
           report.accepted_global_energy_residual<=interface.tolerance) {
            Tfluid=std::move(current_fluid);
            Tsolid=std::move(current_solid);
            irradiation=std::move(current_G);
            emission_minus_absorption=std::move(current_Q);
            report.converged=true;
            return report;
        }
    }
    return report;
}
} // namespace cfdx::physics::v41
