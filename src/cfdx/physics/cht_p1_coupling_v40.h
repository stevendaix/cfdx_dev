#pragma once
// CFDX V40 opt-in experiment: nonconformal monolithic CHT / gray participating
// P1 with face-resolved Marshak exchange. Volume emission is subtracted from
// the fluid ENERGY source; the SAME P1 wall watts are deposited on the solid
// interface overlap owners. All fields are transactional until the final
// nonlinearly re-evaluated fixed-point gates pass. Restricted to a fully closed
// two-region interface with reflecting non-interface P1 walls.
#include "cfdx/physics/cht_nonconformal_v37.h"
#include "cfdx/physics/p1_marshak_v39.h"
#include <algorithm>
#include <cmath>
#include <limits>

#include <vector>

namespace cfdx::physics::v40 {
using namespace cfdx::core;
using namespace cfdx::physics::v37;
struct ChtP1Controls {
    double emissivity=0.85;
    double absorption=0.7;
    double scattering=1.3;
    double coupling_tolerance=1.e-10;
    double outer_relaxation=0.1;
    std::size_t max_outer_iterations=80;
    P1MarshakControls p1;
};
struct ChtP1Report {
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
};
inline ChtP1Report solve_nonconformal_cht_p1(
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
    const ChtP1Controls& controls={},
    const ScalarBoundaryConditions& thermal_fluid={},
    const ScalarBoundaryConditions& thermal_solid={})
{
    if(controls.max_outer_iterations==0 || !(controls.emissivity>0&&controls.emissivity<=1) ||
       !(controls.absorption>0) || !(controls.scattering>=0) ||
       !(controls.coupling_tolerance>0) ||
       !(controls.outer_relaxation>0&&controls.outer_relaxation<=1) ||
       !std::isfinite(controls.emissivity)||!std::isfinite(controls.absorption)||
       !std::isfinite(controls.scattering)||!std::isfinite(controls.coupling_tolerance))
        throw std::invalid_argument("CHT-P1 V40 invalid controls");
    if(Tfluid.size()!=mfluid.n_cells() || Tsolid.size()!=msolid.n_cells() ||
       irradiation.size()!=mfluid.n_cells() ||
       emission_minus_absorption.size()!=mfluid.n_cells() ||
       base_source_fluid.size()!=mfluid.n_cells() ||
       base_source_solid.size()!=msolid.n_cells())
        throw std::invalid_argument("CHT-P1 V40 incompatible fields");
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
    ChtP1Report report;report.interface_overlaps=map.overlaps.size();
    // The region pair is solved once for every staged radiation load, but the
    // original physical time level must remain frozen throughout this loop.
    auto evaluate_radiation=[&](const Field<double,Location::CELL>& f,
                                const Field<double,Location::CELL>& s,
                                Field<double,Location::CELL>& G,
                                Field<double,Location::CELL>& Q,
                                Field<double,Location::CELL>& srcf,
                                Field<double,Location::CELL>& srcs,
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
        long double total_radiation_cell=0.0L, gross_radiation_cell=0.0L;
        for(std::size_t i=0;i<mfluid.n_cells();++i){
            srcf(i)-=Q(i); // physical matter energy source is -emission
            const long double power=static_cast<long double>(Q(i))*gfluid.cell_volumes[i];
            total_radiation_cell+=power;
            gross_radiation_cell+=std::abs(power);
        }
        long double watt_solid=0.0L, gross_wall=0.0L;
        for(const auto& e:map.overlaps){
            const double Af=gfluid.face_area_vectors[e.face1].mag();
            const double frac=e.area/Af;
            const double watt=rad.face_outward_W.at(e.face1)*frac;
            srcs(e.cell2)+=watt/gsolid.cell_volumes[e.cell2];
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
        double watt_old=0,balance_old=0;
        if(!evaluate_radiation(current_fluid,current_solid,Gpred,Qpred,
                               srcf,srcs,watt_old,balance_old))return report;
        auto trial_fluid=current_fluid,trial_solid=current_solid;
        const auto cht=solve_nonconformal_cht(
            mfluid,gfluid,msolid,gsolid,phi_fluid,phi_solid,
            trial_fluid,trial_solid,srcf,srcs,energyfluid,energysolid,
            interface,thermal_fluid,thermal_solid,nullptr,
            &original_fluid,&original_solid);
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
        double watt_new=0,balance_new=0;
        if(!evaluate_radiation(trial_fluid,trial_solid,Gnew,Qnew,
                               srcf_new,srcs_new,watt_new,balance_new))return report;
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
        report.cht_max_cell_residual=cht.maximum_cell_balance;
        report.cht_global_balance=cht.global_balance;
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
           report.cht_global_balance<=interface.tolerance) {
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
} // namespace cfdx::physics::v40
