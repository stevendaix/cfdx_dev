#pragma once
// Read-only near-wall regime audit: y+ from Spalding U and from k,
// omegaWallFunction face-local references, and their mismatch with current
// owner omega. No wall BC is changed. This is NOT OpenFOAM equivalence.
#include "cfdx/physics/sst_openfoam_wall_coupling_v8.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
struct BfsWallSampleV74 {
    std::size_t face=0,owner=0;
    double distance=0.,speed_tangential=0.;
    double molecular_nu=0.,k=0.,omega=0.,nut_owner=0.;
};
struct BfsWallFaceV74 {
    std::size_t face=0,owner=0;
    double yplus_velocity=0.,yplus_k=0.;
    double omega_step=0.,omega_blended=0.;
    double omega_owner=0.,viscous_omega_6=0.,viscous_omega_60=0.;
    double nut_spalding=0.,nut_owner=0.;
    bool wall_resolved=false,buffer_or_transition=false,log_layer=false;
    bool out_of_model_range=false;
};
struct BfsWallRegimeV74 {
    std::vector<BfsWallFaceV74> faces;
    std::size_t wall_resolved=0,buffer_or_transition=0,log_layer=0;
    std::size_t out_of_model_range=0,corner_owner_count=0;
    double max_yplus_velocity=0.,min_yplus_velocity=std::numeric_limits<double>::infinity();
    double max_omega_step_relative_gap=0.,max_omega_blended_relative_gap=0.;
    double max_nut_wall_relative_gap=0.;
    bool same_operator_wall_triplet_certified=false;
    bool physical_bfs_validated=false;
};
inline BfsWallRegimeV74 bfs_wall_regime_v74(const std::vector<BfsWallSampleV74>& samples) {
    if(samples.empty())throw std::invalid_argument("V74 wall audit needs wall faces");
    BfsWallRegimeV74 result;
    std::map<std::size_t,std::size_t> face_count_by_owner;
    for(const auto& s:samples) {
        if(!(s.distance>0.) || !(s.molecular_nu>0.) || s.speed_tangential<0. ||
           s.k<0. || !(s.omega>0.) || s.nut_owner<0. ||
           !std::isfinite(s.distance) || !std::isfinite(s.molecular_nu) ||
           !std::isfinite(s.speed_tangential) || !std::isfinite(s.k) ||
           !std::isfinite(s.omega) || !std::isfinite(s.nut_owner))
            throw std::invalid_argument("V74 invalid wall reference sample");
        const double grad=s.speed_tangential/s.distance;
        const auto step=openfoam_omega_wall_reference(
            s.k,s.molecular_nu,s.distance,grad,ReferenceOmegaWallBlend::Stepwise);
        const auto blended=openfoam_omega_wall_reference(
            s.k,s.molecular_nu,s.distance,grad,ReferenceOmegaWallBlend::LaminarWeighted);
        const double wall_nut=spalding_wall_nut_v8(
            s.speed_tangential,s.distance,s.molecular_nu);
        double yplus=0.;
        if(s.speed_tangential>0.) {
            const double u_tau=std::sqrt(
                (wall_nut+s.molecular_nu)*s.speed_tangential/s.distance);
            yplus=s.distance*u_tau/s.molecular_nu;
        }
        if(!std::isfinite(yplus) || !std::isfinite(wall_nut))
            throw std::runtime_error("V74 nonfinite wall regime reference");
        BfsWallFaceV74 f;
        f.face=s.face;f.owner=s.owner;
        f.yplus_velocity=yplus;f.yplus_k=step.y_plus_from_k;
        f.omega_step=step.reference_omega;f.omega_blended=blended.reference_omega;
        f.omega_owner=s.omega;f.nut_owner=s.nut_owner;
        f.viscous_omega_6=step.viscous_omega;
        f.viscous_omega_60=10.*step.viscous_omega;
        f.nut_spalding=wall_nut;
        f.wall_resolved=yplus<=1.;
        f.buffer_or_transition=yplus>1. && yplus<30.;
        f.log_layer=yplus>=30. && yplus<=300.;
        f.out_of_model_range=yplus>300.;
        result.wall_resolved+=f.wall_resolved;
        result.buffer_or_transition+=f.buffer_or_transition;
        result.log_layer+=f.log_layer;
        result.out_of_model_range+=f.out_of_model_range;
        result.max_yplus_velocity=std::max(result.max_yplus_velocity,yplus);
        result.min_yplus_velocity=std::min(result.min_yplus_velocity,yplus);
        result.max_omega_step_relative_gap=std::max(
            result.max_omega_step_relative_gap,
            std::abs(s.omega-f.omega_step)/std::max(s.omega,f.omega_step));
        result.max_omega_blended_relative_gap=std::max(
            result.max_omega_blended_relative_gap,
            std::abs(s.omega-f.omega_blended)/std::max(s.omega,f.omega_blended));
        result.max_nut_wall_relative_gap=std::max(
            result.max_nut_wall_relative_gap,
            std::abs(s.nut_owner-wall_nut)/
                std::max({s.molecular_nu,s.nut_owner,wall_nut}));
        ++face_count_by_owner[s.owner];
        result.faces.push_back(f);
    }
    for(const auto& [owner,count]:face_count_by_owner)
        if(count>1)++result.corner_owner_count;
    return result;
}
} // namespace cfdx::physics
