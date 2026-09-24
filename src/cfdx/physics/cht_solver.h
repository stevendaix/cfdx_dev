#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/physics/energy_solver.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct ChtInterfaceControls {
    std::string region1_patch;
    std::string region2_patch;
    double conductivity1 = 1.0;
    double conductivity2 = 1.0;
    double tolerance = 1e-8;
    std::size_t max_iterations = 100;
    double relaxation = 0.8;
    double matching_tolerance = 1e-8;
    double area_relative_tolerance = 1e-6;
};

struct ChtInterfaceFacePair {
    std::size_t face1 = 0;
    std::size_t face2 = 0;
    std::size_t cell1 = 0;
    std::size_t cell2 = 0;
};

inline std::vector<ChtInterfaceFacePair> match_cht_interface(
    const cfdx::core::Mesh& mesh1, const FvGeometry& g1,
    const cfdx::core::Mesh& mesh2, const FvGeometry& g2,
    const ChtInterfaceControls& c)
{
    if(c.region1_patch.empty() || c.region2_patch.empty() ||
       c.conductivity1<=0.0 || c.conductivity2<=0.0 ||
       c.matching_tolerance<=0.0 || c.area_relative_tolerance<=0.0)
        throw std::invalid_argument("invalid CHT interface controls");

    std::size_t p1=mesh1.boundary().n_patches();
    std::size_t p2=mesh2.boundary().n_patches();
    for(std::size_t p=0;p<p1;++p)
        if(mesh1.boundary().patch(p).name==c.region1_patch) { p1=p; break; }
    for(std::size_t p=0;p<p2;++p)
        if(mesh2.boundary().patch(p).name==c.region2_patch) { p2=p; break; }
    if(p1>=mesh1.boundary().n_patches() || p2>=mesh2.boundary().n_patches())
        throw std::invalid_argument("CHT interface patch not found");

    const auto& f1=mesh1.boundary().patch(p1).face_ids;
    const auto& f2=mesh2.boundary().patch(p2).face_ids;
    if(f1.size()!=f2.size())
        throw std::invalid_argument("CHT interface face count mismatch");

    std::vector<bool> used(f2.size(),false);
    std::vector<ChtInterfaceFacePair> pairs;
    pairs.reserve(f1.size());

    for(const auto face1:f1) {
        double best=std::numeric_limits<double>::infinity();
        std::size_t bestj=f2.size();
        for(std::size_t j=0;j<f2.size();++j) {
            if(used[j]) continue;
            const double d=(g1.face_centres[face1]-g2.face_centres[f2[j]]).mag();
            if(d<best) { best=d; bestj=j; }
        }
        if(bestj==f2.size() || best>c.matching_tolerance)
            throw std::runtime_error("CHT interface face matching failed");
        const double a1=g1.face_area_vectors[face1].mag();
        const double a2=g2.face_area_vectors[f2[bestj]].mag();
        const double area_scale=std::max({a1,a2,1e-30});
        if(std::abs(a1-a2)/area_scale>c.area_relative_tolerance)
            throw std::runtime_error("CHT interface face area mismatch");
        used[bestj]=true;
        pairs.push_back({
            face1,f2[bestj],
            mesh1.ownership().owner(face1),
            mesh2.ownership().owner(f2[bestj])
        });
    }
    return pairs;
}

struct ChtSolveResult {
    bool converged=false;
    std::size_t iterations=0;
    double interface_imbalance=0.0;
    double interface_temperature_change=0.0;
};

inline ChtSolveResult solve_two_region_cht(
    const cfdx::core::Mesh& mesh1, const FvGeometry& g1,
    const cfdx::core::Mesh& mesh2, const FvGeometry& g2,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& flux1,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& flux2,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& T1,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& T2,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& source1,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& source2,
    const EnergySolverControls& energy1,
    const EnergySolverControls& energy2,
    const ChtInterfaceControls& controls,
    const ScalarBoundaryConditions& bcs1 = {},
    const ScalarBoundaryConditions& bcs2 = {})
{
    const auto pairs=match_cht_interface(mesh1,g1,mesh2,g2,controls);
    ChtSolveResult result;
    std::vector<double> previous_interface_temperature(pairs.size(), std::numeric_limits<double>::quiet_NaN());

    for(std::size_t iter=1;iter<=controls.max_iterations;++iter) {
        ScalarBoundaryFaceValues fv1,fv2;

        // Compute the exact two-layer interface temperature for the current
        // cell-centre values. This enforces equal normal heat flux on both
        // sides in the fixed-point limit.
        for(std::size_t i=0;i<pairs.size();++i) {
            const auto& p=pairs[i];
            const double d1=(g1.face_centres[p.face1]-g1.cell_centres[p.cell1]).mag();
            const double d2=(g2.face_centres[p.face2]-g2.cell_centres[p.cell2]).mag();
            const double h1=controls.conductivity1/d1;
            const double h2=controls.conductivity2/d2;
            const double Tint_new=(h1*T1(p.cell1)+h2*T2(p.cell2))/(h1+h2);
            double Tint=Tint_new;
            if(std::isfinite(previous_interface_temperature[i]))
                Tint=controls.relaxation*Tint_new +
                     (1.0-controls.relaxation)*previous_interface_temperature[i];
            previous_interface_temperature[i]=Tint;
            auto& values1=fv1.values[controls.region1_patch];
            auto& values2=fv2.values[controls.region2_patch];
            if(values1.size()!=mesh1.n_faces()) values1.resize(mesh1.n_faces(),0.0);
            if(values2.size()!=mesh2.n_faces()) values2.resize(mesh2.n_faces(),0.0);
            values1[p.face1]=Tint;
            values2[p.face2]=Tint;
        }

        auto r1=solve_energy(mesh1,g1,flux1,T1,source1,energy1,bcs1,&fv1);
        auto r2=solve_energy(mesh2,g2,flux2,T2,source2,energy2,bcs2,&fv2);
        if(!r1.converged || !r2.converged)
            throw std::runtime_error("CHT region energy solve did not converge");

        double qimb=0.0, qscale=1.0, dtint=0.0;
        for(std::size_t i=0;i<pairs.size();++i) {
            const auto& p=pairs[i];
            const double d1=(g1.face_centres[p.face1]-g1.cell_centres[p.cell1]).mag();
            const double d2=(g2.face_centres[p.face2]-g2.cell_centres[p.cell2]).mag();
            const double Tint=fv1.values.at(controls.region1_patch)[p.face1];
            const double qflux1=controls.conductivity1*
                (T1(p.cell1)-Tint)/d1;
            const double qflux2=controls.conductivity2*
                (Tint-T2(p.cell2))/d2;
            qimb=std::max(qimb,std::abs(qflux1-qflux2));
            qscale=std::max(qscale,std::abs(qflux1));
            qscale=std::max(qscale,std::abs(qflux2));
            const double current_tint=fv1.values.at(controls.region1_patch)[p.face1];
            if(std::isfinite(previous_interface_temperature[i]))
                dtint=std::max(dtint,std::abs(current_tint-previous_interface_temperature[i]));
        }
        result.interface_imbalance=qimb/qscale;
        result.interface_temperature_change=dtint;
        result.iterations=iter;

        if(result.interface_imbalance<=controls.tolerance &&
           result.interface_temperature_change<=controls.matching_tolerance) {
            result.converged=true;
            break;
        }
    }
    return result;
}

} // namespace cfdx::physics
