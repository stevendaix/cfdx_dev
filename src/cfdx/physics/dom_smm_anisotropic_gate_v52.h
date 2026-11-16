#pragma once
// V52 research-only SMM/AMG anisotropic qualification.
// The auxiliary SPD AMG hierarchy is not the true WLS second-moment matrix.
// Two separate gates: strict true GMRES residual AND independently reconstructed
// physical face/cell/global radiation energy; neither can imply the other.
// No production solver registration, no tolerance changes in V49/V50/V51.
#include "cfdx/physics/dom_smm_gray_wall_v49.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {
struct DomSmmAnisotropicV52Options {
    std::size_t max_krylov=300;
    int restart=30;
    double linear_tolerance=1.e-13;
    double physical_tolerance=1.e-10;
    cfdx::core::AMGInterpolationPolicy interpolation=
        cfdx::core::AMGInterpolationPolicy::DirectionalSemiCoarsening;
    cfdx::core::AMGSmootherPolicy smoother=
        cfdx::core::AMGSmootherPolicy::AlgebraicLineGaussSeidel;
};
struct DomSmmAnisotropicV52Result {
    bool amg_setup_ok=false;
    bool linear_strict_pass=false;
    bool physical_pass=false;
    bool qualified=false;
    std::size_t iterations=0,amg_levels=0;
    double true_linear_l2=std::numeric_limits<double>::infinity();
    double true_linear_max=std::numeric_limits<double>::infinity();
    double max_cell_balance=std::numeric_limits<double>::infinity();
    double global_energy_balance=std::numeric_limits<double>::infinity();
    double setup_seconds=0.,total_seconds=0.;
    std::string diagnostic;
    // Never publish unqualified data. For negative tests, metrics alone remain.
    std::vector<double> accepted_J,accepted_wall_watts;
};
inline DomSmmAnisotropicV52Result solve_smm_anisotropic_v52(
        const cfdx::core::Mesh& mesh,const FvGeometry& geo,
        const DomSmmGrayWallV49& system,
        const std::vector<DomSmmSecondMomentV48::Tensor>& E,
        const std::vector<double>& extinction,
        const std::vector<double>& absorption,
        const std::vector<double>& gas_blackbody,
        const DomSmmAnisotropicV52Options& options={}) {
    using namespace cfdx::core;
    const auto nc=mesh.n_cells(),nf=mesh.n_faces();
    if(nc!=system.dimension() || E.size()!=nc || extinction.size()!=nc ||
       absorption.size()!=nc || gas_blackbody.size()!=nc ||
       geo.cell_volumes.size()!=nc || options.max_krylov==0 || options.restart<2 ||
       !(options.linear_tolerance>0.) || !std::isfinite(options.linear_tolerance) ||
       !(options.physical_tolerance>0.) || !std::isfinite(options.physical_tolerance))
        throw std::invalid_argument("V52 SMM: invalid audit dimensions or controls");
    if(!system.frozen_operator_matches_v52(mesh,geo,E,extinction,absorption))
        throw std::invalid_argument("V52 SMM: physical auditor and solved operator differ");
    DomSmmAnisotropicV52Result report;
    const auto started=std::chrono::steady_clock::now();
    const auto seconds=[&](){return std::chrono::duration<double>(
        std::chrono::steady_clock::now()-started).count();};
    const auto bvalues=system.rhs(gas_blackbody);
    Vector b(nc,0.),x(nc,0.);
    for(std::size_t c=0;c<nc;++c){
        b(c)=bvalues[c];x(c)=gas_blackbody[c];
    }
    LinearOperator matvec;
    matvec.size=nc;
    matvec.apply=[&system,nc](const Vector& input,Vector& output){
        std::vector<double> v(nc);
        for(std::size_t c=0;c<nc;++c)v[c]=input(c);
        const auto Av=system.apply(v);
        if(output.size()!=nc)output.resize(nc);
        for(std::size_t c=0;c<nc;++c)output(c)=Av[c];
    };
    FunctionalLinearOperator true_wls(nc,matvec.apply);
    const auto auxiliary=system.auxiliary_matrix();
    MatrixFreeVcyclePreconditioner amg(true_wls,.55,2,2,.25,12,
        options.interpolation,false,false,options.smoother);
    if(!amg.setup(auxiliary)){
        report.diagnostic="V52 SMM: AMG setup failed: "+amg.last_setup_error();
        report.total_seconds=seconds();return report;
    }
    report.amg_setup_ok=true;
    report.amg_levels=amg.hierarchy_level_sizes().size();
    report.setup_seconds=seconds();
    const auto linear=solve_gmres(matvec,b,x,options.restart,options.max_krylov,
                                 options.linear_tolerance,&amg);
    report.iterations=linear.iterations;
    // Audit the ORIGINAL operator, not the AMG auxiliary matrix and not the
    // recursive Hessenberg residual. The candidate may have been rejected
    // by GMRES; this audit must never convert it to a published state.
    std::vector<double> candidate(nc);
    for(std::size_t c=0;c<nc;++c){
        candidate[c]=x(c);
        if(!std::isfinite(candidate[c]) || !(candidate[c]>=0.)){
            report.diagnostic="V52 SMM: invalid candidate J";
            report.total_seconds=seconds();return report;
        }
    }
    const auto Ax=system.apply(candidate);
    long double r2=0.,b2=0.;
    double maxr=0.,maxb=1.e-300;
    for(std::size_t c=0;c<nc;++c){
        const long double defect=static_cast<long double>(Ax[c])-bvalues[c];
        r2+=defect*defect;
        b2+=static_cast<long double>(bvalues[c])*bvalues[c];
        maxr=std::max(maxr,std::abs(static_cast<double>(defect)));
        maxb=std::max(maxb,std::abs(bvalues[c]));
    }
    if(!(b2>0.) || !std::isfinite(static_cast<double>(b2))) {
        report.diagnostic="V52 SMM: no finite RHS norm for relative audit";
        report.total_seconds=seconds();return report;
    }
    report.true_linear_l2=std::sqrt(static_cast<double>(r2/b2));
    report.true_linear_max=maxr/maxb;
    report.linear_strict_pass=linear.status==SolverStatus::CONVERGED &&
        report.true_linear_l2<=options.linear_tolerance;
    std::vector<double> wall_current(nf,0.);
    for(const auto& w:system.walls())
        wall_current[w.face]=w.conductance*(candidate[w.owner]-w.ib);
    const DomSmmSecondMomentV48 smm(mesh,geo);
    const auto physical=smm.evaluate(candidate,E,extinction,wall_current);
    constexpr double fourpi=12.56637061435917295385;
    long double global=0.,scale=1.;
    report.max_cell_balance=0.;
    for(std::size_t c=0;c<nc;++c){
        const double v=fourpi*absorption[c]*geo.cell_volumes[c]*
            (candidate[c]-gas_blackbody[c]);
        const double current=physical.divergence_watts[c];
        report.max_cell_balance=std::max(report.max_cell_balance,
            std::abs(current+v)/std::max({1.,std::abs(current),std::abs(v)}));
        global+=v;scale+=std::abs(v);
    }
    for(const auto& w:system.walls()){
        const auto current=wall_current[w.face];
        if(w.emissivity==0. && current!=0.){
            report.diagnostic="V52 SMM: reflective wall leaking power";
            report.total_seconds=seconds();return report;
        }
        global+=current;scale+=std::abs(current);
    }
    report.global_energy_balance=std::abs(static_cast<double>(global/scale));
    report.physical_pass=report.max_cell_balance<=options.physical_tolerance &&
        report.global_energy_balance<=options.physical_tolerance &&
        report.true_linear_max<=options.physical_tolerance;
    report.qualified=report.physical_pass && report.linear_strict_pass;
    report.diagnostic=report.qualified?"V52 SMM: both independent gates pass":
        report.physical_pass?"V52 SMM: physical pass, algebraic strict gate FAIL":
                             "V52 SMM: physical energy gate FAIL";
    if(report.qualified){
        report.accepted_J=std::move(candidate);
        report.accepted_wall_watts=std::move(wall_current);
    }
    report.total_seconds=seconds();
    return report;
}
} // namespace cfdx::physics
