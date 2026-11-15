#pragma once
// V76 experimental PRIVATE 4x4 SST Newton/block-regularized update.
// Four log unknowns: (k_a, omega_a, k_b, omega_b). A complete physical
// residual callback reassembles eddy viscosity, rAU, pressure and RC flux.
// No native linear certificates and no production commit are issued here.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace cfdx::physics {

struct BfsSstBlockReportV76 {
    bool attempted=false, selected=false, derivative_regular=true;
    std::string reason="not_requested";
    std::array<std::size_t,2> cells{};
    std::array<double,4> residual_before{}, residual_after{}, log_direction{};
    double scaled_local_before=0., scaled_local_after=0.;
    double accepted_alpha=0., max_branch_discrepancy=0.;
    unsigned int residual_evaluations=0, accepted_backtracks=0;
    bool native_sst_linear_certificate=false, native_momentum_linear_certificate=false;
    bool production_commit_authorized=false;
};

inline bool bfs_sst_solve4_v76(
    std::array<std::array<double,4>,4> A, std::array<double,4> rhs,
    std::array<double,4>& answer) {
    // Row scales also cover the k/omega units and avoid a fixed unit pivot.
    for(std::size_t i=0;i<4;++i){
        double scale=0.;
        for(const double x:A[i]){
            if(!std::isfinite(x)) return false;
            scale=std::max(scale,std::abs(x));
        }
        if(!(scale>0.) || !std::isfinite(rhs[i]))return false;
        for(auto& x:A[i])x/=scale;
        rhs[i]/=scale;
    }
    for(std::size_t k=0;k<4;++k){
        std::size_t pivot=k;
        for(std::size_t i=k+1;i<4;++i)
            if(std::abs(A[i][k])>std::abs(A[pivot][k]))pivot=i;
        if(!(std::abs(A[pivot][k])>1e-12))return false;
        std::swap(A[k],A[pivot]);std::swap(rhs[k],rhs[pivot]);
        for(std::size_t i=k+1;i<4;++i){
            const double mult=A[i][k]/A[k][k];
            for(std::size_t j=k;j<4;++j)A[i][j]-=mult*A[k][j];
            rhs[i]-=mult*rhs[k];
        }
    }
    answer.fill(0.);
    for(int i=3;i>=0;--i){
        double value=rhs[static_cast<std::size_t>(i)];
        for(std::size_t j=static_cast<std::size_t>(i)+1;j<4;++j)
            value-=A[static_cast<std::size_t>(i)][j]*answer[j];
        answer[static_cast<std::size_t>(i)]=value/A[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
        if(!std::isfinite(answer[static_cast<std::size_t>(i)]))return false;
    }
    return true;
}

template<class Candidate, class Residual, class Evaluate, class ExternalGuard>
BfsSstBlockReportV76 bfs_private_implicit_sst_block_v76(
    const Candidate& initial, const Residual& initial_f,
    std::size_t cell_a, std::size_t cell_b, double k_min,double omega_min,
    Evaluate&& eval,ExternalGuard&& external_guard,
    Candidate* accepted=nullptr, Residual* accepted_f=nullptr,
    double diagonal_regularization=0.,unsigned int max_halves=9) {
    BfsSstBlockReportV76 out;
    const std::size_t n=initial.cells;
    if(n==0 || cell_a>=n || cell_b>=n || cell_a==cell_b ||
       initial.fields.size()!=6*n || initial_f.signed_raw.size()!=6*n ||
       !(k_min>=0.) || !(omega_min>0.) || !(diagonal_regularization>=0.) ||
       !std::isfinite(k_min) || !std::isfinite(omega_min) ||
       !std::isfinite(diagonal_regularization) || max_halves>12)
        throw std::invalid_argument("V76 invalid private 4x4 SST block domain");
    out.cells={cell_a,cell_b};
    if(accepted)*accepted=initial;
    if(accepted_f)*accepted_f=initial_f;
    const std::array<std::size_t,4> idx{
        4*n+cell_a,5*n+cell_a,4*n+cell_b,5*n+cell_b};
    std::array<double,4> starting{};
    for(std::size_t i=0;i<4;++i){
        const double q=initial.fields[idx[i]];
        if(!(q>0.) || !std::isfinite(q)){
            out.reason="nonpositive_start";return out;
        }
        if(q<(i%2==0?k_min:omega_min)){
            out.reason="below_native_floor";return out;
        }
        starting[i]=std::log(q);
        out.residual_before[i]=initial_f.signed_raw[idx[i]];
        if(!std::isfinite(out.residual_before[i])){
            out.reason="nonfinite_base_residual";return out;
        }
    }
    out.attempted=true;
    // This is a LOCAL correction. The external guard must check the COMPLETE
    // physical residual and the historical V65 growth gates on each trial.
    auto evaluate_log=[&](const std::array<double,4>& z,
                          Candidate& cand,Residual& phys)->bool{
        cand=initial;
        for(std::size_t i=0;i<4;++i){
            if(!std::isfinite(z[i]) || z[i]>700. || z[i]<-700.)return false;
            const double q=std::exp(z[i]);
            if(!std::isfinite(q) || q<(i%2==0?k_min:omega_min))return false;
            cand.fields[idx[i]]=q;
        }
        try {
            phys=eval(cand);
            ++out.residual_evaluations;
            if(phys.signed_raw.size()!=6*n)return false;
            for(std::size_t i=0;i<4;++i)
                if(!std::isfinite(phys.signed_raw[idx[i]]))return false;
            return true;
        }catch(const std::exception&){return false;}
    };
    constexpr double h=1e-4;
    std::array<std::array<double,4>,4> J{};
    for(std::size_t j=0;j<4;++j){
        auto z_plus=starting,z_minus=starting;
        z_plus[j]+=h;z_minus[j]-=h;
        Candidate xp,xm;
        Residual rp,rm;
        if(!evaluate_log(z_plus,xp,rp) || !evaluate_log(z_minus,xm,rm)){
            out.reason="finite_difference_evaluation_failed";return out;
        }
        for(std::size_t i=0;i<4;++i){
            const double fw=(rp.signed_raw[idx[i]]-out.residual_before[i])/h;
            const double bw=(out.residual_before[i]-rm.signed_raw[idx[i]])/h;
            J[i][j]=(rp.signed_raw[idx[i]]-rm.signed_raw[idx[i]])/(2.*h);
            const double scale=std::max({1e-16,std::abs(fw),std::abs(bw)});
            const double mismatch=std::abs(fw-bw)/scale;
            out.max_branch_discrepancy=std::max(out.max_branch_discrepancy,mismatch);
            if(!std::isfinite(J[i][j])){
                out.reason="nonfinite_J";return out;
            }
        }
    }
    if(out.max_branch_discrepancy>0.5){
        out.derivative_regular=false;
        out.reason="nonsmooth_limiter_or_wall_branch";return out;
    }
    // Optional purely algebraic regularization of the local Newton matrix.
    // This is NOT physical pseudo-transient continuation: rho*V/dt and the
    // correct unknown-dependent transformation are still unimplemented.
    for(std::size_t i=0;i<4;++i)J[i][i]+=diagonal_regularization;
    std::array<double,4> rhs{};
    for(std::size_t i=0;i<4;++i)rhs[i]=-out.residual_before[i];
    if(!bfs_sst_solve4_v76(J,rhs,out.log_direction)){
        out.reason="singular_SST_block";return out;
    }
    // Limit the logarithmic step, then use Armijo-like backtracking with a
    // full reassembly. The outer gamma/alpha_p are NEVER changed by this.
    double direction_linf=0.;
    for(const double x:out.log_direction)direction_linf=std::max(direction_linf,std::abs(x));
    const double initial_alpha=direction_linf>0.5?0.5/direction_linf:1.;
    const auto local_norm=[&](const std::array<double,4>& r){
        double max_relative=0.;
        for(std::size_t i=0;i<4;++i){
            const double scale=std::max(1e-12,std::abs(out.residual_before[i]));
            max_relative=std::max(max_relative,std::abs(r[i])/scale);
        }
        return max_relative;
    };
    out.scaled_local_before=local_norm(out.residual_before);
    for(unsigned int back=0;back<=max_halves;++back){
        const double alpha=std::ldexp(initial_alpha,-static_cast<int>(back));
        auto z=starting;
        for(std::size_t i=0;i<4;++i)z[i]+=alpha*out.log_direction[i];
        Candidate cand;
        Residual phys;
        if(!evaluate_log(z,cand,phys))continue;
        std::array<double,4> after{};
        for(std::size_t i=0;i<4;++i)after[i]=phys.signed_raw[idx[i]];
        const double norm=local_norm(after);
        if(!(norm < out.scaled_local_before*(1.-1e-4*alpha)))continue;
        if(!external_guard(initial,initial_f,cand,phys))continue;
        out.selected=true;out.reason="private_sst_descent_only";
        out.residual_after=after;
        out.scaled_local_after=norm;
        out.accepted_alpha=alpha;out.accepted_backtracks=back;
        if(accepted)*accepted=std::move(cand);
        if(accepted_f)*accepted_f=std::move(phys);
        return out;
    }
    out.reason="no_full_physical_descent";
    return out;
}
} // namespace cfdx::physics
