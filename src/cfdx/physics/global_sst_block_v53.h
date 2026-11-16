#pragma once
// CFDX V53 experimental 6-field RIGHT preconditioner.
// Pmat is frozen for one Newton step. A Jv always comes from the independent
// physical candidate function (V51); Pmat never changes that function.
// Explicit scalar finite-volume SST advection/diffusion stencils are reused,
// while the LOCAL SST source/wall and momentum/SST cross blocks are probed at
// the SAME candidate. Nonlocal cross derivatives remain in matrix-free Jv.
// Application is linear and uses two native 4N Schur actions and one 2N ILU0.
#include "cfdx/physics/physical_six_candidate_v51.h"
#include "cfdx/core/linalg/advanced_preconditioners.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::physics {

struct GlobalSSTBlockStatsV53 {
    std::size_t colors=0;
    std::size_t candidate_evaluations=0;
    std::size_t transport_nnz=0;
    double maximum_locality_relative_error=0.0;
};

class GlobalSSTBlockV53 {
public:
    enum class Composition { BlockDiagonal, Lower, Upper, Symmetric };
    using Scales=std::array<double,6>;
    using NativeAction=std::function<std::vector<double>(
        const std::vector<double>&,const Scales&,const Scales&,std::size_t)>;
    using Evaluator=std::function<PhysicalResidualV51(const PhysicalCandidateV51&)>;

    // An adjacency entry must correspond to an actual internal FV face.
    // Color distance <= 3 to separate local source and two-hop gradients.
    static std::vector<std::size_t> distance_three_colors(
        const std::vector<std::vector<std::size_t>>& adjacency) {
        const std::size_t n=adjacency.size();
        if(!n)throw std::invalid_argument("V53 empty mesh adjacency");
        std::vector<std::size_t> colors(n,n);
        for(std::size_t i=0;i<n;++i){
            std::vector<bool> blocked(n,false),seen(n,false);
            std::queue<std::pair<std::size_t,int>> frontier;
            frontier.push({i,0});seen[i]=true;
            while(!frontier.empty()) {
                const auto [node,dist]=frontier.front();frontier.pop();
                if(node!=i && colors[node]<n)blocked[colors[node]]=true;
                if(dist==3)continue;
                for(const auto next:adjacency[node]) {
                    if(next>=n || next==node)
                        throw std::invalid_argument("V53 inconsistent FV adjacency");
                    if(!seen[next]) {seen[next]=true;frontier.push({next,dist+1});}
                }
            }
            std::size_t c=0;
            while(c<n && blocked[c])++c;
            if(c==n)throw std::logic_error("V53 coloring exhausted");
            colors[i]=c;
        }
        return colors;
    }

    GlobalSSTBlockV53(
        const PhysicalCandidateV51& base,const PhysicalResidualV51& residual,
        const Scales& state_scale,const Scales& row_scale,
        const cfdx::core::SparseMatrix& k_fv,
        const cfdx::core::SparseMatrix& omega_fv,
        const std::vector<std::vector<std::size_t>>& adjacency,
        const Evaluator& evaluate,const NativeAction& native,
        bool pressure_gauge,std::size_t pressure_pin,
        double normalized_difference=1e-6,
        Composition composition=Composition::Symmetric)
      : n_(base.cells), composition_(composition),state_scale_(state_scale),row_scale_(row_scale),
        native_(native),pressure_gauge_(pressure_gauge),pressure_pin_(pressure_pin) {
        validate_physical_candidate_v51(base);
        validate_physical_result_v51(base,residual);
        if(!native_ || !evaluate || adjacency.size()!=n_ ||
           !std::isfinite(normalized_difference) ||
           normalized_difference<1e-9 || normalized_difference>1e-3 ||
           (pressure_gauge_ && pressure_pin_>=n_))
            throw std::invalid_argument("V53 missing preconditioner contract");
        for(std::size_t b=0;b<6;++b)
            if(!(state_scale_[b]>0.) || !std::isfinite(state_scale_[b]) ||
               !(row_scale_[b]>0.) || !std::isfinite(row_scale_[b]))
                throw std::invalid_argument("V53 nonfinite frozen scales");
        const auto colors=distance_three_colors(adjacency);
        stats_.colors=*std::max_element(colors.begin(),colors.end())+1;
        std::array<std::vector<std::array<double,8>>,2> cross_data;
        // [0]: local B_4x2, [1]: local C_2x4 stored per-cell.
        cross_data[0].resize(n_);cross_data[1].resize(n_);
        std::vector<std::array<double,4>> diagonal(n_);
        for(std::size_t block=0;block<6;++block) {
            for(std::size_t color=0;color<stats_.colors;++color) {
                PhysicalCandidateV51 shifted=base;
                bool used=false;
                for(std::size_t i=0;i<n_;++i)if(colors[i]==color) {
                    shifted.fields[block*n_+i]+=normalized_difference*state_scale_[block];
                    used=true;
                }
                if(!used)continue;
                const auto perturbed=evaluate(shifted);
                validate_physical_result_v51(shifted,perturbed);
                ++stats_.candidate_evaluations;
                for(std::size_t i=0;i<n_;++i)if(colors[i]==color){
                    for(std::size_t out=0;out<6;++out){
                        if(block<4 && out<4)continue;
                        {
                            double d=(perturbed.signed_raw[out*n_+i]-
                                residual.signed_raw[out*n_+i])/
                                (normalized_difference*row_scale_[out]);
                            if(pressure_gauge_ && out==3 && i==pressure_pin_)d=0.0;
                            if(!std::isfinite(d))
                                throw std::runtime_error("V53 nonfinite local cross probe");
                            if(block>=4 && out>=4)
                                diagonal[i][(out-4)*2+(block-4)]=d;
                            else if(block>=4 && out<4)
                                cross_data[0][i][out*2+(block-4)]=d;
                            else if(block<4 && out>=4)
                                cross_data[1][i][(out-4)*4+block]=d;
                        }
                    }
                }
            }
        }
        B_=std::move(cross_data[0]);C_=std::move(cross_data[1]);
        // Use the complete scalar-FV transport stencil from the SAME physical
        // state. Replace local diagonal by the V51 finite-difference source,
        // wall and viscosity response. Only nonlocal cross couplings are
        // approximated/omitted from Pmat, never from matrix-free Jv.
        cfdx::core::SparseMatrix full(2*n_,2*n_);
        const std::array<const cfdx::core::SparseMatrix*,2> fv{{&k_fv,&omega_fv}};
        for(std::size_t part=0;part<2;++part) {
            const auto& A=*fv[part];
            if(A.n_rows()!=n_ || A.n_cols()!=n_ || !A.is_consistent())
                throw std::invalid_argument("V53 SST FV operator is not valid CSR");
            const auto* ro=A.row_offsets_data();
            const auto* ci=A.columns_data();
            const auto* av=A.values_data();
            for(std::size_t i=0;i<n_;++i) {
                for(auto at=ro[i];at<ro[i+1];++at) {
                    if(ci[at]>=n_ || !std::isfinite(av[at]))
                        throw std::invalid_argument("V53 bad FV stencil");
                    if(ci[at]==i)continue;
                    const double v=av[at]*state_scale_[4+part]/row_scale_[4+part];
                    if(v!=0.)full.push_back(part*n_+i,part*n_+ci[at],v);
                }
                for(std::size_t col=0;col<2;++col)
                    full.push_back(part*n_+i,col*n_+i,diagonal[i][part*2+col]);
            }
        }
        full.finalize();
        stats_.transport_nnz=full.nnz();
        transport_=std::make_unique<cfdx::core::ILU0Preconditioner>();
        if(!transport_->setup(full))
            throw std::runtime_error("V53 SST transport ILU0 setup failed (no fallback)");
        // Sample independent probes: the local row in a color must not
        // change because of another simultaneously perturbed cell. This does
        // not establish global sparsity of the complete Jacobian.
        const std::size_t samples=std::min<std::size_t>(n_,6);
        for(std::size_t sample=0;sample<samples;++sample) {
            const std::size_t cell=sample*(n_-1)/std::max<std::size_t>(1,samples-1);
            const std::size_t block=sample%6;
            PhysicalCandidateV51 x=base;
            x.fields[block*n_+cell]+=normalized_difference*state_scale_[block];
            const auto r=evaluate(x);
            validate_physical_result_v51(x,r);
            ++stats_.candidate_evaluations;
            for(std::size_t out=0;out<6;++out) {
                if(block<4 && out<4)continue;
                const double actual=(r.signed_raw[out*n_+cell]-
                    residual.signed_raw[out*n_+cell])/
                    (normalized_difference*row_scale_[out]);
                double compressed=0.;
                if(block>=4 && out>=4)
                    compressed=diagonal[cell][(out-4)*2+(block-4)];
                else if(block>=4)compressed=B_[cell][out*2+(block-4)];
                else compressed=C_[cell][(out-4)*4+block];
                if(pressure_gauge_ && out==3 && cell==pressure_pin_)continue;
                const double discrepancy=std::abs(actual-compressed)/
                    std::max({1e-10,std::abs(actual),std::abs(compressed)});
                stats_.maximum_locality_relative_error=std::max(
                    stats_.maximum_locality_relative_error,discrepancy);
            }
        }
        if(stats_.maximum_locality_relative_error>1e-2)
            throw std::runtime_error("V53 colored local cross derivatives contaminated");
    }

    const GlobalSSTBlockStatsV53& statistics() const noexcept{return stats_;}
    std::vector<double> apply(const std::vector<double>& rhs) const {
        if(rhs.size()!=6*n_)throw std::invalid_argument("V53 incomplete 6N RHS");
        for(const double x:rhs)if(!std::isfinite(x))
            throw std::invalid_argument("V53 nonfinite 6N RHS");
        const auto schur=[&](const std::vector<double>& r){
            auto z=native_(r,state_scale_,row_scale_,n_);
            if(z.size()!=6*n_)throw std::runtime_error("V53 native Schur wrong size");
            for(double x:z)if(!std::isfinite(x))
                throw std::runtime_error("V53 nonfinite native Schur output");
            return z;
        };
        const auto solve_transport=[&](const std::vector<double>& r) {
            cfdx::core::Vector v(2*n_,0.0),solution(2*n_,0.0);
            for(std::size_t i=0;i<2*n_;++i)v(i)=r[4*n_+i];
            if(!transport_->apply(v,solution)||!solution.is_valid())
                throw std::runtime_error("V53 SST ILU0 solve failed");
            std::vector<double> z(2*n_);
            for(std::size_t i=0;i<2*n_;++i)z[i]=solution(i);
            return z;
        };
        const auto only4=[&](const std::vector<double>& r){
            auto a=r;
            std::fill(a.begin()+static_cast<std::ptrdiff_t>(4*n_),a.end(),0.0);
            return schur(a);
        };
        std::vector<double> z(6*n_,0.0);
        if(composition_==Composition::Upper){
            const auto ss=solve_transport(rhs);
            auto r=rhs;
            for(std::size_t i=0;i<n_;++i)
                for(std::size_t m=0;m<4;++m)
                    for(std::size_t s=0;s<2;++s)
                        r[m*n_+i]-=B_[i][m*2+s]*ss[s*n_+i];
            z=only4(r);
            std::copy(ss.begin(),ss.end(),z.begin()+static_cast<std::ptrdiff_t>(4*n_));
        }else{
            z=only4(rhs);
            auto r=rhs;
            if(composition_!=Composition::BlockDiagonal)
                for(std::size_t i=0;i<n_;++i)
                    for(std::size_t s=0;s<2;++s)
                        for(std::size_t m=0;m<4;++m)
                            r[(4+s)*n_+i]-=C_[i][s*4+m]*z[m*n_+i];
            const auto ss=solve_transport(r);
            std::copy(ss.begin(),ss.end(),z.begin()+static_cast<std::ptrdiff_t>(4*n_));
            if(composition_==Composition::Symmetric){
                std::vector<double> correction(6*n_,0.0);
                for(std::size_t i=0;i<n_;++i)
                    for(std::size_t m=0;m<4;++m)
                        for(std::size_t s=0;s<2;++s)
                            correction[m*n_+i]-=B_[i][m*2+s]*ss[s*n_+i];
                const auto feedback=only4(correction);
                for(std::size_t i=0;i<4*n_;++i)z[i]+=feedback[i];
            }
        }
        for(double x:z)if(!std::isfinite(x))
            throw std::runtime_error("V53 nonfinite global preconditioner action");
        return z;
    }
private:
    std::size_t n_=0;
    Composition composition_;
    Scales state_scale_{};
    Scales row_scale_{};
    NativeAction native_;
    bool pressure_gauge_=false;
    std::size_t pressure_pin_=0;
    GlobalSSTBlockStatsV53 stats_;
    std::vector<std::array<double,8>> B_,C_;
    std::unique_ptr<cfdx::core::ILU0Preconditioner> transport_;
};

} // namespace cfdx::physics
