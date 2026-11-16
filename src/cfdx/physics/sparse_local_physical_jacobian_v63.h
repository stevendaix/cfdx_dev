#pragma once
// V63: opt-in SPARSE colored-Fphys Jacobian experiment, not an analytic
// differentiator nor a certified generalized Newton matrix. It differs from
// V54's O((6N)^2) dense Jacobian: locality comes from an explicit cell graph.
// All coefficients differentiate the complete candidate-owned Fphys with
// viscosity, Rhie-Chow, momentum and SST rebuilt on every probe. Missing
// long-range dependencies, active-set changes and roundoff must FAIL CLOSED.
// Sparse CSR and a fixed color schedule are immutable during one Arnoldi run.

#include "cfdx/physics/physical_six_candidate_v51.h"
#include "cfdx/physics/generalized_anchor_v62.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct SparseLocalJacobianStatsV63 {
    std::size_t colors=0,physical_evaluations=0,nonzeros=0,rejected_columns=0;
    std::size_t largest_color=0;
    double maximum_outside_stencil=0.0,maximum_sided=0.0;
    double anchor_displacement=0.0,anchor_normalized_distance=0.0;
    bool trusted=true;
};
struct SparseLocalDirectionalAuditV63 {
    double fidelity=0.,halving=0.,sided=0.;
    std::array<double,6> block_fidelity{};
    std::array<double,6> block_sided{};
    std::size_t evaluations=0;
};

class SparseLocalPhysicalJacobianV63 {
public:
    using Scale=std::array<double,6>;
    using Graph=std::vector<std::vector<std::size_t>>;
    template<class Evaluate>
    SparseLocalPhysicalJacobianV63(const PhysicalCandidateV51& base,
        const PhysicalResidualV51& original,const Scale& state_scale,
        const Scale& row_scale,const Graph& adjacency,std::size_t radius,
        double step,bool pin_pressure,std::size_t reference_cell,
        double reference_pressure,double consistency_tolerance,
        Evaluate&& evaluate,double k_min=0.,double omega_min=0.,
        bool use_common_anchor=false,double anchor_tolerance=0.05)
        : base_(base),states_(state_scale),rowscale_(row_scale),
          step_(step),reference_pressure_(reference_pressure),
          kmin_(k_min),wmin_(omega_min),pin_(pin_pressure),
          reference_(reference_cell) {
        validate_physical_candidate_v51(base_,kmin_,wmin_);
        validate_physical_result_v51(base_,original);
        n_=base_.cells;size_=6*n_;
        if(adjacency.size()!=n_ || radius<1 || radius>3 ||
            !(step>0.) || !std::isfinite(step) || reference_>=n_ ||
            !(consistency_tolerance>0.) || consistency_tolerance>0.05)
            throw std::invalid_argument("V63 invalid graph, step or strict tolerance");
        if(!std::isfinite(anchor_tolerance) || !(anchor_tolerance>0.) ||
           anchor_tolerance>0.5)
            throw std::invalid_argument("V63 invalid anchored fidelity guard");
        for(std::size_t b=0;b<6;++b)
            if(!(states_[b]>0.) || !std::isfinite(states_[b]) ||
               !(rowscale_[b]>0.) || !std::isfinite(rowscale_[b]))
                throw std::invalid_argument("V63 invalid frozen scales");
        for(std::size_t c=0;c<n_;++c)
            for(auto next:adjacency[c]) {
                if(next>=n_ || next==c ||
                   std::find(adjacency[next].begin(),adjacency[next].end(),c)==
                       adjacency[next].end())
                    throw std::invalid_argument("V63 asymmetric/invalid face graph");
            }
        support_.resize(n_);
        for(std::size_t c=0;c<n_;++c) {
            std::vector<int> distance(n_,-1);
            distance[c]=0;
            std::vector<std::size_t> queue{c};
            for(std::size_t p=0;p<queue.size();++p) {
                const auto u=queue[p];
                if(static_cast<std::size_t>(distance[u])==radius)continue;
                for(auto next:adjacency[u])if(distance[next]<0) {
                    distance[next]=distance[u]+1;queue.push_back(next);
                }
            }
            support_[c]=std::move(queue);
            std::sort(support_[c].begin(),support_[c].end());
        }
        // Greedy distance-2*radius coloring of the exact column supports.
        // Pressure/velocity/SST blocks cannot share a color on the same row.
        std::vector<std::vector<unsigned char>> covered;
        std::vector<std::vector<std::size_t>> groups;
        for(std::size_t block=0;block<6;++block)
            for(std::size_t c=0;c<n_;++c) {
                std::size_t color=0;
                for(;color<groups.size();++color) {
                    bool overlap=false;
                    for(auto row:support_[c])if(covered[color][row]) {
                        overlap=true;break;
                    }
                    if(!overlap)break;
                }
                if(color==groups.size()) {
                    groups.emplace_back();
                    covered.emplace_back(n_,0);
                }
                groups[color].push_back(block*n_+c);
                for(auto row:support_[c])covered[color][row]=1;
            }
        stats_.colors=groups.size();
        for(const auto& group:groups)
            stats_.largest_color=std::max(stats_.largest_color,group.size());
        center_=base_;
        auto central=original;
        if(use_common_anchor) {
            center_=generalized_anchor_candidate_v62(base_,states_,pin_,reference_,32.*step_);
            validate_physical_candidate_v51(center_,kmin_,wmin_);
            central=evaluate_physical_candidate_v51(center_,evaluate,kmin_,wmin_);
            ++stats_.physical_evaluations;
            const auto f0=linear_rows(base_,original);
            const auto fa=linear_rows(center_,central);
            stats_.anchor_displacement=distance_relative(f0,fa);
            long double shift2=0.;
            for(std::size_t i=0;i<size_;++i) {
                const long double d=(center_.fields[i]-base_.fields[i])/
                                    states_[i/n_];
                shift2+=d*d;
            }
            stats_.anchor_normalized_distance=std::sqrt(static_cast<double>(shift2));
            stats_.trusted &=stats_.anchor_displacement<=anchor_tolerance;
        }
        const auto fcenter=linear_rows(center_,central);
        cfdx::core::SparseMatrix matrix(size_,size_);
        for(const auto& group:groups) {
            std::vector<int> row_owner(n_,-1);
            for(auto col:group)
                for(auto cell:support_[col%n_]) {
                    if(row_owner[cell]!=-1)
                        throw std::logic_error("V63 color row collision");
                    row_owner[cell]=static_cast<int>(col);
                }
            auto plus=center_,minus=center_;
            for(auto col:group) {
                const double shift=step_*states_[col/n_];
                plus.fields[col]+=shift;
                minus.fields[col]-=shift;
            }
            // Never fall back to one-sided samples silently at SST bounds.
            validate_physical_candidate_v51(plus,kmin_,wmin_);
            validate_physical_candidate_v51(minus,kmin_,wmin_);
            const auto rp=evaluate_physical_candidate_v51(plus,evaluate,kmin_,wmin_);
            const auto rm=evaluate_physical_candidate_v51(minus,evaluate,kmin_,wmin_);
            stats_.physical_evaluations+=2;
            const auto fp=linear_rows(plus,rp),fm=linear_rows(minus,rm);
            long double outside=0.,total=0.;
            std::vector<long double> sided(group.size(),0.),sided_scale(group.size(),0.);
            // map from column index to group slot; only on cells used by a color
            for(std::size_t block=0;block<6;++block)
                for(std::size_t cell=0;cell<n_;++cell) {
                    const auto row=block*n_+cell;
                    const long double diff=static_cast<long double>(fp[row])-fm[row];
                    total+=diff*diff;
                    const int col=row_owner[cell];
                    if(col<0) {outside+=diff*diff;continue;}
                    const auto delta=(fp[row]-fm[row])/(2.*step_);
                    if(!std::isfinite(delta))
                        throw std::runtime_error("V63 nonfinite local physical derivative");
                    if(delta!=0.)matrix.push_back(row,static_cast<std::size_t>(col),delta);
                    const auto slot=static_cast<std::size_t>(std::distance(
                        group.begin(),std::find(group.begin(),group.end(),
                        static_cast<std::size_t>(col))));
                    const long double pos=(fp[row]-fcenter[row])/step_;
                    const long double neg=(fcenter[row]-fm[row])/step_;
                    sided[slot]+=(pos-neg)*(pos-neg);
                    sided_scale[slot]+=std::max(pos*pos,neg*neg);
                }
            const double lost=std::sqrt(static_cast<double>(outside)) /
                std::max(std::sqrt(static_cast<double>(total)),1e-12);
            stats_.maximum_outside_stencil=std::max(stats_.maximum_outside_stencil,lost);
            if(!(lost<=consistency_tolerance))stats_.trusted=false;
            for(std::size_t idx=0;idx<group.size();++idx) {
                const double disagreement=std::sqrt(static_cast<double>(sided[idx]))/
                    std::max(std::sqrt(static_cast<double>(sided_scale[idx])),1e-12);
                stats_.maximum_sided=std::max(stats_.maximum_sided,disagreement);
                if(!std::isfinite(disagreement) || disagreement>consistency_tolerance) {
                    ++stats_.rejected_columns;
                    stats_.trusted=false;
                }
            }
        }
        matrix.finalize();
        matrix_=std::move(matrix);
        stats_.nonzeros=matrix_.nnz();
    }

    std::vector<double> apply(const std::vector<double>& y) const {
        if(y.size()!=size_ || !std::all_of(y.begin(),y.end(),
                                           [](double v){return std::isfinite(v);}))
            throw std::invalid_argument("V63 nonfinite/invalid Krylov direction");
        return matrix_.matvec(y);
    }

    template<class Evaluate>
    SparseLocalDirectionalAuditV63 audit_direction(const std::vector<double>& v,
                                                    Evaluate&& evaluate) const {
        if(v.size()!=size_)
            throw std::invalid_argument("V63 physical audit direction dimension");
        long double norm2=0.;
        for(double x:v) {
            if(!std::isfinite(x))throw std::invalid_argument("V63 nonfinite audit direction");
            norm2+=static_cast<long double>(x)*x;
        }
        if(!(norm2>0.))throw std::invalid_argument("V63 zero audit direction");
        const double norm=std::sqrt(static_cast<double>(norm2));
        const auto av=apply(v);
        const auto base_eval=evaluate_physical_candidate_v51(base_,evaluate,kmin_,wmin_);
        const auto fbase=linear_rows(base_,base_eval);
        SparseLocalDirectionalAuditV63 result;
        ++result.evaluations;
        std::vector<double> previous;
        for(int level=0;level<2;++level) {
            const double h=std::ldexp(step_,-level)/norm;
            auto plus=base_,minus=base_;
            for(std::size_t i=0;i<size_;++i) {
                plus.fields[i]+=h*v[i]*states_[i/n_];
                minus.fields[i]-=h*v[i]*states_[i/n_];
            }
            validate_physical_candidate_v51(plus,kmin_,wmin_);
            validate_physical_candidate_v51(minus,kmin_,wmin_);
            const auto fp=linear_rows(plus,evaluate_physical_candidate_v51(
                                        plus,evaluate,kmin_,wmin_));
            const auto fm=linear_rows(minus,evaluate_physical_candidate_v51(
                                        minus,evaluate,kmin_,wmin_));
            result.evaluations+=2;
            std::vector<double> derivative(size_);
            long double difference=0.,size=0.,sided=0.;
            std::array<long double,6> blockdiff{},blocksize{},blocksided{};
            for(std::size_t i=0;i<size_;++i) {
                derivative[i]=(fp[i]-fm[i])/(2.*h);
                const long double d=static_cast<long double>(derivative[i])-av[i];
                const long double s=static_cast<long double>(fp[i])-2.*fbase[i]+fm[i];
                const long double slope=s/(h);
                const long double scale=std::max(std::abs(static_cast<long double>(av[i])),
                                                 std::abs(static_cast<long double>(derivative[i])));
                difference+=d*d;size+=scale*scale;sided+=slope*slope;
                const auto block=i/n_;
                blockdiff[block]+=d*d;blocksize[block]+=scale*scale;
                blocksided[block]+=slope*slope;
            }
            const double denominator=std::max(std::sqrt(static_cast<double>(size)),1e-12);
            result.fidelity=std::max(result.fidelity,
                       std::sqrt(static_cast<double>(difference))/denominator);
            result.sided=std::max(result.sided,
                       std::sqrt(static_cast<double>(sided))/denominator);
            for(std::size_t b=0;b<6;++b) {
                const double scale=std::max(std::sqrt(static_cast<double>(blocksize[b])),1e-12);
                result.block_fidelity[b]=std::max(result.block_fidelity[b],
                    std::sqrt(static_cast<double>(blockdiff[b]))/scale);
                result.block_sided[b]=std::max(result.block_sided[b],
                    std::sqrt(static_cast<double>(blocksided[b]))/scale);
            }
            if(!previous.empty()) {
                long double discrepancy=0.;
                for(std::size_t i=0;i<size_;++i) {
                    const long double d=static_cast<long double>(derivative[i])-previous[i];
                    discrepancy+=d*d;
                }
                result.halving=std::sqrt(static_cast<double>(discrepancy))/denominator;
            }
            previous=std::move(derivative);
        }
        return result;
    }
    const SparseLocalJacobianStatsV63& statistics() const noexcept {return stats_;}
    const cfdx::core::SparseMatrix& matrix() const noexcept {return matrix_;}

private:
    std::vector<double> linear_rows(const PhysicalCandidateV51& state,
                                    const PhysicalResidualV51& f) const {
        std::vector<double> rows(size_);
        for(std::size_t i=0;i<size_;++i)
            rows[i]=f.signed_raw[i]/rowscale_[i/n_];
        if(pin_)
            rows[3*n_+reference_]=(state.fields[3*n_+reference_]-
                                    reference_pressure_)/states_[3];
        return rows;
    }
    static double distance_relative(const std::vector<double>& a,
                                    const std::vector<double>& b) {
        long double d=0.,scale=0.;
        for(std::size_t i=0;i<a.size();++i) {
            const long double diff=static_cast<long double>(a[i])-b[i];
            d+=diff*diff;
            scale+=static_cast<long double>(a[i])*a[i];
        }
        return std::sqrt(static_cast<double>(d))/
            std::max(std::sqrt(static_cast<double>(scale)),1e-12);
    }
    PhysicalCandidateV51 base_,center_;
    Scale states_{},rowscale_{};
    double step_=0.,reference_pressure_=0.,kmin_=0.,wmin_=0.;
    bool pin_=false;
    std::size_t reference_=0,n_=0,size_=0;
    std::vector<std::vector<std::size_t>> support_;
    cfdx::core::SparseMatrix matrix_;
    SparseLocalJacobianStatsV63 stats_;
};
} // namespace cfdx::physics