#pragma once
// CFDX V54: bounded, opt-in, explicit frozen finite-difference Jacobian.
// The expensive dense operator is a FORENSIC reference for small problems:
// it is NOT a substitute for an efficient industrial matrix-free Jv.
// Amat is assembled from the candidate-only Fphys; Pmat remains independent.
// The linear action used by ALL Arnoldi products and by the true-residual
// check is exactly the same frozen map. Finite difference fidelity to Fphys
// is audited independently and can refuse the Newton direction.
#include "cfdx/physics/physical_six_candidate_v51.h"
#include "cfdx/physics/convective_face_tangent_v55.h"
#include "cfdx/physics/coupled_newton_krylov.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
#include <utility>

namespace cfdx::physics {

struct FrozenPhysicalJacobianStatsV54 {
    std::size_t rows=0;
    std::size_t evaluations=0;
    std::size_t centered_columns=0;
    std::size_t one_sided_columns=0;
    std::size_t stored_bytes=0;
    double maximum_sample_fidelity=0.0;
    double maximum_sample_halving=0.0;
    double maximum_sample_sided_nonsmoothness=0.0;
    std::array<double,6> maximum_block_fidelity{};
    std::array<double,6> maximum_block_sided{};
    std::size_t maximum_flux_sign_flips=0;
    std::size_t maximum_near_zero_flux_crossings=0;
    double maximum_flux_sided_asymmetry=0.;
};

struct FrozenDerivativeAuditV54 {
    double fidelity=0.0;
    double halving=0.0;
    double sided_nonsmoothness=0.0;
    std::array<double,6> block_fidelity{};
    std::array<double,6> block_sided{};
    std::size_t flux_sign_flips=0;
    std::size_t near_zero_flux_crossings=0;
    std::vector<FaceBranchEventV55> face_branch_events;
    double flux_sided_asymmetry=0.;
};

class FrozenPhysicalJacobianV54 {
public:
    using Scales=std::array<double,6>;
    // Include only a bounded small mesh; memory usage is O((6N)^2).
    static constexpr std::size_t max_unknowns=1024;

    template<class Evaluate>
    FrozenPhysicalJacobianV54(const PhysicalCandidateV51& base,
                              const PhysicalResidualV51& physical_base,
                              const Scales& state_scale,
                              const Scales& row_scale,
                              bool pin_pressure,std::size_t pin_cell,
                              double reference_pressure,double relative_step,
                              Evaluate&& evaluate,double k_min=0.,
                              double omega_min=0.)
      : n_(base.cells),size_(base.fields.size()),x_(base),
        state_scale_(state_scale),row_scale_(row_scale),
        pin_pressure_(pin_pressure),pin_cell_(pin_cell),
        pressure_reference_(reference_pressure),
        k_min_(k_min),omega_min_(omega_min) {
        validate_physical_candidate_v51(base,k_min,omega_min);
        validate_physical_result_v51(base,physical_base);
        if(size_>max_unknowns || pin_cell_>=n_ ||
           !std::isfinite(relative_step) || relative_step<1e-9 ||
           relative_step>1e-3 || !std::isfinite(pressure_reference_))
            throw std::invalid_argument("V54 invalid bounded frozen Jacobian settings");
        for(std::size_t b=0;b<6;++b)
            if(!std::isfinite(state_scale[b]) || !(state_scale[b]>0.) ||
               !std::isfinite(row_scale[b]) || !(row_scale[b]>0.))
                throw std::invalid_argument("V54 invalid frozen scales");
        for(std::size_t i=0;i<size_;++i)
            y_.push_back(x_.fields[i]/state_scale_[i/n_]);
        f0_=rows(x_,physical_base);
        flux0_=physical_base.face_flux;
        // Walker-Pernice-inspired common perturbation magnitude for UNIT
        // directions in normalized state coordinates. Explicit columns
        // are individually differentiated once, then frozen for Krylov.
        const double radius=relative_step*std::sqrt(1.0+coupled_newton_norm2(y_));
        if(!(radius>0.) || !std::isfinite(radius))
            throw std::runtime_error("V54 nonfinite WP differencing radius");
        step_=radius;
        coefficients_.assign(size_*size_,0.0);
        stats_.rows=size_;
        stats_.stored_bytes=coefficients_.size()*sizeof(double);
        std::vector<double> shifted=y_;
        for(std::size_t column=0;column<size_;++column) {
            auto probe=[&](double increment, bool& valid) {
                shifted[column]=y_[column]+increment;
                PhysicalCandidateV51 candidate=x_;
                for(std::size_t i=0;i<size_;++i)
                    candidate.fields[i]=shifted[i]*state_scale_[i/n_];
                shifted[column]=y_[column];
                try {
                    validate_physical_candidate_v51(candidate,k_min_,omega_min_);
                }catch(const std::invalid_argument&) {
                    valid=false;
                    return std::vector<double>{};
                }
                const auto physical=evaluate_physical_candidate_v51(
                    candidate,evaluate,k_min_,omega_min_);
                ++stats_.evaluations;
                valid=true;
                return rows(candidate,physical);
            };
            bool positive=false,negative=false;
            const auto plus=probe(step_,positive);
            const auto minus=probe(-step_,negative);
            if(!positive && !negative)
                throw std::runtime_error("V54 neither finite difference probe admissible");
            if(positive && negative)++stats_.centered_columns;
            else ++stats_.one_sided_columns;
            for(std::size_t row=0;row<size_;++row) {
                const double entry=positive && negative
                    ?(plus[row]-minus[row])/(2.0*step_)
                    :positive?(plus[row]-f0_[row])/step_
                             :(f0_[row]-minus[row])/step_;
                if(!std::isfinite(entry))
                    throw std::runtime_error("V54 nonfinite explicit Jacobian entry");
                coefficients_[row*size_+column]=entry;
            }
        }
    }

    std::vector<double> apply(const std::vector<double>& v) const {
        if(v.size()!=size_)
            throw std::invalid_argument("V54 wrong frozen Jacobian action dimension");
        std::vector<double> result(size_,0.0);
        for(std::size_t row=0;row<size_;++row) {
            // Extended accumulator lowers false residuals from the linear
            // algebra, independently of the physical residual evaluation.
            long double sum=0.;
            for(std::size_t col=0;col<size_;++col) {
                if(!std::isfinite(v[col]))
                    throw std::invalid_argument("V54 nonfinite direction");
                sum+=static_cast<long double>(coefficients_[row*size_+col])*v[col];
            }
            result[row]=static_cast<double>(sum);
            if(!std::isfinite(result[row]))
                throw std::runtime_error("V54 nonfinite frozen matrix product");
        }
        return result;
    }

    // Independent Fphys check on an actual Krylov direction or final step.
    // Both h and h/2 are checked. The same frozen A is NOT used to form the
    // finite difference comparator. This can detect noisy or nonsmooth F.
    template<class Evaluate>
    FrozenDerivativeAuditV54 audit_direction(const std::vector<double>& direction,
                                               Evaluate&& evaluate) {
        const double norm=coupled_newton_norm2(direction);
        if(direction.size()!=size_ || !(norm>0.) || !std::isfinite(norm))
            throw std::invalid_argument("V54 invalid audit direction");
        const auto av=apply(direction);
        double sided=0.0,flux_asymmetry=0.0;
        std::size_t flips=0,near_flips=0;
        std::vector<FaceBranchEventV55> face_events;
        std::array<double,6> sided_by_block{};
        const auto direct=[&](double h) {
            std::vector<double> flux_plus,flux_minus;
            auto probe=[&](double sign,bool& valid) {
                auto trial=x_;
                for(std::size_t i=0;i<size_;++i)
                    trial.fields[i]=(y_[i]+sign*h*direction[i]/norm)*
                                     state_scale_[i/n_];
                try {validate_physical_candidate_v51(trial,k_min_,omega_min_);}
                catch(const std::invalid_argument&) {
                    valid=false;return std::vector<double>{};
                }
                const auto r=evaluate_physical_candidate_v51(
                    trial,evaluate,k_min_,omega_min_);
                ++stats_.evaluations;valid=true;
                if(sign>0.0)flux_plus=r.face_flux;
                else flux_minus=r.face_flux;
                return rows(trial,r);
            };
            bool plus_ok=false,minus_ok=false;
            const auto plus=probe(+1.,plus_ok);
            const auto minus=probe(-1.,minus_ok);
            if(!plus_ok && !minus_ok)
                throw std::runtime_error("V54 no valid physical audit probes");
            std::vector<double> output(size_);
            const double multiplier=norm/h;
            if(plus_ok && minus_ok) {
                if(flux_plus.size()!=flux0_.size() ||
                   flux_minus.size()!=flux0_.size())
                    throw std::runtime_error("V54 changed face ordering in physical probes");
                long double numerator=0.,denominator=0.;
                std::size_t count=0;
                for(std::size_t face=0;face<flux0_.size();++face) {
                    if(flux_plus[face]*flux_minus[face]<0.)++count;
                    const long double asym=static_cast<long double>(flux_plus[face])+
                        flux_minus[face]-2.*flux0_[face];
                    const long double slope=static_cast<long double>(flux_plus[face])-
                        flux_minus[face];
                    numerator+=asym*asym;denominator+=slope*slope;
                }
                flips=std::max(flips,count);
                const auto branches=audit_face_branches_v55(
                    flux0_,flux_plus,flux_minus);
                near_flips=std::max(near_flips,branches.near_zero_crossings);
                // Record deterministic face IDs from the tighter h/2 probe;
                // these are diagnostics only, not derivatives of limiter/RC.
                face_events=std::move(branches.events);
                flux_asymmetry=std::max(flux_asymmetry,
                    std::sqrt(static_cast<double>(numerator))/
                    std::max(std::sqrt(static_cast<double>(denominator)),1e-12));
                std::vector<double> one_side_diff(size_);
                for(std::size_t i=0;i<size_;++i)
                    one_side_diff[i]=(plus[i]-2.*f0_[i]+minus[i])*multiplier;
                const double reference=std::max({coupled_newton_norm2(av),
                    coupled_newton_norm2(plus),coupled_newton_norm2(minus),1e-12});
                sided=std::max(sided,coupled_newton_norm2(one_side_diff)/reference);
                for(std::size_t block=0;block<6;++block) {
                    long double subtotal=0.,blockdenom=0.;
                    for(std::size_t i=block*n_;i<(block+1)*n_;++i) {
                        subtotal+=static_cast<long double>(one_side_diff[i])*one_side_diff[i];
                        blockdenom+=static_cast<long double>(av[i])*av[i];
                    }
                    sided_by_block[block]=std::max(sided_by_block[block],
                        std::sqrt(static_cast<double>(subtotal)) /
                        std::max(std::sqrt(static_cast<double>(blockdenom)),1e-12));
                }
            }
            for(std::size_t i=0;i<size_;++i) {
                output[i]=(plus_ok && minus_ok
                    ?(plus[i]-minus[i])*0.5
                    :plus_ok?(plus[i]-f0_[i]):(f0_[i]-minus[i]))*multiplier;
            }
            return output;
        };
        const auto f1=direct(step_),f2=direct(0.5*step_);
        std::vector<double> fidelity(size_),halving(size_);
        for(std::size_t i=0;i<size_;++i) {
            fidelity[i]=av[i]-f2[i];
            halving[i]=f1[i]-f2[i];
        }
        const double denominator=std::max({coupled_newton_norm2(av),
            coupled_newton_norm2(f1),coupled_newton_norm2(f2),1e-12});
        const double err=coupled_newton_norm2(fidelity)/denominator;
        const double half_err=coupled_newton_norm2(halving)/denominator;
        FrozenDerivativeAuditV54 details;
        details.fidelity=err;details.halving=half_err;
        details.sided_nonsmoothness=sided;
        details.flux_sign_flips=flips;
        details.near_zero_flux_crossings=near_flips;
        details.face_branch_events=std::move(face_events);
        details.flux_sided_asymmetry=flux_asymmetry;
        for(std::size_t block=0;block<6;++block) {
            long double num=0.,denom=0.;
            for(std::size_t i=block*n_;i<(block+1)*n_;++i) {
                num+=static_cast<long double>(fidelity[i])*fidelity[i];
                denom+=static_cast<long double>(av[i])*av[i];
            }
            details.block_fidelity[block]=std::sqrt(static_cast<double>(num)) /
                    std::max(std::sqrt(static_cast<double>(denom)),1e-12);
            details.block_sided[block]=sided_by_block[block];
            stats_.maximum_block_fidelity[block]=std::max(
                stats_.maximum_block_fidelity[block],details.block_fidelity[block]);
            stats_.maximum_block_sided[block]=std::max(
                stats_.maximum_block_sided[block],details.block_sided[block]);
        }
        stats_.maximum_flux_sign_flips=std::max(stats_.maximum_flux_sign_flips,flips);
        stats_.maximum_near_zero_flux_crossings=std::max(
            stats_.maximum_near_zero_flux_crossings,near_flips);
        stats_.maximum_flux_sided_asymmetry=std::max(
            stats_.maximum_flux_sided_asymmetry,flux_asymmetry);
        stats_.maximum_sample_fidelity=std::max(stats_.maximum_sample_fidelity,err);
        stats_.maximum_sample_halving=std::max(stats_.maximum_sample_halving,half_err);
        stats_.maximum_sample_sided_nonsmoothness=std::max(
            stats_.maximum_sample_sided_nonsmoothness,sided);
        return details;
    }

    const FrozenPhysicalJacobianStatsV54& statistics() const noexcept {return stats_;}
    std::size_t size() const noexcept {return size_;}
private:
    std::vector<double> rows(const PhysicalCandidateV51& candidate,
                             const PhysicalResidualV51& r) const {
        std::vector<double> output(size_);
        for(std::size_t i=0;i<size_;++i)
            output[i]=r.signed_raw[i]/row_scale_[i/n_];
        if(pin_pressure_)
            output[3*n_+pin_cell_]=
                (candidate.fields[3*n_+pin_cell_]-pressure_reference_)/
                state_scale_[3];
        return output;
    }
    std::size_t n_=0,size_=0;
    PhysicalCandidateV51 x_;
    Scales state_scale_{},row_scale_{};
    bool pin_pressure_=false;
    std::size_t pin_cell_=0;
    double pressure_reference_=0.,step_=0.;
    double k_min_=0.,omega_min_=0.;
    std::vector<double> y_,f0_,flux0_,coefficients_;
    FrozenPhysicalJacobianStatsV54 stats_;
};

} // namespace cfdx::physics
