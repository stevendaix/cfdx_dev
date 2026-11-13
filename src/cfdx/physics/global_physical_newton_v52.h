#pragma once
// CFDX V52: opt-in, PRIVATE-state six-field nonlinear Newton–Krylov trial.
// The only nonlinear callback is V51's candidate-owned finite-volume F(X).
// This diagnostic does not publish fields into the production solve. Each
// Armijo decision is made from a newly reassembled V51 physical residual.
// The gauge equation replaces one redundant continuity row ONLY in the
// Krylov algebra; the full continuity block remains in the physical merit.
#include "cfdx/physics/coupled_newton_krylov.h"
#include "cfdx/physics/physical_six_candidate_v51.h"
#include "cfdx/physics/frozen_physical_jacobian_v54.h"
#include "cfdx/physics/hybrid_physical_jv_v57.h"
#include "cfdx/physics/physical_remainder_audit_v60.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct GlobalPhysicalNewtonControlsV52 {
    std::size_t maximum_steps=1;
    std::size_t krylov_dimension=6;
    // Larger FGMRES window is permitted only with the opt-in V53 Pmat.
    bool v53_extended_krylov=false;
    bool v53_central_difference=false;
    // V54: bounded explicit 6N Amat, separate from the V53 Schur Pmat.
    // Never used in ordinary production (2 evaluations per Jacobian column).
    bool v54_frozen_jacobian=false;
    // V57 forensic only: verify TRUE Krylov residual and action additivity
    // at every Arnoldi update, rather than trusting the final Hessenberg norm.
    bool v57_true_krylov_each_step=false;
    bool v57_hybrid_convective_jv=false;
    // V58: a separately audited generalized derivative, never the physical F.
    bool v58_symmetric_sou_ties=false;
    // V59: fixed SOU selection at the nonlinear base, not physical certification.
    bool v59_frozen_sou_active_set=false;
    // Read-only audit on two REAL preconditioned Krylov vectors, OFF by default.
    bool v60_remainder_audit=false;
    double v57_convective_fidelity_tolerance=1e-3;
    double v57_krylov_consistency_tolerance=1e-6; // algebraic, not CFD tolerance

    double v54_fidelity_tolerance=0.05;
    std::size_t maximum_backtracks=9;
    double difference_step=1e-6;
    double linear_relative_tolerance=0.3;
    double sufficient_decrease=1e-4;
    // Only a pure-Neumann pressure system gets a replaced algebraic row.
    bool pin_pressure=true;
    std::size_t reference_cell=0;
    double reference_pressure=0.0;
    double k_min=0.0;
    double omega_min=0.0;
    // Experimental two-probe row-sum approximation of the local SST 2x2.
    // NOT the exact local Jacobian; neighbor couplings remain in matrix-free Jv.
    bool sst_two_by_two_preconditioner=false;
};

struct KrylovTruthSnapshotV57 {
    std::size_t direction=0;
    double arnoldi_relative=0.;
    double true_relative=0.;
    double action_nonadditivity_relative=0.;
    double arnoldi_gap=0.;
};

struct GlobalPhysicalNewtonStepV52 {
    bool accepted=false;
    double merit_before=std::numeric_limits<double>::infinity();
    double merit_after=std::numeric_limits<double>::infinity();
    double linear_true_relative=std::numeric_limits<double>::infinity();
    double linear_arnoldi_relative=std::numeric_limits<double>::infinity();
    double jv_linearity_relative=std::numeric_limits<double>::infinity();
    double alpha=0.0;
    std::size_t candidate_evaluations=0;
    std::size_t rejected_trials=0;
    std::string reason;
    std::vector<KrylovTruthSnapshotV57> v57_krylov_trace;
    bool v57_krylov_consistent=true;
    double v57_maximum_arnoldi_gap=0.;
    double v57_maximum_action_nonadditivity=0.;
    double v57_maximum_convective_discrepancy=0.;
    std::size_t v57_nonsmooth_faces=0;
    std::size_t v57_hybrid_jv_applications=0;
    bool v57_convective_fidelity_ok=true;
    bool v58_generalized_sou=false;
    PhysicalRemainderAdditivityV60 v60_remainder;
    bool v60_remainder_consistent=false;
    bool v60_full_additive=false;
    bool sst_two_by_two_ready=false;
    bool v53_global_block_ready=false;
    std::size_t v53_probe_evaluations=0;
    bool v54_frozen_jacobian_ready=false;
    bool v54_krylov_fidelity_ok=true;
    std::size_t v54_jacobian_evaluations=0;
    std::size_t v54_basis_audits=0;
    double v54_maximum_jv_fidelity=0.0;
    double v54_maximum_halving_discrepancy=0.0;
    double v54_maximum_sided_nonsmoothness=0.0;
    std::array<double,6> v54_block_fidelity{};
    std::array<double,6> v54_block_sided{};
    std::size_t v54_centered_columns=0;
    std::size_t v54_one_sided_columns=0;
    std::size_t v54_flux_sign_flips=0;
    double v54_flux_asymmetry=0.0;
    std::size_t v53_transport_nnz=0;
    double v53_maximum_locality_error=0.0;
};

struct GlobalPhysicalNewtonResultV52 {
    // A private sequence: never automatically publish to U,p,phi or SST.
    PhysicalCandidateV51 last_private_state;
    PhysicalResidualV51 last_private_residual;
    std::vector<GlobalPhysicalNewtonStepV52> steps;
    std::size_t accepted_steps=0;
    bool production_committed=false;
    bool physically_qualified=false;
};

inline void validate_global_physical_newton_v52(
    const GlobalPhysicalNewtonControlsV52& c, const PhysicalCandidateV51& x) {
    validate_physical_candidate_v51(x,c.k_min,c.omega_min);
    if (!c.maximum_steps || c.maximum_steps>8 || !c.krylov_dimension ||
        c.krylov_dimension>(c.v53_extended_krylov?96u:32u) ||
        c.maximum_backtracks>20 ||
        !std::isfinite(c.difference_step) || c.difference_step<1e-9 ||
        c.difference_step>1e-3 ||
        !std::isfinite(c.v54_fidelity_tolerance) ||
        c.v54_fidelity_tolerance<=0.0 || c.v54_fidelity_tolerance>0.5 ||
        (c.v54_frozen_jacobian && x.fields.size()>FrozenPhysicalJacobianV54::max_unknowns) ||
        !std::isfinite(c.linear_relative_tolerance) ||
        !(c.linear_relative_tolerance>0.0 && c.linear_relative_tolerance<1.0) ||
        !std::isfinite(c.sufficient_decrease) ||
        !(c.sufficient_decrease>0.0 && c.sufficient_decrease<0.5) ||
        !std::isfinite(c.reference_pressure) ||
        c.reference_cell>=x.cells ||
        !std::isfinite(c.v57_krylov_consistency_tolerance) ||
        c.v57_krylov_consistency_tolerance<=0. ||
        c.v57_krylov_consistency_tolerance>0.1 ||
        !std::isfinite(c.v57_convective_fidelity_tolerance) ||
        c.v57_convective_fidelity_tolerance<=0. ||
        c.v57_convective_fidelity_tolerance>0.05 ||
        (c.v57_hybrid_convective_jv && !c.v57_true_krylov_each_step) ||
        (c.v57_hybrid_convective_jv && c.v54_frozen_jacobian) ||
        (c.v58_symmetric_sou_ties &&
         (!c.v57_hybrid_convective_jv || !c.v57_true_krylov_each_step)) ||
        (c.v59_frozen_sou_active_set &&
         (!c.v57_hybrid_convective_jv || !c.v57_true_krylov_each_step ||
          c.v58_symmetric_sou_ties)) ||
        (c.v60_remainder_audit &&
         (!c.v57_hybrid_convective_jv || !c.v57_true_krylov_each_step)))
        throw std::invalid_argument("V52 invalid JFNK controls or pressure gauge");
}

// The preconditioner is expressed in nondimensional state coordinates and
// must be a FIXED LINEAR map throughout one GMRES solve. If inner pressure
// solves vary, replace GMRES by FGMRES before enabling a production path.
using GlobalPhysicalRightPreconditionerV52 = CoupledRightPreconditioner;
// Gives a production 4N Schur action access to the frozen row/state scales.
// The first 4N residual coordinates are mapped back to physical units before
// the Schur action, and its physical correction is mapped back to y units.
// Two SST blocks are included in the same 6N vector (identity in V52).
using GlobalPhysicalNativeSchurV52 = std::function<std::vector<double>(
    const std::vector<double>&,const std::array<double,6>&,
    const std::array<double,6>&,std::size_t)>;

// Optional Pmat builder invoked ONCE per accepted Newton iterate, with its
// current candidate and both fixed variable/row scales. This is independent
// from the matrix-free Jv. A failed build is a failed private experiment.
struct GlobalPhysicalFrozenPmatV53 {
    GlobalPhysicalRightPreconditionerV52 apply;
    std::size_t probes=0;
    std::size_t transport_nnz=0;
    double locality_error=0.0;
};
using GlobalPhysicalPmatBuilderV53 = std::function<GlobalPhysicalFrozenPmatV53(
    const PhysicalCandidateV51&, const PhysicalResidualV51&,
    const std::array<double,6>&, const std::array<double,6>&)>;

// Flexible right-preconditioned Arnoldi: z_j=P_j^{-1}v_j is retained
// individually. This is essential for PCD's *inexact* nested solves. The
// final check always applies the original matrix-free Jv to sum(y_j*z_j).
template<class Apply,class Recover>
bool global_physical_fgmres_v52(Apply&& jv, Recover&& recover,
                               const std::vector<double>& rhs,
                               std::size_t dimension,double relative_tolerance,
                               std::vector<double>& step,double& true_relative,
                               double* arnoldi_relative=nullptr,
                               double* jv_linearity_relative=nullptr,
                               const std::function<bool(const std::vector<double>&)>&
                                   v54_basis_fidelity={},
                               const std::function<void(const KrylovTruthSnapshotV57&)>&
                                   v57_truth_observer={},
                               const std::function<void(const std::vector<double>&)>&
                                   v60_basis_observer={}) {
    const std::size_t n=rhs.size();
    const double beta=coupled_newton_norm2(rhs);
    step.assign(n,0.0);
    if(beta==0.0) {true_relative=0.0;return true;}
    const std::size_t m=std::min(n,dimension);
    std::vector<std::vector<double>> q,z,images;
    q.reserve(m+1);z.reserve(m);images.reserve(m);
    q.push_back(rhs);
    for(double& v:q.back())v/=beta;
    std::vector<std::vector<double>> h(m+1,std::vector<double>(m,0.0));
    std::vector<double> cs(m),sn(m),g(m+1,0.0);
    g[0]=beta;
    std::size_t used=0;
    for(std::size_t j=0;j<m;++j) {
        auto correction=recover(q[j]);
        if(correction.size()!=n)
            throw std::runtime_error("V52 incomplete flexible Schur action");
        (void)coupled_newton_norm2(correction);
        z.push_back(std::move(correction));
        if(v60_basis_observer) v60_basis_observer(z.back());
        if(v54_basis_fidelity && !v54_basis_fidelity(z.back())) {
            true_relative=std::numeric_limits<double>::infinity();
            return false;
        }
        auto w=jv(z.back());
        if(w.size()!=n) throw std::runtime_error("V52 incomplete directional Jv");
        images.push_back(w);
        for(int pass=0;pass<2;++pass)
            for(std::size_t i=0;i<=j;++i) {
                double dot=0.0;
                for(std::size_t k=0;k<n;++k)dot+=w[k]*q[i][k];
                h[i][j]+=dot;
                for(std::size_t k=0;k<n;++k)w[k]-=dot*q[i][k];
            }
        h[j+1][j]=coupled_newton_norm2(w);
        const bool breakdown=h[j+1][j]<1e-13*std::max(1.,beta);
        if(!breakdown) {
            for(double& value:w)value/=h[j+1][j];
            q.push_back(std::move(w));
        }
        for(std::size_t i=0;i<j;++i) {
            const double a=h[i][j],b=h[i+1][j];
            h[i][j]=cs[i]*a+sn[i]*b;
            h[i+1][j]=-sn[i]*a+cs[i]*b;
        }
        const double length=std::hypot(h[j][j],h[j+1][j]);
        if(length<=1e-18)break;
        cs[j]=h[j][j]/length;sn[j]=h[j+1][j]/length;
        h[j][j]=length;h[j+1][j]=0.0;
        g[j+1]=-sn[j]*g[j];g[j]=cs[j]*g[j];
        used=j+1;
        if(v57_truth_observer){
            // Use the *same* preconditioned z_i and frozen Amat action as
            // Krylov. Reconstruct the partial least-squares step after
            // every Arnoldi update; one independent Amat product per step.
            std::vector<double> y(used,0.);
            for(std::size_t back=0;back<used;++back){
                const auto i=used-1-back;
                double value=g[i];
                for(std::size_t col=i+1;col<used;++col)
                    value-=h[i][col]*y[col];
                if(std::abs(h[i][i])<1e-18)
                    throw std::runtime_error("V57 singular Arnoldi triangular factor");
                y[i]=value/h[i][i];
            }
            std::vector<double> partial(n,0.),predicted(n,0.);
            for(std::size_t col=0;col<used;++col)
                for(std::size_t i=0;i<n;++i){
                    partial[i]+=y[col]*z[col][i];
                    predicted[i]+=y[col]*images[col][i];
                }
            const auto actual=jv(partial);
            if(actual.size()!=n)throw std::runtime_error("V57 invalid partial true action");
            std::vector<double> error(n),nonadd(n);
            for(std::size_t i=0;i<n;++i){
                error[i]=rhs[i]-actual[i];
                nonadd[i]=actual[i]-predicted[i];
            }
            KrylovTruthSnapshotV57 snap;
            snap.direction=used;
            snap.arnoldi_relative=std::abs(g[used])/beta;
            snap.true_relative=coupled_newton_norm2(error)/beta;
            snap.action_nonadditivity_relative=coupled_newton_norm2(nonadd)/beta;
            snap.arnoldi_gap=std::abs(snap.true_relative-snap.arnoldi_relative);
            v57_truth_observer(snap);
        }
        if(std::abs(g[j+1])<=relative_tolerance*beta || breakdown)break;
    }
    if(!used) {true_relative=std::numeric_limits<double>::infinity();return false;}
    if(arnoldi_relative)*arnoldi_relative=std::abs(g[used])/beta;
    std::vector<double> coefficients(used,0.0);
    for(std::size_t back=0;back<used;++back) {
        const std::size_t i=used-1-back;
        double total=g[i];
        for(std::size_t j=i+1;j<used;++j)total-=h[i][j]*coefficients[j];
        if(std::abs(h[i][i])<=1e-18)
            return false;
        coefficients[i]=total/h[i][i];
    }
    for(std::size_t j=0;j<used;++j)
        for(std::size_t i=0;i<n;++i)
            step[i]+=coefficients[j]*z[j][i];
    const auto actual=jv(step);
    if(actual.size()!=n)throw std::runtime_error("V52 invalid verified Jv");
    if(jv_linearity_relative){
        std::vector<double> predicted(n,0.0);
        for(std::size_t j=0;j<used;++j)
            for(std::size_t i=0;i<n;++i)
                predicted[i]+=coefficients[j]*images[j][i];
        for(std::size_t i=0;i<n;++i)predicted[i]-=actual[i];
        *jv_linearity_relative=coupled_newton_norm2(predicted)/beta;
    }
    std::vector<double> error(n);
    for(std::size_t i=0;i<n;++i)error[i]=rhs[i]-actual[i];
    true_relative=coupled_newton_norm2(error)/beta;
    return true_relative<=relative_tolerance;
}

template<class Evaluate>
GlobalPhysicalNewtonResultV52 global_physical_newton_v52(
    const PhysicalCandidateV51& initial, Evaluate&& candidate_operator,
    const GlobalPhysicalNewtonControlsV52& controls,
    const GlobalPhysicalRightPreconditionerV52& right_preconditioner={},
    const GlobalPhysicalNativeSchurV52& native_schur={},
    const GlobalPhysicalPmatBuilderV53& v53_builder={},
    const GlobalPhysicalConvectiveTangentV57& v57_convective={}) {
    validate_global_physical_newton_v52(controls,initial);
    if((controls.v53_extended_krylov ||
        (controls.v53_central_difference && !controls.v57_hybrid_convective_jv)) &&
       !v53_builder)
        throw std::invalid_argument("V53 Krylov/central mode requires an explicit frozen Pmat builder");
    if (controls.v57_hybrid_convective_jv && !v57_convective)
        throw std::invalid_argument("V57 hybrid requires the native V56 convective tangent");
    const std::size_t nc=initial.cells, size=initial.fields.size();
    auto evaluate=[&](const PhysicalCandidateV51& x) {
        return evaluate_physical_candidate_v51(
            x,candidate_operator,controls.k_min,controls.omega_min);
    };
    GlobalPhysicalNewtonResultV52 answer;
    answer.last_private_state=initial;
    answer.last_private_residual=evaluate(initial);
    // Freeze state scales for this Newton *experiment*. These are variable
    // coordinates only, never physical residual scales or convergence gates.
    std::array<double,6> state_scale{};
    for (std::size_t b=0;b<6;++b) {
        double maximum=0.;
        for (std::size_t i=0;i<nc;++i)
            maximum=std::max(maximum,std::abs(initial.fields[b*nc+i]));
        state_scale[b]=std::max(maximum,1.0);
    }
    const auto to_state=[&](const std::vector<double>& y) {
        PhysicalCandidateV51 candidate{nc,y};
        if (y.size()!=size)
            throw std::invalid_argument("V52 incomplete six-field iterate");
        for (std::size_t i=0;i<size;++i)
            candidate.fields[i]*=state_scale[i/nc];
        return candidate;
    };
    auto to_y=[&](const PhysicalCandidateV51& candidate) {
        std::vector<double> y(size);
        for (std::size_t i=0;i<size;++i)
            y[i]=candidate.fields[i]/state_scale[i/nc];
        return y;
    };
    auto is_admissible=[&](const PhysicalCandidateV51& candidate) {
        try {
            validate_physical_candidate_v51(candidate,controls.k_min,
                                            controls.omega_min);
            return true;
        } catch(const std::invalid_argument&) {return false;}
    };
    for (std::size_t outer=0;outer<controls.maximum_steps;++outer) {
        GlobalPhysicalNewtonStepV52 record;
        const auto& current=answer.last_private_state;
        const auto& physical_base=answer.last_private_residual;
        record.merit_before=physical_base.merit;
        // Matrix-free Jacobian of RAW FV rows; fixed row scalings serve only
        // linear algebra. Never difference V51's normalized diagnostic rows.
        std::array<double,6> row_scale{};
        for (std::size_t b=0;b<6;++b) {
            double maximum=0.;
            for (std::size_t i=0;i<nc;++i)
                maximum=std::max(maximum,
                    std::abs(physical_base.signed_raw[b*nc+i]));
            row_scale[b]=std::max(maximum,1e-12);
        }
        GlobalPhysicalFrozenPmatV53 frozen_v53;
        if(v53_builder) {
            if(controls.sst_two_by_two_preconditioner)
                throw std::invalid_argument("V53 block and V52 local SST probes are exclusive");
            frozen_v53=v53_builder(current,physical_base,state_scale,row_scale);
            if(!frozen_v53.apply)
                throw std::runtime_error("V53 preconditioner builder returned no operator");
            record.v53_global_block_ready=true;
            record.v53_probe_evaluations=frozen_v53.probes;
            record.v53_transport_nnz=frozen_v53.transport_nnz;
            record.v53_maximum_locality_error=frozen_v53.locality_error;
        }
        std::vector<std::array<double,4>> sst_inverse;
        if(controls.sst_two_by_two_preconditioner) {
            std::array<PhysicalResidualV51,2> perturbed;
            const std::array<std::size_t,2> blocks{{4,5}};
            std::array<double,2> h{{0.,0.}};
            for(std::size_t slot=0;slot<2;++slot) {
                const std::size_t b=blocks[slot];
                h[slot]=controls.difference_step*state_scale[b];
                auto probe=current;
                for(std::size_t i=0;i<nc;++i)
                    probe.fields[b*nc+i]+=h[slot];
                perturbed[slot]=evaluate(probe);
                ++record.candidate_evaluations;
            }
            sst_inverse.resize(nc);
            bool usable=true;
            for(std::size_t i=0;i<nc;++i) {
                const double a=(perturbed[0].signed_raw[4*nc+i]-
                                physical_base.signed_raw[4*nc+i])/h[0];
                const double c=(perturbed[0].signed_raw[5*nc+i]-
                                physical_base.signed_raw[5*nc+i])/h[0];
                const double b=(perturbed[1].signed_raw[4*nc+i]-
                                physical_base.signed_raw[4*nc+i])/h[1];
                const double d=(perturbed[1].signed_raw[5*nc+i]-
                                physical_base.signed_raw[5*nc+i])/h[1];
                const double det=a*d-b*c;
                if(!std::isfinite(det) ||
                   std::abs(det)<=1e-10*std::max(1e-30,std::abs(a*d)+std::abs(b*c))) {
                    usable=false;break;
                }
                sst_inverse[i]={d/det,-b/det,-c/det,a/det};
                for(double coefficient:sst_inverse[i])
                    if(!std::isfinite(coefficient) || std::abs(coefficient)>1e12)
                        usable=false;
                if(!usable)break;
            }
            if(!usable) {
                record.reason="SST 2x2 probe singular or nonfinite: no silent fallback";
                answer.steps.push_back(std::move(record));
                break;
            }
            record.sst_two_by_two_ready=true;
        }
        const auto make_linear_rows=[&](const PhysicalCandidateV51& state,
                                       const PhysicalResidualV51& physical) {
            std::vector<double> rows(size);
            for (std::size_t i=0;i<size;++i)
                rows[i]=physical.signed_raw[i]/row_scale[i/nc];
            if (controls.pin_pressure)
                rows[3*nc+controls.reference_cell]=
                    (state.fields[3*nc+controls.reference_cell]-
                     controls.reference_pressure)/state_scale[3];
            return rows;
        };
        const auto rhs0=make_linear_rows(current,physical_base);
        std::vector<double> rhs(rhs0.size());
        for (std::size_t i=0;i<size;++i)rhs[i]=-rhs0[i];
        const double rhs_norm=coupled_newton_norm2(rhs);
        if (rhs_norm==0.0) {
            record.reason="zero algebraic residual (physical gates not certified)";
            answer.steps.push_back(std::move(record));
            break;
        }
        const auto y0=to_y(current);
        // Optional fixed preconditioner in normalized coordinates. Never
        // replace a failed action with the identity silently.
        const auto recover=[&](const std::vector<double>& v) {
            auto z=frozen_v53.apply ? frozen_v53.apply(v)
                : native_schur ? native_schur(v,state_scale,row_scale,nc)
                : right_preconditioner ? right_preconditioner(v) : v;
            if(z.size()!=size)
                throw std::runtime_error("V52 preconditioner returned wrong dimension");
            if (controls.sst_two_by_two_preconditioner) {
                for(std::size_t i=0;i<nc;++i) {
                    const auto& inv=sst_inverse[i];
                    const double rk=v[4*nc+i]*row_scale[4];
                    const double rw=v[5*nc+i]*row_scale[5];
                    z[4*nc+i]=(inv[0]*rk+inv[1]*rw)/state_scale[4];
                    z[5*nc+i]=(inv[2]*rk+inv[3]*rw)/state_scale[5];
                }
            }
            (void)coupled_newton_norm2(z);
            return z;
        };
        // V54: Amat is a fixed, explicit linear map for every Arnoldi action,
        // independent of Pmat and of the previous direction history. Only
        // bounded small cases may build this expensive forensic reference.
        std::unique_ptr<FrozenPhysicalJacobianV54> v54_amatrix;
        if(controls.v54_frozen_jacobian) {
            v54_amatrix=std::make_unique<FrozenPhysicalJacobianV54>(
                current,physical_base,state_scale,row_scale,
                controls.pin_pressure,controls.reference_cell,
                controls.reference_pressure,controls.difference_step,
                evaluate,controls.k_min,controls.omega_min);
            record.v54_frozen_jacobian_ready=true;
            record.v54_jacobian_evaluations=v54_amatrix->statistics().evaluations;
            record.v54_centered_columns=v54_amatrix->statistics().centered_columns;
            record.v54_one_sided_columns=v54_amatrix->statistics().one_sided_columns;
            record.candidate_evaluations+=record.v54_jacobian_evaluations;
        }
        const auto apply=[&](const std::vector<double>& step) {
            if(v54_amatrix)return v54_amatrix->apply(step);
            const double norm=coupled_newton_norm2(step);
            if(norm==0.0) return std::vector<double>(size,0.0);
            double h=controls.difference_step*(1.0+coupled_newton_norm2(y0)) /
                     std::max(norm,1e-300);
            PhysicalCandidateV51 trial,trial_minus;
            // Stay strictly inside the admissible k/omega region for Jv.
            bool valid=false;
            for(int attempt=0;attempt<22;++attempt) {
                auto y=y0;
                for(std::size_t i=0;i<size;++i)y[i]+=h*step[i];
                trial=to_state(y);
                const bool positive_valid=is_admissible(trial);
                bool negative_valid=true;
                if(controls.v53_central_difference) {
                    auto ym=y0;
                    for(std::size_t i=0;i<size;++i)ym[i]-=h*step[i];
                    trial_minus=to_state(ym);
                    negative_valid=is_admissible(trial_minus);
                }
                if (positive_valid && negative_valid) {valid=true;break;}
                h*=0.5;
            }
            if(!valid || !(h>1e-15))
                throw std::runtime_error("V52 no admissible directional probe");
            const auto physical=evaluate(trial);
            ++record.candidate_evaluations;
            const auto perturbed=make_linear_rows(trial,physical);
            std::vector<double> result(size);
            if(controls.v53_central_difference) {
                const auto physical_minus=evaluate(trial_minus);
                ++record.candidate_evaluations;
                const auto minus_rows=make_linear_rows(trial_minus,physical_minus);
                if(controls.v57_hybrid_convective_jv){
                    // dX/dy, including frozen nondimensional state scaling.
                    std::vector<double> physical_direction(size);
                    for(std::size_t i=0;i<size;++i)
                        physical_direction[i]=step[i]*state_scale[i/nc];
                    const auto hybrid=combine_physical_jv_v57(
                        current,physical_base,physical,physical_minus,
                        physical_direction,state_scale,row_scale,h,
                        v57_convective,controls.pin_pressure,
                        controls.reference_cell);
                    ++record.v57_hybrid_jv_applications;
                    record.v57_maximum_convective_discrepancy=std::max(
                        record.v57_maximum_convective_discrepancy,
                        hybrid.convective_relative_discrepancy);
                    record.v57_nonsmooth_faces=std::max(
                        record.v57_nonsmooth_faces,hybrid.nonsmooth_faces);
                    // V58 generalized selection is still experimental.
                    // Never weaken the V57 active-set gate or the physical
                    // line-search/true Krylov gates without independent proof.
                    record.v58_generalized_sou=controls.v58_symmetric_sou_ties;
                    record.v57_convective_fidelity_ok &=
                        hybrid.nonsmooth_faces==0 &&
                        hybrid.convective_relative_discrepancy<=
                            controls.v57_convective_fidelity_tolerance;
                    result=hybrid.action;
                } else {
                    for(std::size_t i=0;i<size;++i)
                        result[i]=(perturbed[i]-minus_rows[i])/(2.0*h);
                }
            }else {
                for(std::size_t i=0;i<size;++i)
                    result[i]=(perturbed[i]-rhs0[i])/h;
            }
            (void)coupled_newton_norm2(result);
            return result;
        };
        // Each actual preconditioned Krylov basis vector is compared to
        // independent finite differences Fphys(X+h*z) at h and h/2. Thus
        // freezing Amat cannot silently hide its inaccuracy to the nonlinear F.
        const std::function<bool(const std::vector<double>&)> basis_audit=
            v54_amatrix ? std::function<bool(const std::vector<double>&)>(
                [&](const std::vector<double>& v) {
                    const auto errors=v54_amatrix->audit_direction(v,evaluate);
                    ++record.v54_basis_audits;
                    record.v54_maximum_jv_fidelity=std::max(
                        record.v54_maximum_jv_fidelity,errors.fidelity);
                    record.v54_maximum_halving_discrepancy=std::max(
                        record.v54_maximum_halving_discrepancy,errors.halving);
                    record.v54_maximum_sided_nonsmoothness=std::max(
                        record.v54_maximum_sided_nonsmoothness,
                        errors.sided_nonsmoothness);
                    const bool ok=errors.fidelity<=controls.v54_fidelity_tolerance &&
                        errors.halving<=controls.v54_fidelity_tolerance &&
                        errors.sided_nonsmoothness<=controls.v54_fidelity_tolerance;
                    for(std::size_t block=0;block<6;++block) {
                        record.v54_block_fidelity[block]=std::max(
                            record.v54_block_fidelity[block],errors.block_fidelity[block]);
                        record.v54_block_sided[block]=std::max(
                            record.v54_block_sided[block],errors.block_sided[block]);
                    }
                    record.v54_flux_sign_flips=std::max(
                        record.v54_flux_sign_flips,errors.flux_sign_flips);
                    record.v54_flux_asymmetry=std::max(
                        record.v54_flux_asymmetry,errors.flux_sided_asymmetry);
                    record.v54_krylov_fidelity_ok &=ok;
                    // When V57 explicitly asks for a full Krylov trace,
                    // keep sampling the other basis directions after a
                    // failed fidelity probe; the final Newton decision
                    // remains fail-closed.
                    return ok || controls.v57_true_krylov_each_step;
                }) : std::function<bool(const std::vector<double>&)>{};
        std::vector<double> direction;
        const std::function<void(const KrylovTruthSnapshotV57&)> v57_observer =
            controls.v57_true_krylov_each_step ?
            std::function<void(const KrylovTruthSnapshotV57&)>(
                [&](const KrylovTruthSnapshotV57& snap){
                    record.v57_krylov_trace.push_back(snap);
                    record.v57_maximum_arnoldi_gap=std::max(
                        record.v57_maximum_arnoldi_gap,snap.arnoldi_gap);
                    record.v57_maximum_action_nonadditivity=std::max(
                        record.v57_maximum_action_nonadditivity,
                        snap.action_nonadditivity_relative);
                    record.v57_krylov_consistent &=
                        snap.arnoldi_gap<=controls.v57_krylov_consistency_tolerance &&
                        snap.action_nonadditivity_relative<=
                            controls.v57_krylov_consistency_tolerance;
                }) : std::function<void(const KrylovTruthSnapshotV57&)>{};
        std::vector<std::vector<double>> v60_bases;
        const std::function<void(const std::vector<double>&)> v60_observer=
            controls.v60_remainder_audit ?
            std::function<void(const std::vector<double>&)>(
                [&](const std::vector<double>& z){
                    if(v60_bases.size()<2) v60_bases.push_back(z);
                }) : std::function<void(const std::vector<double>&)>{};
        const bool linear_accepted=global_physical_fgmres_v52(
            apply,recover,rhs,controls.krylov_dimension,
            controls.linear_relative_tolerance,direction,
            record.linear_true_relative,&record.linear_arnoldi_relative,
            &record.jv_linearity_relative,basis_audit,v57_observer,v60_observer);
        if(controls.v60_remainder_audit) {
            if(v60_bases.size()==2) {
                std::vector<double> v(size),w(size),total(size);
                for(std::size_t i=0;i<size;++i) {
                    const double factor=state_scale[i/nc];
                    v[i]=v60_bases[0][i]*factor;
                    w[i]=v60_bases[1][i]*factor;
                    total[i]=v60_bases[0][i]+v60_bases[1][i];
                }
                // Same displacement in normalized coordinates for v,w,v+w.
                // Every sample independently rebuilds Fphys and the V57 tape.
                const double h=controls.difference_step*
                    (1.+coupled_newton_norm2(y0))/
                    std::max({coupled_newton_norm2(v60_bases[0]),
                        coupled_newton_norm2(v60_bases[1]),
                        coupled_newton_norm2(total),1e-12});
                record.v60_remainder=audit_physical_remainder_additivity_v60(
                    current,v,w,h,evaluate,controls.k_min,controls.omega_min);
                record.candidate_evaluations+=record.v60_remainder.evaluations;
                const auto& audit=record.v60_remainder;
                record.v60_remainder_consistent=audit.audited &&
                    audit.global_remainder<=controls.v57_krylov_consistency_tolerance &&
                    std::all_of(audit.remainder_blocks.begin(),audit.remainder_blocks.end(),
                        [&](double defect){return defect<=controls.v57_krylov_consistency_tolerance;});
                // If a cell block or rebuilt RC/convective face term is
                // nonadditive, do not certify a full physical Jacobian.
                record.v60_full_additive=record.v60_remainder_consistent &&
                    audit.rhie_chow_flux<=controls.v57_krylov_consistency_tolerance &&
                    audit.convection<=controls.v57_krylov_consistency_tolerance;
            }
            // A missing second basis is NOT a certificate of additivity.
        }
        if(controls.v60_remainder_audit && !record.v60_full_additive) {
            record.reason="V60 freshly reassembled blocks/RC/convection fail Krylov additivity";
            answer.steps.push_back(std::move(record));
            break;
        }
        if(!record.v57_convective_fidelity_ok) {
            record.reason="V57 active-set convection differs from freshly reassembled FV derivative";
            answer.steps.push_back(std::move(record));
            break;
        }
        if(!record.v57_krylov_consistent) {
            record.reason="V57 true residual / Arnoldi mismatch at an intermediate Krylov direction";
            answer.steps.push_back(std::move(record));
            break;
        }
        if (!record.v54_krylov_fidelity_ok) {
            record.reason="V54 physical derivative fidelity failed for at least one Krylov direction";
            answer.steps.push_back(std::move(record));
            break;
        }
        if (!linear_accepted) {
            record.reason=!record.v54_krylov_fidelity_ok
                ? "V54 frozen Amat fails independent Fphys Krylov basis fidelity"
                : (record.jv_linearity_relative>
                    std::max(0.1,0.5*controls.linear_relative_tolerance))
                  ? "nonadditive matrix-free Jv reconstruction: false Arnoldi contraction"
                  : "unqualified true Jv-Krylov residual";
            answer.steps.push_back(std::move(record));
            break;
        }
        if(basis_audit && !basis_audit(direction)) {
            record.reason="V54 final Newton direction disagrees with independent Fphys derivative";
            answer.steps.push_back(std::move(record));
            break;
        }
        for(std::size_t bt=0;bt<=controls.maximum_backtracks;++bt) {
            const double alpha=std::ldexp(1.0,-static_cast<int>(bt));
            auto y=y0;
            for(std::size_t i=0;i<size;++i)y[i]+=alpha*direction[i];
            auto trial=to_state(y);
            if(!is_admissible(trial)) {++record.rejected_trials;continue;}
            PhysicalResidualV51 candidate;
            try {candidate=evaluate(trial);++record.candidate_evaluations;}
            catch(const std::exception&) {
                ++record.rejected_trials;continue;
            }
            // The only acceptance measure is the *newly reassembled* V51
            // physical merit. A small Krylov residual is never a certificate.
            if (candidate.merit<record.merit_before &&
                candidate.merit<=record.merit_before*
                    (1.0-controls.sufficient_decrease*alpha)) {
                record.accepted=true;
                record.alpha=alpha;
                record.merit_after=candidate.merit;
                record.reason="fresh physical merit contracted (private state)";
                answer.last_private_state=std::move(trial);
                answer.last_private_residual=std::move(candidate);
                ++answer.accepted_steps;
                break;
            }
            ++record.rejected_trials;
        }
        if(!record.accepted)
            record.reason="no candidate satisfies original physical merit decrease";
        answer.steps.push_back(std::move(record));
        if(!answer.steps.back().accepted)break;
    }
    // This runner cannot certify production: no mutation of U,p,phi,SST.
    return answer;
}
} // namespace cfdx::physics
