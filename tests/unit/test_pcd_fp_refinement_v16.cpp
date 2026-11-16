#include "cfdx/physics/pcd_fp_refinement_v16.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#undef assert
#define assert(pred) do { if (!(pred)) throw std::runtime_error("V16 check failed: " #pred); } while (0)
#include <limits>
using cfdx::core::Vector;
using cfdx::core::SolverResult;
using cfdx::core::SolverStatus;
using cfdx::physics::refine_pcd_fp_v16;

// Non-trivial diagonal operator; independent original-operator residual is
// recomputed, rather than trusting the nested solver's CONVERGED status.
int main() {
    const double tolerance=1.e-8;
    Vector x(2,0.0);
    x(0)=1.0; x(1)=2.0+8.e-8;
    SolverResult result;
    result.status=SolverStatus::CONVERGED;
    unsigned requested=0;
    const auto true_residual=[](const Vector& state) {
        return std::hypot(3.0-3.0*state(0),4.0-2.0*state(1))/5.0;
    };
    auto solve=[&](Vector& candidate, double requested_tolerance) {
        ++requested;
        assert(requested_tolerance<tolerance);
        candidate(0)=1.0;
        candidate(1)=2.0+requested_tolerance;
        SolverResult r;
        r.status=SolverStatus::CONVERGED;
        return r;
    };
    auto verify=[&](const SolverResult& r, const Vector& candidate) {
        return r.status==SolverStatus::CONVERGED &&
            true_residual(candidate)<=tolerance;
    };
    const auto certified=refine_pcd_fp_v16(x,result,tolerance,2,solve,verify);
    assert(certified.certified && certified.extra_solves==1 && requested==1);
    assert(true_residual(x)<=tolerance);
    // Status-only false positive must never pass if the returned iterate is
    // still outside the true-residual gate even after tighter requests.
    x(1)=2.0+8.e-8;
    requested=0;
    const auto failed=refine_pcd_fp_v16(x,result,tolerance,2,
        [&](Vector&,double){++requested;SolverResult r;
            r.status=SolverStatus::CONVERGED;return r;},verify);
    assert(!failed.certified && failed.extra_solves==2 && requested==2);
    assert(true_residual(x)>tolerance);
    bool invalid=false;
    try {refine_pcd_fp_v16(x,result,0.,2,solve,verify);}
    catch(const std::invalid_argument&){invalid=true;}
    assert(invalid);
    // A failed solve producing non-finite state must not turn into a pass.
    x(1)=2.0+8.e-8;
    const auto nonfinite=refine_pcd_fp_v16(x,result,tolerance,2,
       [&](Vector& trial,double) { trial(1)=std::numeric_limits<double>::quiet_NaN();
          SolverResult r;r.status=SolverStatus::CONVERGED;return r;},verify);
    assert(!nonfinite.certified && x.is_valid());
    std::cout<<"PCD V16 strict refinement: true residual certified; failures rejected\n";
}
