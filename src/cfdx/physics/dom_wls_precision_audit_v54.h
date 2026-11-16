#pragma once
// V54 reference-only precision audit for <=640 unknowns.  The explicit dense
// 80-bit solve is a DIAGNOSTIC, not a scalable method or a production fallback.
// The authoritative certificate remains the true V48 double WLS residual.
#include "cfdx/physics/dom_stable_wls_v54.h"
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
struct DomWlsPrecisionAuditV54 {
    std::size_t cells=0;
    double extended_l2=std::numeric_limits<double>::infinity();
    double rounded_dense_l2=std::numeric_limits<double>::infinity();
    double split_hi_lo_l2=std::numeric_limits<double>::infinity();
    double rounded_stable_l2=std::numeric_limits<double>::infinity();
    double authoritative_v48_l2=std::numeric_limits<double>::infinity();
    bool extended_strict=false;
    bool original_strict=false;
    // Deliberately no accepted_J: numerical tests MUST NOT publish the 80-bit
    // state as a converged production/CFDX low-order solution.
};

inline DomWlsPrecisionAuditV54 audit_dense_precision_v54(
    const DomStableWlsV54& stable,const DomSmmGrayWallV49& reference,
    const std::vector<double>& rhs,double strict_tolerance=1e-13,
    std::size_t max_reference_dimension=640,
    std::vector<DomHiLoV55>* diagnostic_hilo_state_v55=nullptr) {
    if(!std::isfinite(strict_tolerance) || !(strict_tolerance>0.) ||
       rhs.size()!=reference.dimension() || rhs.size()<8 ||
       rhs.size()>max_reference_dimension || max_reference_dimension>640 ||
       std::numeric_limits<long double>::digits<=
           std::numeric_limits<double>::digits)
        throw std::invalid_argument("V54 dense precision audit: unsupported size/precision");
    const auto n=rhs.size();
    DomWlsPrecisionAuditV54 out;out.cells=n;
    std::vector<long double> matrix(n*n,0.),b(n),e(n,0.);
    for(std::size_t i=0;i<n;++i){
        if(!std::isfinite(rhs[i]))
            throw std::invalid_argument("V54 dense precision audit: nonfinite RHS");
        b[i]=static_cast<long double>(rhs[i]);
    }
    for(std::size_t col=0;col<n;++col) {
        e[col]=1.;
        auto a=stable.apply_extended(e);
        e[col]=0.;
        for(std::size_t row=0;row<n;++row)
            matrix[row*n+col]=a[row];
    }
    // Deliberately copy for LU: the unfactored matrix is also used for an
    // independently accumulated explicit-matrix residual.
    auto lu=matrix;
    for(std::size_t k=0;k<n;++k){
        std::size_t pivot=k;
        for(std::size_t i=k+1;i<n;++i)
            if(std::abs(lu[i*n+k])>std::abs(lu[pivot*n+k]))pivot=i;
        if(!(std::abs(lu[pivot*n+k])>0.) ||
           !std::isfinite(lu[pivot*n+k]))
            throw std::runtime_error("V54 dense precision audit: singular matrix");
        if(pivot!=k){
            for(std::size_t j=k;j<n;++j)
                std::swap(lu[pivot*n+j],lu[k*n+j]);
            std::swap(b[pivot],b[k]);
        }
        const auto diagonal=lu[k*n+k];
        for(std::size_t i=k+1;i<n;++i) {
            const auto factor=lu[i*n+k]/diagonal;
            lu[i*n+k]=0.;
            for(std::size_t j=k+1;j<n;++j)
                lu[i*n+j]-=factor*lu[k*n+j];
            b[i]-=factor*b[k];
        }
    }
    std::vector<long double> x(n);
    for(std::size_t i=n;i-->0;){
        long double sum=b[i];
        for(std::size_t j=i+1;j<n;++j)
            sum-=lu[i*n+j]*x[j];
        x[i]=sum/lu[i*n+i];
    }
    std::vector<double> rounded(n);
    for(std::size_t i=0;i<n;++i){
        if(!(x[i]>=0.) || !std::isfinite(x[i]))
            throw std::runtime_error("V54 dense precision audit: nonphysical reference moment");
        rounded[i]=static_cast<double>(x[i]);
    }
    // V55 explicit OPT-IN: research state is diagnostic and MUST NOT be
    // passed to a production solver or published from this routine.
    if(diagnostic_hilo_state_v55) {
        diagnostic_hilo_state_v55->resize(n);
        for(std::size_t i=0;i<n;++i)
            (*diagnostic_hilo_state_v55)[i]=dom_from_long_double_v55(x[i]);
    }
    // Storage experiment only: two double components preserve the 80-bit
    // reference state, but the evaluator STILL uses 80-bit arithmetic.
    // This does not qualify a GPU double-double kernel or a production solve.
    std::vector<long double> split_state(n);
    for(std::size_t i=0;i<n;++i){
        const double hi=rounded[i];
        const double lo=static_cast<double>(x[i]-static_cast<long double>(hi));
        split_state[i]=static_cast<long double>(hi)+static_cast<long double>(lo);
    }
    const auto exact=stable.apply_extended(x);
    const auto split=stable.apply_extended(split_state);
    const auto rounded_stable=stable.apply(rounded);
    const auto rounded_original=reference.apply(rounded);
    long double r80=0.,rsplit=0.,rdense=0.,rstable=0.,roriginal=0.,den=0.;
    for(std::size_t i=0;i<n;++i){
        long double sum=0.;
        for(std::size_t j=0;j<n;++j)
            sum+=matrix[i*n+j]*static_cast<long double>(rounded[j]);
        const auto delta80=exact[i]-static_cast<long double>(rhs[i]);
        const auto delta_split=split[i]-static_cast<long double>(rhs[i]);
        const auto delta_dense=sum-static_cast<long double>(rhs[i]);
        const auto delta_stable=static_cast<long double>(rounded_stable[i])-rhs[i];
        const auto delta_original=static_cast<long double>(rounded_original[i])-rhs[i];
        r80+=delta80*delta80;
        rsplit+=delta_split*delta_split;
        rdense+=delta_dense*delta_dense;
        rstable+=delta_stable*delta_stable;
        roriginal+=delta_original*delta_original;
        den+=static_cast<long double>(rhs[i])*rhs[i];
    }
    if(!(den>0.))throw std::invalid_argument("V54 dense precision audit: zero RHS norm");
    out.extended_l2=static_cast<double>(std::sqrt(r80/den));
    out.split_hi_lo_l2=static_cast<double>(std::sqrt(rsplit/den));
    out.rounded_dense_l2=static_cast<double>(std::sqrt(rdense/den));
    out.rounded_stable_l2=static_cast<double>(std::sqrt(rstable/den));
    out.authoritative_v48_l2=static_cast<double>(std::sqrt(roriginal/den));
    out.extended_strict=out.extended_l2<=strict_tolerance;
    out.original_strict=out.authoritative_v48_l2<=strict_tolerance;
    return out;
}
} // namespace cfdx::physics
