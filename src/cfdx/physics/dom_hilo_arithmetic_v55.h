#pragma once
#if defined(__FAST_MATH__)
#error "V55 hi-lo arithmetic is invalid with -ffast-math"
#endif
// CFDX V55: deterministic, CPU-only two-double arithmetic for RESEARCH.
// No production, CUDA, MPI, or FP-contract/fast-math qualification is implied.
// The public representation contains the physically stored hi AND lo fields.
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {
struct DomHiLoV55 {
    double hi=0.,lo=0.;
    DomHiLoV55()=default;
    DomHiLoV55(double x) : hi(x),lo(0.) {}
    DomHiLoV55(double h,double l) : hi(h),lo(l) {}
};
static_assert(std::numeric_limits<double>::is_iec559,
              "V55 requires IEEE-754 double arithmetic");

inline DomHiLoV55 dom_two_sum_v55(double a,double b) noexcept {
    const double s=a+b;
    const double bb=s-a;
    return {s,(a-(s-bb))+(b-bb)};
}
inline DomHiLoV55 dom_normalize_v55(double h,double l) noexcept {
    return dom_two_sum_v55(h,l);
}
inline DomHiLoV55 operator+(DomHiLoV55 a,DomHiLoV55 b) noexcept {
    const auto top=dom_two_sum_v55(a.hi,b.hi);
    const auto bottom=dom_two_sum_v55(a.lo,b.lo);
    const auto t=dom_two_sum_v55(top.lo,bottom.hi);
    const auto s=dom_two_sum_v55(top.hi,t.hi);
    return dom_normalize_v55(s.hi,s.lo+t.lo+bottom.lo);
}
inline DomHiLoV55 operator-(DomHiLoV55 a) noexcept {
    return {-a.hi,-a.lo};
}
inline DomHiLoV55 operator-(DomHiLoV55 a,DomHiLoV55 b) noexcept {
    return a+(-b);
}
inline DomHiLoV55 operator*(DomHiLoV55 a,DomHiLoV55 b) noexcept {
    const double p=a.hi*b.hi;
    const double e=std::fma(a.hi,b.hi,-p);
    // The cross terms and retained low*low term protect cancellation.
    const auto t=dom_two_sum_v55(a.hi*b.lo,a.lo*b.hi);
    const auto u=dom_two_sum_v55(e,t.hi);
    return dom_normalize_v55(p,u.hi+(u.lo+t.lo+a.lo*b.lo));
}
inline DomHiLoV55 operator/(DomHiLoV55 a,DomHiLoV55 b) {
    if((b.hi==0.&&b.lo==0.) || !std::isfinite(b.hi) ||
       !std::isfinite(b.lo) || !std::isfinite(a.hi) ||
       !std::isfinite(a.lo))
        throw std::invalid_argument("V55 hi-lo: nonfinite division or zero divisor");
    const double q1=a.hi/b.hi;
    const auto r1=a-b*DomHiLoV55(q1);
    const double q2=r1.hi/b.hi;
    const auto r2=r1-b*DomHiLoV55(q2);
    const double q3=r2.hi/b.hi;
    const auto out=DomHiLoV55(q1)+DomHiLoV55(q2)+DomHiLoV55(q3);
    if(!std::isfinite(out.hi)||!std::isfinite(out.lo))
        throw std::overflow_error("V55 hi-lo: division overflow");
    return out;
}
inline DomHiLoV55& operator+=(DomHiLoV55& a,DomHiLoV55 b) noexcept {
    a=a+b;return a;
}
inline DomHiLoV55& operator-=(DomHiLoV55& a,DomHiLoV55 b) noexcept {
    a=a-b;return a;
}
inline bool dom_finite_v55(DomHiLoV55 a) noexcept {
    return std::isfinite(a.hi)&&std::isfinite(a.lo);
}
inline long double dom_to_long_double_v55(DomHiLoV55 a) noexcept {
    // AUDIT ONLY; neither the operator nor the state update calls long double.
    return static_cast<long double>(a.hi)+static_cast<long double>(a.lo);
}
inline DomHiLoV55 dom_from_long_double_v55(long double v) {
    if(!std::isfinite(v))
        throw std::invalid_argument("V55 hi-lo: nonfinite reference value");
    const double h=static_cast<double>(v);
    const double l=static_cast<double>(v-static_cast<long double>(h));
    if(!std::isfinite(h)||!std::isfinite(l))
        throw std::overflow_error("V55 hi-lo: nonfinite split");
    return {h,l};
}
// Split from the correctly rounded 4*pi, not from a truncated decimal
// literal. The low component is needed when comparing with the V54 80-bit
// extended reference; this is NOT the V48 double constant.
inline DomHiLoV55 dom_four_pi_v55() noexcept {
    return {12.566370614359172,4.8985871965894128e-16};
}
} // namespace cfdx::physics
