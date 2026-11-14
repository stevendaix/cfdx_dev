#pragma once
// Experimental shared pseudo-time schedule for steady U/p/k/omega sweeps.
// Pressure remains an incompressibility constraint: it receives NO pseudo
// mass. One frozen CFL regularizes momentum and the coupled k-omega Newton
// direction. A new CFL is published only after the FULL physical transaction.
#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
namespace cfdx::physics {
struct SSTV31OuterPTC {
    double cfl=0.5;
    bool baseline_valid=false;
    std::array<double,5> last{{0.,0.,0.,0.,0.}}; // U, continuity, k, omega, nut
    std::size_t committed_outer=0;
    static constexpr double minimum_cfl=1./128.;
    static constexpr double maximum_cfl=16.;
    void validate() const {
        if (!(cfl>=minimum_cfl && cfl<=maximum_cfl) || !std::isfinite(cfl))
            throw std::invalid_argument("V31 outer PTC: invalid CFL");
    }
    static double max_merit(const std::array<double,5>& x) {
        double out=0.;
        for (double v:x) {
            if (!(v>=0.) || !std::isfinite(v))
                throw std::invalid_argument("V31 outer PTC: invalid physical residual");
            out=std::max(out,v);
        }
        return out;
    }
    // Called on a COPY before transaction.commit(), applied after commit.
    // On a changed discrete convection operator, never compare unlike PDEs.
    void committed(const std::array<double,5>& now, bool same_operator,
                   bool globally_converged=false) {
        validate();
        const double fresh=max_merit(now);
        const double old=max_merit(last);
        if (!baseline_valid || !same_operator) {
            last=now;
            baseline_valid=true;
            ++committed_outer;
            return;
        }
        // V14 gives an independently certified converged state precedence
        // over a machine-noise merit ratio. Do not expand CFL in that case.
        if (globally_converged) {
            last=now;
            ++committed_outer;
            return;
        }
        if (!(fresh<old) || !(old>0.))
            throw std::invalid_argument("V31 outer PTC: global physical merit did not contract");
        // Prefer shrinking a shared pseudo CFL if even one physics component
        // degrades considerably. Never relax certification thresholds.
        bool component_growth=false;
        for (std::size_t i=0;i<5;++i)
            if (now[i]>1.05*last[i]) component_growth=true;
        const double ratio=fresh/old;
        const double gain=component_growth ? 0.7 : (ratio<0.8 ? 1.2 : 1.05);
        cfl=std::clamp(cfl*gain,minimum_cfl,maximum_cfl);
        last=now;
        ++committed_outer;
    }
};
} // namespace cfdx::physics
