#pragma once
// Safeguarded depth-two nonlinear SST candidate: history is local to the
// current SST transaction. The caller must verify the fully reassembled
// original nonlinear equations and preserve positivity before accepting it.
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
namespace cfdx::physics {
struct SSTAndersonV20Candidate {
    bool available=false;
    double theta=0.;
    std::vector<double> k,omega;
};
class SSTAndersonV20 {
public:
    SSTAndersonV20Candidate propose(const std::vector<double>& xk,
        const std::vector<double>& xw,const std::vector<double>& gk,
        const std::vector<double>& gw,double kmin,double wmin) const {
        const auto n=xk.size();
        if (!n || xw.size()!=n || gk.size()!=n || gw.size()!=n ||
            !(kmin>0.) || !(wmin>0.))
            throw std::invalid_argument("V20 Anderson: inconsistent dimensions or floors");
        SSTAndersonV20Candidate out;
        if (!ready_) return out;
        if (old_xk_.size()!=n) throw std::invalid_argument("V20 Anderson: stale history");
        long double sk=kmin,sw=wmin;
        for (std::size_t i=0;i<n;++i) {
            for (double v : {xk[i],xw[i],gk[i],gw[i],old_xk_[i],old_xw_[i],old_gk_[i],old_gw_[i]})
                if (!std::isfinite(v)) throw std::invalid_argument("V20 Anderson: nonfinite data");
            sk=std::max(sk,(long double)std::max({std::abs(xk[i]),std::abs(gk[i]),std::abs(old_xk_[i]),std::abs(old_gk_[i])}));
            sw=std::max(sw,(long double)std::max({std::abs(xw[i]),std::abs(gw[i]),std::abs(old_xw_[i]),std::abs(old_gw_[i])}));
        }
        long double dot=0.,norm=0.;
        for (std::size_t i=0;i<n;++i) {
            long double fk=((long double)gk[i]-xk[i])/sk;
            long double fw=((long double)gw[i]-xw[i])/sw;
            long double dk=((long double)old_gk_[i]-old_xk_[i])/sk-fk;
            long double dw=((long double)old_gw_[i]-old_xw_[i])/sw-fw;
            dot+=fk*dk+fw*dw;norm+=dk*dk+dw*dw;
        }
        if (!(norm>1e-34L) || !std::isfinite(norm)) return out;
        const long double theta=std::clamp(-dot/norm,0.L,1.L);
        if (!(theta>1e-8L) || !std::isfinite(theta)) return out;
        out.k.resize(n);out.omega.resize(n);
        for (std::size_t i=0;i<n;++i) {
            out.k[i]=(double)((1.L-theta)*gk[i]+theta*old_gk_[i]);
            out.omega[i]=(double)((1.L-theta)*gw[i]+theta*old_gw_[i]);
            if (!std::isfinite(out.k[i]) || !std::isfinite(out.omega[i]) ||
                out.k[i]<kmin || out.omega[i]<wmin) return {};
        }
        out.available=true;out.theta=(double)theta;
        return out;
    }
    void record(const std::vector<double>& xk,const std::vector<double>& xw,
                const std::vector<double>& gk,const std::vector<double>& gw) {
        if (xk.empty() || xk.size()!=xw.size() || xk.size()!=gk.size() || xk.size()!=gw.size())
            throw std::invalid_argument("V20 Anderson: invalid record dimensions");
        for (std::size_t i=0;i<xk.size();++i)
            if (!std::isfinite(xk[i]) || !std::isfinite(xw[i]) ||
                !std::isfinite(gk[i]) || !std::isfinite(gw[i]))
                throw std::invalid_argument("V20 Anderson: nonfinite record");
        old_xk_=xk;old_xw_=xw;old_gk_=gk;old_gw_=gw;ready_=true;
    }
private:
    bool ready_=false;
    std::vector<double> old_xk_,old_xw_,old_gk_,old_gw_;
};
}
