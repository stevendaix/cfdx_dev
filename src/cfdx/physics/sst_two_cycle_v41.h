#pragma once
// V41: detect a persistent two-step nonlinear SST orbit. This is an
// optional NUMERICAL globalization aid, not a turbulence-model change.
// A single bad cell cannot silently certify a solution: this class only
// returns evidence that permits an explicitly selected damping reduction.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::physics {

struct SSTTwoCycleEvidenceV41 {
    bool detected = false;
    std::size_t cell = 0;
    std::uint8_t consecutive = 0;
    double relative_one_step = 0.0;
    double relative_two_step = 0.0;
};

class SSTTwoCycleAuditV41 {
public:
    explicit SSTTwoCycleAuditV41(std::size_t cells)
        : cells_(cells), consecutive_(cells,0)
    {
        if (!cells) throw std::invalid_argument("SST V41 cycle audit: empty field");
    }

    SSTTwoCycleEvidenceV41 observe(const std::vector<double>& k,
                                    const std::vector<double>& omega) {
        if (k.size()!=cells_ || omega.size()!=cells_)
            throw std::invalid_argument("SST V41 cycle audit: field shape mismatch");
        for (std::size_t i=0;i<cells_;++i)
            if (!std::isfinite(k[i]) || k[i]<0.0 ||
                !std::isfinite(omega[i]) || omega[i]<=0.0)
                throw std::invalid_argument("SST V41 cycle audit: invalid k/omega state");
        SSTTwoCycleEvidenceV41 evidence;
        if (samples_>=2) {
            for (std::size_t i=0;i<cells_;++i) {
                const auto relative=[](double a,double b) {
                    return std::abs(a-b)/std::max({1.0,std::abs(a),std::abs(b)});
                };
                const double one=std::max(relative(k[i],last_k_[i]),
                                          relative(omega[i],last_omega_[i]));
                const double two=std::max(relative(k[i],older_k_[i]),
                                          relative(omega[i],older_omega_[i]));
                // Strong alternating motion with a nearly repeated state two
                // steps earlier. Consecutive observations rule out an isolated
                // transition through a non-smooth SST blending branch.
                const bool recurring=one>0.05 && two<1e-3 && two<0.01*one;
                consecutive_[i]=recurring
                    ? static_cast<std::uint8_t>(std::min<unsigned>(
                          255u,static_cast<unsigned>(consecutive_[i])+1u)) : 0;
                if (consecutive_[i]>=3 &&
                    (!evidence.detected || one>evidence.relative_one_step)) {
                    evidence={true,i,consecutive_[i],one,two};
                }
            }
        }
        // Swap the two retained states rather than reallocating O(N) host
        // buffers on every iteration. After the first two observations both
        // buffers keep their capacity (important for the 20,540-cell BFS).
        if (samples_>=1) {
            older_k_.swap(last_k_);
            older_omega_.swap(last_omega_);
        }
        last_k_.assign(k.begin(),k.end());
        last_omega_.assign(omega.begin(),omega.end());
        ++samples_;
        return evidence;
    }
private:
    std::size_t cells_;
    std::size_t samples_ = 0;
    std::vector<std::uint8_t> consecutive_;
    std::vector<double> older_k_,older_omega_,last_k_,last_omega_;
};

} // namespace cfdx::physics
