#pragma once
// V43 opt-in, read-only diagnostic. Contrast COMPLETE committed outer states
// n, n-1 and n-4 instead of inferring a field cycle from a periodic scalar
// merit. Does not participate in any nonlinear decision or acceptance gate.
// Memory is allocated once (five history snapshots plus a staged buffer); no allocation
// for this feature when disabled. Restricted to a bounded-size CPU audit.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct OuterFieldPeriodReportV43 {
    bool has_four_back = false;
    std::array<double, 8> relative_four_back{};
    std::array<double, 8> relative_previous{};
    double maximum_four_back = 0.0;
    double maximum_previous = 0.0;
    std::size_t worst_channel_four_back = 0;
    std::size_t worst_index_four_back = 0;
    std::size_t consecutive_four_back = 0;
    bool period_four_candidate = false;
};

class OuterFieldPeriodProbeV43 {
public:
    // Channels: Ux, Uy, Uz, gauge-invariant p, conservative mass flux,
    // k, omega, nu_t. Face flux has n_faces values, the others n_cells.
    static constexpr std::size_t channels = 8;
    explicit OuterFieldPeriodProbeV43(std::size_t n_cells,
                                     std::size_t n_faces)
        : n_cells_(n_cells), n_faces_(n_faces) {
        if (!n_cells || !n_faces || n_cells > 500000 || n_faces > 3000000)
            throw std::invalid_argument("V43 field cycle: invalid diagnostic size");
        for (auto& snapshot : snapshots_)
            for (std::size_t ch = 0; ch < channels; ++ch)
                snapshot[ch].resize(length(ch));
        for (std::size_t ch=0;ch<channels;++ch)
            pending_[ch].resize(length(ch));
    }

    // Getter(ch, index) gives the COMMITTED state, not a pending trial.
    // No dynamic allocation on subsequent samples and no hidden I/O.
    template<class Getter>
    OuterFieldPeriodReportV43 record(Getter&& get) {
        // Strong exception guarantee: stage the complete sample separately
        // from the retained orbit. If ANY getter throws/returns nonfinite,
        // not even the ring slot about to be overwritten is modified. Once
        // all values are validated, exchange preallocated vector buffers
        // without allocation. Cost O(N) CPU and one extra snapshot only when
        // the diagnostic is explicitly enabled.
        for (std::size_t ch=0; ch<channels; ++ch)
            for (std::size_t i=0; i<length(ch); ++i) {
                const double value = get(ch,i);
                if (!std::isfinite(value))
                    throw std::invalid_argument("V43 field cycle: nonfinite committed state");
                pending_[ch][i]=value;
            }
        const auto slot = count_ % snapshots_.size();
        auto& current = snapshots_[slot];
        for (std::size_t ch=0; ch<channels; ++ch)
            current[ch].swap(pending_[ch]);
        OuterFieldPeriodReportV43 report;
        if (count_ >= 1) {
            const auto& previous=snapshots_[(slot+snapshots_.size()-1)%snapshots_.size()];
            compare(current, previous, report.relative_previous,
                    report.maximum_previous, nullptr, nullptr);
        }
        if (count_ >= 4) {
            report.has_four_back=true;
            const auto& old=snapshots_[(slot+1)%snapshots_.size()];
            compare(current,old,report.relative_four_back,
                    report.maximum_four_back,&report.worst_channel_four_back,
                    &report.worst_index_four_back);
            // The four-step return must be small for EVERY active channel,
            // while successive states differ substantially; one coincidental
            // near-return is insufficient. This is a diagnostic hypothesis,
            // not a mathematical proof of a periodic nonlinear orbit.
            const bool near=report.maximum_four_back <= 1e-3 &&
                report.maximum_previous >=
                    std::max(1e-4, 10.0*report.maximum_four_back);
            consecutive_=near ? consecutive_+1 : 0;
            report.consecutive_four_back=consecutive_;
            report.period_four_candidate=consecutive_ >= 8;
        }
        ++count_;
        period_candidate_=report.period_four_candidate;
        return report;
    }

    std::size_t observed() const noexcept { return count_; }
    bool period_four_candidate() const noexcept { return period_candidate_; }
    // Read-only four-state Cesaro average, for a separate experimental
    // predictor. The caller must reassemble U/p/phi and SST on the candidate
    // and reject it transactionally if physical gates are not met.
    double average_last_four(std::size_t channel, std::size_t index) const {
        if (count_<4 || channel>=channels || index>=length(channel))
            throw std::invalid_argument("V43 four-state average unavailable");
        double sum=0.0;
        for (std::size_t j=0;j<4;++j)
            sum+=snapshots_[(count_-1-j)%snapshots_.size()][channel][index];
        return sum*0.25;
    }

private:
    using Snapshot=std::array<std::vector<double>, channels>;
    std::size_t length(std::size_t channel) const noexcept {
        return channel==4 ? n_faces_ : n_cells_;
    }
    void compare(const Snapshot& now, const Snapshot& earlier,
                 std::array<double,channels>& result, double& max_relative,
                 std::size_t* worst_channel, std::size_t* worst_index) const {
        max_relative=0.0;
        for (std::size_t ch=0; ch<channels; ++ch) {
            double abs_max=0.0, magnitude=0.0;
            std::size_t argmax=0;
            for (std::size_t i=0; i<length(ch); ++i) {
                const double distance=std::abs(now[ch][i]-earlier[ch][i]);
                magnitude=std::max({magnitude,std::abs(now[ch][i]),
                                   std::abs(earlier[ch][i])});
                if (distance>abs_max) { abs_max=distance;argmax=i; }
            }
            // Remove only components that are identically zero to machine
            // precision (e.g., Uz on a truly 2D extrusion).
            result[ch] = magnitude <= 1e-25 ? 0.0 :
                abs_max/std::max(magnitude,1e-25);
            if (result[ch]>max_relative) {
                max_relative=result[ch];
                if (worst_channel) *worst_channel=ch;
                if (worst_index) *worst_index=argmax;
            }
        }
    }

    std::size_t n_cells_,n_faces_;
    std::array<Snapshot,5> snapshots_;
    Snapshot pending_; // transactionally staged, preallocated once
    std::size_t count_=0;
    std::size_t consecutive_=0;
    bool period_candidate_=false;
};

} // namespace cfdx::physics
