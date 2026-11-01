#ifdef NDEBUG
#undef NDEBUG
#endif

#include "cfdx/physics/n256_continuation.h"

#include <cassert>
#include <cmath>
#include <vector>

using namespace cfdx_n256;

int main() {
    Config c;
    c.good_needed = 2;
    c.full_nonorth_needed = 2;
    c.final_needed = 3;
    Continuation s(c);

    assert(s.controls().phase == Phase::robust_start);
    double r = 100.0;
    s.observe(r, true);
    for (int i = 0; i < 8; ++i) {
        r *= 0.5;
        const auto x = s.observe(r, true);
        assert(x.nonorth_factor >= 0.0 && x.nonorth_factor <= 1.0);
    }

    auto x = s.controls();
    assert(x.phase == Phase::full_nonorth && x.nonorth_factor == 1.0);
    for (int i = 0; i < 8; ++i) {
        r *= 0.5;
        x = s.observe(r, true);
    }
    assert(x.phase == Phase::second_order && x.second_order_weight == 1.0);
    assert(!x.qualification_ready);

    for (int i = 0; i < 3; ++i) {
        r *= 0.5;
        x = s.observe(r, true, i == 2);
    }
    assert(x.qualification_ready);

    x = s.observe(r * 0.5, true, false);
    assert(!x.qualification_ready);

    const double old_cfl = x.cfl;
    x = s.observe(r * 0.5, false);
    assert(!x.qualification_ready && x.cfl < old_cfl);
    assert(x.phase == Phase::second_order && x.second_order_weight == 1.0);

    assert(neighbor_volume_ratio(0.00543, 1.0) > 0.0054);
    const std::vector<double> h{1.0, 1.2, 1.4};
    assert(std::abs(max_adjacent_growth(h.begin(), h.end()) - 1.2) < 1e-12);
}
