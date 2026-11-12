#include "cfdx/physics/outer_field_period_v43.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace cfdx::physics;
static void require(bool condition,const char* message) {
    if(!condition) throw std::runtime_error(message);
}
int main() {
    try {
        OuterFieldPeriodProbeV43 p(3,5);
        OuterFieldPeriodReportV43 report;
        for (int k=0;k<13;++k) {
            const double phase=static_cast<double>(k%4+1);
            report=p.record([&](std::size_t channel,std::size_t i) {
                if (channel==3) return phase + 1000.*k; // deliberately invalid gauge drift!
                return phase*(channel+1)*(i+1);
            });
        }
        // A drifting pressure offset is NOT gauge-invariant. The caller must
        // subtract the pressure reference value before recording it.
        require(!report.period_four_candidate,"changing uncorrected pressure not detected");
        OuterFieldPeriodProbeV43 gauge(3,5);
        for (int k=0;k<13;++k) {
            const double phase=static_cast<double>(k%4+1);
            report=gauge.record([&](std::size_t channel,std::size_t i) {
                (void)i;
                // The production caller provides p(cell)-p(reference) here.
                return phase*(channel+1)*(i+1);
            });
        }
        require(report.has_four_back && report.period_four_candidate,
                "period-four physical fields not detected");
        require(report.maximum_four_back==0.0 &&
                report.maximum_previous>0.0,
                "field four-step contrast corrupted");
        // On the four periodic phases [1,2,3,4], the mean is phase 2.5,
        // independent of the last phase (nontrivial affine predictor).
        require(std::abs(gauge.average_last_four(0,0)-2.5)<1e-15,
                "four-state velocity predictor wrong");
        require(std::abs(gauge.average_last_four(3,1)-20.0)<1e-14,
                "four-state gauge-pressure predictor wrong");
        // Nonfinite at the LAST entry of a snapshot must not corrupt any
        // earlier ring slot. After rejection, the next valid observation
        // continues from the same committed count and retains the orbit.
        bool observer_failed=false;
        const std::size_t before=gauge.observed();
        try {
            gauge.record([](std::size_t ch,std::size_t i)->double {
                if (ch==7 && i==2) return std::nan("");
                return -123.0;
            });
        } catch(const std::invalid_argument&) { observer_failed=true; }
        require(observer_failed && gauge.observed()==before,
                "rejected field sample advanced the committed record");
        report=gauge.record([&](std::size_t ch,std::size_t i) {
            const double phase=static_cast<double>(before%4+1);
            return phase*static_cast<double>(ch+1)*static_cast<double>(i+1);
        });
        require(report.period_four_candidate && report.maximum_four_back==0.0,
                "failed field observation corrupted the four-phase ring");
        OuterFieldPeriodProbeV43 contract(2,3);
        bool bad=false;
        try { (void)contract.record([](std::size_t,std::size_t) {return std::nan("");}); }
        catch (const std::invalid_argument&) {bad=true;}
        require(bad,"nonfinite physical state accepted");
        // Nonperiodic sequence that has a periodic scalar merit (fixed Ux
        // maximum, evolving pressure/nu_t) must NOT be classified as a cycle.
        OuterFieldPeriodProbeV43 nonperiodic(2,3);
        for (int k=0;k<16;++k) {
            report=nonperiodic.record([&](std::size_t ch,std::size_t i) {
                return ch==0 ? 4.0 : 1.0+i+ch + 0.03*k;
            });
        }
        require(!report.period_four_candidate && report.maximum_four_back>0.,
                "nonperiodic fields wrongly classified from a scalar norm");
        std::cout<<"OUTER_FIELD_PERIOD_V43_PASS gauge=1 four_states=1 no_false_cycle=1\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<"OUTER_FIELD_PERIOD_V43_FAIL "<<error.what()<<'\n';
        return 1;
    }
}
