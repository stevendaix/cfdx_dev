#include "cfdx/physics/sst_two_cycle_v41.h"
#include "cfdx/physics/turbulence_transport.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
using cfdx::physics::SSTTwoCycleAuditV41;
using cfdx::physics::TurbulenceTransportControls;
using cfdx::physics::validate_turbulence_controls;
static void require(bool ok) {
    if (!ok) throw std::runtime_error("SST V41 cycle qualification failed");
}
int main() {
    // Release builds must execute checks, not compile them out via NDEBUG.
    SSTTwoCycleAuditV41 cycle(2);
    bool detected=false;
    for (int n=0;n<12;++n) {
        const double s=(n%2==0)?1.0:0.0;
        const auto r=cycle.observe({1.0,8.0+3.0*s}, {2.0,3500.0+3000.0*s});
        if (r.detected) {
            require(n>=4 && r.cell==1 && r.consecutive>=3);
            require(r.relative_one_step>0.3);
            require(r.relative_two_step==0.0);
            detected=true;
        }
    }
    require(detected);
    // Smooth contraction must not be classified as a two-cycle.
    SSTTwoCycleAuditV41 smooth(1);
    for (int n=0;n<50;++n)
        require(!smooth.observe({1.0+std::exp(-0.1*n)},
                                {10.0+std::exp(-0.1*n)}).detected);
    // Irregular alternating sequence should not qualify if the n-2 match
    // remains above the guard's 1e-3 threshold.
    SSTTwoCycleAuditV41 irregular(1);
    for (int n=0;n<30;++n) {
        const double drift=0.03*n;
        const double s=(n%2==0)?1.0:0.0;
        require(!irregular.observe({10.0+s*3.0+drift},
                                   {3500.0+s*3000.0+100.0*drift}).detected);
    }
    bool threw=false;
    try {smooth.observe({1.},{});} catch(const std::invalid_argument&){threw=true;}
    require(threw);
    threw=false;
    try {smooth.observe({1.},{std::numeric_limits<double>::quiet_NaN()});}
    catch(const std::invalid_argument&){threw=true;}
    require(threw);
    TurbulenceTransportControls controls;
    require(!controls.sst_two_cycle_rescue_v41_experimental);
    validate_turbulence_controls(controls);
    controls.sst_two_cycle_rescue_v41_experimental=true;
    validate_turbulence_controls(controls);
    controls.sst_nonlinear_line_search=true;
    threw=false;
    try {validate_turbulence_controls(controls);}
    catch(const std::invalid_argument&){threw=true;}
    require(threw);
}
