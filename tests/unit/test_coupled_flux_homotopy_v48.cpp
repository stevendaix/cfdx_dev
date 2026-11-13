#include "cfdx/physics/coupled_flux_homotopy_v48.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void check(bool b) { if (!b) throw std::runtime_error("V48 homotopy test failed"); }
int main() {
    CoupledInverseMomentumDiagonalV48 old{{{1.,2.},{3.,4.},{5.,6.}}};
    CoupledInverseMomentumDiagonalV48 now{{{3.,4.},{7.,8.},{9.,10.}}};
    auto a=blend_inverse_momentum_v48(old,now,0.0);
    auto b=blend_inverse_momentum_v48(old,now,1.0);
    auto m=blend_inverse_momentum_v48(old,now,0.25);
    check(a==old && b==now && m[0][0]==1.5 && m[1][1]==5.);
    check(!coupled_flux_homotopy_final_gate_v48(true,0.25));
    check(coupled_flux_homotopy_final_gate_v48(true,1.0));
    check(coupled_flux_homotopy_final_gate_v48(false,0.0));
    bool failed=false;
    try { validate_coupled_flux_homotopy_v48(true,false); }
    catch (const std::invalid_argument&) {failed=true;}
    check(failed);
    validate_coupled_flux_homotopy_v48(true,true);
    const std::vector<double> previous{0.01,-0.02,0.04};
    const std::vector<double> matching=previous;
    auto audit=check_coupled_flux_identity_v48(previous.size(),
        [&](std::size_t f){return previous[f];},
        [&](std::size_t f){return matching[f];});
    check(audit.relative_l1==0.0 && audit.absolute_linf==0.0);
    const std::vector<double> drift{0.01,-0.019,0.04};
    audit=check_coupled_flux_identity_v48(previous.size(),
        [&](std::size_t f){return previous[f];},
        [&](std::size_t f){return drift[f];});
    check(audit.worst_face==1 && audit.relative_l1>0.0);
    failed=false;
    try { blend_inverse_momentum_v48(old,now,1.1); }
    catch (const std::invalid_argument&) { failed=true; }
    check(failed);
    failed=false;
    try { auto corrupt=old; corrupt[0][0]=0.;
          blend_inverse_momentum_v48(corrupt,now,0.5); }
    catch (const std::invalid_argument&) { failed=true; }
    check(failed);
    failed=false;
    try { auto corrupt=old; corrupt[0].pop_back();
          blend_inverse_momentum_v48(corrupt,now,0.5); }
    catch (const std::invalid_argument&) { failed=true; }
    check(failed);
    failed=false;
    try { coupled_flux_homotopy_final_gate_v48(true,
        std::numeric_limits<double>::quiet_NaN()); }
    catch (const std::invalid_argument&) { failed=true; }
    check(failed);
}
