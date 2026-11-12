#include "cfdx/physics/coupled_phase_escape_v45.h"
#include "cfdx/physics/coupled_outer_step_v43.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>
using namespace cfdx::physics;
static void check(bool value,const char* why) {
    if (!value) throw std::runtime_error(why);
}
int main() {
    try {
        CoupledPhaseEscapeLedgerV45 ledger;
        // Disabled is a strict no-op even if a four-cycle is reported.
        auto t=ledger.trial(59,true,0.0,12,1.0);
        check(!t.active && t.beta==1.0 && ledger.active_through()==0,
              "disabled V45 mutated nonlinear state");
        // A trigger detected from committed full fields is a trial only.
        t=ledger.trial(59,true,0.5,4,1.0);
        check(t.active && t.triggered && t.beta==0.5 &&
              t.remaining_including_this==4, "fresh V45 trigger wrong");
        auto retry=ledger.trial(59,true,0.5,4,0.5);
        check(retry.beta==0.25 && ledger.active_through()==0,
              "failed trial advanced V45 window or did not shrink");
        // Commit the original accepted trial. Only accepted states advance it.
        ledger.commit(t,59,4);
        check(ledger.active_through()==62, "V45 commit window end wrong");
        for (std::size_t n=60;n<=62;++n) {
            t=ledger.trial(n,false,0.5,4,1.0);
            check(t.active && !t.triggered &&
                  t.remaining_including_this==63-n,
                  "V45 continuation shortened incorrectly");
            check(coupled_outer_blend_v43(2.0,6.0,t.beta)==4.0,
                  "V45 carrier interpolation wrong");
            ledger.commit(t,n,4);
            check(ledger.active_through()==62,"ordinary commit reset V45 window");
        }
        t=ledger.trial(63,false,0.5,4,1.0);
        check(!t.active && t.beta==1.0,"V45 damping persisted beyond bound");
        t=ledger.trial(71,true,0.75,3,0.5);
        check(t.active && t.beta==0.375,"V45 new trigger bad");
        ledger.commit(t,71,3);
        check(ledger.active_through()==73,"V45 re-arm failed");
        bool threw=false;
        try { (void)ledger.trial(74,false,0.5,129,1.0); }
        catch(const std::invalid_argument&){threw=true;}
        check(threw,"V45 accepts unbounded steps");
        threw=false;
        try {(void)ledger.trial(74,false,0.5,3,
                               std::numeric_limits<double>::quiet_NaN());}
        catch(const std::invalid_argument&){threw=true;}
        check(threw,"V45 accepts nonfinite retry");
        threw=false;
        try {(void)ledger.trial(75,true,-0.2,4,1.0);}
        catch(const std::invalid_argument&){threw=true;}
        check(threw,"V45 accepts negative beta");
        threw=false;
        try {(void)ledger.trial(std::numeric_limits<std::size_t>::max(),
                               true,0.5,4,1.0);}
        catch(const std::invalid_argument&){threw=true;}
        check(threw,"V45 accepts overflowing iteration window");
        std::cout<<"COUPLED_PHASE_ESCAPE_V45_PASS rollback=1 bounded=1 "
                    "damping=1 no_default=1\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"COUPLED_PHASE_ESCAPE_V45_FAIL "<<e.what()<<'\n';
        return 1;
    }
}
