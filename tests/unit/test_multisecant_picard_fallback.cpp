#include "cfdx/physics/nonlinear_multisecant_acceleration.h"
#include "cfdx/physics/timestep_control.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
static void demand(bool ok, const char* reason) {
    if (!ok) throw std::runtime_error(reason);
}
int main() {
    NonlinearMultisecantControls multi;
    multi.method=NonlinearMultisecantMethod::ANDERSON;
    NonlinearMultisecantHistory history(multi);
    history.append({{1.0},{0.5},{0.1}},{1.0});
    history.append({{0.5},{0.8},{0.2}},{1.0});
    Field<double,Location::CELL> U(1,"U","m/s",3),p(1,"p","Pa",1);
    U.fill(0.0);p.fill(0.0);
    U(0,0)=0.2;U(0,1)=-0.1;p(0)=0.04;
    std::vector<double> phi{0.3,-0.3};
    const std::vector<double> accepted_flux=phi;
    const double accepted_u=U(0,0), accepted_v=U(0,1),accepted_p=p(0);
    NonlinearRetryController retry;
    NonlinearStateRollback state(U,p);
    state.begin();
    U(0,0)=8.0;U(0,1)=4.0;p(0)=11.0;
    phi={19.0,-18.0};
    // Same rollback path invoked for a rejected Anderson candidate.
    state.reject();phi=accepted_flux;history.reset();
    retry.retry_without_relaxation();
    demand(retry.retries()==1,"fallback retry is not counted");
    demand(retry.alpha_u(0.6)==0.6 && retry.alpha_p(0.25)==0.25,
           "fallback incorrectly changed momentum/pressure relaxation");
    demand(U(0,0)==accepted_u && U(0,1)==accepted_v && p(0)==accepted_p &&
           phi==accepted_flux,"fallback did not restore U,p,phi");
    demand(history.size()==0,"rejected multisecant history was retained");
    // A real Picard failure after fallback may now shrink alpha normally.
    state.begin();U(0,0)=7.0;p(0)=4.0;phi={2.0,2.0};
    state.reject();phi=accepted_flux;
    retry.reject();
    demand(retry.alpha_u(0.6)==0.3 && retry.alpha_p(0.25)==0.125,
           "ordinary Picard retry no longer shrinks relaxation");
    std::cout<<"MULTISECANT_PICARD_FALLBACK_PASS\n";
}
