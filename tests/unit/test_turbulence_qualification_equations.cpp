#include "cfdx/physics/turbulence_models.h"
#include "cfdx/physics/turbulence_transport.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void expect_close(double a,double b,double tol=1e-12) {
    if(std::abs(a-b)>tol*std::max({1.0,std::abs(a),std::abs(b)}))
        throw std::runtime_error("turbulence equation reference mismatch");
}
}

int main() {
    using namespace cfdx::physics;
    TurbulenceTransportControls c;
    const double k=2.0, eps=4.0, omega=5.0, S=3.0, delta=0.2, volume=delta*delta*delta, y=0.05;

    c.model=TurbulenceModel::LAMINAR;
    expect_close(turbulence_nu_t(k,eps,S,y,c),0.0);

    c.model=TurbulenceModel::KEPSILON;
    expect_close(turbulence_nu_t(k,eps,S,y,c),0.09);

    c.model=TurbulenceModel::RNG_KEPSILON;
    expect_close(turbulence_nu_t(k,eps,S,y,c),0.0845);

    c.model=TurbulenceModel::REALIZABLE_KEPSILON;
    expect_close(turbulence_nu_t(k,eps,S,y,c,volume,1.0,0.0,0.0,0.0),0.25);

    c.model=TurbulenceModel::KOMEGA;
    expect_close(turbulence_nu_t(k,omega,S,y,c),k/omega);

    c.model=TurbulenceModel::SST;
    expect_close(turbulence_nu_t(k,omega,S,y,c,volume,1.0),0.31*k/std::max(0.31*omega,S));

    c.model=TurbulenceModel::SPALART_ALLMARAS;
    expect_close(turbulence_nu_t(0.0,0.01,S,y,c),0.01*std::pow(10.0,3.0)/(std::pow(10.0,3.0)+std::pow(c.sa_cv1,3.0)),1e-10);

    c.model=TurbulenceModel::SMAGORINSKY;
    expect_close(turbulence_nu_t(k,eps,S,y,c,volume),(c.smagorinsky_Cs*delta)*(c.smagorinsky_Cs*delta)*S);

    c.model=TurbulenceModel::WALE;
    const double S2=2.0, Sd2=0.5;
    const double expected_wale=std::pow(c.wale_Cw*delta,2.0)*std::pow(Sd2,1.5)/
        (std::pow(S2,2.5)+std::pow(Sd2,1.25));
    expect_close(turbulence_nu_t(k,eps,S,y,c,volume,1.0,0.0,0.0,0.0,1.0,S2,Sd2),expected_wale);

    c.model=TurbulenceModel::DYNAMIC_KEQN;
    expect_close(turbulence_nu_t(k,0.25,S,y,c,volume),0.1*delta*std::sqrt(0.25));

    c.model=TurbulenceModel::DES;
    expect_close(turbulence_nu_t(k,eps,S,y,c,volume),0.17*0.17*y*y*S);

    c.model=TurbulenceModel::DDES;
    const double rd=0.5;
    const double fd=1.0-std::tanh(std::pow(8.0*rd,3.0));
    const double lddes=y-fd*std::max(0.0,y-c.des_Cdes*delta);
    expect_close(turbulence_nu_t(k,eps,S,y,c,volume,1.0,0.0,0.0,rd),0.17*0.17*lddes*lddes*S);

    c.model=TurbulenceModel::IDDES;
    expect_close(turbulence_nu_t(k,eps,S,y,c,volume,1.0,0.0,0.0,rd,0.25),
                 0.17*0.17*lddes*lddes*S);

    std::cout << "turbulence qualification equation references: PASS\n";
    return 0;
}
