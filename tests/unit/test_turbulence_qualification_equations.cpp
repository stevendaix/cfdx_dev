#include "cfdx/physics/turbulence_models.h"
#include "cfdx/physics/turbulence_transport.h"
#include <algorithm>
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
    {
        const double rho=1.0, nu=1e-3, Cmu=0.0845, C1=1.42, C2=1.68;
        const double sigk=0.71942, sige=0.71942, eta0=4.38, beta=0.012;
        const double nut=Cmu*k*k/eps, P=rho*nut*S*S, eta=S*k/eps;
        const double C1star=C1-eta*(1.0-eta/eta0)/(1.0+beta*eta*eta*eta);
        const auto q=rng_kepsilon_cell_sources(k,eps,S,c);
        expect_close(q.nut,nut); expect_close(q.sk,P);
        expect_close(q.spk,-rho*eps/k);
        expect_close(q.se,C1star*P*eps/k);
        expect_close(q.spe,-rho*C2*eps/k);
        expect_close(q.gamma_k,rho*(nu+nut/sigk));
        expect_close(q.gamma_epsilon,rho*(nu+nut/sige));
    }

    c.model=TurbulenceModel::REALIZABLE_KEPSILON;
    {
        const RealizableKEpsilonInvariants inv{S,1.2,0.15};
        const auto q=realizable_kepsilon_cell_sources(k,eps,inv,c);
        const double W=std::clamp(inv.third_invariant,-1.0/std::sqrt(6.0),1.0/std::sqrt(6.0));
        const double phi=std::acos(std::sqrt(6.0)*W)/3.0;
        const double As=std::sqrt(6.0)*std::cos(phi);
        const double Ustar=std::sqrt(S*S+inv.rotation_magnitude*inv.rotation_magnitude);
        const double cmu=1.0/(c.realizable_A0+As*Ustar*k/eps);
        const double nut=cmu*k*k/eps, P=c.density*nut*S*S;
        const double c1=std::max(0.43,(S*k/eps)/((S*k/eps)+5.0));
        expect_close(q.cmu,cmu); expect_close(q.nut,nut);
        expect_close(q.sk,P); expect_close(q.spk,-c.density*eps/k);
        expect_close(q.se,c1*c.density*S*eps);
        expect_close(q.spe,-c.density*c.realizable_C2*eps/(k+std::sqrt(c.molecular_viscosity*eps)));
        expect_close(q.gamma_k,c.density*(c.molecular_viscosity+nut/c.realizable_sigma_k));
        expect_close(q.gamma_epsilon,c.density*(c.molecular_viscosity+nut/c.realizable_sigma_epsilon));
    }

    c.model=TurbulenceModel::KOMEGA;
    {
        const double grad_cross=-0.7;
        const auto q=komega2006_cell_sources(k,omega,S,grad_cross,c);
        const double nut=k/omega, P=c.density*nut*S*S;
        expect_close(q.nut,nut); expect_close(q.sk,P);
        expect_close(q.spk,-c.density*c.beta_star*omega);
        expect_close(q.sw,c.komega_alpha*(omega/k)*P);
        expect_close(q.spw,-c.density*c.komega_beta0*omega);
        expect_close(q.gamma_k,c.density*(c.molecular_viscosity+c.komega_sigma_k*nut));
        expect_close(q.gamma_w,c.density*(c.molecular_viscosity+c.komega_sigma_w*nut));
        const auto qp=komega2006_cell_sources(k,omega,S,0.7,c);
        expect_close(qp.sw,c.komega_alpha*(omega/k)*P+c.density*c.komega_sigma_d0/omega*0.7);
    }

    c.model=TurbulenceModel::SST;
    {
        const double F1=0.3,F2=0.8,cross=-0.4;
        const auto q=sst_cell_sources(k,omega,S,F1,F2,c);
        const double sk=std::sqrt(k)/(c.beta_star*omega*y);
        const double visc=500.0*c.molecular_viscosity/(y*y*omega);
        const double cdkw=std::max(2.0*c.density*c.sst_sigma_w2/omega*cross,1e-10);
        const double arg1=std::min(std::max(sk,visc),4.0*c.density*c.sst_sigma_w2*k/(cdkw*y*y));
        const double arg2=std::max(2.0*std::sqrt(k)/(c.beta_star*omega*y),visc);
        const auto blend_ref=compute_sst_blending(k,omega,y,c.molecular_viscosity,c.beta_star,
                                                   cross,c.density,c.sst_sigma_w2);
        expect_close(blend_ref.first,std::tanh(std::pow(arg1,4.0)));
        expect_close(blend_ref.second,std::tanh(arg2*arg2));
        const auto blend=[F1](double a,double b){return F1*a+(1.0-F1)*b;};
        const double nut=c.a1*k/std::max(c.a1*omega,S*F2);
        const double P=std::min(c.density*nut*S*S,c.sst_production_limiter*c.beta_star*c.density*k*omega);
        const double cross_term=2.0*(1.0-F1)*c.density*c.sst_sigma_w2/omega*cross;
        expect_close(q.nut,nut); expect_close(q.sk,P);
        expect_close(q.spk,-c.density*c.beta_star*omega);
        expect_close(q.sw,blend(c.gamma1,c.gamma2)*c.density*S*S+std::max(cross_term,0.0));
        expect_close(q.spw,-c.density*blend(c.beta1,c.beta2)*omega+std::min(cross_term,0.0)/omega);
        expect_close(q.gamma_k,c.density*(c.molecular_viscosity+blend(c.sst_sigma_k1,c.sst_sigma_k2)*nut));
        expect_close(q.gamma_w,c.density*(c.molecular_viscosity+blend(c.sst_sigma_w1,c.sst_sigma_w2)*nut));
    }

    c.model=TurbulenceModel::SPALART_ALLMARAS;
    {
        SpalartAllmarasModel sa;
        const double rho=c.density, nu=c.molecular_viscosity, wt=0.01, vort=2.0, d=y;
        const double chi=wt/nu;
        const double chi3=chi*chi*chi;
        const double fv1=chi3/(chi3+std::pow(sa.cv1,3.0));
        const double fv2=1.0-chi/(1.0+chi*fv1);
        const double sbar=wt*fv2/(sa.kappa*sa.kappa*d*d);
        constexpr double cv2=0.7, cv3=0.9;
        const double st=sbar>=-cv2*vort ? vort+sbar :
            vort+vort*(cv2*cv2*vort+cv3*sbar)/((cv3-2.0*cv2)*vort-sbar);
        const double r=std::min(wt/(st*sa.kappa*sa.kappa*d*d),10.0);
        const double g=r+sa.cw2*(std::pow(r,6.0)-r);
        const double fw=g*std::pow((1.0+std::pow(sa.cw3,6.0))/(std::pow(g,6.0)+std::pow(sa.cw3,6.0)),1.0/6.0);
        const double ft2=sa.ct3*std::exp(-sa.ct4*chi*chi);
        const double prod=sa.cb1*(1.0-ft2)*st;
        const double destr=(sa.cw1*fw-sa.cb1*ft2/(sa.kappa*sa.kappa))*wt/(d*d);
        const double chi3_ref=std::pow(chi,3.0);
        const double expected_nut=wt*chi3_ref/(chi3_ref+std::pow(sa.cv1,3.0));
        expect_close(turbulence_nu_t(0.0,wt,S,d,c),expected_nut,1e-10);
        const double grad2=0.13;
        const auto src=sa_cell_sources(wt,vort,d,grad2,c,sa);
        expect_close(src.production,prod*wt);
        expect_close(src.destruction,destr);
        expect_close(src.gradient_source,sa.cb2/sa.sigma*grad2);
        expect_close(src.gamma,rho*(nu+wt)/sa.sigma);
    }

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
    const double rd=0.05;
    const double fd=1.0-std::tanh(std::pow(8.0*rd,3.0));
    const double lddes=y-fd*std::max(0.0,y-c.des_Cdes*delta);
    expect_close(turbulence_nu_t(k,eps,S,y,c,volume,1.0,0.0,0.0,rd),0.17*0.17*lddes*lddes*S);

    c.model=TurbulenceModel::IDDES;
    const double iddes_shield=(1.0-0.25)+0.25*fd;
    const double liddes=y-iddes_shield*std::max(0.0,y-c.des_Cdes*delta);
    expect_close(turbulence_nu_t(k,eps,S,y,c,volume,1.0,0.0,0.0,rd,0.25),
                 0.17*0.17*liddes*liddes*S);

    const double dynamic_Ck=0.13;
    const double dynamic_expected=dynamic_Ck*delta*std::sqrt(0.25);
    expect_close(dynamic_one_equation_eddy_viscosity(0.25,delta,dynamic_Ck),dynamic_expected);
    expect_close(dynamic_les_coefficient(0.5,0.2,2.0),0.15);

    expect_close(ddes_shielding(0.0),1.0);
    expect_close(ddes_shielding(1.0),1.0-std::tanh(std::pow(8.0,3.0)));
    expect_close(iddes_shielding(0.0,0.0),1.0);
    expect_close(iddes_shielding(0.0,1.0),1.0);
    expect_close(iddes_shielding(1.0,1.0),ddes_shielding(1.0));
    expect_close(dynamic_les_coefficient(2.0,1.0,2.0),0.5);
    expect_close(dynamic_les_coefficient(0.0,2.0,2.0),0.0);
    expect_close(dynamic_les_coefficient(3.0,0.0,2.0,0.0,0.23),0.23);

    bool rejected=false;
    try { (void)dynamic_les_coefficient(1.0,0.0,0.0); } catch(const std::invalid_argument&) { rejected=true; }
    if(!rejected) throw std::runtime_error("dynamic LES denominator guard failed");
    rejected=false;
    try { (void)dynamic_one_equation_eddy_viscosity(-1.0,delta); } catch(const std::invalid_argument&) { rejected=true; }
    if(!rejected) throw std::runtime_error("dynamic one-equation input guard failed");

    std::cout << "turbulence qualification equation references: PASS\n";
    return 0;
}
