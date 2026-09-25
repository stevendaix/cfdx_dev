#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/turbulence_transport.h"
#include "cfdx/physics/spalart_allmaras.h"
#include "cfdx/core/numerics/gradient.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::physics {

struct TurbulenceTransportResult {
    bool converged=false;
    std::size_t iterations=0;
    double k_residual=0.0;
    double second_residual=0.0;
};

inline void compute_kepsilon_production_field(
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& epsilon,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& production,
    const TurbulenceTransportControls& c)
{
    for(std::size_t i=0;i<mesh.n_cells();++i) {
        const double nut=c.C_mu*std::max(k(i),c.k_min)*std::max(k(i),c.k_min)/
                         std::max(epsilon(i),c.epsilon_min);
        production(i)=c.density*nut*strain_rate(i)*strain_rate(i);
    }
}

inline TurbulenceTransportResult solve_kepsilon_transport(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& epsilon,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const TurbulenceTransportControls& controls,
    const ScalarBoundaryConditions& k_bcs = {},
    const ScalarBoundaryConditions& epsilon_bcs = {},
    std::size_t max_iterations = 100,
    double tolerance = 1e-8)
{
    validate_turbulence_controls(controls);
    if(strain_rate.size()!=mesh.n_cells() || k.size()!=mesh.n_cells() ||
       epsilon.size()!=mesh.n_cells() || mass_flux.size()!=mesh.n_faces())
        throw std::invalid_argument("k-epsilon field size mismatch");
    if(controls.model!=TurbulenceModel::KEPSILON)
        throw std::invalid_argument("solve_kepsilon_transport requires KEPSILON model");

    const std::size_t n=mesh.n_cells();
    cfdx::core::Field<double,cfdx::core::Location::CELL> P(
        n,"Pk","W/m3",1);
    TurbulenceTransportResult result;

    for(std::size_t iter=1;iter<=max_iterations;++iter) {
        auto oldk=k;
        auto olde=epsilon;
        compute_kepsilon_production_field(mesh,strain_rate,k,epsilon,P,controls);

        cfdx::core::Field<double,cfdx::core::Location::CELL> sk(
            n,"Sk","W/m3",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> se(
            n,"Se","W/m3",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> spk(
            n,"Spk","kg/m3/s",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> spe(
            n,"Spe","kg/m3/s",1);

        for(std::size_t i=0;i<n;++i) {
            const double ki=std::max(k(i),controls.k_min);
            const double ei=std::max(epsilon(i),controls.epsilon_min);
            sk(i)=P(i);
            spk(i)=-controls.density*ei/ki;
            se(i)=controls.C1*P(i)*ei/ki;
            spe(i)=-controls.density*controls.C2*ei/ki;
        }

        std::vector<double> gamma_k(n), gamma_e(n);
        for(std::size_t i=0;i<n;++i) {
            const double ki=std::max(k(i),controls.k_min);
            const double ei=std::max(epsilon(i),controls.epsilon_min);
            const double nut=controls.C_mu*ki*ki/ei;
            gamma_k[i]=controls.density*
                (controls.molecular_viscosity+nut/controls.sigma_k);
            gamma_e[i]=controls.density*
                (controls.molecular_viscosity+nut/controls.sigma_epsilon);
        }

        auto eqk=assemble_scalar_equation(
            mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,
            nullptr,nullptr,nullptr,&gamma_k);
        auto eqe=assemble_scalar_equation(
            mesh,geometry,mass_flux,0.0,se,spe,epsilon_bcs,true,
            nullptr,nullptr,nullptr,&gamma_e);

        ScalarSolveControls sc;
        sc.max_iterations=2000;
        sc.tolerance=tolerance;
        sc.relaxation=0.7;
        cfdx::core::Vector k_solution(n,0.0);
        cfdx::core::Vector epsilon_solution(n,0.0);
        for(std::size_t i=0;i<n;++i) {
            k_solution(i)=k(i);
            epsilon_solution(i)=epsilon(i);
        }
        auto rk=solve_scalar_equation(eqk,k_solution,sc);
        auto re=solve_scalar_equation(eqe,epsilon_solution,sc);
        for(std::size_t i=0;i<n;++i) {
            k(i)=k_solution(i);
            epsilon(i)=epsilon_solution(i);
        }
        enforce_turbulence_bounds(k,epsilon,controls);

        result.k_residual=scalar_equation_residual_inf(eqk,k_solution);
        result.second_residual=scalar_equation_residual_inf(eqe,epsilon_solution);
        result.iterations=iter;

        double dk=0.0,de=0.0;
        for(std::size_t i=0;i<n;++i) {
            dk=std::max(dk,std::abs(k(i)-oldk(i)));
            de=std::max(de,std::abs(epsilon(i)-olde(i)));
        }
        if(rk.status==cfdx::core::SolverStatus::CONVERGED &&
           re.status==cfdx::core::SolverStatus::CONVERGED &&
           std::max(dk,de)<=tolerance) {
            result.converged=true;
            break;
        }
    }
    return result;
}


inline TurbulenceTransportResult solve_rng_kepsilon_transport(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& epsilon,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const TurbulenceTransportControls& controls,
    const ScalarBoundaryConditions& k_bcs = {},
    const ScalarBoundaryConditions& epsilon_bcs = {},
    std::size_t max_iterations = 100, double tolerance = 1e-8)
{
    validate_turbulence_controls(controls);
    if (controls.model != TurbulenceModel::RNG_KEPSILON)
        throw std::invalid_argument("solve_rng_kepsilon_transport requires RNG_KEPSILON model");
    const std::size_t n = mesh.n_cells();
    if (k.size()!=n || epsilon.size()!=n || strain_rate.size()!=n ||
        mass_flux.size()!=mesh.n_faces())
        throw std::invalid_argument("RNG k-epsilon field size mismatch");

    TurbulenceTransportResult result;
    for (std::size_t iter=1; iter<=max_iterations; ++iter) {
        auto oldk=k, olde=epsilon;
        cfdx::core::Field<double,cfdx::core::Location::CELL> sk(n,"Sk","W/m3",1), se(n,"Se","W/m3",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> spk(n,"Spk","kg/m3/s",1), spe(n,"Spe","kg/m3/s",1);
        std::vector<double> gamma_k(n), gamma_e(n);
        for (std::size_t i=0;i<n;++i) {
            const double ki=std::max(k(i),controls.k_min), ei=std::max(epsilon(i),controls.epsilon_min);
            const double nut=controls.rng_C_mu*ki*ki/ei;
            const double S=std::max(strain_rate(i),0.0);
            const double P=controls.density*nut*S*S;
            const double eta=S*ki/ei;
            const double C1star=controls.rng_C1-
                eta*(1.0-eta/controls.rng_eta0)/(1.0+controls.rng_beta*eta*eta*eta);
            sk(i)=P;
            spk(i)=-controls.density*ei/ki;
            se(i)=C1star*P*ei/ki;
            spe(i)=-controls.density*controls.rng_C2*ei/ki;
            gamma_k[i]=controls.density*(controls.molecular_viscosity+nut/controls.rng_sigma_k);
            gamma_e[i]=controls.density*(controls.molecular_viscosity+nut/controls.rng_sigma_epsilon);
        }
        auto eqk=assemble_scalar_equation(mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,nullptr,nullptr,&gamma_k);
        auto eqe=assemble_scalar_equation(mesh,geometry,mass_flux,0.0,se,spe,epsilon_bcs,true,nullptr,nullptr,&gamma_e);
        ScalarSolveControls sc{2000,tolerance,0.7};
        cfdx::core::Vector ks(n,0.0), es(n,0.0);
        for(std::size_t i=0;i<n;++i){ks(i)=k(i);es(i)=epsilon(i);}
        const auto rk=solve_scalar_equation(eqk,ks,sc);
        const auto re=solve_scalar_equation(eqe,es,sc);
        for(std::size_t i=0;i<n;++i){k(i)=ks(i);epsilon(i)=es(i);}
        enforce_turbulence_bounds(k,epsilon,controls);
        double dk=0.0,de=0.0;
        for(std::size_t i=0;i<n;++i){dk=std::max(dk,std::abs(k(i)-oldk(i)));de=std::max(de,std::abs(epsilon(i)-olde(i)));}
        result.k_residual=scalar_equation_residual_inf(eqk,ks);
        result.second_residual=scalar_equation_residual_inf(eqe,es);
        result.iterations=iter;
        if(rk.status==cfdx::core::SolverStatus::CONVERGED &&
           re.status==cfdx::core::SolverStatus::CONVERGED &&
           std::max(dk,de)<=tolerance){result.converged=true;break;}
    }
    return result;
}


struct KOmegaCellSources {
    double nut=0.0, sk=0.0, spk=0.0, sw=0.0, spw=0.0, gamma_k=0.0, gamma_w=0.0;
};

inline KOmegaCellSources komega2006_cell_sources(
    double k, double omega, double strain, double grad_k_dot_grad_omega,
    const TurbulenceTransportControls& c)
{
    if(!std::isfinite(k)||!std::isfinite(omega)||!std::isfinite(strain)||
       !std::isfinite(grad_k_dot_grad_omega))
        throw std::invalid_argument("komega2006_cell_sources: non-finite input");
    const double ki=std::max(k,c.k_min), wi=std::max(omega,c.omega_min);
    const double S=std::max(strain,0.0);
    const double omega_t=std::max(wi,c.komega_clim*S/std::sqrt(c.beta_star));
    KOmegaCellSources q;
    q.nut=ki/omega_t;
    const double Pk=c.density*q.nut*S*S;
    q.sk=Pk; q.spk=-c.density*c.beta_star*wi;
    const double sigma_d=grad_k_dot_grad_omega>0.0?c.komega_sigma_d0:0.0;
    q.sw=c.komega_alpha*(wi/ki)*Pk+c.density*sigma_d/wi*grad_k_dot_grad_omega;
    q.spw=-c.density*c.komega_beta0*wi;
    q.gamma_k=c.density*(c.molecular_viscosity+c.komega_sigma_k*ki/wi);
    q.gamma_w=c.density*(c.molecular_viscosity+c.komega_sigma_w*ki/wi);
    return q;
}

inline std::vector<double> cell_gradient_dot(
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& a,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& b,
    const cfdx::core::Mesh& mesh)
{
    const auto ga=cfdx::core::compute_gradient_gauss(a,mesh);
    const auto gb=cfdx::core::compute_gradient_gauss(b,mesh);
    std::vector<double> dot(mesh.n_cells());
    for(std::size_t i=0;i<mesh.n_cells();++i)
        dot[i]=ga.component_data(0)[i]*gb.component_data(0)[i]+
               ga.component_data(1)[i]*gb.component_data(1)[i]+
               ga.component_data(2)[i]*gb.component_data(2)[i];
    return dot;
}

inline TurbulenceTransportResult solve_komega_transport(
    const cfdx::core::Mesh& mesh,const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const TurbulenceTransportControls& controls,const ScalarBoundaryConditions& k_bcs={},
    const ScalarBoundaryConditions& omega_bcs={},std::size_t max_iterations=100,double tolerance=1e-8)
{
    validate_turbulence_controls(controls);
    if(controls.model!=TurbulenceModel::KOMEGA) throw std::invalid_argument("solve_komega_transport requires KOMEGA model");
    const std::size_t n=mesh.n_cells();
    if(k.size()!=n||omega.size()!=n||strain_rate.size()!=n||mass_flux.size()!=mesh.n_faces())
        throw std::invalid_argument("k-omega field size mismatch");
    TurbulenceTransportResult result;
    for(std::size_t iter=1;iter<=max_iterations;++iter){
        auto old_k=k,old_w=omega;
        cfdx::core::Field<double,cfdx::core::Location::CELL> sk(n,"Sk","W/m3",1),sw(n,"Sw","W/m3",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> spk(n,"Spk","kg/m3/s",1),spw(n,"Spw","kg/m3/s",1);
        std::vector<double> gk(n),gw(n);
        const auto cross=cell_gradient_dot(k,omega,mesh);
        for(std::size_t i=0;i<n;++i){
            const auto q=komega2006_cell_sources(k(i),omega(i),strain_rate(i),cross[i],controls);
            sk(i)=q.sk; spk(i)=q.spk; sw(i)=q.sw; spw(i)=q.spw;
            gk[i]=q.gamma_k; gw[i]=q.gamma_w;
        }
        auto eqk=assemble_scalar_equation(mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,nullptr,nullptr,nullptr,&gk);
        auto eqw=assemble_scalar_equation(mesh,geometry,mass_flux,0.0,sw,spw,omega_bcs,true,nullptr,nullptr,nullptr,&gw);
        ScalarSolveControls sc{2000,tolerance,0.7}; cfdx::core::Vector ks(n,0.0),ws(n,0.0);
        for(std::size_t i=0;i<n;++i){ks(i)=k(i);ws(i)=omega(i);}
        const auto rk=solve_scalar_equation(eqk,ks,sc), rw=solve_scalar_equation(eqw,ws,sc);
        for(std::size_t i=0;i<n;++i){k(i)=ks(i);omega(i)=ws(i);}
        enforce_turbulence_bounds(k,omega,controls);
        double dk=0,dw=0; for(std::size_t i=0;i<n;++i){dk=std::max(dk,std::abs(k(i)-old_k(i)));dw=std::max(dw,std::abs(omega(i)-old_w(i)));}
        result.k_residual=scalar_equation_residual_inf(eqk,ks); result.second_residual=scalar_equation_residual_inf(eqw,ws); result.iterations=iter;
        if(rk.status==cfdx::core::SolverStatus::CONVERGED&&rw.status==cfdx::core::SolverStatus::CONVERGED&&std::max(dk,dw)<=tolerance){result.converged=true;break;}
    }
    return result;
}

inline double sa_modified_vorticity(
    double vorticity, double nu_tilde, double chi, double wall_distance,
    const SpalartAllmarasModel& sa)
{
    if(!(wall_distance>0.0) || !std::isfinite(vorticity) || vorticity<0.0 ||
       !std::isfinite(nu_tilde) || nu_tilde<0.0)
        throw std::invalid_argument("sa_modified_vorticity: invalid inputs");
    constexpr double cv2=0.7, cv3=0.9;
    const double sbar=nu_tilde*sa.fv2(chi)/(sa.kappa*sa.kappa*wall_distance*wall_distance);
    if(sbar>=-cv2*vorticity) return vorticity+sbar;
    return vorticity+vorticity*(cv2*cv2*vorticity+cv3*sbar)/
                     ((cv3-2.0*cv2)*vorticity-sbar);
}

inline TurbulenceTransportResult solve_spalart_allmaras_transport(
    const cfdx::core::Mesh& mesh,const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& nu_tilde,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& vorticity,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& wall_distance,
    const TurbulenceTransportControls& controls,const ScalarBoundaryConditions& bcs={},
    std::size_t max_iterations=100,double tolerance=1e-8)
{
    validate_turbulence_controls(controls);
    if(controls.model!=TurbulenceModel::SPALART_ALLMARAS) throw std::invalid_argument("solve_spalart_allmaras_transport requires SPALART_ALLMARAS model");
    const std::size_t n=mesh.n_cells();
    if(nu_tilde.size()!=n||vorticity.size()!=n||wall_distance.size()!=n||mass_flux.size()!=mesh.n_faces())
        throw std::invalid_argument("Spalart-Allmaras field size mismatch");
    SpalartAllmarasModel sa; sa.cb1=controls.sa_cb1;sa.cb2=controls.sa_cb2;sa.sigma=controls.sa_sigma;sa.kappa=controls.sa_kappa;sa.cw2=controls.sa_cw2;sa.cw3=controls.sa_cw3;sa.cv1=controls.sa_cv1;sa.ct3=controls.sa_ct3;sa.ct4=controls.sa_ct4;
    sa.cw1=sa.cb1/(sa.kappa*sa.kappa)+(1.0+sa.cb2)/sa.sigma;
    TurbulenceTransportResult result;
    for(std::size_t iter=1;iter<=max_iterations;++iter){
        auto old=nu_tilde;
        cfdx::core::Field<double,cfdx::core::Location::CELL> su(n,"Snu","m2/s3",1),spu(n,"Spnu","1/s",1);
        std::vector<double> gamma(n);
        const auto grad_nu_tilde = cfdx::core::compute_gradient_gauss(nu_tilde,mesh);
        const double* gx = grad_nu_tilde.component_data(0);
        const double* gy = grad_nu_tilde.component_data(1);
        const double* gz = grad_nu_tilde.component_data(2);
        for(std::size_t i=0;i<n;++i){
            const double wt=std::max(nu_tilde(i),0.0), nu=controls.molecular_viscosity, d=wall_distance(i), vort=std::max(vorticity(i),1e-20);
            if(!(d>0)||!std::isfinite(d)||!std::isfinite(vort)) throw std::invalid_argument("invalid SA wall/vorticity field");
            const double chi=wt/nu;
            const double st=std::max(sa_modified_vorticity(vort,wt,chi,d,sa),1e-20);
            const double r=std::min(wt/(st*sa.kappa*sa.kappa*d*d),10.0);
            const double fw=sa.destruction_coefficient(r), ft2=sa.ft2(chi);
            const double prod=sa.cb1*(1-ft2)*st;
            const double destr=std::max(sa.cw1*fw-sa.cb1*ft2/(sa.kappa*sa.kappa),0.0)*wt/(d*d);
            const double grad2=gx[i]*gx[i]+gy[i]*gy[i]+gz[i]*gz[i];
            su(i)=controls.density*(prod*wt + sa.cb2/sa.sigma*grad2); spu(i)=-controls.density*destr;
            gamma[i]=controls.density*(nu+wt)/sa.sigma;
        }
        auto eq=assemble_scalar_equation(mesh,geometry,mass_flux,0.0,su,spu,bcs,true,nullptr,nullptr,nullptr,&gamma);
        ScalarSolveControls sc{2000,tolerance,0.7}; cfdx::core::Vector sol(n,0.0); for(std::size_t i=0;i<n;++i)sol(i)=nu_tilde(i);
        const auto rr=solve_scalar_equation(eq,sol,sc); for(std::size_t i=0;i<n;++i)nu_tilde(i)=std::max(sol(i),0.0);
        double du=0;for(std::size_t i=0;i<n;++i)du=std::max(du,std::abs(nu_tilde(i)-old(i)));
        result.k_residual=scalar_equation_residual_inf(eq,sol);result.second_residual=result.k_residual;result.iterations=iter;
        if(rr.status==cfdx::core::SolverStatus::CONVERGED&&du<=tolerance){result.converged=true;break;}
    }
    return result;
}

inline ScalarBoundaryCondition wall_k_fixed_value(double k_value)
{
    if(!std::isfinite(k_value) || k_value<0.0)
        throw std::invalid_argument("invalid wall k value");
    return {ScalarBoundaryType::FIXED_VALUE,k_value,0.0};
}

inline ScalarBoundaryCondition wall_k_epsilon(double k_value, double epsilon_value)
{
    if(epsilon_value<=0.0) throw std::invalid_argument("invalid wall turbulence values");
    return wall_k_fixed_value(k_value);
}

inline double wall_epsilon_from_k(double k, double y, double Cmu=0.09, double kappa=0.41)
{
    if(k<0.0 || y<=0.0 || Cmu<=0.0 || kappa<=0.0)
        throw std::invalid_argument("invalid wall epsilon input");
    return std::pow(Cmu,0.75)*std::pow(k,1.5)/(kappa*y);
}

inline double wall_omega_from_k(double k, double y, double betaStar=0.09, double kappa=0.41)
{
    if(k<0.0 || y<=0.0 || betaStar<=0.0 || kappa<=0.0)
        throw std::invalid_argument("invalid wall omega input");
    return std::sqrt(k)/(std::pow(betaStar,0.25)*kappa*y);
}

} // namespace cfdx::physics
