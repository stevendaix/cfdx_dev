#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {
enum class ExtrapolationPolicy { CLAMP, LINEAR, REJECT };
struct TabulatedProperty {
    std::vector<double> temperature, value;
    ExtrapolationPolicy extrapolation=ExtrapolationPolicy::CLAMP;
    void validate(const char* name) const {
        if(temperature.size()<2 || temperature.size()!=value.size()) throw std::invalid_argument(std::string(name)+": at least two pairs required");
        for(std::size_t i=0;i<temperature.size();++i) {
            if(!std::isfinite(temperature[i]) || !std::isfinite(value[i])) throw std::invalid_argument(std::string(name)+": non-finite table value");
            if(i && !(temperature[i]>temperature[i-1])) throw std::invalid_argument(std::string(name)+": temperatures must increase");
        }
    }
    double evaluate(double T,const char* name="table") const {
        validate(name); if(!std::isfinite(T)) throw std::invalid_argument(std::string(name)+": non-finite temperature");
        if(T<temperature.front() || T>temperature.back()) {
            if(extrapolation==ExtrapolationPolicy::REJECT) throw std::out_of_range(std::string(name)+": temperature outside range");
            if(extrapolation==ExtrapolationPolicy::CLAMP) return T<temperature.front()?value.front():value.back();
            const std::size_t i=T<temperature.front()?0:temperature.size()-2;
            const double f=(T-temperature[i])/(temperature[i+1]-temperature[i]);
            return value[i]+f*(value[i+1]-value[i]);
        }
        const auto it=std::upper_bound(temperature.begin(),temperature.end(),T);
        const std::size_t i=it==temperature.begin()?0:static_cast<std::size_t>(it-temperature.begin()-1);
        if(i>=temperature.size()-1) return value.back();
        const double f=(T-temperature[i])/(temperature[i+1]-temperature[i]);
        return value[i]+f*(value[i+1]-value[i]);
    }
};
enum class ScalarPropertyModel { CONSTANT, LINEAR, POLYNOMIAL, POWER_LAW, EXPONENTIAL, ARRHENIUS, TABLE };
struct ScalarPropertyControls {
    ScalarPropertyModel model=ScalarPropertyModel::CONSTANT;
    double reference_temperature=300.0, reference_value=1.0;
    std::vector<double> coefficients;
    double exponent=0.0, activation_temperature=0.0, minimum=0.0, maximum=0.0;
    bool enforce_bounds=false;
    TabulatedProperty table;
};
inline double evaluate_scalar_property(const ScalarPropertyControls& c,double T) {
    if(!std::isfinite(T)||T<=0||!std::isfinite(c.reference_temperature)||c.reference_temperature<=0||
       !std::isfinite(c.reference_value)||c.reference_value<=0) throw std::invalid_argument("scalar property: invalid reference");
    const double x=T-c.reference_temperature; double v=c.reference_value;
    switch(c.model) {
    case ScalarPropertyModel::CONSTANT: break;
    case ScalarPropertyModel::LINEAR: v+= (c.coefficients.empty()?0.0:c.coefficients[0])*x; break;
    case ScalarPropertyModel::POLYNOMIAL: {
        if(c.coefficients.empty()) throw std::invalid_argument("polynomial property requires coefficients");
        v=0; double p=1; for(double a:c.coefficients){ if(!std::isfinite(a)) throw std::invalid_argument("polynomial coefficient is non-finite"); v+=a*p; p*=x; } break;
    }
    case ScalarPropertyModel::POWER_LAW: v=c.reference_value*std::pow(T/c.reference_temperature,c.exponent); break;
    case ScalarPropertyModel::EXPONENTIAL: v=c.reference_value*std::exp(c.exponent*x/c.reference_temperature); break;
    case ScalarPropertyModel::ARRHENIUS: v=c.reference_value*std::exp(c.activation_temperature*(1/c.reference_temperature-1/T)); break;
    case ScalarPropertyModel::TABLE: v=c.table.evaluate(T,"scalar property table"); break;
    }
    if(!std::isfinite(v)||v<=0) throw std::domain_error("scalar property produced invalid value");
    if(c.enforce_bounds) {
        if(!std::isfinite(c.minimum)||!std::isfinite(c.maximum)||c.minimum<=0||c.maximum<c.minimum) throw std::invalid_argument("invalid property bounds");
        v=std::clamp(v,c.minimum,c.maximum);
    }
    return v;
}
enum class DensityModel { CONSTANT, LINEAR_TEMPERATURE, BOUSSINESQ, IDEAL_GAS, TABLE };
struct DensityControls {
    DensityModel model=DensityModel::CONSTANT;
    double rho0=1.0, reference_temperature=300.0, beta=0.0, pressure_reference=101325.0, gas_constant=287.05;
    double max_boussinesq_relative_change=0.1;
    bool enforce_boussinesq_limit=true;
    TabulatedProperty table;
};
inline double evaluate_density(const DensityControls& c,double T,double p=101325.0) {
    if(!std::isfinite(T)||T<=0||!std::isfinite(c.rho0)||c.rho0<=0) throw std::invalid_argument("density: invalid inputs");
    double rho=c.rho0;
    switch(c.model) {
    case DensityModel::CONSTANT: break;
    case DensityModel::LINEAR_TEMPERATURE:
        rho=c.rho0*(1-c.beta*(T-c.reference_temperature)); break;
    case DensityModel::BOUSSINESQ: {
        if(!std::isfinite(c.reference_temperature)||c.reference_temperature<=0||!std::isfinite(c.beta))
            throw std::invalid_argument("Boussinesq density: invalid reference temperature or expansion coefficient");
        const double relative_change=std::abs(c.beta*(T-c.reference_temperature));
        if(c.enforce_boussinesq_limit &&
           (!std::isfinite(c.max_boussinesq_relative_change) || c.max_boussinesq_relative_change<=0 ||
            relative_change>c.max_boussinesq_relative_change))
            throw std::domain_error("Boussinesq density: |beta*DeltaT| exceeds validity limit");
        rho=c.rho0*(1-c.beta*(T-c.reference_temperature)); break;
    }
    case DensityModel::IDEAL_GAS:
        if(!std::isfinite(p)||p<=0||!std::isfinite(c.gas_constant)||c.gas_constant<=0) throw std::invalid_argument("ideal gas density: invalid inputs");
        rho=p/(c.gas_constant*T); break;
    case DensityModel::TABLE: rho=c.table.evaluate(T,"density table"); break;
    }
    if(!std::isfinite(rho)||rho<=0) throw std::domain_error("density model produced invalid density");
    return rho;
}
struct ThermophysicalProperties {
    ScalarPropertyControls conductivity, heat_capacity, viscosity, turbulent_prandtl;
    DensityControls density;
    double k(double T) const{return evaluate_scalar_property(conductivity,T);}
    double cp(double T) const{return evaluate_scalar_property(heat_capacity,T);}
    double rho(double T,double p=101325) const{return evaluate_density(density,T,p);}
    double mu(double T) const{return evaluate_scalar_property(viscosity,T);}
    double pr_t(double T) const{return evaluate_scalar_property(turbulent_prandtl,T);}
    double pr_t(double T,double y_plus) const {
        if(!std::isfinite(y_plus)||y_plus<0) throw std::invalid_argument("turbulent Prandtl: invalid y+");
        return evaluate_scalar_property(turbulent_prandtl,T);
    }
    double enthalpy(double T,double Tref=300,double intervals=64) const {
        if(!std::isfinite(T)||!std::isfinite(Tref)||T<=0||Tref<=0||intervals<=0)
            throw std::invalid_argument("enthalpy: invalid inputs");
        if (!std::isfinite(intervals) || intervals<1.0 || intervals>1000000.0)
            throw std::invalid_argument("enthalpy: invalid quadrature interval count");
        const double a=std::min(T,Tref), b=std::max(T,Tref);
        (void)cp(a); (void)cp(b); // Validate property and table extrapolation.
        if (a==b) return 0.0;
        double integral=0.0;
        const auto model=heat_capacity.model;
        if (model==ScalarPropertyModel::CONSTANT) {
            integral=cp(a)*(b-a);
        } else if (model==ScalarPropertyModel::TABLE || model==ScalarPropertyModel::LINEAR) {
            // Every segment is affine. Insert table knots and each bound
            // crossing; trapezoids then integrate the ACTUAL cp exactly.
            std::vector<double> knots{a,b};
            if(model==ScalarPropertyModel::TABLE) {
                heat_capacity.table.validate("heat capacity table");
                for(double x:heat_capacity.table.temperature)
                    if(x>a && x<b) knots.push_back(x);
            }
            std::sort(knots.begin(),knots.end());
            const auto raw=[&](double x) {
                if(model==ScalarPropertyModel::TABLE)
                    return heat_capacity.table.evaluate(x,"heat capacity table");
                return heat_capacity.reference_value+
                    (heat_capacity.coefficients.empty()?0.0:heat_capacity.coefficients[0])*
                    (x-heat_capacity.reference_temperature);
            };
            for(std::size_t i=1;i<knots.size();++i) {
                const double l=knots[i-1],r=knots[i],vl=raw(l),vr=raw(r);
                std::vector<double> cuts{l,r};
                if(heat_capacity.enforce_bounds && vr!=vl) {
                    for(double bound:{heat_capacity.minimum,heat_capacity.maximum}) {
                        const double w=(bound-vl)/(vr-vl);
                        if(w>0.0 && w<1.0) cuts.push_back(l+w*(r-l));
                    }
                }
                std::sort(cuts.begin(),cuts.end());
                for(std::size_t j=1;j<cuts.size();++j)
                    integral+=0.5*(cp(cuts[j-1])+cp(cuts[j]))*(cuts[j]-cuts[j-1]);
            }
        } else if(model==ScalarPropertyModel::POLYNOMIAL && !heat_capacity.enforce_bounds) {
            for(std::size_t i=0;i<heat_capacity.coefficients.size();++i) {
                const double power=static_cast<double>(i+1);
                integral+=heat_capacity.coefficients[i]*
                    (std::pow(b-heat_capacity.reference_temperature,power)-
                     std::pow(a-heat_capacity.reference_temperature,power))/power;
            }
        } else {
            // Nonlinear/bounded laws: adaptive Simpson integration of cp(),
            // not of the unbounded raw polynomial. Exhaustion is an error.
            std::function<double(double,double,double,double,double,double,double,int)> adapt;
            adapt=[&](double l,double r,double fl,double fm,double fr,double whole,double tol,int depth) {
                const double m=(l+r)*0.5;
                const double f1=cp((l+m)*0.5), f2=cp((m+r)*0.5);
                const double left=(m-l)*(fl+4*f1+fm)/6;
                const double right=(r-m)*(fm+4*f2+fr)/6;
                const double delta=left+right-whole;
                if(std::abs(delta)<=15*tol) return left+right+delta/15;
                if(depth==0) throw std::runtime_error("enthalpy quadrature did not converge");
                return adapt(l,m,fl,f1,fm,left,tol/2,depth-1)+
                       adapt(m,r,fm,f2,fr,right,tol/2,depth-1);
            };
            const std::size_t n=static_cast<std::size_t>(intervals);
            for(std::size_t i=0;i<n;++i) {
                const double l=a+(b-a)*i/n,r=a+(b-a)*(i+1)/n;
                const double fl=cp(l),fm=cp((l+r)*0.5),fr=cp(r);
                const double whole=(r-l)*(fl+4*fm+fr)/6;
                integral+=adapt(l,r,fl,fm,fr,whole,1e-11*std::max(1.0,std::abs(whole)),20);
            }
        }
        if (!std::isfinite(integral)) throw std::runtime_error("enthalpy integral non-finite");
        return T>=Tref?integral:-integral;
    }
};
inline double turbulent_thermal_diffusivity(double mu_t,double rho,double Pr_t) {
    if(!std::isfinite(mu_t)||mu_t<0||!std::isfinite(rho)||rho<=0||!std::isfinite(Pr_t)||Pr_t<=0) throw std::invalid_argument("invalid turbulent thermal diffusivity inputs");
    return mu_t/(rho*Pr_t);
}
}
