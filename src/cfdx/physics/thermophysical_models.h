#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
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
        const double a=std::min(T,Tref), b=std::max(T,Tref);
        double integral=0.0;
        if (heat_capacity.model == ScalarPropertyModel::CONSTANT) {
            integral=heat_capacity.reference_value*(b-a);
        } else if (heat_capacity.model == ScalarPropertyModel::LINEAR) {
            const double c0=heat_capacity.reference_value;
            const double c1=heat_capacity.coefficients.empty()?0.0:heat_capacity.coefficients[0];
            const double xa=a-heat_capacity.reference_temperature, xb=b-heat_capacity.reference_temperature;
            integral=c0*(b-a)+0.5*c1*(xb*xb-xa*xa);
        } else if (heat_capacity.model == ScalarPropertyModel::POLYNOMIAL) {
            if (heat_capacity.coefficients.empty())
                throw std::invalid_argument("enthalpy: polynomial heat capacity requires coefficients");
                integral=0.0;
            for (std::size_t i=0;i<heat_capacity.coefficients.size();++i) {
                const double ai=heat_capacity.coefficients[i];
                if(!std::isfinite(ai)) throw std::invalid_argument("enthalpy: non-finite heat-capacity coefficient");
                const double power=static_cast<double>(i+1);
                integral += ai*(std::pow(b-heat_capacity.reference_temperature,power)
                              - std::pow(a-heat_capacity.reference_temperature,power))/power;
            }
        } else if (heat_capacity.model == ScalarPropertyModel::TABLE) {
            heat_capacity.table.validate("heat capacity table");
            const auto& tab=heat_capacity.table;
            auto integrate_in_range = [&](double lo,double hi) {
                if (hi<=lo) return 0.0;
                double result=0.0;
                double x0=lo, y0=tab.evaluate(lo,"heat capacity table");
                for (std::size_t i=0;i+1<tab.temperature.size();++i) {
                    const double l=std::max(lo,tab.temperature[i]);
                    const double r=std::min(hi,tab.temperature[i+1]);
                    if (r<=l) continue;
                    const double yl=tab.value[i] + (l-tab.temperature[i]) *
                        (tab.value[i+1]-tab.value[i])/(tab.temperature[i+1]-tab.temperature[i]);
                    const double yr=tab.value[i] + (r-tab.temperature[i]) *
                        (tab.value[i+1]-tab.value[i])/(tab.temperature[i+1]-tab.temperature[i]);
                    result += 0.5*(yl+yr)*(r-l);
                    x0=r; y0=yr;
                }
                return result;
            };
            if (a < tab.temperature.front() || b > tab.temperature.back()) {
                if (tab.extrapolation == ExtrapolationPolicy::REJECT)
                    throw std::out_of_range("enthalpy: heat-capacity table outside range");
                const double lo=tab.temperature.front(), hi=tab.temperature.back();
                integral=integrate_in_range(std::max(a,lo),std::min(b,hi));
                if (a<lo) {
                    if (tab.extrapolation == ExtrapolationPolicy::CLAMP)
                        integral += tab.value.front()*(lo-a);
                    else {
                        const double slope=(tab.value[1]-tab.value[0])/(tab.temperature[1]-tab.temperature[0]);
                        integral += tab.value.front()*(lo-a) + 0.5*slope*((lo-a)*(lo-a));
                    }
                }
                if (b>hi) {
                    if (tab.extrapolation == ExtrapolationPolicy::CLAMP)
                        integral += tab.value.back()*(b-hi);
                    else {
                        const double slope=(tab.value.back()-tab.value[tab.value.size()-2])/
                                           (tab.temperature.back()-tab.temperature[tab.temperature.size()-2]);
                        const double d=b-hi;
                        integral += tab.value.back()*d + 0.5*slope*d*d;
                    }
                }
            } else {
                integral=integrate_in_range(a,b);
            }
        } else {
            const std::size_t n=static_cast<std::size_t>(intervals);
            const double h=(b-a)/n;
            double s=cp(a)+cp(b);
            for(std::size_t i=1;i<n;++i) s+=2*cp(a+h*i);
            integral=0.5*h*s;
        }
        return T>=Tref?integral:-integral;
    }
};
inline double turbulent_thermal_diffusivity(double mu_t,double rho,double Pr_t) {
    if(!std::isfinite(mu_t)||mu_t<0||!std::isfinite(rho)||rho<=0||!std::isfinite(Pr_t)||Pr_t<=0) throw std::invalid_argument("invalid turbulent thermal diffusivity inputs");
    return mu_t/(rho*Pr_t);
}
}