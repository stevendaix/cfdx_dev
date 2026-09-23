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
    TabulatedProperty table;
};
inline double evaluate_density(const DensityControls& c,double T,double p=101325.0) {
    if(!std::isfinite(T)||T<=0||!std::isfinite(c.rho0)||c.rho0<=0) throw std::invalid_argument("density: invalid inputs");
    double rho=c.rho0;
    switch(c.model) {
    case DensityModel::CONSTANT: break;
    case DensityModel::LINEAR_TEMPERATURE:
    case DensityModel::BOUSSINESQ: rho=c.rho0*(1-c.beta*(T-c.reference_temperature)); break;
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
            if (a < heat_capacity.table.temperature.front() ||
                b > heat_capacity.table.temperature.back()) {
                if (heat_capacity.table.extrapolation == ExtrapolationPolicy::REJECT)
                    throw std::out_of_range("enthalpy: heat-capacity table outside range");
            }
            const double xa=std::max(a,heat_capacity.table.temperature.front());
            const double xb=std::min(b,heat_capacity.table.temperature.back());
            if (xb>xa) {
                double prevT=xa, prevV=heat_capacity.table.evaluate(xa,"heat capacity table");
                for (std::size_t i=0;i<heat_capacity.table.temperature.size();++i) {
                    const double t=heat_capacity.table.temperature[i];
                    if(t<=xa || t>=xb) continue;
                    const double v=heat_capacity.table.evaluate(t,"heat capacity table");
                    integral += 0.5*(prevV+v)*(t-prevT); prevT=t; prevV=v;
                }
                const double endV=heat_capacity.table.evaluate(xb,"heat capacity table");
                integral += 0.5*(prevV+endV)*(xb-prevT);
            }
            if (a < xa || b > xb) {
                if (heat_capacity.table.extrapolation == ExtrapolationPolicy::CLAMP) {
                    if(a<xa) integral += heat_capacity.table.value.front()*(xa-a);
                    if(b>xb) integral += heat_capacity.table.value.back()*(b-xb);
                } else if (heat_capacity.table.extrapolation == ExtrapolationPolicy::LINEAR) {
                    const double va=heat_capacity.table.evaluate(a,"heat capacity table");
                    const double vb=heat_capacity.table.evaluate(b,"heat capacity table");
                    const double inside=(xb>xa)?(0.5*(heat_capacity.table.evaluate(xa,"heat capacity table")+heat_capacity.table.evaluate(xb,"heat capacity table))*(xb-xa)):0.0;
                    integral = inside + 0.5*(va+heat_capacity.table.evaluate(xa,"heat capacity table"))*(xa-a)
                                     + 0.5*(heat_capacity.table.evaluate(xb,"heat capacity table")+vb)*(b-xb);
                }
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