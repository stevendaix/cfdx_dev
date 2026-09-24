#pragma once

#include "cfdx/physics/radiation.h"
#include "cfdx/physics/radiation_advanced.h"
#include "cfdx/physics/radiation_models.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct S2SControls {
    double tolerance = 1e-10;
    std::size_t max_iterations = 1;
    double ambient_irradiation = 0.0;
};

struct S2SResult {
    bool converged = false;
    std::vector<double> radiosity;
    std::vector<double> irradiation;
    std::vector<double> net_flux;
    double total_power = 0.0;
    double energy_balance_error = 0.0;
};

inline void validate_s2s_view_factors(
    const std::vector<double>& F,
    const std::vector<double>& areas,
    double tolerance = 1e-10)
{
    const std::size_t n = areas.size();
    if (n == 0 || F.size() != n*n)
        throw std::invalid_argument("S2S view-factor matrix size mismatch");
    if (!(tolerance > 0.0) || !std::isfinite(tolerance))
        throw std::invalid_argument("invalid S2S tolerance");

    for (double A : areas)
        if (!std::isfinite(A) || A <= 0.0)
            throw std::invalid_argument("S2S surface areas must be positive");

    for (std::size_t i=0; i<n; ++i) {
        double row_sum=0.0;
        for (std::size_t j=0; j<n; ++j) {
            const double fij=F[i*n+j];
            if (!std::isfinite(fij) || fij < -tolerance || fij > 1.0+tolerance)
                throw std::invalid_argument("invalid S2S view factor");
            if (i==j && std::abs(fij)>tolerance)
                throw std::invalid_argument("S2S self-view factor must be zero");
            row_sum += fij;
        }
        if (row_sum > 1.0+tolerance)
            throw std::invalid_argument("S2S view-factor row is not conservative");
    }

    for (std::size_t i=0; i<n; ++i)
        for (std::size_t j=i+1; j<n; ++j) {
            const double lhs=areas[i]*F[i*n+j];
            const double rhs=areas[j]*F[j*n+i];
            if (std::abs(lhs-rhs) > tolerance*std::max(1.0,std::max(std::abs(lhs),std::abs(rhs))))
                throw std::invalid_argument("S2S view factors violate reciprocity");
        }
}

inline S2SResult solve_s2s_radiosity(
    const std::vector<double>& areas,
    const std::vector<double>& emissivities,
    const std::vector<double>& temperatures,
    const std::vector<double>& view_factors,
    const S2SControls& controls = {},
    const std::vector<double>& external_irradiation = {})
{
    const std::size_t n=areas.size();
    if (n==0 || emissivities.size()!=n || temperatures.size()!=n)
        throw std::invalid_argument("S2S surface vector size mismatch");
    if (!external_irradiation.empty() && external_irradiation.size()!=n)
        throw std::invalid_argument("S2S external irradiation size mismatch");
    if (!(controls.tolerance>0.0) || controls.max_iterations==0 ||
        !std::isfinite(controls.ambient_irradiation) || controls.ambient_irradiation<0.0)
        throw std::invalid_argument("invalid S2S controls");

    validate_s2s_view_factors(view_factors,areas,controls.tolerance);

    std::vector<double> b(n,0.0);
    std::vector<double> M(n*n,0.0);
    for (std::size_t i=0;i<n;++i) {
        if (!std::isfinite(emissivities[i]) || emissivities[i]<=0.0 || emissivities[i]>1.0 ||
            !std::isfinite(temperatures[i]) || temperatures[i]<=0.0)
            throw std::invalid_argument("invalid S2S surface state");
        const double E=STEFAN_BOLTZMANN*std::pow(temperatures[i],4);
        const double Gext=external_irradiation.empty()
            ? controls.ambient_irradiation : external_irradiation[i];
        if (!std::isfinite(Gext) || Gext<0.0)
            throw std::invalid_argument("invalid S2S external irradiation");
        const double rho=1.0-emissivities[i];
        b[i]=E+rho*Gext;
        M[i*n+i]=1.0;
        for (std::size_t j=0;j<n;++j) M[i*n+j]-=rho*view_factors[i*n+j];
    }

    // Dense pivoted Gaussian elimination is intentional here: the S2S
    // radiosity system is a surface-space system, normally much smaller than
    // the volume discretization. It also provides a deterministic oracle.
    std::vector<double> A=M, x=b;
    for (std::size_t k=0;k<n;++k) {
        std::size_t p=k;
        for (std::size_t i=k+1;i<n;++i)
            if (std::abs(A[i*n+k])>std::abs(A[p*n+k])) p=i;
        if (std::abs(A[p*n+k]) <= std::numeric_limits<double>::epsilon())
            throw std::runtime_error("singular S2S radiosity system");
        if (p!=k) {
            for (std::size_t j=k;j<n;++j) std::swap(A[k*n+j],A[p*n+j]);
            std::swap(x[k],x[p]);
        }
        for (std::size_t i=k+1;i<n;++i) {
            const double q=A[i*n+k]/A[k*n+k];
            A[i*n+k]=0.0;
            for (std::size_t j=k+1;j<n;++j) A[i*n+j]-=q*A[k*n+j];
            x[i]-=q*x[k];
        }
    }
    for (std::size_t ii=0;ii<n;++ii) {
        const std::size_t i=n-1-ii;
        double s=x[i];
        for (std::size_t j=i+1;j<n;++j) s-=A[i*n+j]*x[j];
        x[i]=s/A[i*n+i];
    }

    S2SResult r;
    r.converged=true;
    r.radiosity=x;
    r.irradiation.resize(n,0.0);
    r.net_flux.resize(n,0.0);
    double absorbed=0.0, emitted=0.0;
    for (std::size_t i=0;i<n;++i) {
        const double Gext=external_irradiation.empty()
            ? controls.ambient_irradiation : external_irradiation[i];
        for (std::size_t j=0;j<n;++j) r.irradiation[i]+=view_factors[i*n+j]*x[j];
        r.irradiation[i]+=Gext;
        r.net_flux[i]=x[i]-r.irradiation[i];
        r.total_power+=areas[i]*r.net_flux[i];
        emitted+=areas[i]*emissivities[i]*STEFAN_BOLTZMANN*std::pow(temperatures[i],4);
        absorbed+=areas[i]*emissivities[i]*r.irradiation[i];
    }
    r.energy_balance_error=std::abs(r.total_power)/
        std::max(1.0,std::max(emitted,absorbed));
    return r;
}

struct S2SMonteCarloControls {
    std::size_t samples = 16384;
};

inline std::vector<double> estimate_s2s_view_factors_monte_carlo(
    const std::vector<std::vector<RadiationTriangle>>& surfaces,
    const std::vector<RadiationTriangle>& blockers = {},
    const S2SMonteCarloControls& controls = {})
{
    if (surfaces.empty() || controls.samples==0)
        throw std::invalid_argument("S2S Monte Carlo requires surfaces and samples");
    const std::size_t n=surfaces.size();
    std::vector<double> F(n*n,0.0);
    for (std::size_t i=0;i<n;++i) {
        for (std::size_t j=0;j<n;++j) {
            if (i==j) continue;
            F[i*n+j]=estimate_view_factor_ray_traced(
                surfaces[i], surfaces[j], blockers, controls.samples);
        }
    }
    return F;
}

} // namespace cfdx::physics
