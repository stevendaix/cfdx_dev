#pragma once
#include "cfdx/physics/wall_distance.h"
#include "cfdx/physics/turbulence_models.h"
#include "cfdx/physics/turbulence_solver.h"
#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/spalart_allmaras.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
namespace cfdx::physics {
struct WallDistanceTurbulenceDiagnostics { double min_distance=0,max_distance=0,min_y_plus=0,max_y_plus=0,des_les_fraction=0,ddes_les_fraction=0,iddes_les_fraction=0; std::size_t samples=0; };
inline void validate_wall_distance_field(const std::vector<double>& d){for(double v:d)if(!std::isfinite(v)||v<=0)throw std::invalid_argument("wall-distance field must contain finite positive values");}
inline std::vector<double> wall_distance_y_plus(
    const std::vector<double>& d, const std::vector<double>& ut, double nu) {
    if (d.size() != ut.size() || d.empty() || nu <= 0.0 || !std::isfinite(nu))
        throw std::invalid_argument("invalid wall-distance/y+ inputs");
    validate_wall_distance_field(d);
    std::vector<double> y(d.size());
    for (std::size_t i = 0; i < d.size(); ++i) {
        if (!std::isfinite(ut[i]) || ut[i] < 0.0)
            throw std::invalid_argument("friction velocity must be finite and non-negative");
        double value = d[i] * ut[i] / nu;
        if (!std::isfinite(value) || (value == 0.0 && ut[i] > 0.0)) {
            // Recover from intermediate overflow *or* underflow without
            // penalizing the normal path. Use an alternative double ordering
            // first (portable where long double has double precision).
            const double alternative = (d[i] / nu) * ut[i];
            if (std::isfinite(alternative) && alternative > 0.0) {
                value = alternative;
            } else {
                const long double extended =
                    (static_cast<long double>(d[i]) * ut[i]) / nu;
                if (!std::isfinite(extended) ||
                    extended > std::numeric_limits<double>::max())
                    throw std::overflow_error("wall-distance y+ is non-representable");
                value = static_cast<double>(extended);
            }
        }
        if (!std::isfinite(value) || value < 0.0)
            throw std::overflow_error("wall-distance y+ is non-representable");
        y[i] = value;
    }
    return y;
}
inline double des_wall_transition_ratio(double d, double delta, double Cdes=.65) {
    if (!std::isfinite(d) || !std::isfinite(delta) || !std::isfinite(Cdes) ||
        d<=0.0 || delta<=0.0 || Cdes<=0.0)
        throw std::invalid_argument("invalid DES wall-distance inputs");
    const double ratio = d/(Cdes*delta);
    if (!std::isfinite(ratio))
        throw std::overflow_error("non-representable DES wall-distance ratio");
    return ratio;
}
inline double ddes_shielding_from_wall_distance(double nu_tilde,double s_hat,double d,double nu,double kappa=.41){
 if(!std::isfinite(nu_tilde)||!std::isfinite(s_hat)||!std::isfinite(d)||!std::isfinite(nu)||!std::isfinite(kappa)||nu_tilde<0||s_hat<=0||d<=0||nu<=0||kappa<=0)throw std::invalid_argument("invalid DDES shielding inputs");
 const double chi=nu_tilde/nu;
 if (!std::isfinite(chi))
     throw std::overflow_error("non-representable DDES viscosity ratio");
 SpalartAllmarasModel sa; sa.kappa=kappa;
 const double st=std::max(sa_modified_vorticity(s_hat,nu_tilde,chi,d,sa),1e-30);
 const double rd=std::min(nu_tilde/(st*kappa*kappa*d*d),10.0); return ddes_shielding(rd);
}
inline WallDistanceTurbulenceDiagnostics audit_wall_distance_turbulence_coupling(const std::vector<double>& d,const std::vector<double>& ut,double nu,double delta,double Cdes=.65){
 if(d.size()!=ut.size()||d.empty())throw std::invalid_argument("invalid wall-distance diagnostic fields");
 const auto yp=wall_distance_y_plus(d,ut,nu); WallDistanceTurbulenceDiagnostics o; o.samples=d.size();
 o.min_distance=*std::min_element(d.begin(),d.end()); o.max_distance=*std::max_element(d.begin(),d.end());
 o.min_y_plus=*std::min_element(yp.begin(),yp.end()); o.max_y_plus=*std::max_element(yp.begin(),yp.end());
 std::size_t a=0,b=0,c=0; for(std::size_t i=0;i<d.size();++i){if(des_wall_transition_ratio(d[i],delta,Cdes)>1)++a; const double fd=ddes_shielding_from_wall_distance(1e-4,1.0,d[i],nu); if(fd<.5)++b; if(iddes_shielding(0.5*fd,0.5)<.5)++c;}
 const double sample_count = static_cast<double>(o.samples); o.des_les_fraction=static_cast<double>(a)/sample_count;o.ddes_les_fraction=static_cast<double>(b)/sample_count;o.iddes_les_fraction=static_cast<double>(c)/sample_count;return o;
}
inline void require_wall_distance_resolves_surface(
    const WallDistanceGrid& g, const std::vector<double>& d) {
    if (d.size()!=g.points.size() || g.solid.size()!=g.points.size())
        throw std::invalid_argument("wall-distance field/grid mask size mismatch");
    for (std::size_t i=0; i<d.size(); ++i)
        if (!g.solid[i] && (!std::isfinite(d[i]) || d[i]<=0.0))
            throw std::domain_error("invalid fluid wall distance");
}
} // namespace cfdx::physics
