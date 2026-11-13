#pragma once
// V39 experimental angular refinement candidate. Tensor Gauss-Legendre in
// cos(theta), periodic midpoint trapezoid in phi, and antipodal symmetry.
// Does NOT certify ray effects or volumetric transport for any geometry.
#include "cfdx/physics/radiation_solver.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
inline std::vector<DiscreteDirection> make_dom_tensor_gauss_legendre_v39(
    std::size_t polar_order,std::size_t azimuth_count)
{
    if(polar_order<2 || polar_order>48 || polar_order%2!=0 ||
       azimuth_count<8 || azimuth_count>192 || azimuth_count%4!=0)
        throw std::invalid_argument("DOM tensor GL: even polar 2..48, azimuth 8..192 divisible by 4");
    std::vector<DiscreteDirection> directions;
    directions.reserve(polar_order*azimuth_count);
    std::vector<double> nodes(polar_order),weights(polar_order);
    const double pi=std::acos(-1.0);
    for(std::size_t i=0;i<polar_order/2;++i){
        double x=std::cos(pi*(double(i)+.75)/(double(polar_order)+.5));
        double derivative=0.;
        bool converged=false;
        for(int sweep=0;sweep<40;++sweep){
            double pm2=1.,pm1=x;
            for(std::size_t n=2;n<=polar_order;++n){
                const double p=(double(2*n-1)*x*pm1-double(n-1)*pm2)/double(n);
                pm2=pm1;pm1=p;
            }
            derivative=double(polar_order)*(x*pm1-pm2)/(x*x-1.);
            if(!std::isfinite(derivative)||derivative==0.)
                throw std::runtime_error("DOM tensor GL: invalid Legendre derivative");
            const double dx=pm1/derivative;
            x-=dx;
            if(std::abs(dx)<4.e-16){converged=true;break;}
        }
        if(!converged || !(x>0.) || !(x<1.))
            throw std::runtime_error("DOM tensor GL: root did not converge");
        // Re-evaluate the derivative at the final root, without using the
        // pre-Newton value from the previous iteration.
        double pm2=1.,pm1=x;
        for(std::size_t n=2;n<=polar_order;++n){
            const double p=(double(2*n-1)*x*pm1-double(n-1)*pm2)/double(n);
            pm2=pm1;pm1=p;
        }
        derivative=double(polar_order)*(x*pm1-pm2)/(x*x-1.);
        const double w=2./((1.-x*x)*derivative*derivative);
        if(!(w>0.) || !std::isfinite(w))
            throw std::runtime_error("DOM tensor GL: invalid weight");
        nodes[i]=-x;nodes[polar_order-1-i]=x;
        weights[i]=w;weights[polar_order-1-i]=w;
    }
    for(std::size_t i=0;i<polar_order;++i){
        const double z=nodes[i], r=std::sqrt(std::max(0.,1.-z*z));
        for(std::size_t j=0;j<azimuth_count;++j){
            const double phi=2.*pi*(double(j)+.5)/double(azimuth_count);
            const double w=weights[i]*2.*pi/double(azimuth_count);
            directions.push_back({r*std::cos(phi),r*std::sin(phi),z,w});
        }
    }
    validate_discrete_directions(directions);
    return directions;
}
}
