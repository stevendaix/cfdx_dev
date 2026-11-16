#include "cfdx/physics/convective_face_tangent_v55.h"
#include <cmath>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void check(bool b){if(!b)throw std::runtime_error("V55 face tangent regression");}
int main(){
    for(double phi:{-2.0,-1.e-4,1.e-4,2.0}){
        const double a=3.,b=8.,dp=.23,da=.1,db=-.15,h=1.e-7;
        const auto t=donor_face_frozen_tangent_v55(phi,a,b,dp,da,db);
        const auto fd=(donor_face_flux_v55(phi+h*dp,a+h*da,b+h*db)-
                       donor_face_flux_v55(phi-h*dp,a-h*da,b-h*db))/(2*h);
        check(std::abs(t-fd)<1.e-7);
        std::vector<double> r(2);
        scatter_conservative_face_v55(r,0,1,t);
        check(std::abs(r[0]+r[1])<1.e-15);
    }
    // At zero flux, the function has different directional derivatives;
    // never advertise a unique smooth Jacobian at a branch change.
    const double h=1.e-6;
    const auto right=(donor_face_flux_v55(h,3.,8.)-
                      donor_face_flux_v55(0.,3.,8.))/h;
    const auto left=(donor_face_flux_v55(0.,3.,8.)-
                     donor_face_flux_v55(-h,3.,8.))/h;
    check(std::abs(right-3.)<1.e-12 && std::abs(left-8.)<1.e-12);
    const auto audit=audit_face_branches_v55(
        {0.,0.1,2.,-0.000001},{0.5,0.2,2.,0.1},
        {-0.5,-0.1,1.,-0.1},1.e-4);
    check(audit.crossings==3 && audit.near_zero_crossings==2);
    check(audit.events.size()==3 && audit.events[0].face==0 &&
          audit.events[1].face==1 && audit.events[2].face==3);
    bool rejected=false;
    try{(void)audit_face_branches_v55({0.},{1.},{});}
    catch(const std::invalid_argument&){rejected=true;}
    check(rejected);
}
