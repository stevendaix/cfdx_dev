// V55 finite arithmetic and lost-in-double-component negative controls.
#include "cfdx/physics/dom_hilo_arithmetic_v55.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void require(bool v,const char* text){if(!v)throw std::runtime_error(text);}
int main(){
    try{
        const auto recovered=(DomHiLoV55(1.e16)+DomHiLoV55(1.))-DomHiLoV55(1.e16);
        require(recovered.hi==1. && recovered.lo==0.,
                "V55 split sum lost the canceled unit");
        const DomHiLoV55 x(1.0000000000000002,2.4e-17);
        const DomHiLoV55 y(1.0000000000000002,1.3e-17);
        const auto product=x*y;
        const auto quotient=x/y;
        const auto expected_product=dom_to_long_double_v55(x)*dom_to_long_double_v55(y);
        const auto expected_quotient=dom_to_long_double_v55(x)/dom_to_long_double_v55(y);
        require(std::abs(dom_to_long_double_v55(product)-expected_product)<1e-18L,
                "V55 multiplication differs from 80-bit oracle");
        require(std::abs(dom_to_long_double_v55(quotient)-expected_quotient)<1e-18L,
                "V55 division differs from 80-bit oracle");
        const auto one_third=DomHiLoV55(1.)/DomHiLoV55(3.);
        require(std::abs(dom_to_long_double_v55(one_third)-1.L/3.L)<1e-19L,
                "V55 division did not preserve low component");
        bool rejected=false;
        try{(void)(DomHiLoV55(1.)/DomHiLoV55(0.));}
        catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"V55 divisor-zero gate accepted");
        std::cout<<"V55_HILO_ARITHMETIC_PASS recovered="<<recovered.hi
                 <<" multiply_error="<<static_cast<double>(
                         std::abs(dom_to_long_double_v55(product)-expected_product))
                 <<" quotient_error="<<static_cast<double>(
                         std::abs(dom_to_long_double_v55(quotient)-expected_quotient))
                 <<"\n";
        return 0;
    }catch(const std::exception& e){
        std::cerr<<"V55_HILO_ARITHMETIC_FAIL "<<e.what()<<'\n';return 1;
    }
}
