#include "cfdx/physics/coupled_fingerprint_v40.h"
#include <cassert>
#include <limits>
#include <stdexcept>
#include <iostream>
using cfdx::core::SparseMatrix;
using cfdx::core::Vector;
using cfdx::physics::fingerprint_coupled_v40;
int main() {
    SparseMatrix A(2,2); A.push_back(0,0,2.0);A.push_back(1,1,3.0);A.finalize();
    Vector b(2,1.0),x(2,0.0);
    const auto h=fingerprint_coupled_v40(A,b,x);
    const auto same=fingerprint_coupled_v40(A,b,x);
    assert(h.A==same.A && h.b==same.b && h.x0==same.x0);
    assert(h.A.size()==16 && h.b.size()==16 && h.x0.size()==16);
    b(0)=2.0;
    const auto changed_b=fingerprint_coupled_v40(A,b,x);
    assert(h.A==changed_b.A && h.b!=changed_b.b && h.x0==changed_b.x0);
    x(1)=1.0;
    assert(fingerprint_coupled_v40(A,b,x).x0!=h.x0);
    SparseMatrix B(2,2);B.push_back(0,0,2.0);B.push_back(1,1,4.0);B.finalize();
    assert(fingerprint_coupled_v40(B,b,x).A!=h.A);
    b(0)=std::numeric_limits<double>::quiet_NaN();
    try {(void)fingerprint_coupled_v40(A,b,x);assert(false);}
    catch(const std::invalid_argument&){}
    std::cout<<"V40_COUPLED_FINGERPRINT 6/6 PASS\n";
}
