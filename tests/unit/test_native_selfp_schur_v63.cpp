#include "cfdx/core/linalg/native_selfp_schur_v63.h"
#include <stdexcept>
#include <cmath>
#include <iostream>
using namespace cfdx::core;
#define CHECK(expr) do { if (!(expr)) { \
    std::cerr << "CHECK FAILED: " << #expr << "\n"; return 1; \
} } while(false)
int main() {
    const std::size_t n=4;
    SparseMatrix A(4*n,4*n);
    for(std::size_t i=0;i<3*n;++i) {
        A.push_back(i,i,3.0+0.1*i);
        if(i%4!=3)A.push_back(i,i+1,-0.15);
        if(i%4!=0)A.push_back(i,i-1,-0.15);
        const auto c=i%n;
        // Pressure gauge row 0: no D coupling, C(0,0)=1.
        A.push_back(i,3*n+c,0.2);
    }
    A.push_back(3*n,3*n,1.0);
    for(std::size_t c=1;c<n;++c) {
        A.push_back(3*n+c,3*n+c,0.25);
        A.push_back(3*n+c,c,-0.3);
        A.push_back(3*n+c,n+c,-0.15);
    }
    A.finalize();
    NativeSelfpILUSchurV63 pc(n);
    CHECK(pc.setup(A));
    CHECK(pc.schur_matrix().n_rows()==n);
    CHECK(std::abs(pc.schur_matrix()(0,0)-1.)<1e-14);
    CHECK(pc.schur_matrix()(1,1)>0.25);
    Vector r(4*n),z(4*n),image(4*n);
    for(std::size_t i=0;i<4*n;++i)r(i)=0.1+0.03*i;
    CHECK(pc.apply(r,z));
    CHECK(z.is_valid());
    A.matvec_into(z,image);
    double rr=0,nn=0;
    for(std::size_t i=0;i<4*n;++i) {
        rr+=(r(i)-image(i))*(r(i)-image(i));nn+=r(i)*r(i);
    }
    if(!(std::sqrt(rr/nn)<0.15))return 1;
    // A missing Auu diagonal must fail rather than change physical gauge.
    SparseMatrix invalid(4*n,4*n);
    for(std::size_t row=0;row<4*n;++row) {
        if(row!=3)invalid.push_back(row,row,1.0);
    }
    invalid.finalize();
    CHECK(!pc.setup(invalid));
    CHECK(!pc.apply(r,z));
    std::cout<<"V63_SELFP_NATIVE_PASS relative="<<std::sqrt(rr/nn)<<"\n";
}