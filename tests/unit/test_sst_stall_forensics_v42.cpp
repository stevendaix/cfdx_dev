#include "cfdx/physics/sst_stall_forensics_v42.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#define CHECK(x) do { if(!(x)) throw std::runtime_error("V42 SST failed: " #x); } while(false)
int main() {
    using namespace cfdx::physics;
    using namespace cfdx::core;
    // Same frozen finite-volume rows: 2x-1y=2, -x+3y=4 at (x,y)=(2,2).
    ScalarEquation eq;
    eq.matrix=SparseMatrix(2,2);
    eq.matrix.push_back(0,0,2.);
    eq.matrix.push_back(0,1,-1.);
    eq.matrix.push_back(1,0,-1.);
    eq.matrix.push_back(1,1,3.);
    eq.matrix.finalize();
    eq.rhs=Vector(2);
    eq.rhs(0)=2.;eq.rhs(1)=4.;
    Vector solution(2);solution(0)=2.;solution(1)=2.;
    auto exact=sst_stall_overview_v42(eq,solution,1e-9);
    CHECK(exact.native_row_scaled_worst.native_row_scaled==0.0);
    CHECK(exact.compensated_local_worst.compensated_local_relative==0.0);
    solution(1)+=0.1;
    auto row0=sst_stall_row_v42(eq,solution,0,1e-9);
    CHECK(std::abs(row0.compensated_signed+0.1)<1e-14);
    CHECK(std::abs(row0.native_row_scaled-0.1/3.0)<1e-14);
    auto rows=sst_stall_overview_v42(eq,solution,1e-9);
    CHECK(rows.native_row_scaled_worst.cell==1);
    CHECK(rows.native_local_worst.cell==1);
    CHECK(rows.compensated_row_scaled_worst.cell==1);
    CHECK(rows.native_row_scaled_worst.native_row_scaled>0.0);
    // The arithmetic diagnostic must not replace or relax the gate.
    CHECK(rows.native_row_scaled_worst.native_row_scaled>1e-8);
    bool bad=false;
    try {(void)sst_stall_row_v42(eq,solution,2,1e-9);}
    catch(const std::invalid_argument&) {bad=true;}
    CHECK(bad);
    bad=false;
    try {(void)sst_stall_row_v42(eq,solution,0,0.0);}
    catch(const std::invalid_argument&) {bad=true;}
    CHECK(bad);
    std::cout<<"SST_STALL_FORENSICS_V42_PASS exact_rows=1 no_qualifier=1\n";
}
