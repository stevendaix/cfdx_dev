#include "cfdx/physics/post_sst_momentum_hotspot_v45.h"
#include "cfdx/physics/finite_volume_transport.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace cfdx::physics;
using namespace cfdx::core;
static void need(bool c,const char* e) { if (!c) throw std::runtime_error(e); }
int main() {
    try {
        ScalarEquation eq;
        eq.matrix=SparseMatrix(3,3);
        eq.matrix.push_back(0,0,2);eq.matrix.push_back(0,1,-1);
        eq.matrix.push_back(1,0,-1);eq.matrix.push_back(1,1,3);
        eq.matrix.push_back(2,2,4);
        eq.matrix.finalize();
        eq.rhs=Vector(3);
        eq.rhs(0)=2;eq.rhs(1)=4;eq.rhs(2)=8;
        Vector x(3);x(0)=2;x(1)=2;x(2)=2;
        auto d=post_sst_momentum_hotspot_v45(eq,x,1,4.0);
        need(d.original_abs==0 && d.cell==0 && d.original_rhs_relative==0,
             "exact momentum FV state failed");
        x(1)=2.1;
        d=post_sst_momentum_hotspot_v45(eq,x,1,4.0);
        need(d.cell==1,"momentum worst cell wrong");
        need(d.original_abs==scalar_equation_residual_inf(eq,x),
             "post-SST operator not identical to production residual");
        need(std::abs(d.original_rhs_relative-0.3/4.0)<1e-14,
             "post-SST normalized physical residual changed");
        need(std::abs(d.independent_row.compensated_signed-0.3)<1e-14,
             "compensated momentum cross-check failed");
        bool rejected=false;
        try { (void)post_sst_momentum_hotspot_v45(eq,x,3,4.0); }
        catch(const std::invalid_argument&) {rejected=true;}
        need(rejected,"bad component not rejected");
        rejected=false;
        try { (void)post_sst_momentum_hotspot_v45(
            eq,x,1,std::numeric_limits<double>::infinity()); }
        catch(const std::invalid_argument&) {rejected=true;}
        need(rejected,"nonfinite scale accepted");
        std::cout<<"POST_SST_MOMENTUM_HOTSPOT_V45_PASS true_CSR=1 independent=1 "
                   "no_gate_change=1\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<"POST_SST_MOMENTUM_HOTSPOT_V45_FAIL "<<error.what()<<'\n';
        return 1;
    }
}
