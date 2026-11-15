#include "cfdx/physics/bfs_sst_budget_forensics_v75.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cfdx::physics;
static void need(bool cond,const char* text){if(!cond)throw std::runtime_error(text);}
template<class F>static bool rejected(F fn){try{fn();}catch(const std::exception&){return true;}return false;}
int main(){
    BfsSstCellSnapshotV75 a;
    a.cell=12;a.k=.2;a.omega=2.;a.nut=.03;a.raux=500.;a.phi3=3e-4;
    a.equation[0]={2.,-1.,-3.,1.,.5,-.5};a.raw[0]=-.5;
    a.equation[1]={.1,.2,-.4,-.2,.1,-.2};a.raw[1]=-.2;
    auto b=a;
    b.k=.21;b.omega=2.01;b.nut=.028;b.raux=510.;b.phi3=2.9e-4;
    b.equation[0]={1.5,-.5,-2.,.5,.7,.2};b.raw[0]=.2;
    b.equation[1]={.2,.1,-.3,-.15,.15,0.};b.raw[1]=0.;
    const auto d=bfs_compare_sst_budget_v75(a,b);
    need(d.independently_closed && d.cell==12,"V75 common-state closure failed");
    need(std::abs(d.equation[0].sum-.7)<1e-13 &&
         std::abs(d.equation[1].sum-.2)<1e-13,"V75 scalar decomposition is wrong");
    need(d.equation[0].largest_absolute_contribution==std::string("production") &&
         std::abs(d.phi3_change+1e-5)<1e-14,"V75 signal attribution failed");
    need(rejected([&]{auto f=b;f.raw[0]+=.01;(void)bfs_compare_sst_budget_v75(a,f);}),
         "V75 missed a raw Fphys vs FV row mismatch");
    need(rejected([&]{auto f=b;f.equation[1].production+=.03;
        (void)bfs_compare_sst_budget_v75(a,f);}),
         "V75 missed a nonconservative term difference");
    need(rejected([&]{auto f=b;f.omega=std::numeric_limits<double>::quiet_NaN();
        (void)bfs_compare_sst_budget_v75(a,f);}),"V75 allowed NaN");
    need(rejected([&]{auto f=b;f.cell=0;(void)bfs_compare_sst_budget_v75(a,f);}),
         "V75 accepted a cross-cell comparison");
    need(rejected([&]{auto f=b;f.raux=-1.;(void)bfs_compare_sst_budget_v75(a,f);}),
         "V75 accepted a nonphysical momentum diagonal");
    // A tiny-scale defect must not pass due to the old max(1, |residual|) floor.
    auto tiny_a=a,tiny_b=a;
    tiny_a.equation[0]={0.,0.,0.,0.,0.,0.};tiny_a.raw[0]=0.;
    tiny_b.equation[0]={1e-7,0.,0.,0.,0.,1e-7};tiny_b.raw[0]=1e-7;
    tiny_b.equation[0].diffusion=1e-9;
    need(rejected([&]{(void)bfs_compare_sst_budget_v75(tiny_a,tiny_b);}),
         "V75 masked small-scale SST imbalance");
    std::cout<<"V75_SST_BUDGET_PASS independent_decomposition=2 counterfeit_rejected=6 committed=0\n";
}
