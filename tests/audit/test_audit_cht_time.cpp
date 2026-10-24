// Two independent thermal capacities coupled by an analytic series resistance.
// Vary the physical step and interface relaxation without changing the oracle.
#include "chain_fixture.h"
#include "cfdx/physics/cht_solver.h"
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace cfdx::core;
using namespace cfdx::physics;
int main() {
    int checks=0, failures=0;
    for (const double dt : {0.01, 0.1, 1.0}) {
        for (const double relaxation : {0.5, 1.0}) {
            auto a=make_1d_chain(1), b=make_1d_chain(1);
            for(std::size_t i=0;i<b.n_points();++i)
                b.points().set(i,b.points().x(i)+1,b.points().y(i),b.points().z(i));
            auto ga=build_fv_geometry(a), gb=build_fv_geometry(b);
            Field<double,Location::FACE> fa(a.n_faces(),"F","",1),fb(b.n_faces(),"F","",1);
            fa.fill(0); fb.fill(0);
            Field<double,Location::CELL> ta(1,"T","K",1),tb(1,"T","K",1),qa(1,"Q","",1),qb(1,"Q","",1);
            ta.fill(300);tb.fill(400);qa.fill(10);qb.fill(-5);
            EnergySolverControls ea,eb;
            ea.density=2;eb.density=3;ea.cp=eb.cp=1;
            ea.conductivity=1;eb.conductivity=4;ea.dt=eb.dt=dt;
            ea.relaxation=eb.relaxation=1;ea.tolerance=eb.tolerance=1e-12;
            ChtInterfaceControls c;
            c.region1_patch="outlet";c.region2_patch="inlet";
            c.conductivity1=1;c.conductivity2=4;c.relaxation=relaxation;
            c.tolerance=1e-11;c.temperature_tolerance=1e-11;c.max_iterations=500;
            ScalarBoundaryConditions ba,bb;
            for(const auto* n : {"inlet","outlet","wall"}) {
                ba[n]={ScalarBoundaryType::ZERO_GRADIENT,0,0};bb[n]=ba[n];
            }
            ba["outlet"]={ScalarBoundaryType::FIXED_VALUE,350,0};bb["inlet"]=ba["outlet"];
            const double g=1.0/(0.5/1.0+0.5/4.0), m1=2/dt, m2=3/dt;
            const double d=(m1+g)*(m2+g)-g*g, b1=m1*300+10, b2=m2*400-5;
            const double exact1=((m2+g)*b1+g*b2)/d, exact2=(g*b1+(m1+g)*b2)/d;
            auto check=[&](const char* name,bool ok,double value){++checks;failures+=!ok;
                std::cout<<name<<"\tdt="<<dt<<"\trelax="<<relaxation<<'\t'<<(ok?"PASS":"FAIL")<<'\t'<<std::setprecision(17)<<value<<'\n';};
            try {
                auto r=solve_two_region_cht(a,ga,b,gb,fa,fb,ta,tb,qa,qb,ea,eb,c,ba,bb);
                check("cht/status",r.converged,r.iterations);
                check("cht/temperature1",std::isfinite(ta(0))&&std::abs(ta(0)-exact1)<1e-8,ta(0)-exact1);
                check("cht/temperature2",std::isfinite(tb(0))&&std::abs(tb(0)-exact2)<1e-8,tb(0)-exact2);
                const double balance=2*(ta(0)-300)+3*(tb(0)-400)-dt*5;
                check("cht/energy",std::isfinite(balance)&&std::abs(balance)<1e-7,balance);
            } catch(const std::exception& e) {
                ++checks;++failures;std::cerr<<"cht/exception "<<e.what()<<'\n';
            }
        }
    }
    std::cout<<"SUMMARY\tchecks="<<checks<<"\tfailures="<<failures<<'\n';
    return failures?1:0;
}
