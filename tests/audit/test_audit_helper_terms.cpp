// F12/F13 audit 2026-10-24: evaluate_momentum_terms / evaluate_energy_terms
// must return homogeneous per-unit-volume terms. Independent oracles are
// hand-derived on a two-hex stack (unit and uniformly scaled by 2) plus a
// single non-unit cube; the rescaling checks catch an integrated (flux-sum)
// convection term, which scales as s^2 instead of 1/s.
#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/physics/incompressible.h"
#include "cfdx/physics/thermal.h"
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

Mesh single_cube(double L)
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                          {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for(std::size_t i=0;i<8;++i)
        m.points().set(i,L*p[i][0],L*p[i][1],L*p[i][2]);
    m.faces().push_face({0,3,2,1});m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for(std::size_t f=0;f<6;++f){
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    return m;
}

// Two unit hexes stacked along y, uniformly scaled by s. Points 8..15 are
// points 0..7 shifted by +1 in y (scaled). The shared face f3 = {3,7,6,2}.
Mesh two_hex_stack(double s)
{
    Mesh m;
    m.points().resize(16);
    const double p[16][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                           {0,0,1},{1,0,1},{1,1,1},{0,1,1},
                           {0,1,0},{1,1,0},{1,2,0},{0,2,0},
                           {0,1,1},{1,1,1},{1,2,1},{0,2,1}};
    for(std::size_t i=0;i<16;++i)
        m.points().set(i,s*p[i][0],s*p[i][1],s*p[i][2]);
    // cell 0
    m.faces().push_face({0,3,2,1});   // f0 z-
    m.faces().push_face({4,5,6,7});   // f1 z+
    m.faces().push_face({0,1,5,4});   // f2 y- (boundary)
    m.faces().push_face({3,7,6,2});   // f3 y+ shared with cell 1
    m.faces().push_face({0,4,7,3});   // f4 x-
    m.faces().push_face({1,2,6,5});   // f5 x+
    // cell 1
    m.faces().push_face({8,11,10,9}); // f6 z-
    m.faces().push_face({12,13,14,15}); // f7 z+
    m.faces().push_face({11,15,14,10}); // f8 y+ (boundary)
    m.faces().push_face({8,12,15,11});  // f9 x-
    m.faces().push_face({9,10,14,13});  // f10 x+
    m.ownership().resize(11);
    for(std::size_t f=0;f<11;++f){
        const bool shared = (f==3);
        m.ownership().set_owner(f,shared?0u:(f<6?0u:1u));
        m.ownership().set_neighbour(f,shared?static_cast<std::int64_t>(1):FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    m.cells().push_cell({3,6,7,8,9,10});
    return m;
}

void set_vec(Field<double,Location::CELL>& u,std::size_t c,double x,double y,double z)
{
    u.component_data(0)[c]=x;u.component_data(1)[c]=y;u.component_data(2)[c]=z;
}

int checks=0,failures=0;
void check(const char* name,double expected,double actual,double tol=1e-12)
{
    const bool ok=std::isfinite(actual)&&std::abs(actual-expected)<=tol;
    ++checks;failures+=!ok;
    std::cout<<(ok?"PASS":"FAIL")<<'\t'<<name
              <<"\texpected="<<expected
              <<"\tvalue="<<std::setprecision(17)<<actual<<'\n';
}

} // namespace

int main()
{
    // ---- F12: momentum helper returns per-unit-volume terms ----
    // Unit stack, U0=(0.5,1,0), U1=(1.5,1,0), nu=2.
    // Cell 0: every face value of Ux is the owner value 0.5 and the net
    // volumetric flux sums to zero, so conv_x(0)=0. Cell 1: outflow faces
    // carry U1x=1.5 except the shared inflow face (upwind owner value 0.5):
    // +1.5 (y+) - 2.25 (x-) + 2.25 (x+) - 0.5 (shared) = 1.0.
    {
        auto m=two_hex_stack(1.0);
        Field<double,Location::CELL> u(2,"U","m/s",3),p(2,"p","Pa",1);
        set_vec(u,0,0.5,1.0,0.0);set_vec(u,1,1.5,1.0,0.0);
        p.fill(0.0);
        const auto t=evaluate_momentum_terms(u,p,m,2.0);
        check("F12/unit/convection_x_cell0",0.0,t.convection.component_data(0)[0]);
        check("F12/unit/convection_x_cell1",1.0,t.convection.component_data(0)[1]);
        check("F12/unit/convection_y_cell0",0.0,t.convection.component_data(1)[0]);
        // diffusion: nu * alpha*(U1x-U0x)/V = 2*1*1/1 = 2 (cell 0), -2 (cell 1)
        check("F12/unit/diffusion_x_cell0",2.0,t.diffusion.component_data(0)[0]);
        check("F12/unit/diffusion_x_cell1",-2.0,t.diffusion.component_data(0)[1]);
    }
    // Same field on the stack scaled by 2: per-volume convection scales 1/s
    // (an integrated flux-sum would scale s^2), diffusion scales 1/s^2.
    {
        auto m=two_hex_stack(2.0);
        Field<double,Location::CELL> u(2,"U","m/s",3),p(2,"p","Pa",1);
        set_vec(u,0,0.5,1.0,0.0);set_vec(u,1,1.5,1.0,0.0);
        p.fill(0.0);
        const auto t=evaluate_momentum_terms(u,p,m,2.0);
        check("F12/scaled2/convection_x_cell1",0.5,t.convection.component_data(0)[1]);
        check("F12/scaled2/convection_x_cell0",0.0,t.convection.component_data(0)[0]);
        check("F12/scaled2/diffusion_x_cell0",0.5,t.diffusion.component_data(0)[0]);
        check("F12/scaled2/diffusion_x_cell1",-0.5,t.diffusion.component_data(0)[1]);
    }

    // ---- F13: energy helper returns K/s terms ----
    // Single cube of side L, uniform T=400 K, mass flux 2 kg/s leaving
    // through face 0 only, rho=5, cp=2, k=1:
    //   convection = 2*400/(rho*V) [K/s], conduction = 0.
    // The old contract returned cp*sum(mdot*T) = 3200 W on BOTH sizes.
    for(const double L:{1.0,2.0}) {
        auto m=single_cube(L);
        Field<double,Location::CELL> T(1,"T","K",1);T.fill(400.0);
        Field<double,Location::FACE> mdot(6,"mdot","kg/s",1);mdot.fill(0.0);
        mdot(0)=2.0;
        const auto e=evaluate_energy_terms(T,mdot,m,5.0,2.0,1.0);
        const double expected=2.0*400.0/(5.0*L*L*L);
        check("F13/cube/convection_K_per_s",expected,e.convection(0));
        check("F13/cube/conduction_zero",0.0,e.conduction(0));
    }
    // Stack with T0=300, T1=400, zero mass flux: two-point orthogonal
    // conduction alpha=1 (Sf=(0,1,0), d=(0,1,0)), V=1 per cell:
    //   conduction(0) = k/(rho*cp) * 1*100/1 = 10 K/s, conduction(1) = -10.
    // On the x2 stack: alpha=2, V=8 -> 2.5 / -2.5. The old k/rho scaling
    // gave 20 (unit) which is dimensionally inconsistent with K/s convection.
    for(const double s:{1.0,2.0}) {
        auto m=two_hex_stack(s);
        Field<double,Location::CELL> T(2,"T","K",1);
        T(0)=300.0;T(1)=400.0;
        Field<double,Location::FACE> mdot(11,"mdot","kg/s",1);mdot.fill(0.0);
        const auto e=evaluate_energy_terms(T,mdot,m,5.0,2.0,1.0);
        const double expected=10.0/(s*s);
        check("F13/stack/conduction_cell0",expected,e.conduction(0));
        check("F13/stack/conduction_cell1",-expected,e.conduction(1));
        check("F13/stack/convection_zero",0.0,e.convection(0));
        check("F13/stack/convection_zero_cell1",0.0,e.convection(1));
    }

    std::cout<<"SUMMARY checks="<<checks<<" failures="<<failures<<'\n';
    return failures?1:0;
}
