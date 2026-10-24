// H05: a one-cell transient balance with an independent normal-distance oracle.
#include "chain_fixture.h"
#include "cfdx/physics/energy_solver.h"
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace cfdx::core;
using namespace cfdx::physics;
int main() {
    int checks=0,failures=0;
    for(const double shear:{0.,1.,3.}) for(const bool override_value:{false,true}) {
        auto m=make_1d_chain(1);
        for(std::size_t i=0;i<m.n_points();++i)
            m.points().set(i,m.points().x(i),m.points().y(i)+shear*(m.points().x(i)-1),m.points().z(i));
        const auto geometry=build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"F","",1);flux.fill(0);
        Field<double,Location::CELL> t(1,"T","K",1),old(1,"Told","K",1),source(1,"Q","",1);
        old.fill(300);source.fill(0);t.fill((10.*300+4.*400)/14.);
        EnergySolverControls c;c.cp=1;c.density=1;c.conductivity=2;c.dt=.1;
        ScalarBoundaryConditions bc;
        ScalarBoundaryFaceValues fv;
        if(override_value) fv.values["outlet"].assign(m.n_faces(),400.);
        else bc["outlet"]={ScalarBoundaryType::FIXED_VALUE,400,0};
        const auto* values=override_value?&fv:nullptr;
        const auto eq=assemble_energy_equation(m,geometry,flux,source,old,c,bc,values);
        Vector x(1,t(0));
        const double matrix=energy_balance_relative_from_equation(eq,x);
        const double diagnostic=energy_balance_relative(m,geometry,flux,t,old,source,c,bc,values);
        for(const auto& pair: {std::make_pair("matrix",matrix),std::make_pair("diagnostic",diagnostic)}) {
            const bool ok=std::isfinite(pair.second)&&pair.second<1e-12;
            ++checks;failures+=!ok;
            std::cout<<(ok?"PASS":"FAIL")<<'\t'<<pair.first<<"\tshear="<<shear<<"\toverride="<<override_value
                     <<"\tvalue="<<std::setprecision(17)<<pair.second<<'\n';
        }
    }
    std::cout<<"SUMMARY checks="<<checks<<" failures="<<failures<<'\n';
    return failures?1:0;
}
