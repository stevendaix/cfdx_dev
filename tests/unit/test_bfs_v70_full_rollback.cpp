@@ -0,0 +1,85 @@
// V70: exercise BOTH the atomic SST transaction helper and the actual
// production rollback pattern (NonlinearStateRollback + phi/SST snapshots).
// A projected state must either publish all five components or none.
#include "cfdx/physics/sst_production_coupling.h"
#include "cfdx/physics/timestep_control.h"
#include <iostream>
#include <stdexcept>
using namespace cfdx::physics;
using namespace cfdx::core;
static void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
int main(){
    Field<double,Location::CELL> u(2,"u","m/s",3),p(2,"p","Pa",1);
    Field<double,Location::FACE> phi(3,"phi","kg/s",1);
    SSTProductionState sst(2);
    for(std::size_t i=0;i<2;++i) {
        u.component_data(0)[i]=10.+static_cast<double>(i);
        u.component_data(1)[i]=2.;u.component_data(2)[i]=0.;
        p(i)=3.+static_cast<double>(i);
        sst.k(i)=.1+static_cast<double>(i)*.01;
        sst.omega(i)=20.+static_cast<double>(i);
        sst.nu_t(i)=.005;sst.wall_distance(i)=.01;
        sst.strain_rate(i)=3.;
    }
    phi(0)=-1.;phi(1)=1.;phi(2)=1.;
    const auto old_u=u; const auto old_p=p; const auto old_phi=phi;
    const auto old_sst=snapshot_sst_production_state(sst);
    SSTCoupledStateTransaction transaction(u,p,phi,sst);
    transaction.begin();
    u.component_data(0)[0]+=12.;p(0)+=100.;phi(0)+=.25;
    sst.k(0)+=.1;sst.omega(0)+=3.;sst.nu_t(0)+=.02;
    // Simulate an independent F_phys rejection and require FULL rollback.
    transaction.rollback();
    require(!transaction.active(),"rejected transaction remained active");
    for(std::size_t i=0;i<2;++i) {
        for(std::size_t d=0;d<3;++d)
            require(u.component_data(d)[i]==old_u.component_data(d)[i],"U leaked");
        require(p(i)==old_p(i),"pressure leaked");
        require(sst.k(i)==old_sst.k(i),"k leaked");
        require(sst.omega(i)==old_sst.omega(i),"omega leaked");
        require(sst.nu_t(i)==old_sst.nu_t(i),"nu_t leaked");
        require(sst.wall_distance(i)==old_sst.wall_distance(i),"wall distance leaked");
        require(sst.strain_rate(i)==old_sst.strain_rate(i),"strain leaked");
    }
    for(std::size_t f=0;f<3;++f)
        require(phi(f)==old_phi(f),"mass flux leaked");
    // Positive counterpart: an atomic commit publishes all fields together.
    transaction.begin();
    u.component_data(0)[0]+=1.;p(0)+=2.;phi(0)+=.1;
    sst.k(0)+=.01;sst.omega(0)+=.02;sst.nu_t(0)+=.03;
    transaction.commit();
    require(!transaction.active() && u.component_data(0)[0]==11. &&
        p(0)==5. && std::abs(phi(0)+.9)<1e-14 &&
        std::abs(sst.k(0)-.11)<1e-14 &&
        std::abs(sst.omega(0)-20.02)<1e-14 &&
        std::abs(sst.nu_t(0)-.035)<1e-14,
        "atomic full state commit failed");
    // Match the actual steady COUPLED loop: only U/p use
    // NonlinearStateRollback; carrier phi and SST are explicitly restored
    // from the SAME outer-iteration snapshots in the catch/retry paths.
    const auto accepted_u=u; const auto accepted_p=p;
    const auto accepted_flux=phi;
    const auto accepted_sst=snapshot_sst_production_state(sst);
    NonlinearStateRollback actual_outer(u,p);
    actual_outer.begin();
    const auto flux_snapshot=phi;
    const auto sst_snapshot=snapshot_sst_production_state(sst);
    u.component_data(0)[1]+=5.;p(1)+=10.;phi(2)+=.5;
    sst.k(1)+=.12;sst.omega(1)+=.33;sst.nu_t(1)+=.01;
    actual_outer.reject();
    phi=flux_snapshot;
    restore_sst_production_state(sst,sst_snapshot);
    require(!actual_outer.active(),"production rollback remained active");
    for(std::size_t i=0;i<2;++i){
        for(std::size_t d=0;d<3;++d)
            require(u.component_data(d)[i]==accepted_u.component_data(d)[i],
                    "production rollback U leak");
        require(p(i)==accepted_p(i),"production rollback p leak");
        require(sst.k(i)==accepted_sst.k(i),"production rollback k leak");
        require(sst.omega(i)==accepted_sst.omega(i),"production rollback omega leak");
        require(sst.nu_t(i)==accepted_sst.nu_t(i),"production rollback nut leak");
    }
    for(std::size_t f=0;f<3;++f)
        require(phi(f)==accepted_flux(f),"production rollback phi leak");
    std::cout<<"V70_ATOMIC_CARRIER_SST_ROLLBACK_PASS rejected=2 committed=1\n";
}
