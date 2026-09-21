#include "cfdx/core/linalg/linear_operator.h"
#include "cfdx/core/memory/reuse_pool.h"
#include "cfdx/core/memory/field_lifetime.h"
#include "cfdx/core/linalg/mpi_overlap.h"
#include "cfdx/core/mesh/sfc_ordering.h"
#include "cfdx/core/linalg/krylov_controls.h"
#include "cfdx/core/linalg/communication_avoiding.h"
#include "cfdx/core/linalg/chebyshev_smoother.h"
#include "cfdx/physics/advanced_convergence.h"
#include "cfdx/physics/local_time_stepping.h"
#include "cfdx/physics/low_mach.h"
#include "cfdx/physics/radiation_models.h"
#include "cfdx/thermodynamics/thermo_cache.h"
#include "cfdx/core/field/field.h"
#include <cmath>
#include <iostream>
int main() {
    using namespace cfdx;
    core::FunctionalLinearOperator op(3, [](const core::Vector& x, core::Vector& y) {
        y.resize(3);
        y(0)=2.0*x(0); y(1)=3.0*x(1); y(2)=4.0*x(2);
    });
    core::Vector rhs(3,1.0), x(3,0.0);
    core::ChebyshevSmoother smoother({2,1.0,4.0,1.0,4});
    smoother.apply(op,rhs,x);
    if (!(x.norm2()>0.0)) { std::cerr<<"performance check 1 failed\n"; return 1; }
    auto tol=physics::eisenstat_walker_tolerance(1e-4,1e-2);
    if (!(tol>0.0 && tol<1.0)) { std::cerr<<"performance check 2 failed tol="<<tol<<"\n"; return 2; }
    if (physics::choose_pressure_correctors(1e-2)<2) { std::cerr<<"performance check 3 failed\n"; return 3; }
    physics::RadiationModelSelector selector;
    if (selector.select(2.0,0.0,2.0)!=physics::RadiationApproximation::Rosseland) { std::cerr<<"performance check 4 failed\n"; return 4; }
    if (!(physics::rosseland_conductivity(1000.0,1.0)>0.0)) { std::cerr<<"performance check 5 failed\n"; return 5; }
    core::Vector b(3,2.0);
    auto red=core::fused_reduction(rhs,b);
    if (std::abs(red.dot-6.0)>1e-12) { std::cerr<<"performance check 6 failed dot="<<red.dot<<"\n"; return 6; }
    auto restart=core::choose_gmres_restart(30,0.5);
    if (restart>=30) { std::cerr<<"performance check 7 failed restart="<<restart<<"\n"; return 7; }
    thermodynamics::ThermoCache cache;
    thermodynamics::IdealGasThermoModel gas;
    cache.update(gas,{101325.0,101325.0},{300.0,310.0});
    if (!cache.valid || cache.state.size()!=2) { std::cerr<<"performance check 8 failed\n"; return 8; }
    core::memory::ReusePool pool;
    const auto aoff = pool.acquire(128);
    const auto boff = pool.acquire(64);
    (void)boff;
    pool.release(aoff);
    const auto coff = pool.acquire(32);
    if (coff != aoff || pool.allocated_bytes() != 96) { std::cerr<<"performance check 9 failed coff="<<coff<<" aoff="<<aoff<<" active="<<pool.allocated_bytes()<<"\n"; return 9; }

    core::memory::FieldLifetime fa{"a",64,0,2,core::memory::Residency::Ephemeral};
    core::memory::FieldLifetime fb{"b",64,2,4,core::memory::Residency::Ephemeral};
    if (!core::memory::reusable(fa,fb)) { std::cerr<<"performance check 10 failed\n"; return 10; }

    bool overlap = false;
    core::HaloOverlap schedule([&](){ overlap = true; }, [](){}, [](){}, [](){});
    schedule.execute();
    if (!overlap) { std::cerr<<"performance check 11 failed\n"; return 11; }

    std::cout << "performance runtime: PASS\\n";
    return 0;
}
