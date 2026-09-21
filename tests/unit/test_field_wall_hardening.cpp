#include "cfdx/core/field/field.h"
#include "cfdx/core/field/storage.h"
#include "cfdx/physics/wall_functions.h"
#include "common/test_harness.h"
#include <cstdint>
#include <stdexcept>

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("field_storage_is_single_aligned_block", [] {
        Field<double,Location::CELL> f(16,"p","Pa",3);
        const auto p0=reinterpret_cast<std::uintptr_t>(f.component_data(0));
        const auto p1=reinterpret_cast<std::uintptr_t>(f.component_data(1));
        EXPECT_TRUE(p0%64==0);
        EXPECT_TRUE(p1%64==0);
        EXPECT_TRUE(p1-p0==16*sizeof(double));
    });
    run_case("vec3_field_has_real_soa_components", [] {
        Vec3CellField u(4,"U","m/s");
        u.set(0,1.0,2.0,3.0);
        EXPECT_NEAR(u.component_data(0)[0],1.0,1e-12);
        EXPECT_NEAR(u.component_data(1)[0],2.0,1e-12);
        EXPECT_NEAR(u.component_data(2)[0],3.0,1e-12);
        EXPECT_TRUE(reinterpret_cast<std::uintptr_t>(u.component_data(0))%64==0);
    });
    run_case("storage_core_does_not_silently_claim_gpu_sync", [] {
        StorageHandle h;
        EXPECT_THROW(h.sync_host_to_device(),std::runtime_error);
        EXPECT_THROW(h.sync_device_to_host(),std::runtime_error);
    });
    run_case("wall_function_reference_relations", [] {
        const double nu=1e-5, utau=0.5, y=0.001;
        const double yp=cfdx::cfdx::physics::wall::y_plus(y,utau,nu);
        EXPECT_NEAR(yp,50.0,1e-12);
        EXPECT_TRUE(cfdx::physics::wall::u_plus_log(yp)>0.0);
        EXPECT_TRUE(cfdx::physics::wall::turbulent_viscosity_log(y,utau,nu)>0.0);
        EXPECT_TRUE(cfdx::physics::wall::omega_log(y,utau)>0.0);
    });
    return run_all();
}
