#include "cfdx/physics/production_buoyancy_stage.h"
#include "common/test_harness.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace cfdx::physics;
using namespace cfdx::core;
using namespace cfdx::testing;

int main()
{
    const ProductionBuoyancyCoupling valid{0.002, 300.0, {0.0, 0.0, -9.81}};
    run_case("buoyancy_late_invalid_cell_preserves_existing_acceleration", [&] {
        Field<double,Location::CELL> T(2,"T","K",1);
        T(0)=310.0;
        T(1)=std::numeric_limits<double>::quiet_NaN();
        std::vector<Vec3> acc{{11, 12, 13}, {21, 22, 23}, {31, 32, 33}};
        EXPECT_THROW_WITH(build_cell_buoyancy_acceleration(acc,T,valid),
                          std::invalid_argument, "temperature must be finite");
        EXPECT_TRUE(acc.size()==3);
        EXPECT_NEAR(acc[0].x,11.0,0.0);
        EXPECT_NEAR(acc[1].z,23.0,0.0);
        EXPECT_NEAR(acc[2].z,33.0,0.0);
    });
    run_case("buoyancy_rejects_vector_temperature_field", [&] {
        Field<double,Location::CELL> T(1,"T","K",3);
        T.fill(310.0);
        std::vector<Vec3> acc{{12, 13, 14}};
        EXPECT_THROW_WITH(build_cell_buoyancy_acceleration(acc,T,valid),
                          std::invalid_argument, "temperature must be a scalar");
        EXPECT_NEAR(acc[0].x,12.0,0.0);
    });
    run_case("buoyancy_rejects_finite_inputs_with_overflowing_force", [&] {
        Field<double,Location::CELL> T(1,"T","K",1);
        T(0)=1.0e300;
        const ProductionBuoyancyCoupling extreme{1.0e100, 300.0, {0, 0, -9.81}};
        std::vector<Vec3> acc{{12, 13, 14}};
        EXPECT_THROW_WITH(build_cell_buoyancy_acceleration(acc,T,extreme),
                          std::runtime_error, "non-representable density perturbation");
        EXPECT_NEAR(acc[0].x,12.0,0.0);
    });
    run_case("buoyancy_rejects_nonfinite_gravity", [&] {
        Field<double,Location::CELL> T(1,"T","K",1);
        T(0)=310.0;
        auto bad=valid;
        bad.gravity.z=std::numeric_limits<double>::quiet_NaN();
        std::vector<Vec3> acc{{12, 13, 14}};
        EXPECT_THROW_WITH(build_cell_buoyancy_acceleration(acc,T,bad),
                          std::invalid_argument, "invalid Boussinesq coupling");
        EXPECT_NEAR(acc[0].x,12.0,0.0);
    });
    run_case("buoyancy_valid_hot_cold_symmetry_and_reference_zero", [&] {
        Field<double,Location::CELL> T(3,"T","K",1);
        T(0)=310.0;T(1)=290.0;T(2)=300.0;
        std::vector<Vec3> acc;
        build_cell_buoyancy_acceleration(acc,T,valid);
        EXPECT_TRUE(acc.size()==3);
        EXPECT_NEAR(acc[0].z,0.1962,1e-12);
        EXPECT_NEAR(acc[1].z,-0.1962,1e-12);
        EXPECT_NEAR(acc[2].z,0.0,0.0);
    });
    return run_all();
}
