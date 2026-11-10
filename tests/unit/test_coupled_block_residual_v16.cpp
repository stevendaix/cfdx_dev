#include "cfdx/physics/coupled_block_residual_v16.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

#define REQUIRE(v) do { if (!(v)) throw std::runtime_error("V16 block certificate failed: " #v); } while (false)

int main() {
    using namespace cfdx::core;
    using cfdx::physics::independent_coupled_block_residual_v16;
    SparseMatrix A(4, 4);
    for (std::size_t i = 0; i < 4; ++i) A.push_back(i, i, 1.0);
    A.finalize();
    Vector b(4), x(4);
    b(0)=1.; b(1)=2.; b(2)=3.; b(3)=4.;
    x(0)=1.; x(1)=2.; x(2)=0.; x(3)=4.;
    const auto pressure_bad=independent_coupled_block_residual_v16(A,b,x,2);
    REQUIRE(pressure_bad.momentum.absolute_l2 == 0.);
    REQUIRE(pressure_bad.continuity_and_gauge.absolute_l2 == 3.);
    REQUIRE(pressure_bad.continuity_and_gauge.worst_row == 2);
    REQUIRE(std::abs(pressure_bad.relative_to_global_rhs()-3./std::sqrt(30.))<1.e-14);
    REQUIRE(std::abs(pressure_bad.continuity_and_gauge.relative_to_block_rhs()-0.6L)<1.e-14);
    REQUIRE(pressure_bad.momentum_fraction_of_squared_defect() == 0.);
    x(2)=3.; x(1)=1.;
    const auto momentum_bad=independent_coupled_block_residual_v16(A,b,x,2);
    REQUIRE(momentum_bad.continuity_and_gauge.absolute_l2 == 0.);
    REQUIRE(momentum_bad.momentum.absolute_l2 == 1.);
    REQUIRE(momentum_bad.momentum.worst_row == 1);
    REQUIRE(momentum_bad.momentum_fraction_of_squared_defect() == 1.);
    x=b;
    REQUIRE(independent_coupled_block_residual_v16(A,b,x,2).global_l2 == 0.);
    bool rejected=false;
    try {independent_coupled_block_residual_v16(A,b,x,4);} catch(const std::invalid_argument&){rejected=true;}
    REQUIRE(rejected);
    x(0)=std::numeric_limits<double>::quiet_NaN();
    rejected=false;
    try {independent_coupled_block_residual_v16(A,b,x,2);} catch(const std::runtime_error&){rejected=true;}
    REQUIRE(rejected);
    std::cout << "V16 true CSR residual split: independent momentum and continuity PASS\n";
}
