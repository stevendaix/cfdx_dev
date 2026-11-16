#include "cfdx/physics/sst_production_coupling.h"
#include "cfdx/physics/sst_case_contract.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

static void require(bool b,const char* m) {if(!b) throw std::runtime_error(m);}
int main() {
    using namespace cfdx;
    io::CaseSetup setup;
    io::BoundarySpec wall;
    wall.patch_name="downstreamWall";
    wall.type=io::BCType::WALL;
    wall.sst_scalar_conditions={
        {"k",io::SSTScalarBoundaryType::WALL_KQ_R,0.0},
        {"omega",io::SSTScalarBoundaryType::WALL_OMEGA_SST,0.0}
    };
    setup.boundary_conditions={wall};
    physics::validate_sst_scalar_boundary_contract(setup);
    const auto runtime=physics::make_sst_runtime_boundary_contract(setup);
    require(runtime.k.at("downstreamWall").type==physics::ScalarBoundaryType::ZERO_GRADIENT,
        "kqR must be a zero-gradient boundary");
    require(runtime.wall_omega_sst_patches.size()==1,
        "omega wall obligation must remain explicit");
    // The old resolved-wall k=0 choice has different discrete physics.
    setup.boundary_conditions[0].sst_scalar_conditions[0].type=
        io::SSTScalarBoundaryType::WALL_K_ZERO;
    physics::validate_sst_scalar_boundary_contract(setup);
    require(physics::make_sst_runtime_boundary_contract(setup)
            .k.at("downstreamWall").type==physics::ScalarBoundaryType::FIXED_VALUE,
        "WALL_K_ZERO unexpectedly changed");
    setup.boundary_conditions[0].sst_scalar_conditions[0].type=
        io::SSTScalarBoundaryType::ZERO_GRADIENT;
    bool refused=false;
    try { physics::validate_sst_scalar_boundary_contract(setup); }
    catch(const std::invalid_argument&){refused=true;}
    require(refused,"generic zero gradient must not masquerade as kqR");
    std::cout<<"SST V16 kqR zero-gradient distinguished from resolved k=0\n";
}
