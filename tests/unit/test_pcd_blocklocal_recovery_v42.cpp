#include "cfdx/physics/pcd_blocklocal_recovery_v42.h"
#include "cfdx/physics/pressure_velocity_algorithms.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#define CHECK(x) do { if(!(x)) throw std::runtime_error("V42 PCD failed: " #x); } while(false)
int main() {
    using namespace cfdx::core;
    using namespace cfdx::physics;
    SparseMatrix A(4,4);
    for(std::size_t i=0;i<4;++i) A.push_back(i,i,1.0);
    A.finalize();
    Vector b(4),candidate(4);
    for(std::size_t i=0;i<4;++i) b(i)=candidate(i)=static_cast<double>(i+1);
    auto d=certify_pcd_blocklocal_recovery_v42(A,b,candidate,3,1e-12,true);
    CHECK(d.adopt && d.true_relative==0.0);
    d=certify_pcd_blocklocal_recovery_v42(A,b,candidate,3,1e-12,false);
    CHECK(!d.adopt && d.true_relative==0.0);
    candidate(3)+=1e-7;
    d=certify_pcd_blocklocal_recovery_v42(A,b,candidate,3,1e-12,true);
    CHECK(!d.adopt && d.true_relative>1e-12 &&
          d.blocks.continuity_and_gauge.absolute_l2>0.0);
    CouplingControls controls;
    controls.schur_model=CoupledSchurModel::PCD;
    CHECK(!controls.pcd_blocklocal_recovery_v42_experimental);
    controls.pcd_blocklocal_recovery_v42_experimental=true;
    validate_coupling_controls(controls);
    controls.coupled_monolithic_ilut_fallback_v16_experimental=true;
    bool rejected=false;
    try {validate_coupling_controls(controls);}
    catch(const std::invalid_argument&) { rejected=true; }
    CHECK(rejected);
    controls.coupled_monolithic_ilut_fallback_v16_experimental=false;
    controls.schur_model=CoupledSchurModel::BFBT;
    rejected=false;
    try {validate_coupling_controls(controls);}
    catch(const std::invalid_argument&) { rejected=true; }
    CHECK(rejected);
    controls.schur_model=CoupledSchurModel::PCD;
    controls.pcd_blocklocal_recovery_v42_experimental=false;
    controls.pcd_blocklocal_preflight_iterations_v42=10;
    rejected=false;
    try {validate_coupling_controls(controls);}
    catch(const std::invalid_argument&) {rejected=true;}
    CHECK(rejected);
    controls.pcd_blocklocal_recovery_v42_experimental=true;
    controls.pcd_blocklocal_preflight_iterations_v42=controls.coupled_max_iterations+1;
    rejected=false;
    try {validate_coupling_controls(controls);}
    catch(const std::invalid_argument&) {rejected=true;}
    CHECK(rejected);
    controls.pcd_blocklocal_preflight_iterations_v42=40;
    validate_coupling_controls(controls);
    std::cout<<"PCD_BLOCKLOCAL_RECOVERY_V42_PASS independent_original_A_b=1 fallback_off_default=1\n";
}
