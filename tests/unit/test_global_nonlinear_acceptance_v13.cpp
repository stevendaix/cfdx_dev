#include "cfdx/physics/global_nonlinear_acceptance_v13.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace cfdx::physics;

static GlobalNonlinearResidualV8 certified(double momentum=1e-3, double continuity=1e-4,
                                          double k=1e-5, double omega=1e-4,
                                          double closure=1e-7) {
    GlobalNonlinearResidualV8 x;
    x.algorithm=PressureVelocityAlgorithm::COUPLED;
    x.schur=CoupledSchurModel::PCD;
    x.momentum.max_active_component_relative=momentum;
    x.momentum.global_relative=momentum;
    x.continuity=continuity;
    x.k=k;
    x.omega=omega;
    x.eddy_viscosity_closure=closure;
    x.all_finite=true;
    x.sst=true;
    x.bounded=true;
    x.sst_transport_committed=true;
    return x;
}

static void require(bool yes, const char* msg) {
    if (!yes) throw std::runtime_error(msg);
}

int main() {
    try {
        ConvergenceCriteria gates;
        gates.relative_tolerance=1e-9;
        gates.continuity_tolerance=1e-9;
        const auto old=certified();
        auto next=certified(0.6e-3,0.7e-4,0.8e-5,0.9e-4,1e-7);
        require(evaluate_global_nonlinear_trial_v13(old,next,gates).allowed,"decrease rejected");
        next=certified(1e-1,0.7e-4,0.8e-5,0.9e-4,1e-7);
        auto r=evaluate_global_nonlinear_trial_v13(old,next,gates);
        require(!r.allowed && r.failed_equation=="U_component","momentum growth not caught");
        next=certified(0.6e-3,0.7e-4,0.8e-5,1e-1,1e-7);
        r=evaluate_global_nonlinear_trial_v13(old,next,gates);
        require(!r.allowed && r.failed_equation=="omega","omega growth not caught");
        next=certified(0.6e-3,0.7e-4,1e-1,0.9e-4,1e-7);
        r=evaluate_global_nonlinear_trial_v13(old,next,gates);
        require(!r.allowed && r.failed_equation=="k","k growth not caught");
        next=certified(0.6e-3,1e-1,0.8e-5,0.9e-4,1e-7);
        r=evaluate_global_nonlinear_trial_v13(old,next,gates);
        require(!r.allowed && r.failed_equation=="p_continuity","continuity growth not caught");
        next=certified(0.6e-3,0.7e-4,0.8e-5,0.9e-4,1e-2);
        r=evaluate_global_nonlinear_trial_v13(old,next,gates);
        require(!r.allowed && r.failed_equation=="nu_t_closure","eddy-viscosity mismatch not caught");
        next=certified();next.sst_transport_committed=false;
        bool threw=false;
        try {(void)evaluate_global_nonlinear_trial_v13(old,next,gates);}
        catch (const std::invalid_argument&) {threw=true;}
        require(threw,"uncommitted SST accepted");
        next=certified();next.sst=false;
        threw=false;
        try {(void)evaluate_global_nonlinear_trial_v13(old,next,gates);}
        catch (const std::invalid_argument&) {threw=true;}
        require(threw,"missing SST accepted");
        std::cout << "GLOBAL_NONLINEAR_ACCEPTANCE_V13_PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "GLOBAL_NONLINEAR_ACCEPTANCE_V13_FAIL: " << error.what() << "\n";
        return 1;
    }
}
