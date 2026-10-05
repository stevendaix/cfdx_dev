#include "cfdx/physics/pressure_velocity_algorithms.h"

#include <iostream>
#include <stdexcept>
#include <string>

using namespace cfdx::physics;

int main() {
    try {
        CoupledSchurModel model = CoupledSchurModel::PCD;
        if (!parse_coupled_schur_model("block_local", model) ||
            model != CoupledSchurModel::BlockLocal)
            throw std::runtime_error("block_local parsing failed");
        if (!parse_coupled_schur_model("pcd", model) ||
            model != CoupledSchurModel::PCD)
            throw std::runtime_error("pcd parsing failed");
        if (!parse_coupled_schur_model("lsc", model) ||
            model != CoupledSchurModel::LSC)
            throw std::runtime_error("lsc parsing failed");
        if (!parse_coupled_schur_model("bfbt", model) ||
            model != CoupledSchurModel::BFBT)
            throw std::runtime_error("bfbt parsing failed");
        if (!parse_coupled_schur_model("simple", model) ||
            model != CoupledSchurModel::SIMPLE)
            throw std::runtime_error("simple parsing failed");
        if (!parse_coupled_schur_model("simplec", model) ||
            model != CoupledSchurModel::SIMPLEC)
            throw std::runtime_error("simplec parsing failed");
        if (parse_coupled_schur_model("unknown", model))
            throw std::runtime_error("unknown Schur model was accepted");

        CouplingControls controls;
        if (controls.schur_model != CoupledSchurModel::PCD)
            throw std::runtime_error("existing PCD default changed unexpectedly");
        controls.schur_model = CoupledSchurModel::PCD;
        if (std::string(to_string(controls.schur_model)) != "pcd")
            throw std::runtime_error("PCD model string is not stable");
        controls.schur_model = CoupledSchurModel::LSC;
        if (std::string(to_string(controls.schur_model)) != "lsc")
            throw std::runtime_error("LSC model string is not stable");
        controls.schur_model = CoupledSchurModel::BFBT;
        if (std::string(to_string(controls.schur_model)) != "bfbt")
            throw std::runtime_error("BFBt model string is not stable");
        controls.schur_model = CoupledSchurModel::SIMPLE;
        if (std::string(to_string(controls.schur_model)) != "simple")
            throw std::runtime_error("SIMPLE model string is not stable");
        controls.schur_model = CoupledSchurModel::SIMPLEC;
        if (std::string(to_string(controls.schur_model)) != "simplec")
            throw std::runtime_error("SIMPLEC model string is not stable");

        std::cout << "coupled Schur selection contract PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "coupled Schur selection contract FAIL: " << e.what() << "\n";
        return 1;
    }
}
