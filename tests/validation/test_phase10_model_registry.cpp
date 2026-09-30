#include "cfdx/physics/turbulence_models.h"
#include "cfdx/physics/turbulence_transport.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::physics;

int main() {
    try {
        const std::vector<AdvancedTurbulenceModel> models = {
            AdvancedTurbulenceModel::LAMINAR,
            AdvancedTurbulenceModel::KEPSILON,
            AdvancedTurbulenceModel::RNG_KEPSILON,
            AdvancedTurbulenceModel::REALIZABLE_KEPSILON,
            AdvancedTurbulenceModel::KOMEGA,
            AdvancedTurbulenceModel::SST,
            AdvancedTurbulenceModel::SPALART_ALLMARAS,
            AdvancedTurbulenceModel::SMAGORINSKY,
            AdvancedTurbulenceModel::WALE,
            AdvancedTurbulenceModel::DYNAMIC_KEQN,
            AdvancedTurbulenceModel::DES,
            AdvancedTurbulenceModel::DDES,
            AdvancedTurbulenceModel::IDDES
        };

        TurbulenceModelCoefficients coeff;
        validate_turbulence_model_coefficients(coeff);
        {
            TurbulenceCorrectionControls corr;
            TurbulenceInvariants inv{2.0,1.0,0.5,0.2,0.01};
            corr.rotation_curvature=true;
            corr.compressibility=true;
            corr.kato_launder=true;
            corr.production_limiter=true;
            corr.production_limit=100.0;
            if(!(rotation_curvature_factor(inv,corr) > 1.0)) throw std::runtime_error("rotation correction failed");
            if(!(compressibility_factor(inv,corr,10.0) < 1.0)) throw std::runtime_error("compressibility correction failed");
            if(!(kato_launder_production_factor(inv,corr) > 0.0)) throw std::runtime_error("Kato-Launder correction failed");
            if(!(corrected_turbulence_production(1.0,inv,corr,10.0) > 0.0)) throw std::runtime_error("production correction failed");
            if(std::abs(roughness_factor(0.001,0.01,corr)-1.0)>1e-12) throw std::runtime_error("roughness correction failed");

            // Composition test: independent correction layers must both remain
            // active when applied through the common production hook.
            TurbulenceCorrectionControls composed;
            composed.production_limiter=false;
            composed.kato_launder=true;
            composed.kato_coefficient=1.0;
            composed.roughness=true;
            composed.roughness_height=0.01;
            composed.roughness_coefficient=2.0;
            const double base=2.0;
            const double kato_only=corrected_turbulence_production(
                base, inv, composed, 10.0);
            const double roughness_only_factor=roughness_factor(
                composed.roughness_height, inv.wall_distance, composed);
            const double expected=(base*roughness_only_factor) +
                kato_launder_production_factor(inv, composed);
            if(std::abs(kato_only-expected)>1e-12)
                throw std::runtime_error("correction composition failed");
            corr.qcr=true;
            if(!(qcr_stress_factor(2.0,1.0,corr) > 1.0)) throw std::runtime_error("QCR correction failed");
            if(!(wale_eddy_viscosity({4.0,4.0,0.1},0.325)>0.0)) throw std::runtime_error("WALE kernel failed");
            if(!(dynamic_les_coefficient(2.0,1.0,2.0) >= 0.0)) throw std::runtime_error("dynamic LES kernel failed");
        }

        for (const auto model : models) {
            const auto descriptor = turbulence_model_descriptor(model);
            const auto& capability = turbulence_capability(model);
            if (descriptor.model != model || descriptor.implementation != implementation_kind(model))
                throw std::runtime_error("model descriptor mismatch");
            if (descriptor.status != capability.status)
                throw std::runtime_error("model capability status mismatch");
            if (capability.key == nullptr || capability.label == nullptr ||
                capability.required_fields == nullptr || capability.missing == nullptr)
                throw std::runtime_error("model capability metadata is incomplete");
            if (capability.status == TurbulenceImplementationStatus::SOLVER_READY &&
                capability.model != AdvancedTurbulenceModel::LAMINAR &&
                std::string(capability.missing).empty())
                throw std::runtime_error("solver-ready model has no explicit remaining qualification gap");
        }

        TurbulenceTransportControls c;
        for (const auto model : {
            TurbulenceModel::LAMINAR, TurbulenceModel::KEPSILON,
            TurbulenceModel::RNG_KEPSILON, TurbulenceModel::REALIZABLE_KEPSILON,
            TurbulenceModel::KOMEGA,
            TurbulenceModel::SST, TurbulenceModel::SPALART_ALLMARAS,
            TurbulenceModel::SMAGORINSKY, TurbulenceModel::WALE,
            TurbulenceModel::DYNAMIC_KEQN, TurbulenceModel::DES,
            TurbulenceModel::DDES, TurbulenceModel::IDDES}) {
            c.model = model;
            validate_turbulence_controls(c);
            double nut = 0.0;
            if (model == TurbulenceModel::WALE) {
                // WALE requires the tensor invariants of the resolved velocity
                // gradient; the qualification test must provide them explicitly.
                nut = turbulence_nu_t(0.1, 0.02, 10.0, 0.01, c, 1.0,
                                      1.0, 0.0, 0.0, 1.0, 2.0, 0.5, 0.5);
            } else {
                nut = turbulence_nu_t(0.1, 0.02, 10.0, 0.01, c, 1.0);
            }
            if (!std::isfinite(nut) || nut < 0.0)
                throw std::runtime_error("non-physical turbulent viscosity");
        }

        if (!turbulence_solver_ready(AdvancedTurbulenceModel::SST) ||
            !turbulence_solver_ready(AdvancedTurbulenceModel::SPALART_ALLMARAS))
            throw std::runtime_error("solver-ready RANS models are not registered");
        if (turbulence_solver_ready(AdvancedTurbulenceModel::DES) ||
            turbulence_capability(AdvancedTurbulenceModel::DES).status != TurbulenceImplementationStatus::KERNEL_ONLY)
            throw std::runtime_error("DES capability status is incorrect");
        if (turbulence_capability(AdvancedTurbulenceModel::DYNAMIC_KEQN).status != TurbulenceImplementationStatus::KERNEL_ONLY)
            throw std::runtime_error("dynamic LES capability status is incorrect");
        if (!has_transport_equation(AdvancedTurbulenceModel::REALIZABLE_KEPSILON))
            throw std::runtime_error("realizable k-epsilon transport is not registered");

        // The low-level closure also accepts an explicit Cmu for specialized callers.
        constexpr double realizable_cmu = 0.09;
        const double realizable = realizable_kepsilon_eddy_viscosity(0.1, 0.02, realizable_cmu);
        c.model = TurbulenceModel::RNG_KEPSILON;
        const double rng = turbulence_nu_t(0.1, 0.02, 10.0, 0.01, c);
        if (!(std::isfinite(realizable) && std::isfinite(rng) &&
              realizable > 0.0 && rng > 0.0))
            throw std::runtime_error("RNG/Realizable regression failed");

        c.model = TurbulenceModel::SST;
        if (turbulence_nu_t(0.1, 10.0, 5.0, 0.01, c) <= 0.0)
            throw std::runtime_error("SST closure returned non-positive viscosity");

        c.model = TurbulenceModel::SPALART_ALLMARAS;
        if (turbulence_nu_t(0.1, 0.01, 5.0, 0.01, c) <= 0.0)
            throw std::runtime_error("SA closure returned non-positive viscosity");

        std::cout << "PHASE10_MODEL_REGISTRY: PASS (" << models.size() << " models)\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "PHASE10_MODEL_REGISTRY: FAIL: " << e.what() << "\n";
        return 1;
    }
}
