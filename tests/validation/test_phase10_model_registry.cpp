#include "cfdx/physics/turbulence_transport.h"
#include "cfdx/physics/turbulence_solver.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace cfdx::physics;

int main() {
    try {
        const std::vector<TurbulenceModel> models = {
            TurbulenceModel::LAMINAR,
            TurbulenceModel::KEPSILON,
            TurbulenceModel::RNG_KEPSILON,
            TurbulenceModel::REALIZABLE_KEPSILON,
            TurbulenceModel::KOMEGA,
            TurbulenceModel::SST,
            TurbulenceModel::SPALART_ALLMARAS,
            TurbulenceModel::SMAGORINSKY,
            TurbulenceModel::WALE,
            TurbulenceModel::DES,
            TurbulenceModel::DDES,
            TurbulenceModel::IDDES
        };

        TurbulenceTransportControls c;
        for (const auto model : models) {
            c.model=model;
            validate_turbulence_controls(c);
            const double nut=turbulence_nu_t(0.1,0.02,10.0,0.01,c);
            if (!std::isfinite(nut) || nut < 0.0)
                throw std::runtime_error("non-physical turbulent viscosity");
        }

        c.model=TurbulenceModel::REALIZABLE_KEPSILON;
        const double realizable=turbulence_nu_t(0.1,0.02,10.0,0.01,c);
        c.model=TurbulenceModel::RNG_KEPSILON;
        const double rng=turbulence_nu_t(0.1,0.02,10.0,0.01,c);
        if (!(std::isfinite(realizable) && std::isfinite(rng) &&
              realizable > 0.0 && rng > 0.0))
            throw std::runtime_error("RNG/Realizable regression failed");

        c.model=TurbulenceModel::SST;
        if (turbulence_nu_t(0.1,10.0,5.0,0.01,c) <= 0.0)
            throw std::runtime_error("SST closure returned non-positive viscosity");
        c.model=TurbulenceModel::SPALART_ALLMARAS;
        if (turbulence_nu_t(0.1,0.01,5.0,0.01,c) <= 0.0)
            throw std::runtime_error("SA closure returned non-positive viscosity");

        std::cout << "PHASE10_MODEL_REGISTRY: PASS (" << models.size() << " models)\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "PHASE10_MODEL_REGISTRY: FAIL: " << e.what() << "\n";
        return 1;
    }
}
