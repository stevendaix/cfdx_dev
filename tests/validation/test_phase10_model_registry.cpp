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
        cfdx::core::Mesh mesh;
        mesh.points().resize(8);
        mesh.points().set(0,0,0,0); mesh.points().set(1,1,0,0);
        mesh.points().set(2,1,1,0); mesh.points().set(3,0,1,0);
        mesh.points().set(4,0,0,1); mesh.points().set(5,1,0,1);
        mesh.points().set(6,1,1,1); mesh.points().set(7,0,1,1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> k(1,"k","m2/s2",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> w(1,"omega","1/s",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> s(1,"strain","1/s",1);
        k(0)=0.1; w(0)=10.0; s(0)=5.0;
        cfdx::core::Field<double,cfdx::core::Location::CELL> second(1,"second","",1);
        second(0)=10.0;
        auto algebraic=solve_algebraic_turbulence(mesh,k,second,s,c);
        if (!algebraic.converged || k(0)<=0.0 || second(0)<=0.0)
            throw std::runtime_error("algebraic turbulence driver failed");

        std::cout << "PHASE10_MODEL_REGISTRY: PASS (" << models.size() << " models)\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "PHASE10_MODEL_REGISTRY: FAIL: " << e.what() << "\n";
        return 1;
    }
}
