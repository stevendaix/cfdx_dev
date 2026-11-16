#include "cfdx/physics/state_only_rhie_chow_v49.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void require(bool ok) {
    if (!ok) throw std::runtime_error("V49 state-only Rhie-Chow regression");
}
int main() {
    // This synthetic face/diagonal model has an EXACT fixed point:
    // A(phi)=(2+phi), rAU=1/A, RC(phi)=0.2+0.25*rAU.
    // It tests nonlinear dependence of A on convective face flux.
    auto diag=[](const std::vector<double>& phi) {
        std::array<std::vector<double>,3> result;
        for (auto& d:result) d={1.0/(2.0+phi[0]),1.0/(3.0+phi[1])};
        return result;
    };
    auto reconstruct=[](const std::array<std::vector<double>,3>& d) {
        return std::vector<double>{0.2+0.25*d[0][0],-0.15+0.4*d[1][1]};
    };
    const std::vector<double> seed{0.2,-0.15};
    const auto a=evaluate_state_only_rhie_chow_v49(seed,2,diag,reconstruct);
    const auto b=evaluate_state_only_rhie_chow_v49(seed,2,diag,reconstruct);
    require(a.converged && b.converged && a.iterations==b.iterations);
    require(a.flux==b.flux && a.inverse_diagonal==b.inverse_diagonal);
    require(a.relative_linf<=1e-10);
    require(normalized_flux_fixed_point_defect_v49(
                a.flux,reconstruct(diag(a.flux)))<=1e-10);
    // The closure metric must be relative to the actual face mass flow,
    // not an artificial 1 kg/s scale that hides small BFS face fluxes.
    require(normalized_flux_fixed_point_defect_v49(
                {0.001},{0.0011})>0.09);
    // No previous accepted state or old rAU argument exists in the API.
    // A genuine non-fixed-point mapping must not receive a certificate.
    auto grow=[](const std::array<std::vector<double>,3>& d) {
        return std::vector<double>{2.0/d[0][0],2.0/d[0][1]};
    };
    const auto failure=evaluate_state_only_rhie_chow_v49(
        seed,2,diag,grow,5,1e-12);
    require(!failure.converged && failure.flux.empty() &&
            failure.iterations==5 && std::isfinite(failure.relative_linf));
    bool threw=false;
    try {evaluate_state_only_rhie_chow_v49(seed,2,diag,reconstruct,0);}
    catch(const std::invalid_argument&) {threw=true;}
    require(threw);
    threw=false;
    try {evaluate_state_only_rhie_chow_v49(seed,2,diag,reconstruct,5,1e-10,0.0);}
    catch(const std::invalid_argument&) {threw=true;}
    require(threw);
    threw=false;
    try {std::vector<double> nan=seed;
         nan[0]=std::numeric_limits<double>::quiet_NaN();
         evaluate_state_only_rhie_chow_v49(nan,2,diag,reconstruct);}
    catch(const std::invalid_argument&) {threw=true;}
    require(threw);
    threw=false;
    try {evaluate_state_only_rhie_chow_v49(seed,3,diag,reconstruct);}
    catch(const std::invalid_argument&) {threw=true;}
    require(threw);
    threw=false;
    auto invalid_diag=[](const std::vector<double>&) {
        std::array<std::vector<double>,3> out;
        for(auto& d:out) d={0.,1.};
        return out;
    };
    try {evaluate_state_only_rhie_chow_v49(seed,2,invalid_diag,reconstruct);}
    catch(const std::invalid_argument&) {threw=true;}
    require(threw);
}
