#include "cfdx/physics/sst_implicit_spatial_rescue_v82.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace cfdx::physics;
static void require(bool c,const char* why) {
    if (!c) throw std::runtime_error(why);
}
int main() {
    try {
        // Synthetic two-cell, two-field *coupled* positive FV-like oracle.
        // Rk_i = k_i - 1 + 0.15*(omega_i - 2),
        // Rw_i = omega_i - 2 + 0.2*(k_i - 1).
        SSTV82Evaluator eval=[](const auto& s) {
            SSTV82PhysicalTrial t;t.state=s;
            t.k_local=t.omega_local=0.;
            for(std::size_t i=0;i<2;++i) {
                t.k_local=std::max(t.k_local,std::abs(s[i]-1.+0.15*(s[2+i]-2.)));
                t.omega_local=std::max(t.omega_local,std::abs(s[2+i]-2.+0.2*(s[i]-1.)));
            }
            t.k_row=t.k_local;t.omega_row=t.omega_local;
            return t;
        };
        const auto baseline=eval({1.5,1.5,3.,3.});
        auto r=sst_implicit_spatial_search_v82(baseline,
                    {-0.3,-0.3,-0.3,-0.3}, {},2,1e-8,1e-8,eval);
        require(r.accepted && r.merit_after<r.merit_before,"coupled descent");
        require(r.k_after<=r.k_before && r.omega_after<=r.omega_before,
                "individual physical gates");
        require(r.state.size()==4 && r.state[0]>0 && r.state[2]>0,"positivity");
        auto rejected=sst_implicit_spatial_search_v82(baseline,
                            {1.,1.,1.,1.},{},2,1e-8,1e-8,eval);
        require(!rejected.accepted && rejected.state.empty(),"reject ascent");
        // An ascent in Picard is allowed only as a private continuation,
        // never silently relabeled as an accepted physical descent.
        auto fall=sst_implicit_spatial_search_v82(baseline,
                    {1.,1.,1.,1.},{},2,1e-8,1e-8,eval,
                    {1.8,1.8,3.3,3.3});
        require(!fall.accepted && fall.legacy_fallback_available &&
                fall.legacy_fallback_merit>fall.merit_before,
                "fallback stays explicitly uncertified");
        // The global max-norm cannot worsen, even if one component falls.
        SSTV82Evaluator skewed=[](const auto& s) {
            SSTV82PhysicalTrial t;
            t.state=s;
            t.k_local=s[0];t.omega_local=2.-s[0];
            t.k_row=t.k_local;t.omega_row=t.omega_local;
            return t;
        };
        auto dual=sst_implicit_spatial_search_v82(
                    skewed({1.,1.,2.,2.}),{-0.2,0.,0.,0.},{},2,1e-6,1e-6,skewed);
        require(!dual.accepted,"reject worsening global max norm");
        // Unlike an over-restrictive per-equation Pareto filter, permit
        // a physical global merit reduction with a temporary k tradeoff.
        SSTV82Evaluator stiffness=[](const auto& s) {
            SSTV82PhysicalTrial t;t.state=s;
            t.k_local=s[0];t.omega_local=10.0*s[2];
            t.k_row=t.k_local;t.omega_row=t.omega_local;
            return t;
        };
        const auto stiff_start=stiffness({0.1,0.1,1.,1.});
        auto tradeoff=sst_implicit_spatial_search_v82(stiff_start,{}, {},2,
            1e-8,1e-8,stiffness,{0.2,0.1,0.3,1.});
        require(tradeoff.accepted && tradeoff.selected=="picard_full" &&
                tradeoff.k_after>tradeoff.k_before &&
                tradeoff.omega_after<tradeoff.omega_before &&
                tradeoff.merit_after<tradeoff.merit_before,
                "physical global merit must allow safe component tradeoff");
        // Even if a callback asserts tiny residuals, a NaN or a negative
        // post-wall-projection field makes it inadmissible.
        SSTV82Evaluator bad=[](const auto& s) {
            SSTV82PhysicalTrial t;t.state=s;t.state[2]=-1.;
            t.k_local=t.omega_local=t.k_row=t.omega_row=0.;return t;
        };
        auto invalid=sst_implicit_spatial_search_v82(baseline,
                    {-0.3,-0.3,-0.3,-0.3},{},2,1e-8,1e-8,bad);
        require(!invalid.accepted,"reject invalid wall state");
        SSTV82Evaluator throws=[](const auto&) -> SSTV82PhysicalTrial {
            throw std::runtime_error("nonlinear assembly failed");
        };
        auto thrown=sst_implicit_spatial_search_v82(baseline,
                    {-0.3,-0.3,-0.3,-0.3},{},2,1e-8,1e-8,throws);
        require(!thrown.accepted,"reject assembly exceptions");
        bool bad_input=false;
        try {(void)sst_implicit_spatial_search_v82(baseline,{1.,2.},{},2,0.,1e-8,eval);}
        catch (const std::invalid_argument&) {bad_input=true;}
        require(bad_input,"bad floor fail closed");
        // Invalid direction ignored, with no fallback to false convergence.
        auto bad_direction=sst_implicit_spatial_search_v82(baseline,
             {std::numeric_limits<double>::quiet_NaN(),0.,0.,0.},{},2,1e-8,1e-8,eval);
        require(!bad_direction.accepted,"bad direction fail closed");
        std::cout<<"test_sst_implicit_spatial_v82 PASS\n";
        return EXIT_SUCCESS;
    } catch(const std::exception& e) {
        std::cerr<<"test_sst_implicit_spatial_v82 FAIL: "<<e.what()<<'\n';
        return EXIT_FAILURE;
    }
}
