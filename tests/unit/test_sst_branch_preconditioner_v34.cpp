#include "cfdx/physics/sst_branch_preconditioner_v34.h"
#include "cfdx/physics/turbulence_transport.h"
#include "common/test_harness.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
using namespace cfdx::testing;

static SSTV33SourceBranches smooth_branch() {
    SSTV33SourceBranches s;
    s.frozen_source_derivative_valid=true;
    return s;
}

int main() {
    run_case("v34_two_residual_fields_choose_distinct_worst_cells",[] {
        const std::vector<double> residual{1.,20.,3., 5.,4.,60.};
        const std::vector<double> scale(6,1.);
        const auto ids=sst_v34_hotspot_cells(residual,scale,3);
        EXPECT_TRUE(ids.size()==2 && ids[0]==1 && ids[1]==2);
        const auto tie=sst_v34_hotspot_cells({2.,2.,2.,3.,2.,1.},scale,3);
        EXPECT_TRUE(tie.size()==1 && tie[0]==0);
        EXPECT_TRUE(sst_v34_hotspot_cells(residual,scale,3,0).empty());
        EXPECT_THROW(sst_v34_hotspot_cells(residual,scale,4),std::invalid_argument);
        EXPECT_THROW(sst_v34_hotspot_cells({1.,2.,3.,4.,5.,NAN},scale,3),std::invalid_argument);
    });
    run_case("v34_replaces_two_local_blocks_using_true_full_residual",[] {
        const std::vector<SSTBlockPreconditionerV23::Block> raw(3,{3.,0.,0.,5.});
        const SSTBlockPreconditionerV23 approx(raw);
        const std::vector<double> x{2.,3.,4.,5.,6.,7.}, scales(6,1.);
        const SSTPhysicalResidualV22 exact=[](const auto& y) {
            std::vector<double> f(6);
            for (std::size_t i=0;i<3;++i) {
                const double k=y[i],w=y[3+i];
                f[i]=k*k+0.3*k*w;
                f[3+i]=w*w+0.2*k;
            }
            return f;
        };
        const auto correction=sst_v34_correct_hotspot_blocks(
            approx,{1,2},{smooth_branch(),smooth_branch()},x,scales,scales,exact);
        EXPECT_TRUE(correction.sampled_cells.size()==2);
        EXPECT_TRUE(correction.corrected_cells.size()==2);
        EXPECT_TRUE(correction.residual_evaluations==10);
        EXPECT_TRUE(correction.skipped_nonsmooth==0);
        const auto& b=correction.preconditioner.raw_blocks();
        EXPECT_NEAR(b[0][0],3.,0.);
        for (std::size_t i : {std::size_t{1},std::size_t{2}}) {
            const double k=x[i],w=x[3+i];
            EXPECT_NEAR(b[i][0],2*k*k+0.3*k*w,1e-7);
            EXPECT_NEAR(b[i][1],0.3*k*w,1e-7);
            EXPECT_NEAR(b[i][2],0.2*k,1e-7);
            EXPECT_NEAR(b[i][3],2*w*w,1e-6);
        }
        EXPECT_TRUE(correction.preconditioner.diagonal_fallbacks()==0);
        EXPECT_TRUE(approx.raw_blocks()[1][1]==0.);
    });
    run_case("v34_nonsmooth_kink_is_rejected_not_differentiated_through",[] {
        const SSTBlockPreconditionerV23 approx({{4.,0.,0.,4.}});
        const auto residual=[](const auto& x) {
            return std::vector<double>{std::max(0.,x[0]-2.),x[1]-1.};
        };
        const auto result=sst_v34_correct_hotspot_blocks(approx,{0},{smooth_branch()},
            {2.,3.},{1.,1.},{1.,1.},residual);
        EXPECT_TRUE(result.corrected_cells.empty());
        EXPECT_TRUE(result.skipped_nonsmooth==1);
        EXPECT_TRUE(result.maximum_one_sided_disagreement>0.8);
        EXPECT_TRUE(result.preconditioner.raw_blocks()[0]==approx.raw_blocks()[0]);
    });
    run_case("v34_known_branch_switch_is_skipped_without_evaluation",[] {
        const SSTBlockPreconditionerV23 approx({{2.,0.,0.,3.}});
        auto branch=smooth_branch();branch.close_to_production_switch=true;
        std::size_t evaluations=0;
        const auto residual=[&](const auto& x) {++evaluations;return x;};
        const auto result=sst_v34_correct_hotspot_blocks(approx,{0},{branch},
            {2.,3.},{1.,1.},{1.,1.},residual);
        EXPECT_TRUE(result.skipped_near_switch==1 && evaluations==0);
    });
    run_case("v34_singular_local_jacobian_keeps_previous_inverse",[] {
        const SSTBlockPreconditionerV23 approx({{2.,0.,0.,3.}});
        const auto residual=[](const auto& x) {
            return std::vector<double>{x[0]+x[1],2*x[0]+2*x[1]};
        };
        const auto result=sst_v34_correct_hotspot_blocks(approx,{0},{smooth_branch()},
            {2.,3.},{1.,1.},{1.,1.},residual);
        EXPECT_TRUE(result.skipped_singular==1);
        EXPECT_TRUE(result.corrected_cells.empty());
        EXPECT_TRUE(result.preconditioner.raw_blocks()[0]==approx.raw_blocks()[0]);
    });
    run_case("v34_nonfinite_malformed_and_duplicate_trials_fail_closed",[] {
        const SSTBlockPreconditionerV23 approx({{2.,0.,0.,3.}});
        const auto residual=[](const auto& x) {return x;};
        EXPECT_THROW(sst_v34_correct_hotspot_blocks(approx,{0,0},
            {smooth_branch(),smooth_branch()},{2.,3.},{1.,1.},{1.,1.},residual),
            std::invalid_argument);
        EXPECT_THROW(sst_v34_correct_hotspot_blocks(approx,{1},{smooth_branch()},
            {2.,3.},{1.,1.},{1.,1.},residual),std::invalid_argument);
        EXPECT_THROW(sst_v34_correct_hotspot_blocks(approx,{0},{smooth_branch()},
            {2.,3.},{1.,1.},{1.,1.},[](const auto&) {return std::vector<double>{0.};}),
            std::invalid_argument);
        EXPECT_THROW(sst_v34_correct_hotspot_blocks(approx,{0},{smooth_branch()},
            {2.,NAN},{1.,1.},{1.,1.},residual),std::invalid_argument);
    });
    run_case("v34_cross_cell_audit_captures_missing_nonlocal_k_omega_coupling",[] {
        const std::vector<double> x{2.,3.,5.,7.}, unit(4,1.);
        const SSTPhysicalResidualV22 residual=[](const auto& z) {
            return std::vector<double>{2*z[1]+3*z[3],z[1],
                                       4*z[1]+5*z[3],z[3]};
        };
        const auto audit=sst_v34_cross_cell_jacobian_audit(
            0,1,x,unit,unit,residual,2.,5.);
        EXPECT_TRUE(audit.residual_evaluations==4);
        EXPECT_NEAR(audit.approximate[0],6.,1e-14);
        EXPECT_NEAR(audit.approximate[1],0.,0.);
        EXPECT_NEAR(audit.actual[0],6.,1e-8);
        EXPECT_NEAR(audit.actual[1],21.,1e-8);
        EXPECT_NEAR(audit.actual[2],12.,1e-8);
        EXPECT_NEAR(audit.actual[3],35.,1e-7);
        EXPECT_TRUE(audit.maximum_relative_difference>0.5);
    });
    run_case("v34_cross_cell_diagnostic_parses_and_uses_frozen_CSR",[] {
        auto ids=sst_v34_parse_cross_pair("341:425",712);
        EXPECT_TRUE(ids[0]==341 && ids[1]==425);
        for (const auto* bad:{"341,425","341:341","1:2:3","-1:4","1:900"})
            EXPECT_THROW(sst_v34_parse_cross_pair(bad,712),std::invalid_argument);
        cfdx::core::SparseMatrix A(2,2);
        A.push_back(0,0,4.);A.push_back(0,1,-2.);
        A.push_back(1,0,-3.);A.push_back(1,1,5.);A.finalize();
        EXPECT_NEAR(sst_v34_frozen_offcell_coefficient(A,0,1),-2.,0.);
        EXPECT_NEAR(sst_v34_frozen_offcell_coefficient(A,1,0),-3.,0.);
        EXPECT_THROW(sst_v34_frozen_offcell_coefficient(A,0,0),std::invalid_argument);
    });
    run_case("v34_control_is_opt_in_and_requires_true_jvp_and_spatial_P",[] {
        TurbulenceTransportControls c;c.model=TurbulenceModel::SST;
        EXPECT_TRUE(!c.sst_branch_preconditioner_v34_experimental);
        c.sst_branch_preconditioner_v34_experimental=true;
        EXPECT_THROW(validate_turbulence_controls(c),std::invalid_argument);
        c.sst_jfnk_v22_experimental=true;
        EXPECT_THROW(validate_turbulence_controls(c),std::invalid_argument);
        c.sst_spatial_block_preconditioner_v24_experimental=true;
        validate_turbulence_controls(c);
    });
    return run_all();
}
