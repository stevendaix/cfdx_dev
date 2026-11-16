#include "cfdx/physics/sst_spatial_block_preconditioner_v24.h"
#include "cfdx/physics/sst_jfnk_v22.h"
#include "cfdx/physics/turbulence_transport.h"
#include "common/test_harness.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>

using namespace cfdx::physics;
using namespace cfdx::testing;

static cfdx::core::SparseMatrix transport(std::size_t n,double d,double off) {
    cfdx::core::SparseMatrix A(n,n);
    for (std::size_t i=0;i<n;++i) {
        A.push_back(i,i,d);
        if (i>0) A.push_back(i,i-1,off);
        if (i+1<n) A.push_back(i,i+1,off);
    }
    A.finalize();
    return A;
}

static std::vector<double> frozen_operator(const std::vector<double>& z) {
    const std::size_t n=z.size()/2;
    std::vector<double> rhs(2*n);
    for (std::size_t i=0;i<n;++i) {
        rhs[i]=4.*z[i]+0.2*z[n+i];
        rhs[n+i]=0.1*z[i]+3.*z[n+i];
        if (i>0) {rhs[i]-=z[i-1];rhs[n+i]-=0.7*z[n+i-1];}
        if (i+1<n) {rhs[i]-=z[i+1];rhs[n+i]-=0.7*z[n+i+1];}
    }
    return rhs;
}
static double error_l2(const std::vector<double>& a,const std::vector<double>& b) {
    double s=0.;
    for(std::size_t i=0;i<a.size();++i) s+=(a[i]-b[i])*(a[i]-b[i]);
    return std::sqrt(s);
}

int main() {
    run_case("v24_offcell_transport_reduces_true_frozen_linear_error",[] {
        const std::size_t n=12;
        auto ak=transport(n,4.,-1.),aw=transport(n,3.,-0.7);
        std::vector<SSTBlockPreconditionerV23::Block> d(n,{4.,0.2,0.1,3.});
        const SSTBlockPreconditionerV23 local(d);
        const std::vector<double> unit(2*n,1.);
        const SSTSpatialBlockPreconditionerV24 spatial(ak,aw,unit,unit,unit,local);
        EXPECT_TRUE(spatial.off_cell_coefficients()==4*(n-1));
        std::vector<double> exact(2*n,1.);
        const auto rhs=frozen_operator(exact);
        const auto xl=local.apply(rhs), xs=spatial.apply(rhs);
        const double el=error_l2(frozen_operator(xl),rhs);
        const double es=error_l2(frozen_operator(xs),rhs);
        std::cout << "V24 spatial toy residual: local=" << el << " spatial=" << es << "\n";
        EXPECT_TRUE(es<el*0.2);
    });
    run_case("v24_right_action_is_fixed_linear_and_deterministic",[] {
        auto ak=transport(3,4.,-1.),aw=transport(3,3.,-0.7);
        const SSTBlockPreconditionerV23 local(
            {{4.,0.2,0.1,3.},{4.,0.2,0.1,3.},{4.,0.2,0.1,3.}});
        const SSTSpatialBlockPreconditionerV24 p(
            ak,aw,{1.,2.,3.,4.,5.,6.},{2.,3.,4.,5.,6.,7.},
            {1.,2.,1.,2.,1.,2.},local);
        const std::vector<double> a{1.,-2.,3.,4.,-5.,6.},b{-1.,3.,2.,1.,1.,-7.};
        const auto pa=p.apply(a),pb=p.apply(b),again=p.apply(a);
        std::vector<double> sum(6);
        for(std::size_t i=0;i<6;++i)sum[i]=a[i]+b[i];
        const auto ps=p.apply(sum);
        for(std::size_t i=0;i<6;++i) {
            EXPECT_NEAR(ps[i],pa[i]+pb[i],1e-10);
            EXPECT_NEAR(pa[i],again[i],0.0);
        }
    });
    run_case("v24_zero_offcell_is_identical_to_v23_local",[] {
        auto ak=transport(1,4.,0.),aw=transport(1,3.,0.);
        const SSTBlockPreconditionerV23 local({{4.,0.2,0.1,3.}});
        const SSTSpatialBlockPreconditionerV24 p(ak,aw,{1.,1.},{1.,1.},{1.,1.},local);
        const auto x=p.apply({1.,-2.}),ref=local.apply({1.,-2.});
        EXPECT_NEAR(x[0],ref[0],1e-15);
        EXPECT_NEAR(x[1],ref[1],1e-15);
        EXPECT_TRUE(p.off_cell_coefficients()==0);
    });
    run_case("v24_malformed_operator_and_rhs_fail_closed",[] {
        auto a=transport(2,4.,-1.),w=transport(2,3.,-0.7);
        const SSTBlockPreconditionerV23 local({{4.,0.2,0.1,3.},{4.,0.2,0.1,3.}});
        EXPECT_THROW(SSTSpatialBlockPreconditionerV24(a,w,{1.,1.,1.},
            {1.,1.,1.,1.},{1.,1.,1.,1.},local),std::invalid_argument);
        EXPECT_THROW(SSTSpatialBlockPreconditionerV24(a,w,{1.,1.,1.,1.},
            {1.,1.,1.,1.},{1.,0.,1.,1.},local),std::invalid_argument);
        const SSTSpatialBlockPreconditionerV24 p(a,w,{1.,1.,1.,1.},
            {1.,1.,1.,1.},{1.,1.,1.,1.},local);
        EXPECT_THROW(p.apply({1.,2.}),std::invalid_argument);
        EXPECT_THROW(p.apply({1.,2.,3.,std::numeric_limits<double>::infinity()}),
                     std::invalid_argument);
    });
    run_case("v24_jfnk_keeps_true_residual_gate_and_nonnegative_fields",[] {
        const auto ak=transport(3,4.,-1.),aw=transport(3,3.,-0.7);
        const SSTBlockPreconditionerV23 local(
            {{4.,0.2,0.1,3.},{4.,0.2,0.1,3.},{4.,0.2,0.1,3.}});
        const std::vector<double> unit(6,1.);
        const SSTSpatialBlockPreconditionerV24 spatial(ak,aw,unit,unit,unit,local);
        const SSTPhysicalResidualV22 F=[](const auto& field) {
            auto residual=frozen_operator(field);
            const auto root=frozen_operator({2.,2.,2.,3.,3.,3.});
            for(std::size_t i=0;i<6;++i)residual[i]-=root[i];
            return residual;
        };
        CoupledNewtonKrylovControls ctrl;
        ctrl.krylov_dimension=6;
        ctrl.linear_relative_tolerance=0.25;
        ctrl.finite_difference_step=1e-7;
        const auto candidate=sst_jfnk_propose_v22(unit,unit,3,1e-12,1e-12,F,
            ctrl,{},[&spatial](const auto& rhs){return spatial.apply(rhs);});
        EXPECT_TRUE(candidate.available);
        EXPECT_TRUE(candidate.linear_true_relative<=ctrl.linear_relative_tolerance);
        EXPECT_TRUE(candidate.frozen_merit_after<candidate.frozen_merit_before);
        for(double x:candidate.fields)EXPECT_TRUE(x>0.);
    });
    run_case("v24_local_jacobian_audit_finds_missing_cross_derivatives",[] {
        // Three local physics rows: the frozen approximation intentionally
        // omits a k-to-omega derivative. The true callback detects it.
        const SSTBlockPreconditionerV23 blocks({{4.,0.,0.,3.}});
        const std::vector<double> state{2.,3.};
        const auto F=[](const std::vector<double>& x) {
            return std::vector<double>{2.*x[0]+0.5*x[1],
                                       4.*x[0]+x[1]};
        };
        const auto audit=sst_local_jacobian_audit_v24(
            0,state,{1.,1.},{1.,1.},blocks,F);
        EXPECT_TRUE(audit.residual_evaluations==4);
        EXPECT_NEAR(audit.actual[0],4.,1e-8);
        EXPECT_NEAR(audit.actual[1],1.5,1e-8);
        EXPECT_NEAR(audit.actual[2],8.,1e-8);
        EXPECT_NEAR(audit.actual[3],3.,1e-8);
        EXPECT_TRUE(audit.maximum_relative_difference>0.95);
        EXPECT_THROW(sst_local_jacobian_audit_v24(
            1,state,{1.,1.},{1.,1.},blocks,F),std::invalid_argument);
    });
    run_case("v24_spatial_preconditioner_is_opt_in",[] {
        TurbulenceTransportControls c;
        EXPECT_TRUE(!c.sst_spatial_block_preconditioner_v24_experimental);
        c.sst_spatial_block_preconditioner_v24_experimental=true;
        EXPECT_THROW(validate_turbulence_controls(c),std::invalid_argument);
        c.sst_jfnk_v22_experimental=true;
        validate_turbulence_controls(c);
    });
    return run_all();
}
