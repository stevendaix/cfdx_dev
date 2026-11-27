#include "cfdx/core/linalg/globalized_jfnk_v8.h"
#include "cfdx/physics/coupled_jfnk_driver.h"
#include "cfdx/physics/pressure_velocity_algorithms.h"
#include "common/test_harness.h"
#include <array>
#include <cmath>
#include <iostream>
#include <string>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

// A six-field nonlinear oracle with coupling between Ux, Uy, Uz, p, k,
// omega. This is NOT the finite-volume BFS operator. Its only purpose is to
// certify the matrix-free Krylov/globalization interface for all model tags.
class SixFieldBlockPreconditioner final : public Preconditioner {
public:
    explicit SixFieldBlockPreconditioner(CoupledSchurModel model) : model_(model) {}
    bool setup(const SparseMatrix& a) override {
        return a.n_rows() == 6 && a.n_cols() == 6;
    }
    bool apply(const Vector& r, Vector& z) const override {
        if (r.size() != 6 || z.size() != 6) return false;
        // Distinct fixed linear Schur approximations are tested in the
        // production coupled suite. Here use a stable diagonal block action.
        for (std::size_t i = 0; i < 6; ++i) z(i) = r(i) / (i==4 ? 3. : 1.);
        return true;
    }
    const char* name() const override { return to_string(model_); }
private:
    CoupledSchurModel model_;
};

static bool six_field_residual(const Vector& x, Vector& f) {
    if (x.size() != 6 || f.size() != 6 ||
        x(4) <= 0. || x(5) <= 0.) return false;
    const double u0=x(0)-1., u1=x(1)+0.5, u2=x(2)-0.25;
    const double p=x(3)-0.3, k=x(4)-0.8, w=x(5)-1.2;
    f(0) = u0 + 0.15*p + 0.07*k;
    f(1) = u1 - 0.18*p + 0.05*w;
    f(2) = u2 + 0.12*p;
    f(3) = p + 0.17*u0 - 0.06*u1;
    f(4) = 2.*k + k*k + 0.07*u0;
    f(5) = w + 0.1*w*w + 0.11*k;
    return true;
}

int main() {
    const std::array<CoupledSchurModel,6> family{{
        CoupledSchurModel::BlockLocal, CoupledSchurModel::PCD,
        CoupledSchurModel::LSC, CoupledSchurModel::BFBT,
        CoupledSchurModel::SIMPLE, CoupledSchurModel::SIMPLEC}};
    for (auto model : family) {
        run_case(std::string("v8_jfnk_uses_true_f_for_")+to_string(model), [model] {
            Vector x(6,0.);
            x(4)=0.25; x(5)=0.4;
            SixFieldBlockPreconditioner pre(model);
            GlobalizedJfnkControlsV8 controls;
            controls.nonlinear_tolerance=1e-9;
            controls.linear_tolerance=1e-5;
            const auto report=globalized_jfnk_v8(
                x, six_field_residual, &pre, controls);
            Vector f(6,0.);
            EXPECT_TRUE(six_field_residual(x,f));
            std::cout << "jfnk_result=" << report.reason << " before=" << report.initial_merit << " after=" << report.final_merit << " accepted=" << report.newton_steps << " linear=" << report.krylov_steps << " rejected=" << report.rejected_steps << "\n";
            EXPECT_TRUE(report.converged);
            EXPECT_TRUE(f.norm_inf()<1e-9);
            EXPECT_TRUE(x(4)>0. && x(5)>0.);
            std::cout << to_string(model) << " merit=" << f.norm_inf()
                      << " steps=" << report.newton_steps << '\n';
        });
    }
    run_case("v8_globalization_rejects_unattainable_step", [] {
        Vector x(2,0.);
        x(0)=1.;x(1)=1.;
        auto no_solution=[](const Vector& x,Vector& out) {
            out(0)=x(0)*x(0)+1.;out(1)=x(1)*x(1)+1.;
            return true;
        };
        const auto report=globalized_jfnk_v8(x,no_solution,nullptr);
        EXPECT_TRUE(!report.converged);
        EXPECT_TRUE(x.is_valid());
    });
    run_case("v8_rejects_invalid_start_no_fake_convergence", [] {
        Vector x(2,0.);
        auto invalid=[](const Vector&,Vector&){return false;};
        const auto report=globalized_jfnk_v8(x,invalid,nullptr);
        EXPECT_TRUE(!report.converged);
    });
    run_case("coupled_jfnk_block_scales_accumulate_global_maximum", [] {
        std::vector<double> forward(6, 1.0);
        accumulate_coupled_jfnk_block_scale(forward, 0, 0.5);
        accumulate_coupled_jfnk_block_scale(forward, 0, 12.0);
        accumulate_coupled_jfnk_block_scale(forward, 1, -8.0);
        accumulate_coupled_jfnk_block_scale(forward, 2, 0.25);
        accumulate_coupled_jfnk_block_scale(forward, 2, -3.0);
        EXPECT_NEAR(forward[0], 12.0, 0.0);
        EXPECT_NEAR(forward[1], 8.0, 0.0);
        EXPECT_NEAR(forward[2], 3.0, 0.0);
        EXPECT_NEAR(forward[3], 1.0, 0.0);
        std::vector<double> reverse(6, 1.0);
        accumulate_coupled_jfnk_block_scale(reverse, 0, 12.0);
        accumulate_coupled_jfnk_block_scale(reverse, 0, 0.5);
        accumulate_coupled_jfnk_block_scale(reverse, 1, -8.0);
        accumulate_coupled_jfnk_block_scale(reverse, 2, -3.0);
        accumulate_coupled_jfnk_block_scale(reverse, 2, 0.25);
        EXPECT_TRUE(forward == reverse);
        EXPECT_THROW(accumulate_coupled_jfnk_block_scale(forward, 6, 1.0),
                     std::invalid_argument);
        EXPECT_THROW(accumulate_coupled_jfnk_block_scale(
                         forward, 0, std::numeric_limits<double>::quiet_NaN()),
                     std::invalid_argument);
    });
    return run_all();
}
