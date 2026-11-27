#include "cfdx/core/linalg/globalized_jfnk_v8.h"
#include "cfdx/physics/coupled_jfnk_driver.h"
#include "cfdx/physics/pressure_velocity_algorithms.h"
#include "common/test_harness.h"
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>
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
    run_case("coupled_jfnk_builds_block_variable_scales", [] {
        const std::vector<double> state{
            1.0, -2.0, 0.0, 0.0, 0.0, 0.0,
            1.0e5, 1.0e5, 0.25, 0.5, 3.0, 4.0};
        const auto scales = coupled_jfnk_variable_scales(state, 2);
        EXPECT_TRUE(scales.size() == state.size());
        EXPECT_NEAR(scales[0], 2.0, 0.0);
        EXPECT_NEAR(scales[1], 2.0, 0.0);
        EXPECT_NEAR(scales[2], 1.0, 0.0);
        EXPECT_NEAR(scales[6], 1.0e5, 0.0);
        EXPECT_NEAR(scales[8], 0.5, 0.0);
        EXPECT_THROW(coupled_jfnk_variable_scales(state, 3),
                     std::invalid_argument);
    });
    run_case("coupled_jfnk_scaled_jv_is_not_set_by_pressure_offset", [] {
        const std::vector<double> state{1.0, 1.0e5};
        const std::vector<double> direction{1.0, 0.0};
        const std::vector<double> scales{1.0, 1.0e5};
        std::vector<double> perturbations;
        const CoupledResidual identity = [&](const std::vector<double>& x) {
            perturbations.push_back(std::abs(x[0] - state[0]));
            return x;
        };
        const auto jv = coupled_newton_directional_jvp_scaled(
            identity, state, direction, scales, 1e-6);
        EXPECT_TRUE(perturbations.size() == 2);
        EXPECT_TRUE(perturbations[0] < 3e-6);
        EXPECT_TRUE(perturbations[1] < 3e-6);
        EXPECT_NEAR(jv[0], 1.0, 1e-9);
        EXPECT_NEAR(jv[1], 0.0, 1e-9);
    });
    run_case("coupled_jfnk_does_not_commit_after_one_nonconverged_step", [] {
        CoupledNewtonKrylovResult step;
        step.accepted = true;
        step.physical_merit_after = 2e-3;
        EXPECT_TRUE(!coupled_jfnk_nonlinear_converged(step, 1e-8));
        step.physical_merit_after = 5e-9;
        EXPECT_TRUE(coupled_jfnk_nonlinear_converged(step, 1e-8));
        step.accepted = false;
        EXPECT_TRUE(!coupled_jfnk_nonlinear_converged(step, 1e-8));
        step.accepted = true;
        step.physical_merit_after = std::numeric_limits<double>::quiet_NaN();
        EXPECT_TRUE(!coupled_jfnk_nonlinear_converged(step, 1e-8));
    });
    run_case("coupled_jfnk_residual_scale_handles_exact_zero_without_tiny_floor", [] {
        EXPECT_NEAR(jfnk_residual_scale_or_exact_zero(0.0, 0.0, "momentum"), 1.0, 0.0);
        EXPECT_NEAR(jfnk_residual_scale_or_exact_zero(2.0, 5.0, "momentum"), 5.0, 0.0);
        EXPECT_THROW(jfnk_residual_scale_or_exact_zero(1e-30, 0.0, "continuity"),
                     std::runtime_error);
        EXPECT_THROW(jfnk_residual_scale_or_exact_zero(0.0, -1.0, "SST"),
                     std::runtime_error);
        EXPECT_THROW(jfnk_residual_scale_or_exact_zero(
                         std::numeric_limits<double>::quiet_NaN(), 1.0, "momentum"),
                     std::runtime_error);
    });
    run_case("coupled_jfnk_uses_equation_specific_fixed_residual_scales", [] {
        PhysicalResidualV51 initial;
        initial.signed_raw.assign(12, 0.0);
        initial.normalized.assign(12, 0.0);
        initial.normalization_scales = {
            2.0, 3.0, 2.5, 4.0, 3.5, 4.5,
            5.0, 6.0, 5.5, 7.0, 6.5, 7.5};
        const auto scales = coupled_jfnk_fixed_scales(initial, 2);
        EXPECT_TRUE(scales == initial.normalization_scales);
        EXPECT_THROW(coupled_jfnk_fixed_scales(initial, 3),
                     std::invalid_argument);
        initial.normalization_scales[4] = 0.0;
        EXPECT_THROW(coupled_jfnk_fixed_scales(initial, 2),
                     std::invalid_argument);
    });
    run_case("coupled_jfnk_requires_pressure_reference", [] {
        ScalarBoundaryConditions pressure_bcs;
        EXPECT_TRUE(!jfnk_pressure_has_dirichlet_reference(pressure_bcs));
        pressure_bcs["outlet"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        EXPECT_TRUE(!jfnk_pressure_has_dirichlet_reference(pressure_bcs));
        pressure_bcs["outlet"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
        EXPECT_TRUE(jfnk_pressure_has_dirichlet_reference(pressure_bcs));
    });
    run_case("coupled_jfnk_detects_frozen_face_valued_pressure_snapshot", [] {
        ScalarBoundaryConditions pressure_bcs;
        EXPECT_TRUE(!jfnk_has_face_valued_pressure_snapshot(pressure_bcs));
        auto snapshot = std::make_shared<ScalarBoundaryFaceValues>();
        snapshot->generation = 7;
        pressure_bcs["outlet"].face_values = snapshot;
        EXPECT_TRUE(jfnk_has_face_valued_pressure_snapshot(pressure_bcs));
    });
    run_case("coupled_jfnk_rejects_unimplemented_wall_and_ptc_options", [] {
        const CoupledJfnkWallOperator no_wall_operator{};
        const double disabled_ptc = std::numeric_limits<double>::infinity();
        EXPECT_TRUE(!coupled_jfnk_has_unsupported_extensions(
            no_wall_operator, nullptr, disabled_ptc));
        EXPECT_TRUE(coupled_jfnk_has_unsupported_extensions(
            no_wall_operator, nullptr, 0.25));
        EXPECT_TRUE(coupled_jfnk_has_unsupported_extensions(
            no_wall_operator, nullptr,
            std::numeric_limits<double>::quiet_NaN()));
    });
    run_case("coupled_jfnk_velocity_bc_rejects_mixed_component_types", [] {
        ScalarBoundaryConditions ux, uy, uz;
        ux["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 1.0, 0.0};
        uy["inlet"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        uz["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
        EXPECT_THROW(
            make_jfnk_velocity_boundary_conditions(ux, uy, uz),
            std::invalid_argument);
    });
    run_case("coupled_jfnk_velocity_bc_rejects_unrepresentable_scalar_type", [] {
        ScalarBoundaryConditions ux, uy, uz;
        ux["wall"] = {ScalarBoundaryType::FIXED_GRADIENT, 0.0, 0.25};
        uy["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        uz["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        EXPECT_THROW(
            make_jfnk_velocity_boundary_conditions(ux, uy, uz),
            std::invalid_argument);
    });
    run_case("coupled_jfnk_velocity_bc_preserves_vector_fixed_values", [] {
        ScalarBoundaryConditions ux, uy, uz;
        ux["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 1.25, 0.0};
        uy["inlet"] = {ScalarBoundaryType::FIXED_VALUE, -0.5, 0.0};
        uz["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 0.75, 0.0};
        const auto velocity_bcs =
            make_jfnk_velocity_boundary_conditions(ux, uy, uz);
        const auto it = velocity_bcs.find("inlet");
        EXPECT_TRUE(it != velocity_bcs.end());
        EXPECT_TRUE(it->second.type == VelocityBoundaryCondition::Type::FIXED_VALUE);
        EXPECT_NEAR(it->second.value.x, 1.25, 0.0);
        EXPECT_NEAR(it->second.value.y, -0.5, 0.0);
        EXPECT_NEAR(it->second.value.z, 0.75, 0.0);
    });
    run_case("coupled_jfnk_velocity_bc_rejects_nonfinite_fixed_value", [] {
        ScalarBoundaryConditions ux, uy, uz;
        ux["inlet"] = {
            ScalarBoundaryType::FIXED_VALUE,
            std::numeric_limits<double>::quiet_NaN(), 0.0};
        uy["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
        uz["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
        EXPECT_THROW(
            make_jfnk_velocity_boundary_conditions(ux, uy, uz),
            std::invalid_argument);
    });
    return run_all();
}
