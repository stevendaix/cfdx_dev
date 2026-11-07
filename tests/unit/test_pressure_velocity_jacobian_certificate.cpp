#include "cfdx/physics/pressure_velocity_jacobian_certificate.h"
#include "common/test_harness.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <tuple>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {
SparseMatrix sparse(std::size_t nrows, std::size_t ncols,
                    const std::vector<std::tuple<std::size_t, std::size_t, double>>& entries) {
    SparseMatrix result(nrows, ncols);
    for (const auto& [row, col, val] : entries)
        result.push_back(row, col, val);
    result.finalize();
    return result;
}
} // namespace

int main() {
    // Non-symmetric momentum, non-transpose D/G and a pinned gauge equation.
    const auto M = sparse(2, 2, {{0,0,4.0},{0,1,-0.2},{1,0,0.5},{1,1,3.0}});
    const auto G = sparse(2, 2, {{0,0,-1.2},{0,1,0.3},{1,0,0.4},{1,1,-0.1}});
    const auto D = sparse(2, 2, {{0,0,1.1},{0,1,-0.6}});
    const auto C = sparse(2, 2, {{0,0,0.7},{0,1,-0.25},{1,1,1.0}});
    const BlockOperator blocks(M,G,D,C);
    Vector state(4, 0.0);
    state(0)=0.7;
    state(1)=-0.4;
    state(2)=1.3;
    state(3)=0.2;

    auto actual_residual = [&](const Vector& x, Vector& r) -> bool {
        if (x.size() != 4) return false;
        r = Vector(4, 0.0);
        for (std::size_t row = 0; row < 2; ++row) {
            for (std::size_t col = 0; col < 2; ++col) {
                r(row) += M(row,col) * x(col) + G(row,col) * x(col+2);
                r(row+2) += D(row,col) * x(col) + C(row,col) * x(col+2);
            }
        }
        // A body force/source does not change the Jacobian.
        r(0) -= 0.25;
        r(2) += 0.8;
        return true;
    };

    run_case("pressure_velocity_jacobian_matches_independent_residual", [&] {
        const auto certificate = certify_pressure_velocity_jacobian(
            blocks, state, actual_residual);
        EXPECT_TRUE(certificate.passed);
        EXPECT_TRUE(certificate.columns_checked == 4);
        EXPECT_TRUE(certificate.max_tolerance_ratio <= 1.0);
        std::cout << "PV_JACOBIAN max_absolute_error="
                  << certificate.max_absolute_error << '\n';
    });

    run_case("pressure_velocity_jacobian_catches_missing_Rhie_Chow_like_C_term", [&] {
        const auto wrong_C = sparse(2,2,{{0,0,0.7},{1,1,1.0}});
        const BlockOperator wrong(M,G,D,wrong_C);
        const auto certificate = certify_pressure_velocity_jacobian(
            wrong, state, actual_residual);
        EXPECT_TRUE(!certificate.passed);
        EXPECT_TRUE(certificate.columns_checked == 4);
        EXPECT_TRUE(certificate.max_tolerance_ratio > 1.0);
        EXPECT_TRUE(certificate.worst_row >= 2);
        std::cout << "PV_JACOBIAN_MISSING_C tolerance_ratio="
                  << certificate.max_tolerance_ratio << '\n';
    });

    run_case("pressure_velocity_jacobian_catches_wrong_pressure_gradient_sign", [&] {
        const auto wrong_G = sparse(2,2,{{0,0,1.2},{0,1,0.3},
                                          {1,0,0.4},{1,1,-0.1}});
        const BlockOperator wrong(M,wrong_G,D,C);
        const auto certificate = certify_pressure_velocity_jacobian(
            wrong, state, actual_residual);
        EXPECT_TRUE(!certificate.passed);
        EXPECT_TRUE(certificate.max_tolerance_ratio > 1.0);
        EXPECT_TRUE(certificate.worst_row < 2);
    });

    run_case("pressure_velocity_jacobian_rejects_nonfinite_and_failed_evidence", [&] {
        auto bad_residual = [&](const Vector& x, Vector& y) {
            (void)x;
            y = Vector(4, 0.0);
            y(2) = std::numeric_limits<double>::quiet_NaN();
            return true;
        };
        const auto bad = certify_pressure_velocity_jacobian(
            blocks, state, bad_residual);
        EXPECT_TRUE(!bad.passed);
        EXPECT_TRUE(bad.columns_checked == 0);
        auto failed = [](const Vector&, Vector&) { return false; };
        EXPECT_TRUE(!certify_pressure_velocity_jacobian(
            blocks, state, failed).passed);
    });

    run_case("pressure_velocity_jacobian_sampling_is_reported_not_hidden", [&] {
        PressureVelocityJacobianProbeControls c;
        c.max_columns = 2;
        const auto certificate = certify_pressure_velocity_jacobian(
            blocks, state, actual_residual, c);
        EXPECT_TRUE(certificate.passed);
        EXPECT_TRUE(certificate.columns_checked == 2);
    });

    run_case("pressure_velocity_directional_screening_detects_C_and_G_errors", [&] {
        PressureVelocityDirectionalControls controls;
        controls.directions = 12;
        const auto correct = certify_pressure_velocity_jacobian_directions(
            blocks, state, actual_residual, controls);
        EXPECT_TRUE(correct.passed);
        EXPECT_TRUE(correct.directions_checked == 12);
        EXPECT_TRUE(correct.residual_evaluations == 24);
        const auto wrong_C = sparse(2,2,{{0,0,0.7},{1,1,1.0}});
        const BlockOperator missing_C(M,G,D,wrong_C);
        const auto missing = certify_pressure_velocity_jacobian_directions(
            missing_C, state, actual_residual, controls);
        EXPECT_TRUE(!missing.passed);
        EXPECT_TRUE(missing.worst_row >= 2);
        const auto wrong_G = sparse(2,2,{{0,0,1.2},{0,1,0.3},
                                          {1,0,0.4},{1,1,-0.1}});
        const BlockOperator sign_error(M,wrong_G,D,C);
        EXPECT_TRUE(!certify_pressure_velocity_jacobian_directions(
            sign_error, state, actual_residual, controls).passed);
        std::cout << "PV_DIRECTIONAL missing_C_ratio="
                  << missing.max_tolerance_ratio << '\n';
    });

    run_case("pressure_velocity_directional_deterministic_and_fail_closed", [&] {
        PressureVelocityDirectionalControls controls;
        controls.directions = 6;
        const auto a = certify_pressure_velocity_jacobian_directions(
            blocks, state, actual_residual, controls);
        const auto b = certify_pressure_velocity_jacobian_directions(
            blocks, state, actual_residual, controls);
        EXPECT_TRUE(a.passed && b.passed);
        EXPECT_NEAR(a.max_absolute_error, b.max_absolute_error, 0.0);
        auto failed = [](const Vector&, Vector&) { return false; };
        EXPECT_TRUE(!certify_pressure_velocity_jacobian_directions(
            blocks, state, failed, controls).passed);
        EXPECT_TRUE(certify_pressure_velocity_jacobian_directions(
            blocks, state, failed, controls).directions_checked == 0);
        auto nan = [](const Vector& x, Vector& r) {
            r = Vector(x.size(), std::numeric_limits<double>::quiet_NaN());
            return true;
        };
        EXPECT_TRUE(!certify_pressure_velocity_jacobian_directions(
            blocks, state, nan, controls).passed);
    });

    return run_all();
}
