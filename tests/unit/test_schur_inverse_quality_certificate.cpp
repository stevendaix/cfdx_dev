#include "cfdx/core/linalg/schur_inverse_quality_certificate.h"
#include "common/test_harness.h"

#include <iostream>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
SparseMatrix diag(double a, double b) {
    SparseMatrix A(2,2);
    A.push_back(0,0,a); A.push_back(1,1,b); A.finalize();
    return A;
}

class MockSchurInverse final : public SchurApproximation {
public:
    explicit MockSchurInverse(bool wrong = false) : wrong_(wrong) {}
    const char* name() const noexcept override { return "mock_exact"; }
    SchurAction action() const noexcept override {return SchurAction::InverseOperator;}
    bool setup(const BlockOperator& A) override {
        return A.is_valid() && A.pressure_size() == 2;
    }
    bool update_values(const BlockOperator& A) override {return setup(A);}
    bool apply(const Vector& rhs, Vector& output) const override {
        if (rhs.size() != 2) return false;
        output = Vector(2,0.0);
        output(0) = (wrong_ ? 4.0 : -4.0) * rhs(0);
        output(1) = -2.0 * rhs(1);
        return true;
    }
private:
    bool wrong_;
};
}
int main() {
    const auto Auu = diag(4.0, 2.0);
    const auto G = diag(1.0, 1.0);
    const auto D = diag(1.0, 1.0);
    SparseMatrix C(2, 2); C.finalize();
    const BlockOperator blocks(Auu,G,D,C);
    const auto solve = [](const Vector& r, Vector& x) {
        x = Vector(2, 0.0); x(0) = r(0)/4.0; x(1) = r(1)/2.0;
        return true;
    };
    MockSchurInverse exact(false), wrong(true);
    EXPECT_TRUE(exact.setup(blocks));
    EXPECT_TRUE(wrong.setup(blocks));
    auto probes = schur_inverse_quality_basis(2);
    run_case("schur_inverse_oracle_accepts_exact_action", [&] {
        auto report = certify_schur_inverse_action(blocks, exact, solve, probes);
        EXPECT_TRUE(report.complete);
        EXPECT_TRUE(report.passed);
        EXPECT_TRUE(report.probes_checked == 3);
        EXPECT_NEAR(report.worst_schur_defect, 0.0, 0.0);
    });
    run_case("schur_inverse_oracle_detects_wrong_sign", [&] {
        auto report = certify_schur_inverse_action(blocks, wrong, solve, probes);
        EXPECT_TRUE(report.complete);
        EXPECT_TRUE(!report.passed);
        EXPECT_TRUE(report.worst_relative_defect > 1.5);
        std::cout << "SCHUR_ORACLE wrong_sign_relative_defect="
                  << report.worst_relative_defect << '\n';
    });
    run_case("schur_inverse_oracle_rejects_inaccurate_reference_solve", [&] {
        auto stale = [](const Vector& r, Vector& x) {
            x = r; return true;
        };
        const auto report = certify_schur_inverse_action(blocks, exact, stale, probes);
        EXPECT_TRUE(!report.complete);
        EXPECT_TRUE(!report.passed);
    });
    run_case("schur_inverse_oracle_fail_closed_on_invalid_rhs_or_callback", [&] {
        const std::vector<Vector> bad_rhs{Vector(1,1.0)};
        const auto bad = certify_schur_inverse_action(blocks, exact, solve, bad_rhs);
        EXPECT_TRUE(!bad.complete);
        auto fails = [](const Vector&, Vector&) { return false; };
        const auto failed = certify_schur_inverse_action(blocks, exact, fails, probes);
        EXPECT_TRUE(!failed.passed && !failed.complete);
    });
    return run_all();
}
