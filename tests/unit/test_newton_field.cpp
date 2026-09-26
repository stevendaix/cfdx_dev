#include "cfdx/core/numerics/newton.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::core;
using namespace cfdx::testing;

int main()
{
    run_case("generic_field_newton_solves_componentwise_nonlinear_system", [] {
        Field<double, Location::CELL> phi(2, "phi", "1", 2);
        phi(0) = 1.0; phi(1) = 2.0;
        phi.component_data(1)[0] = 1.0;
        phi.component_data(1)[1] = 3.0;

        NewtonFieldControls controls;
        controls.max_iterations = 12;
        controls.tolerance = 1e-11;

        const auto result = solve_newton_field(
            phi,
            [](const auto& x, auto& r, auto& J) {
                const std::size_t n = x.size();
                J = SparseMatrix(2 * n, 2 * n);
                for (std::size_t d = 0; d < x.dimension(); ++d) {
                    for (std::size_t c = 0; c < n; ++c) {
                        const double v = x.component_data(d)[c];
                        r.component_data(d)[c] = v * v - 2.0;
                        J.push_back(d * n + c, d * n + c, 2.0 * v);
                    }
                }
            },
            controls);

        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.iterations > 0);
        for (std::size_t d = 0; d < phi.dimension(); ++d)
            for (std::size_t c = 0; c < phi.size(); ++c)
                EXPECT_NEAR(phi.component_data(d)[c], std::sqrt(2.0), 1e-9);
        EXPECT_TRUE(result.residual < 1e-10);
    });

    return run_all();
}
