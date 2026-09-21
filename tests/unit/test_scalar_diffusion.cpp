#include "cfdx/core/solvers/scalar_diffusion.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>

using namespace cfdx::core;
using namespace cfdx::testing;

static Mesh cube() {
    Mesh m;
    m.points().resize(8);
    m.points().set(0,0,0,0); m.points().set(1,1,0,0);
    m.points().set(2,1,1,0); m.points().set(3,0,1,0);
    m.points().set(4,0,0,1); m.points().set(5,1,0,1);
    m.points().set(6,1,1,1); m.points().set(7,0,1,1);
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for(std::size_t f=0; f<6; ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    return m;
}

int main() {
    run_case("poisson_one_cell_dirichlet", [] {
        Mesh m = cube();
        DirichletBoundary bc;
        bc.face_values.assign(m.n_faces(), std::numeric_limits<double>::quiet_NaN());
        // x=0 -> 0, x=1 -> 1; remaining faces are zero-gradient.
        bc.face_values[4] = 0.0;
        bc.face_values[5] = 1.0;
        auto result = solve_poisson_dirichlet(m, bc, {0.0});
        EXPECT_TRUE(result.linear_result.status == SolverStatus::CONVERGED);
        EXPECT_NEAR(result.solution(0), 0.5, 1e-12);
    });

    run_case("laplace_zero_source", [] {
        Mesh m = cube();
        DirichletBoundary bc;
        bc.face_values.assign(m.n_faces(), 0.0);
        bc.face_values[5] = 1.0;
        auto result = solve_laplace_dirichlet(m, bc);
        EXPECT_TRUE(result.linear_result.status == SolverStatus::CONVERGED);
        EXPECT_NEAR(result.solution(0), 1.0/6.0, 1e-12);
    });

    run_case("scalar_diffusion_rejects_bad_inputs", [] {
        Mesh m = cube();
        DirichletBoundary bc;
        bc.face_values.assign(m.n_faces(), 0.0);
        EXPECT_THROW(solve_poisson_dirichlet(m, bc, {}, ScalarDiffusionConfig{}), std::invalid_argument);
        ScalarDiffusionConfig cfg;
        cfg.diffusivity = 0.0;
        EXPECT_THROW(solve_poisson_dirichlet(m, bc, {0.0}, cfg), std::invalid_argument);
    });

    return run_all();
}
