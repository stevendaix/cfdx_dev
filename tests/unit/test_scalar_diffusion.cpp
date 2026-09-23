#include "cfdx/core/solvers/scalar_diffusion.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>

using namespace cfdx::core;
using namespace cfdx::testing;

static Mesh cube(double lx = 1.0) {
    Mesh m;
    m.points().resize(8);
    m.points().set(0,0,0,0); m.points().set(1,lx,0,0);
    m.points().set(2,lx,1,0); m.points().set(3,0,1,0);
    m.points().set(4,0,0,1); m.points().set(5,lx,0,1);
    m.points().set(6,lx,1,1); m.points().set(7,0,1,1);
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
        PoissonBoundaryCondition bc =
            PoissonBoundaryCondition::dirichlet(m.n_faces());
        // x=0 -> 0, x=1 -> 1; remaining faces are zero-gradient.
        bc.face_values[4] = 0.0;
        bc.face_values[5] = 1.0;
        auto result = solve_poisson_dirichlet(m, bc, {0.0});
        EXPECT_TRUE(result.linear_result.status == SolverStatus::CONVERGED);
        EXPECT_NEAR(result.solution(0), 0.5, 1e-12);
    });

    run_case("laplace_zero_source", [] {
        Mesh m = cube();
        auto bc = PoissonBoundaryCondition::dirichlet(m.n_faces());
        bc.face_values[5] = 1.0;
        auto result = solve_laplace_dirichlet(m, bc);
        EXPECT_TRUE(result.linear_result.status == SolverStatus::CONVERGED);
        EXPECT_NEAR(result.solution(0), 1.0, 1e-12);
    });

    run_case("volumetric_source_is_integrated_over_cell_volume", [] {
        Mesh m = cube(2.0);
        auto bc = PoissonBoundaryCondition::dirichlet(m.n_faces());
        Vector rhs;
        const auto geometry = make_geometry_cache(m);
        const auto A = assemble_cell_diffusion_matrix(
            m, bc, 1.0, {3.0}, rhs, geometry);
        (void)A;
        EXPECT_NEAR(geometry.cell_volumes[0], 2.0, 1e-12);
        EXPECT_NEAR(rhs(0), 6.0, 1e-12);
    });

    run_case("mixed_dirichlet_neumann_boundary", [] {\n        Mesh m = cube();\n        auto bc = PoissonBoundaryCondition::mixed(m.n_faces());\n        bc.face_types[4] = PoissonBoundaryType::DIRICHLET;\n        bc.face_types[5] = PoissonBoundaryType::NEUMANN;\n        bc.face_values[4] = 0.0;\n        bc.face_values[5] = 1.0;\n        auto result = solve_poisson_mixed(m, bc, {0.0});\n        EXPECT_TRUE(result.linear_result.status == SolverStatus::CONVERGED);\n        EXPECT_NEAR(result.solution(0), 0.5, 1e-12);\n    });\n\n    run_case("neumann_flux_has_canonical_outward_sign", [] {\n        Mesh m = cube();\n        auto bc = PoissonBoundaryCondition::neumann(m.n_faces());\n        bc.face_values.assign(m.n_faces(), 0.0);\n        bc.face_values[5] = 2.0;\n        Vector rhs;\n        const auto geometry = make_geometry_cache(m);\n        const auto A = assemble_cell_diffusion_matrix(\n            m, bc, 1.0, {0.0}, rhs, geometry);\n        (void)A;\n        EXPECT_NEAR(rhs(0), 2.0, 1e-12);\n    });\n\n    run_case("poisson_contract_rejects_wrong_boundary_type", [] {
        Mesh m = cube();
        auto bc = PoissonBoundaryCondition::neumann(m.n_faces());
        bc.face_values.assign(m.n_faces(), 0.0);
        EXPECT_THROW(
            solve_poisson_dirichlet(m, bc, {0.0}),
            std::invalid_argument);
    });

    run_case("scalar_diffusion_rejects_bad_inputs", [] {
        Mesh m = cube();
        auto bc = PoissonBoundaryCondition::dirichlet(m.n_faces());
        bc.face_values.assign(m.n_faces(), 0.0);
        EXPECT_THROW(
            solve_poisson_dirichlet(m, bc, {}),
            std::invalid_argument);
        ScalarDiffusionConfig cfg;
        cfg.diffusivity = 0.0;
        EXPECT_THROW(
            solve_poisson_dirichlet(m, bc, {0.0}, cfg),
            std::invalid_argument);
    });

    return run_all();
}
