#include "cfdx/core/solvers/scalar_diffusion.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>

using namespace cfdx::core;
using namespace cfdx::testing;


static Mesh make_1d_chain(std::size_t n) {
    if (n == 0) throw std::invalid_argument("make_1d_chain: n must be positive");
    Mesh m;
    m.points().resize(8 * n);
    const double dx = 1.0 / static_cast<double>(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double x0 = dx * static_cast<double>(i);
        const double x1 = dx * static_cast<double>(i + 1);
        const std::size_t b = 8 * i;
        const double p[8][3] = {
            {x0,0,0}, {x1,0,0}, {x1,1,0}, {x0,1,0},
            {x0,0,1}, {x1,0,1}, {x1,1,1}, {x0,1,1}
        };
        for (std::size_t q = 0; q < 8; ++q)
            m.points().set(b + q, p[q][0], p[q][1], p[q][2]);
    }

    std::vector<std::size_t> left, right, y0, y1, z0, z1, internal;
    auto add_face = [&](std::initializer_list<std::size_t> v) {
        const std::size_t id = m.faces().n_faces();
        m.faces().push_face(std::vector<FaceIndex>(v.begin(), v.end()));
        return id;
    };
    left.push_back(add_face({0,4,7,3}));
    right.push_back(add_face({8*(n-1)+1,8*(n-1)+2,8*(n-1)+6,8*(n-1)+5}));
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t b = 8 * i;
        y0.push_back(add_face({b,b+1,b+5,b+4}));
        y1.push_back(add_face({b+3,b+7,b+6,b+2}));
        z0.push_back(add_face({b,b+3,b+2,b+1}));
        z1.push_back(add_face({b+4,b+5,b+6,b+7}));
    }
    for (std::size_t i = 0; i + 1 < n; ++i) {
        const std::size_t b = 8 * i;
        internal.push_back(add_face({b+1,b+2,b+6,b+5}));
    }

    m.ownership().resize(m.n_faces());
    m.ownership().set_owner(left[0], 0);
    m.ownership().set_neighbour(left[0], FaceOwnership::BOUNDARY);
    m.ownership().set_owner(right[0], n - 1);
    m.ownership().set_neighbour(right[0], FaceOwnership::BOUNDARY);
    for (std::size_t i = 0; i < n; ++i) {
        for (const auto f : {y0[i], y1[i], z0[i], z1[i]}) {
            m.ownership().set_owner(f, i);
            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        }
    }
    for (std::size_t i = 0; i + 1 < n; ++i) {
        m.ownership().set_owner(internal[i], i);
        m.ownership().set_neighbour(internal[i], static_cast<std::int64_t>(i + 1));
    }
    for (std::size_t i = 0; i < n; ++i) {
        m.cells().push_cell({
            i == 0 ? left[0] : internal[i-1],
            i + 1 == n ? right[0] : internal[i],
            y0[i], y1[i], z0[i], z1[i]
        });
    }
    return m;
}

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

    run_case("mixed_dirichlet_neumann_boundary", [] {
        Mesh m = cube();
        auto bc = PoissonBoundaryCondition::mixed(m.n_faces());
        bc.face_types[4] = PoissonBoundaryType::DIRICHLET;
        bc.face_types[5] = PoissonBoundaryType::NEUMANN;
        bc.face_values[4] = 0.0;
        bc.face_values[5] = 1.0;
        auto result = solve_poisson_mixed(m, bc, {0.0});
        EXPECT_TRUE(result.linear_result.status == SolverStatus::CONVERGED);
        EXPECT_NEAR(result.solution(0), 0.5, 1e-12);
    });

    run_case("neumann_flux_has_canonical_outward_sign", [] {
        Mesh m = cube();
        auto bc = PoissonBoundaryCondition::neumann(m.n_faces());
        bc.face_values.assign(m.n_faces(), 0.0);
        bc.face_values[5] = 2.0;
        Vector rhs;
        const auto geometry = make_geometry_cache(m);
        const auto A = assemble_cell_diffusion_matrix(
            m, bc, 1.0, {0.0}, rhs, geometry);
        (void)A;
        EXPECT_NEAR(rhs(0), 2.0, 1e-12);
    });


    run_case("multicell_dirichlet_matches_quadratic_manufactured_solution", [] {
        const std::size_t n = 8;
        Mesh m = make_1d_chain(n);
        auto bc = PoissonBoundaryCondition::dirichlet(m.n_faces());
        bc.face_values[0] = 0.0;
        bc.face_values[1] = 1.0;
        const std::vector<double> source(n, -2.0);
        ScalarDiffusionConfig cfg;
        cfg.tolerance = 1e-12;
        auto result = solve_poisson_dirichlet(m, bc, source, cfg);
        EXPECT_TRUE(result.linear_result.status == SolverStatus::CONVERGED);
        const auto geometry = make_geometry_cache(m);
        for (std::size_t c = 0; c < n; ++c) {
            const double x = geometry.cell_centres[c].x;
            EXPECT_NEAR(result.solution(c), x * x, 1e-12);
        }
        const double balance =
            poisson_conservation_balance(m, bc, source, result.solution,
                                         geometry, cfg.diffusivity);
        EXPECT_NEAR(balance, 0.0, 1e-12);
    });

    run_case("multicell_mixed_dirichlet_neumann_matches_quadratic", [] {
        const std::size_t n = 8;
        Mesh m = make_1d_chain(n);
        auto bc = PoissonBoundaryCondition::mixed(m.n_faces());
        bc.face_types[0] = PoissonBoundaryType::DIRICHLET;
        bc.face_types[1] = PoissonBoundaryType::NEUMANN;
        bc.face_values[0] = 0.0;
        bc.face_values[1] = 2.0;
        const std::vector<double> source(n, -2.0);
        auto result = solve_poisson_mixed(m, bc, source);
        EXPECT_TRUE(result.linear_result.status == SolverStatus::CONVERGED);
        const auto geometry = make_geometry_cache(m);
        for (std::size_t c = 0; c < n; ++c) {
            const double x = geometry.cell_centres[c].x;
            EXPECT_NEAR(result.solution(c), x * x, 1e-12);
        }
        const double balance =
            poisson_conservation_balance(m, bc, source, result.solution,
                                         geometry, 1.0);
        EXPECT_NEAR(balance, 0.0, 1e-12);
    });

    run_case("multicell_pure_neumann_gauge_and_compatibility", [] {
        const std::size_t n = 8;
        Mesh m = make_1d_chain(n);
        auto bc = PoissonBoundaryCondition::neumann(m.n_faces());
        bc.face_values[0] = 0.0;
        bc.face_values[1] = 2.0;
        ScalarDiffusionConfig cfg;
        const auto geometry = make_geometry_cache(m);
        const double reference = geometry.cell_centres[0].x *
                                 geometry.cell_centres[0].x + 3.0;
        cfg.gauge.reference_cell = 0;
        cfg.gauge.reference_value = reference;
        const std::vector<double> source(n, -2.0);
        auto result = solve_poisson_mixed(m, bc, source, cfg);
        EXPECT_TRUE(result.linear_result.status == SolverStatus::CONVERGED);
        for (std::size_t c = 0; c < n; ++c) {
            const double x = geometry.cell_centres[c].x;
            EXPECT_NEAR(result.solution(c), x * x + 3.0, 1e-12);
        }
        const double balance =
            poisson_conservation_balance(m, bc, source, result.solution,
                                         geometry, 1.0);
        EXPECT_NEAR(balance, 0.0, 1e-12);
    });

    run_case("poisson_conservation_balance_is_zero", [] {
        Mesh m = cube();
        auto bc = PoissonBoundaryCondition::mixed(m.n_faces());
        bc.face_types[4] = PoissonBoundaryType::DIRICHLET;
        bc.face_types[5] = PoissonBoundaryType::DIRICHLET;
        bc.face_values[4] = 0.0;
        bc.face_values[5] = 1.0;
        ScalarDiffusionConfig cfg;
        auto result = solve_poisson_mixed(m, bc, {0.0}, cfg);
        EXPECT_TRUE(result.linear_result.status == SolverStatus::CONVERGED);
        const auto geometry = make_geometry_cache(m);
        const double balance = poisson_conservation_balance(
            m, bc, {0.0}, result.solution, geometry, cfg.diffusivity);
        EXPECT_NEAR(balance, 0.0, 1e-12);
    });

    run_case("poisson_contract_rejects_wrong_boundary_type", [] {
        Mesh m = cube();
        auto bc = PoissonBoundaryCondition::neumann(m.n_faces());
        bc.face_values.assign(m.n_faces(), 0.0);
        EXPECT_THROW(
            solve_poisson_dirichlet(m, bc, {0.0}),
            std::invalid_argument);
    });

    run_case("pure_neumann_rejects_incompatible_data", [] {
        Mesh m = cube();
        auto bc = PoissonBoundaryCondition::neumann(m.n_faces());
        bc.face_values.assign(m.n_faces(), 0.0);
        EXPECT_THROW(
            solve_poisson_mixed(m, bc, {1.0}),
            std::invalid_argument);
    });

    run_case("pure_neumann_compatible_problem_uses_reference_gauge", [] {
        Mesh m = cube();
        auto bc = PoissonBoundaryCondition::neumann(m.n_faces());
        bc.face_values.assign(m.n_faces(), 0.0);
        ScalarDiffusionConfig cfg;
        cfg.gauge.reference_cell = 0;
        cfg.gauge.reference_value = 2.5;
        auto result = solve_poisson_mixed(m, bc, {0.0}, cfg);
        EXPECT_TRUE(result.linear_result.status == SolverStatus::CONVERGED);
        EXPECT_NEAR(result.solution(0), 2.5, 1e-12);
    });

    run_case("pure_neumann_nonzero_source_and_flux_is_compatible", [] {
        Mesh m = cube();
        auto bc = PoissonBoundaryCondition::neumann(m.n_faces());
        bc.face_values.assign(m.n_faces(), 0.0);
        bc.face_values[5] = 1.0;
        ScalarDiffusionConfig cfg;
        cfg.gauge.reference_value = 3.0;
        auto result = solve_poisson_mixed(m, bc, {-1.0}, cfg);
        EXPECT_TRUE(result.linear_result.status == SolverStatus::CONVERGED);
        EXPECT_NEAR(result.solution(0), 3.0, 1e-12);
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
