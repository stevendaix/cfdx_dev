#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_two_cell_channel()
{
    Mesh m;
    m.points().resize(12);
    const double xyz[12][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i, xyz[i][0], xyz[i][1], xyz[i][2]);

    m.faces().push_face({0,3,7,4});       // x = 0
    m.faces().push_face({1,2,6,5});       // internal
    m.faces().push_face({8,9,11,10});     // x = 2
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({1,8,10,5});
    m.faces().push_face({2,6,11,9});
    m.faces().push_face({1,2,9,8});
    m.faces().push_face({5,10,11,6});

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 11; ++f) {
        const auto owner = (f == 2 || f >= 7) ? 1 : 0;
        m.ownership().set_owner(f, owner);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);

    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});

    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,2,3,4,5,6,7,8,9,10};
    m.boundary().add_patch(wall);
    return m;
}

static std::array<std::vector<double>, 3> unit_rAU(std::size_t n)
{
    return {
        std::vector<double>(n, 1.0),
        std::vector<double>(n, 1.0),
        std::vector<double>(n, 1.0)};
}

// Stretched two-cell channel: cell 0 spans x in [0, 0.5], cell 1 spans
// x in [0.5, 2] (stretch ratio 3). Same topology and face winding as
// make_two_cell_channel(); only the x-coordinates of the internal and
// right-boundary points change.
static Mesh make_two_cell_channel_stretched()
{
    Mesh m;
    m.points().resize(12);
    const double xyz[12][3] = {
        {0,0,0},{0.5,0,0},{0.5,1,0},{0,1,0},
        {0,0,1},{0.5,0,1},{0.5,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i, xyz[i][0], xyz[i][1], xyz[i][2]);

    m.faces().push_face({0,3,7,4});       // x = 0
    m.faces().push_face({1,2,6,5});       // internal, x = 0.5
    m.faces().push_face({8,9,11,10});     // x = 2
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({1,8,10,5});
    m.faces().push_face({2,6,11,9});
    m.faces().push_face({1,2,9,8});
    m.faces().push_face({5,10,11,6});

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 11; ++f) {
        const auto owner = (f == 2 || f >= 7) ? 1 : 0;
        m.ownership().set_owner(f, owner);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);

    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});

    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,2,3,4,5,6,7,8,9,10};
    m.boundary().add_patch(wall);
    return m;
}

static std::array<std::vector<double>, 3> scale_rAU(
    const std::array<std::vector<double>, 3>& rAU, double s)
{
    std::array<std::vector<double>, 3> out{};
    for (std::size_t k = 0; k < 3; ++k) {
        out[k].resize(rAU[k].size());
        for (std::size_t i = 0; i < rAU[k].size(); ++i)
            out[k][i] = s * rAU[k][i];
    }
    return out;
}

// The Rhie-Chow pressure correction of an internal face, extracted as
// phi_interpolated - phi_RC. Both fluxes come from the production paths
// (make_mass_flux / make_rhie_chow_mass_flux), so no formula is duplicated
// here: only the difference is taken.
static double rc_correction(std::size_t face,
    const Mesh& mesh,
    const FvGeometry& geometry,
    const Field<double, Location::CELL>& U,
    const Field<double, Location::CELL>& p,
    const std::array<std::vector<double>, 3>& rAU,
    double rho,
    const VelocityBoundaryConditions& ubc,
    const ScalarBoundaryConditions& pbc)
{
    const auto phi_interp = make_mass_flux(mesh, geometry, U, rho, ubc);
    const auto phi_rc = make_rhie_chow_mass_flux(
        mesh, geometry, U, p, rAU, rho, ubc, pbc);
    return phi_interp(face) - phi_rc(face);
}

int main()
{
    run_case("n9_s9_rhie_chow_linear_pressure_consistency", [] {
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        Field<double, Location::CELL> grad_p(2, "grad_p", "Pa/m", 3);
        U.set(0, 2.0, 0.0, 0.0);
        U.set(1, 2.0, 0.0, 0.0);
        p(0) = 1.0;
        p(1) = 3.0;

        // Exact linear pressure field: dp/dx = 2 Pa/m.
        grad_p.set(0, 2.0, 0.0, 0.0);
        grad_p.set(1, 2.0, 0.0, 0.0);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const auto phi = make_rhie_chow_mass_flux(
            mesh, geometry, U, p, rAU, 1.0, ubc, pbc);
        const auto phi_linear = make_mass_flux(mesh, geometry, U, 1.0, ubc);

        // The Rhie-Chow pressure correction must vanish for a pressure field
        // whose cell gradient is exactly the same as its face pressure jump.
        EXPECT_NEAR(phi(1), phi_linear(1), 1e-12);
    });

    run_case("n9_s9_rhie_chow_pressure_offset_invariance", [] {
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p0(2, "p", "Pa", 1);
        Field<double, Location::CELL> p1(2, "p", "Pa", 1);
        U.fill(0.0);
        p0(0) = -2.0;
        p0(1) = 4.0;
        p1(0) = 1002.0;
        p1(1) = 1008.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const auto phi0 = make_rhie_chow_mass_flux(
            mesh, geometry, U, p0, rAU, 1.0, ubc, pbc);
        const auto phi1 = make_rhie_chow_mass_flux(
            mesh, geometry, U, p1, rAU, 1.0, ubc, pbc);

        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
            EXPECT_NEAR(phi0(f), phi1(f), 1e-12);
    });

    run_case("n9_s9_rhie_chow_alternating_pressure_mode_is_coupled", [] {
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p(0) = 1.0;
        p(1) = -1.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const auto phi = make_rhie_chow_mass_flux(
            mesh, geometry, U, p, rAU, 1.0, ubc, pbc);

        // A collocated cell-centred velocity interpolation alone would produce
        // zero mass flux. Rhie-Chow must couple this alternating pressure mode
        // to the face flux instead of allowing a checkerboard pressure mode to
        // become invisible to continuity.
        EXPECT_TRUE(std::abs(phi(1)) > 1e-12);
        EXPECT_NEAR(phi(1), 2.0, 1e-12);
    });


    run_case("n9_s9_rhie_chow_correction_scales_linearly_with_rAU_majumdar", [] {
        // D^RC = V/a_P^r = alpha_u * V/a_P = alpha_u * D^pc: the production
        // RC term weights the face pressure flux with the RELAXED inverse
        // momentum diagonal. Scaling rAU by alpha_u must therefore scale the
        // RC correction by exactly alpha_u. This is the operator-level
        // statement of the Majumdar relation documented in
        // N9_S9_RHIE_CHOW_CONSISTENCY.md.
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p(0) = 1.0;
        p(1) = -1.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const double c1 = rc_correction(
            1, mesh, geometry, U, p, rAU, 1.0, ubc, pbc);
        EXPECT_TRUE(std::abs(c1) > 1e-12);
        for (const double alpha : {0.65, 0.8}) {
            const auto rAU_a = scale_rAU(rAU, alpha);
            const double ca = rc_correction(
                1, mesh, geometry, U, p, rAU_a, 1.0, ubc, pbc);
            EXPECT_NEAR(ca, alpha * c1, 1e-12);
        }
    });

    run_case("n9_s9_rhie_chow_correction_is_linear_in_pressure", [] {
        // The internal RC operator is linear in p (the Green-Gauss gradient
        // and the owner-to-neighbour difference are linear, and the
        // zero-gradient wall treatment preserves linearity):
        // corr(p_a + p_b) = corr(p_a) + corr(p_b).
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> pa(2, "pa", "Pa", 1);
        Field<double, Location::CELL> pb(2, "pb", "Pa", 1);
        Field<double, Location::CELL> psum(2, "psum", "Pa", 1);
        U.fill(0.0);
        pa(0) = 1.0;
        pa(1) = -1.0;      // alternating mode
        pb(0) = 1.0;
        pb(1) = 3.0;       // linear mode
        psum(0) = pa(0) + pb(0);
        psum(1) = pa(1) + pb(1);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const double ca = rc_correction(
            1, mesh, geometry, U, pa, rAU, 1.0, ubc, pbc);
        const double cb = rc_correction(
            1, mesh, geometry, U, pb, rAU, 1.0, ubc, pbc);
        const double cs = rc_correction(
            1, mesh, geometry, U, psum, rAU, 1.0, ubc, pbc);
        EXPECT_TRUE(std::abs(ca) > 1e-12);
        EXPECT_TRUE(std::abs(cb) > 1e-12);
        EXPECT_NEAR(cs, ca + cb, 1e-12);
    });

    run_case("n9_s9_rhie_chow_directional_rAU_face_normal_projection", [] {
        // rfn = cx*ex^2 + cy*ey^2 + cz*ez^2: only the face-normal projection
        // of the directional inverse diagonals enters the correction. The
        // internal face normal is +x, so anisotropic rAU in y/z must not
        // change the correction, while scaling the x component scales it
        // exactly.
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p(0) = 1.0;
        p(1) = -1.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU_iso = unit_rAU(mesh.n_cells());
        std::array<std::vector<double>, 3> rAU_aniso{};
        rAU_aniso[0] = std::vector<double>(mesh.n_cells(), 1.0);
        rAU_aniso[1] = std::vector<double>(mesh.n_cells(), 2.0);
        rAU_aniso[2] = std::vector<double>(mesh.n_cells(), 2.0);
        std::array<std::vector<double>, 3> rAU_x2{};
        rAU_x2[0] = std::vector<double>(mesh.n_cells(), 2.0);
        rAU_x2[1] = std::vector<double>(mesh.n_cells(), 1.0);
        rAU_x2[2] = std::vector<double>(mesh.n_cells(), 1.0);

        const double c_iso = rc_correction(
            1, mesh, geometry, U, p, rAU_iso, 1.0, ubc, pbc);
        const double c_aniso = rc_correction(
            1, mesh, geometry, U, p, rAU_aniso, 1.0, ubc, pbc);
        const double c_x2 = rc_correction(
            1, mesh, geometry, U, p, rAU_x2, 1.0, ubc, pbc);
        EXPECT_TRUE(std::abs(c_iso) > 1e-12);
        EXPECT_NEAR(c_aniso, c_iso, 1e-12);
        EXPECT_NEAR(c_x2, 2.0 * c_iso, 1e-12);
    });

    run_case("n9_s9_rhie_chow_stretched_mesh_volume_weighted_dAU", [] {
        // On a stretched mesh with unequal cell volumes and per-cell inverse
        // diagonals, the face coefficient must be the volume-weighted dAU
        // interpolation dAU_f = 0.5*(V_P*rAU_P + V_N*rAU_N) applied along the
        // owner-to-neighbour line:
        // corr = rho * dAU_f * (p_N - p_P)/d * (S_f . e).
        // The expected value is rebuilt from the same FvGeometry the
        // production path uses - no second implementation of the operator.
        const Mesh mesh = make_two_cell_channel_stretched();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p(0) = 1.0;
        p(1) = -1.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const double rAU0 = 0.5;
        const double rAU1 = 1.7;
        std::array<std::vector<double>, 3> rAU{};
        for (std::size_t k = 0; k < 3; ++k) {
            rAU[k] = std::vector<double>(mesh.n_cells(), 1.0);
            rAU[k][0] = rAU0;
            rAU[k][1] = rAU1;
        }

        const double corr = rc_correction(
            1, mesh, geometry, U, p, rAU, 1.0, ubc, pbc);

        const Vec3 c0 = geometry.cell_centres[0];
        const Vec3 c1 = geometry.cell_centres[1];
        const double dx = c1.x - c0.x;
        const double dy = c1.y - c0.y;
        const double dz = c1.z - c0.z;
        const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        const double ex = dx / d;
        const double ey = dy / d;
        const double ez = dz / d;
        const Vec3 Sf = geometry.face_area_vectors[1];
        const double Sf_dot_e = Sf.x * ex + Sf.y * ey + Sf.z * ez;
        // The mesh is orthogonal: S_non = 0, so the deferred non-orthogonal
        // part vanishes and the expected correction is exact.
        const double dAU_f = 0.5 * (geometry.cell_volumes[0] * rAU0 +
                                   geometry.cell_volumes[1] * rAU1);
        const double expected =
            1.0 * dAU_f * ((p(1) - p(0)) / d) * Sf_dot_e;

        EXPECT_TRUE(std::abs(geometry.cell_volumes[0] - 0.5) < 1e-12);
        EXPECT_TRUE(std::abs(geometry.cell_volumes[1] - 1.5) < 1e-12);
        EXPECT_NEAR(corr, expected, 1e-12);
    });

    run_case("n9_s9_rhie_chow_uniform_mesh_limit_exact_face_operator", [] {
        // On the uniform unit two-cell channel the production correction
        // must collapse to the analytic face operator
        // rho * V * rAU * (p_N - p_P)/d * (S_f . e) with V0 = V1, and the
        // geometry must be the unit construction the analytic constants
        // assume (V = d = S_f.e = 1).
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p(0) = 1.0;
        p(1) = 3.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const double corr = rc_correction(
            1, mesh, geometry, U, p, rAU, 1.0, ubc, pbc);

        const Vec3 c0 = geometry.cell_centres[0];
        const Vec3 c1 = geometry.cell_centres[1];
        const double dx = c1.x - c0.x;
        const double dy = c1.y - c0.y;
        const double dz = c1.z - c0.z;
        const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        const double ex = dx / d;
        const double ey = dy / d;
        const double ez = dz / d;
        const Vec3 Sf = geometry.face_area_vectors[1];
        const double Sf_dot_e = Sf.x * ex + Sf.y * ey + Sf.z * ez;

        EXPECT_NEAR(geometry.cell_volumes[0], 1.0, 1e-12);
        EXPECT_NEAR(geometry.cell_volumes[1], 1.0, 1e-12);
        EXPECT_NEAR(d, 1.0, 1e-12);
        EXPECT_NEAR(Sf_dot_e, 1.0, 1e-12);

        const double dAU_f = 0.5 * (geometry.cell_volumes[0] + geometry.cell_volumes[1]);
        const double expected =
            1.0 * dAU_f * ((p(1) - p(0)) / d) * Sf_dot_e;
        EXPECT_NEAR(corr, expected, 1e-12);
    });


    run_case("n9_s9_pressure_matrix_coefficient_matches_rhie_chow_on_skew_face", [] {
        // The implicit pressure matrix must use the same centre-line
        // projection as rhie_chow_pressure_flux_internal(). On this skew face
        // d points along x while Sf has a tangential y component. The
        // directional inverse diagonal is intentionally anisotropic so that
        // projecting it on Sf (rather than d) is also detected.
        const Vec3 dvec{1.0, 0.0, 0.0};
        const Vec3 Sf{1.0, 0.4, 0.0};
        const Vec3 dAU_diagonal{1.0, 4.0, 9.0};

        const double coefficient =
            pressure_face_orthogonal_coefficient(Sf, dvec, dAU_diagonal);
        EXPECT_NEAR(coefficient, 1.0, 1e-12);

        const double area = Sf.mag();
        const Vec3 face_normal{Sf.x / area, Sf.y / area, Sf.z / area};
        const double old_face_normal_rfn =
            dAU_diagonal.x * face_normal.x * face_normal.x +
            dAU_diagonal.y * face_normal.y * face_normal.y +
            dAU_diagonal.z * face_normal.z * face_normal.z;
        const double old_area_over_distance = old_face_normal_rfn * area / dvec.mag();

        // This is deliberately not a loose numerical tolerance: the old
        // area/d plus face-normal projection differs by more than 0.5 on this
        // analytic geometry, while the corrected coefficient is exactly 1.
        EXPECT_TRUE(std::abs(old_area_over_distance - coefficient) > 0.5);

        // Orthogonal limit: retain the previous coefficient exactly even for
        // anisotropic momentum response.
        const double orthogonal = pressure_face_orthogonal_coefficient(
            Vec3{1.0, 0.0, 0.0}, dvec, dAU_diagonal);
        EXPECT_NEAR(orthogonal, 1.0, 1e-12);
    });

    return 0;
}
