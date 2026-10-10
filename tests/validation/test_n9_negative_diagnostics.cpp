// N9.7 — Negative diagnostics campaign (issue #601, item N9.7).
//
// Each diagnostic deliberately injects a known failure mode into the
// production numerical machinery and asserts that the machinery DETECTS
// it: through an explicit rejection (exception with a pinned diagnostic),
// a descriptive classification (matrix pathologies, solver failure
// classes), or a conservative non-convergence verdict. No gate is ever
// relaxed to make a case pass: a failure that would previously have
// surfaced as a silent false green must surface here as a red.
//
// The ten failure modes are the ones listed in #601 N9.7:
//   1. singular pressure system (all-Neumann)
//   2. wrong claimed null space
//   3. incompatible right-hand side
//   4. non-positive momentum coefficients / invalid controls
//   5. broken face-flux antisymmetry
//   6. inconsistent gradient/divergence operators
//   7. perturbed Rhie-Chow pressure reconstruction
//   8. ill-conditioned Schur complement
//   9. inner linear solve failure
//  10. nonlinear stagnation
// plus a final gate proving the same machinery still solves the reference
// Poiseuille case to the unchanged physical gates (no false green).

#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/matrix_diagnostics.h"
#include "cfdx/core/linalg/null_space.h"
#include "cfdx/core/linalg/simplerc_schur.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {

// ---------------------------------------------------------------------------
// Shared mesh/BC helpers (exact mirror of test_n9_physical_matrix.cpp so the
// negative campaign and the physical matrix exercise the same discretisation).
// ---------------------------------------------------------------------------

Mesh make_channel_mesh(std::size_t nx, std::size_t ny, double skew)
{
    Mesh mesh;
    const std::size_t plane = (nx + 1) * (ny + 1);
    mesh.points().resize(2 * plane);
    const auto id = [nx](std::size_t i, std::size_t j, std::size_t k) {
        return (j * (nx + 1) + i) * 2 + k;
    };

    for (std::size_t j = 0; j <= ny; ++j) {
        const double y = static_cast<double>(j) / static_cast<double>(ny);
        const double sx = skew * y * (1.0 - y);
        for (std::size_t i = 0; i <= nx; ++i) {
            const double x = static_cast<double>(i) / static_cast<double>(nx) + sx;
            mesh.points().set(id(i,j,0), x, y, 0.0);
            mesh.points().set(id(i,j,1), x, y, 1.0);
        }
    }

    std::map<std::vector<std::size_t>, std::size_t> face_map;
    std::vector<std::vector<std::size_t>> cell_faces(nx * ny);

    auto add_face = [&](std::initializer_list<std::size_t> vertices,
                        std::size_t cell) {
        std::vector<std::size_t> key(vertices);
        std::sort(key.begin(), key.end());
        const auto it = face_map.find(key);
        if (it != face_map.end()) {
            mesh.ownership().set_neighbour(it->second, static_cast<int>(cell));
            return it->second;
        }
        const std::size_t f = mesh.faces().n_faces();
        mesh.faces().push_face(vertices);
        mesh.ownership().resize(mesh.faces().n_faces());
        mesh.ownership().set_owner(f, cell);
        mesh.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        face_map.emplace(std::move(key), f);
        return f;
    };

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const std::size_t c = j * nx + i;
            const auto a=id(i,j,0), b=id(i+1,j,0), c0=id(i+1,j+1,0), d=id(i,j+1,0);
            const auto e=id(i,j,1), f=id(i+1,j,1), g=id(i+1,j+1,1), h=id(i,j+1,1);
            cell_faces[c] = {
                add_face({a,d,c0,b},c), add_face({e,f,g,h},c),
                add_face({a,b,f,e},c), add_face({d,h,g,c0},c),
                add_face({a,e,h,d},c), add_face({b,c0,g,f},c)};
        }
    }
    for (const auto& faces : cell_faces) mesh.cells().push_cell(faces);

    Patch inlet{"inlet", PatchType::INLET, {}};
    Patch outlet{"outlet", PatchType::OUTLET, {}};
    Patch bottom{"bottom", PatchType::WALL, {}};
    Patch top{"top", PatchType::WALL, {}};
    Patch front{"front", PatchType::EMPTY, {}};
    Patch back{"back", PatchType::EMPTY, {}};

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f) >= 0) continue;
        const auto& v = mesh.faces().vertices();
        const auto begin = v.begin() + static_cast<std::ptrdiff_t>(mesh.faces().face_offset(f));
        const auto end = begin + static_cast<std::ptrdiff_t>(mesh.faces().face_size(f));
        double x=0.0,y=0.0,z=0.0;
        for (auto it=begin; it!=end; ++it) {
            x += mesh.points().x(*it); y += mesh.points().y(*it); z += mesh.points().z(*it);
        }
        const double n = static_cast<double>(mesh.faces().face_size(f));
        x/=n; y/=n; z/=n;
        constexpr double eps=1e-12;
        if (std::abs(y) < eps) bottom.face_ids.push_back(f);
        else if (std::abs(y-1.0) < eps) top.face_ids.push_back(f);
        else if (std::abs(z) < eps) front.face_ids.push_back(f);
        else if (std::abs(z-1.0) < eps) back.face_ids.push_back(f);
        else {
            bool at_left=true, at_right=true;
            for (auto it=begin; it!=end; ++it) {
                const double xx=mesh.points().x(*it);
                const double yy=mesh.points().y(*it);
                const double logical = xx - skew*yy*(1.0-yy);
                at_left = at_left && std::abs(logical) < eps;
                at_right = at_right && std::abs(logical-1.0) < eps;
            }
            if (at_left) inlet.face_ids.push_back(f);
            else if (at_right) outlet.face_ids.push_back(f);
            else throw std::runtime_error("unclassified channel boundary face");
        }
    }

    mesh.boundary().add_patch(inlet);
    mesh.boundary().add_patch(outlet);
    mesh.boundary().add_patch(bottom);
    mesh.boundary().add_patch(top);
    mesh.boundary().add_patch(front);
    mesh.boundary().add_patch(back);
    return mesh;
}

struct Run {
    Field<double, Location::CELL> U;
    Field<double, Location::CELL> p;
    IncompressibleSolveResult result;
};

IncompressibleSolverControls controls_for(
    PressureVelocityAlgorithm algorithm,
    CoupledSchurModel schur_model = CoupledSchurModel::PCD,
    double viscosity = 0.1)
{
    IncompressibleSolverControls c;
    c.algorithm = algorithm;
    c.density = 1.0;
    c.kinematic_viscosity = viscosity;
    c.coupling.alpha_u = 0.7;
    // SIMPLEC's consistent pressure correction is normally unrelaxed; SIMPLE
    // retains the 0.3 pressure under-relaxation.
    c.coupling.alpha_p = algorithm == PressureVelocityAlgorithm::SIMPLEC ? 1.0 : 0.3;
    c.coupling.n_pressure_correctors =
        (algorithm == PressureVelocityAlgorithm::SIMPLEC ||
         algorithm == PressureVelocityAlgorithm::PISO ||
         algorithm == PressureVelocityAlgorithm::PIMPLE) ? 2 : 1;
    c.coupling.n_outer_correctors =
        algorithm == PressureVelocityAlgorithm::PIMPLE ? 2 : 1;
    c.coupling.n_fractional_steps =
        algorithm == PressureVelocityAlgorithm::FRACTIONAL_STEP ? 2 : 1;
    c.coupling.coupled_max_iterations = 1000;
    c.coupling.schur_model = schur_model;
    c.coupling.coupled_linear_tolerance = 1e-12;
    c.convergence.max_iterations = 1500;
    c.convergence.relative_tolerance = 1e-8;
    c.convergence.continuity_tolerance = 1e-8;
    c.linear_max_iterations = 1000;
    c.linear_tolerance = 1e-9;
    c.pressure_reference_cell = 0;
    c.pressure_reference_value = 0.0;
    c.use_bounded_convection = true;
    c.convection_scheme = ConvectionScheme::UPWIND;
    return c;
}

VelocityBoundaryConditions channel_velocity_bc(double top_u)
{
    VelocityBoundaryConditions bc;
    bc["inlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    bc["outlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    bc["bottom"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    bc["top"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{top_u,0,0}};
    bc["front"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    bc["back"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    return bc;
}

ScalarBoundaryConditions channel_pressure_bc()
{
    ScalarBoundaryConditions bc;
    for (const char* n : {"inlet","outlet","bottom","top","front","back"})
        bc[n] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
    return bc;
}

Run solve_case(
    Mesh mesh,
    PressureVelocityAlgorithm algorithm,
    const VelocityBoundaryConditions& ubc,
    const ScalarBoundaryConditions& pbc,
    double body_force_x,
    double viscosity = 0.1)
{
    Field<double, Location::CELL> U(mesh.n_cells(),"U","m/s",3);
    Field<double, Location::CELL> p(mesh.n_cells(),"p","Pa",1);
    U.fill(0.0);
    p.fill(17.0);

    auto c = controls_for(algorithm, CoupledSchurModel::PCD, viscosity);
    c.body_force = {body_force_x,0.0,0.0};
    const auto result = solve_steady_incompressible(mesh,U,p,ubc,pbc,c);
    return {std::move(U),std::move(p),std::move(result)};
}

double poiseuille_l2(const Run& r, std::size_t nx, std::size_t ny, double G, double nu)
{
    double e2=0.0;
    for (std::size_t c=0;c<r.U.size();++c) {
        const std::size_t j=c/nx;
        const double y=(static_cast<double>(j)+0.5)/static_cast<double>(ny);
        const double exact=G*y*(1.0-y)/(2.0*nu);
        const double e=r.U.component_data(0)[c]-exact;
        e2 += e*e;
    }
    return std::sqrt(e2/static_cast<double>(r.U.size()));
}

// ---------------------------------------------------------------------------
// Negative-campaign helpers.
// ---------------------------------------------------------------------------

SparseMatrix make_sparse(
    std::size_t rows, std::size_t cols,
    const std::initializer_list<std::tuple<std::size_t, std::size_t, double>>& e) {
    SparseMatrix A(rows, cols);
    for (const auto& [i, j, v] : e) A.push_back(i, j, v);
    A.finalize();
    return A;
}

// 1-D finite-volume pressure operator with purely Neumann boundaries: the
// exact campaign-scale shape of an ungauged pressure system. Every row sums
// to zero, so the constant vector spans its null space exactly.
SparseMatrix make_all_neumann_pressure_operator(std::size_t n)
{
    SparseMatrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        const double diagonal = (i == 0 || i + 1 == n) ? 1.0 : 2.0;
        A.push_back(i, i, diagonal);
        if (i > 0) A.push_back(i, i - 1, -1.0);
        if (i + 1 < n) A.push_back(i, i + 1, -1.0);
    }
    A.finalize();
    return A;
}

// Max-norm true residual ||b - A x||_inf computed from raw CSR storage: the
// campaign never trusts a solver-reported residual as its own evidence.
double true_residual(const SparseMatrix& A, const Vector& b, const Vector& x)
{
    const auto* row = A.row_offsets_data();
    const auto* col = A.columns_data();
    const auto* val = A.values_data();
    double worst = 0.0;
    for (std::size_t i = 0; i < A.n_rows(); ++i) {
        double r = -b(i);
        for (auto k = row[i]; k < row[i + 1]; ++k)
            r += val[k] * x(col[k]);
        worst = std::max(worst, std::abs(r));
    }
    return worst;
}

// Per-cell face-flux balance (discrete divergence). With antisymmetric=true
// the face flux leaves the owner with a + sign and enters the neighbour with
// a - sign, which is the only assembly a conservative operator may use.
// antisymmetric=false emulates the classic sign fault where both cells
// accumulate the flux with the same sign. flipped_face negates one internal
// face's flux before the balance is formed.
double cell_mass_balance(
    const Mesh& mesh,
    const Field<double, Location::FACE>& flux,
    bool antisymmetric,
    std::size_t flipped_face = static_cast<std::size_t>(-1))
{
    const auto& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* offsets = cells.offsets_data();
    const auto& own = mesh.ownership();
    double worst = 0.0;
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        double balance = 0.0;
        for (auto k = offsets[c]; k < offsets[c + 1]; ++k) {
            const std::size_t f = cell_faces[k];
            double phi = flux(f);
            if (f == flipped_face) phi = -phi;
            if (own.owner(f) == c)
                balance += phi;
            else
                balance += antisymmetric ? -phi : phi;
        }
        worst = std::max(worst, std::abs(balance));
    }
    return worst;
}

// The reference Poiseuille case solved once and shared by item 5 and the
// final gate: the negative campaign must prove that the machinery that
// rejects the corrupted cases still produces the accepted reference.
struct ReferenceCase {
    Mesh mesh;
    Run run;
};

const ReferenceCase& reference_case()
{
    static const ReferenceCase reference = [] {
        Mesh mesh = make_channel_mesh(12, 16, 0.0);
        Run run = solve_case(
            mesh, PressureVelocityAlgorithm::SIMPLE,
            channel_velocity_bc(0.0), channel_pressure_bc(), 1.0);
        return ReferenceCase{std::move(mesh), std::move(run)};
    }();
    return reference;
}

} // namespace

int main()
{
    // ------------------------------------------------------------------
    // Item 1 — singular pressure system: the null space is detected, the
    // compatible system is solved on the projected space, and the solution
    // carries no null-space component.
    // ------------------------------------------------------------------
    run_case("n9_7_item1_singular_pressure_system", [] {
        const auto A = make_all_neumann_pressure_operator(8);
        const auto constant = NullSpaceProjector::constant(8);
        EXPECT_TRUE(constant.is_null_space(A));

        Vector x_true(8, 0.0);
        for (std::size_t i = 0; i < 8; ++i)
            x_true(i) = static_cast<double>(i) - 3.5;
        Vector b(8, 0.0);
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        for (std::size_t i = 0; i < A.n_rows(); ++i)
            for (auto k = row[i]; k < row[i + 1]; ++k)
                b(i) += val[k] * x_true(col[k]);
        EXPECT_TRUE(constant.is_compatible(b));

        Vector x(8, 0.0);
        const auto result = solve_cg(A, b, x, constant, 200, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(true_residual(A, b, x) < 1e-9);
        EXPECT_TRUE(constant.component_norm(x) < 1e-12);
        std::cout << "N9_7_ITEM1 true_residual=" << true_residual(A, b, x)
                  << " null_component=" << constant.component_norm(x) << "\n";
    });

    // ------------------------------------------------------------------
    // Item 2 — wrongly claimed null space: a checkerboard mode is NOT in
    // the null space of the Neumann operator; the solver must refuse the
    // projected solve instead of silently projecting onto a wrong space.
    // ------------------------------------------------------------------
    run_case("n9_7_item2_wrong_null_space_rejected", [] {
        const auto A = make_all_neumann_pressure_operator(8);
        Vector checkerboard(8, 0.0);
        for (std::size_t i = 0; i < 8; ++i)
            checkerboard(i) = (i % 2 == 0) ? 1.0 : -1.0;
        const NullSpaceProjector wrong({checkerboard});
        EXPECT_FALSE(wrong.is_null_space(A));
        EXPECT_TRUE(wrong.operator_residual(A) > 1e-6);

        Vector b(8, 1.0);
        Vector x(8, 0.0);
        const auto result = solve_cg(A, b, x, wrong, 100, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
        std::cout << "N9_7_ITEM2 wrong_null_operator_residual="
                  << wrong.operator_residual(A)
                  << " solve_status=" << to_string(result.status) << "\n";
    });

    // ------------------------------------------------------------------
    // Item 3 — incompatible right-hand side: a non-zero-mean rhs cannot be
    // solved on the projected space; the projected solver refuses it and the
    // failure classifier attributes it to the rhs, not to the matrix.
    // ------------------------------------------------------------------
    run_case("n9_7_item3_incompatible_rhs_classified", [] {
        const auto A = make_all_neumann_pressure_operator(8);
        const auto constant = NullSpaceProjector::constant(8);
        Vector b(8, 1.0);
        EXPECT_FALSE(constant.is_compatible(b));

        Vector x(8, 0.0);
        const auto result = solve_cg(A, b, x, constant, 100, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);

        const auto diagnostics = diagnose_matrix(A);
        EXPECT_TRUE(classify_solver_failure(
                        result, diagnostics, /*rhs_compatible=*/false) ==
                    SolverFailureClass::IncompatibleRhs);
        std::cout << "N9_7_ITEM3 rhs_component_norm="
                  << constant.component_norm(b)
                  << " classified="
                  << to_string(classify_solver_failure(
                        result, diagnostics, /*rhs_compatible=*/false)) << "\n";
    });

    // ------------------------------------------------------------------
    // Item 4 — non-positive momentum coefficients and invalid controls are
    // rejected with pinned diagnostics, never silently repaired.
    // ------------------------------------------------------------------
    run_case("n9_7_item4_invalid_coefficients_rejected", [] {
        auto bad_density = controls_for(PressureVelocityAlgorithm::SIMPLE);
        bad_density.density = -1.0;
        EXPECT_THROW_WITH(
            validate_incompressible_controls(bad_density, 16),
            std::invalid_argument, "invalid incompressible material properties");

        auto bad_viscosity = controls_for(PressureVelocityAlgorithm::SIMPLE);
        bad_viscosity.kinematic_viscosity = -0.1;
        EXPECT_THROW_WITH(
            validate_incompressible_controls(bad_viscosity, 16),
            std::invalid_argument, "invalid incompressible material properties");

        auto bad_relaxation = controls_for(PressureVelocityAlgorithm::SIMPLE);
        bad_relaxation.coupling.alpha_u = 0.0;
        EXPECT_THROW_WITH(
            validate_incompressible_controls(bad_relaxation, 16),
            std::invalid_argument, "invalid relaxation/corrector count");

        EXPECT_THROW_WITH(
            simplec_consistent_diagonal(-1.0, 0.0),
            std::invalid_argument, "invalid momentum coefficients");
        EXPECT_THROW_WITH(
            simplec_consistent_diagonal(1.0, -2.0),
            std::invalid_argument, "non-positive/non-finite denominator");
        EXPECT_THROW_WITH(
            piso_correction_gain(0.0, 0.0),
            std::invalid_argument, "diagonal must be positive");

        const Mesh mesh = make_channel_mesh(2, 2, 0.0);
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        U.fill(0.0);
        p.fill(0.0);
        std::array<std::vector<double>, 3> rAU;
        for (auto& component : rAU) component.assign(mesh.n_cells(), 1.0);
        rAU[2].assign(mesh.n_cells(), 0.0);
        EXPECT_THROW_WITH(
            make_rhie_chow_mass_flux(mesh, geometry, U, p, rAU, 1.0, {}),
            std::invalid_argument,
            "inverse momentum diagonal must be finite and positive");
        std::cout << "N9_7_ITEM4 all pinned rejections thrown\n";
    });

    // ------------------------------------------------------------------
    // Item 5 — face-flux antisymmetry: the authoritative converged flux is
    // balanced; a non-antisymmetric assembly or one flipped internal face
    // destroys the balance and must be measurable, not invisible.
    // ------------------------------------------------------------------
    run_case("n9_7_item5_flux_antisymmetry_violation_detected", [] {
        const auto& reference = reference_case();
        const auto& flux = reference.run.result.authoritative_mass_flux;
        const double clean = cell_mass_balance(reference.mesh, flux, true);
        EXPECT_TRUE(clean <= 1e-5);

        const double sign_fault = cell_mass_balance(reference.mesh, flux, false);
        EXPECT_TRUE(sign_fault > 1e-2);

        const auto& own = reference.mesh.ownership();
        std::size_t worst_face = 0;
        double worst_flux = 0.0;
        for (std::size_t f = 0; f < reference.mesh.n_faces(); ++f) {
            if (own.neighbour(f) < 0) continue;
            if (std::abs(flux(f)) > worst_flux) {
                worst_flux = std::abs(flux(f));
                worst_face = f;
            }
        }
        const double flipped =
            cell_mass_balance(reference.mesh, flux, true, worst_face);
        EXPECT_TRUE(flipped > 1e-2);
        EXPECT_TRUE(flipped > 100.0 * clean);
        std::cout << "N9_7_ITEM5 clean=" << clean
                  << " sign_fault=" << sign_fault
                  << " flipped=" << flipped
                  << " flipped_face=" << worst_face << "\n";
    });

    // ------------------------------------------------------------------
    // Item 6 — inconsistent D/G operators: the production Gauss gradient is
    // linear-exact on interior cells and the production divergence is
    // conservative for a constant field; an interpolation that collapses to
    // the owner value and drops the owner/neighbour orientation, or a
    // divergence assembly that forgets the antisymmetric sign, produces a
    // measurable inconsistency.
    // ------------------------------------------------------------------
    run_case("n9_7_item6_gradient_divergence_consistency", [] {
        constexpr std::size_t nx = 12, ny = 16;
        const Mesh mesh = make_channel_mesh(nx, ny, 0.0);
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            p(c) = geometry.cell_centres[c].x;
        const auto grad = compute_gradient_gauss(p, mesh);

        double exact_error = 0.0;
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t i = 1; i + 1 < nx; ++i) {
                const std::size_t c = j * nx + i;
                exact_error = std::max(
                    exact_error, std::abs(grad.component_data(0)[c] - 1.0));
                exact_error = std::max(
                    exact_error, std::abs(grad.component_data(1)[c]));
                exact_error = std::max(
                    exact_error, std::abs(grad.component_data(2)[c]));
            }
        EXPECT_TRUE(exact_error <= 1e-12);

        // Collapsed operator: face value = owner value and every face area
        // vector used with the stored (owner) orientation for both cells.
        const auto& cells = mesh.cells();
        const auto* cell_faces = cells.faces_data();
        const auto* offsets = cells.offsets_data();
        const auto& own = mesh.ownership();
        double collapsed_error = 0.0;
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t i = 1; i + 1 < nx; ++i) {
                const std::size_t c = j * nx + i;
                double sum_x = 0.0;
                for (auto k = offsets[c]; k < offsets[c + 1]; ++k) {
                    const std::size_t f = cell_faces[k];
                    const std::size_t o = own.owner(f);
                    sum_x += geometry.face_area_vectors[f].x * p(o);
                }
                const double gx = sum_x / geometry.cell_volumes[c];
                collapsed_error = std::max(collapsed_error, std::abs(gx - 1.0));
            }
        EXPECT_TRUE(collapsed_error > 1e-3);

        // Production divergence of a constant field: exactly conservative.
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        U.fill(0.0);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            U.set(c, 0.37, 0.0, 0.0);
        const auto flux = make_mass_flux(mesh, geometry, U, 1.0, VelocityBoundaryConditions{});
        const double clean = cell_mass_balance(mesh, flux, true);
        EXPECT_TRUE(clean < 1e-12);
        const double sign_fault = cell_mass_balance(mesh, flux, false);
        EXPECT_TRUE(sign_fault > 1e-3);
        std::cout << "N9_7_ITEM6 exact_gradient_error=" << exact_error
                  << " collapsed_gradient_error=" << collapsed_error
                  << " divergence_clean=" << clean
                  << " divergence_sign_fault=" << sign_fault << "\n";
    });

    // ------------------------------------------------------------------
    // Item 7 — perturbed Rhie-Chow reconstruction: the face flux is invariant
    // under a constant pressure offset (checkerboard suppression intact)
    // but responds measurably to a localized pressure perturbation.
    // ------------------------------------------------------------------
    run_case("n9_7_item7_rhie_chow_reconstruction_perturbed", [] {
        const Mesh mesh = make_channel_mesh(12, 16, 0.0);
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        U.fill(0.0);
        Field<double, Location::CELL> p0(mesh.n_cells(), "p", "Pa", 1);
        p0.fill(17.0);
        std::array<std::vector<double>, 3> rAU;
        for (auto& component : rAU) component.assign(mesh.n_cells(), 1.0);
        const auto ubc = channel_velocity_bc(0.0);
        const auto pbc = channel_pressure_bc();

        const auto phi0 = make_rhie_chow_mass_flux(mesh, geometry, U, p0, rAU, 1.0, ubc, pbc);

        Field<double, Location::CELL> p1(mesh.n_cells(), "p", "Pa", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            p1(c) = p0(c) + 1000.0;
        const auto phi1 = make_rhie_chow_mass_flux(mesh, geometry, U, p1, rAU, 1.0, ubc, pbc);
        double offset_delta = 0.0;
        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
            offset_delta = std::max(offset_delta, std::abs(phi0(f) - phi1(f)));
        EXPECT_TRUE(offset_delta <= 1e-12);

        Field<double, Location::CELL> p2(mesh.n_cells(), "p", "Pa", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            p2(c) = p0(c);
        p2(3 * 8 + 3) += 1.0;
        const auto phi2 = make_rhie_chow_mass_flux(mesh, geometry, U, p2, rAU, 1.0, ubc, pbc);
        double bump_delta = 0.0;
        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
            bump_delta = std::max(bump_delta, std::abs(phi0(f) - phi2(f)));
        EXPECT_TRUE(bump_delta > 1e-3);
        std::cout << "N9_7_ITEM7 offset_delta=" << offset_delta
                  << " bump_delta=" << bump_delta << "\n";
    });

    // ------------------------------------------------------------------
    // Item 8 — Schur complement conditioning: healthy blocks set up under
    // both SIMPLE and SIMPLEC; a momentum row whose SIMPLEC consistent
    // denominator collapses to zero is rejected at setup; an ill-conditioned
    // momentum matrix is classified as a matrix pathology, not solved.
    // ------------------------------------------------------------------
    run_case("n9_7_item8_schur_conditioning_and_pathology", [] {
        const auto Auu = make_sparse(2, 2, {{0,0,4.0},{0,1,-1.0},{1,0,-1.0},{1,1,3.0}});
        const auto G = make_sparse(2, 2, {{0,0,1.0},{1,1,1.0}});
        const auto D = make_sparse(2, 2, {{0,0,1.0},{1,1,1.0}});
        const auto C = make_sparse(2, 2, {{0,0,2.0},{1,1,3.0}});
        const BlockOperator blocks(Auu, G, D, C);
        blocks.validate();

        SimplerSchurApproximation simple(SimplerSchurMode::SIMPLE);
        EXPECT_TRUE(simple.setup(blocks));
        SimplerSchurApproximation simplec(SimplerSchurMode::SIMPLEC);
        EXPECT_TRUE(simplec.setup(blocks));

        const auto Auu_degenerate =
            make_sparse(2, 2, {{0,0,1.0},{0,1,-1.0},{1,0,-1.0},{1,1,1.0}});
        const BlockOperator degenerate(Auu_degenerate, G, D, C);
        SimplerSchurApproximation simple_on_degenerate(SimplerSchurMode::SIMPLE);
        EXPECT_TRUE(simple_on_degenerate.setup(degenerate));
        SimplerSchurApproximation simplec_on_degenerate(SimplerSchurMode::SIMPLEC);
        EXPECT_FALSE(simplec_on_degenerate.setup(degenerate));

        const auto Auu_ill = make_sparse(2, 2, {{0,0,4.0},{1,1,4e-15}});
        const auto diagnostics = diagnose_matrix(Auu_ill, 1e-14);
        EXPECT_TRUE(diagnostics.near_zero_diagonal >= 1);
        EXPECT_TRUE(diagnostics.diagonal_dynamic_range > 1e12);
        const auto pathologies = classify_matrix_pathologies(diagnostics);
        EXPECT_FALSE(pathologies.size() == 1 &&
                     pathologies.front() == MatrixPathology::None);

        SolverResult failed;
        failed.status = SolverStatus::MAX_ITER_REACHED;
        failed.residual = 1.0;
        failed.residual_relative = 1.0;
        EXPECT_TRUE(classify_solver_failure(failed, diagnostics) ==
                    SolverFailureClass::MatrixPathology);
        std::cout << "N9_7_ITEM8 near_zero_diagonal="
                  << diagnostics.near_zero_diagonal
                  << " diagonal_dynamic_range="
                  << diagnostics.diagonal_dynamic_range << "\n";
    });

    // ------------------------------------------------------------------
    // Item 9 — inner linear failure: an unconverged momentum Krylov solve
    // must abort the nonlinear solve with the pinned diagnostic instead of
    // injecting an unconverged predictor. The LU/fallback recovery paths
    // are disabled explicitly so the failure is isolated.
    // ------------------------------------------------------------------
    run_case("n9_7_item9_linear_failure_propagates", [] {
        auto c = controls_for(PressureVelocityAlgorithm::SIMPLE);
        c.linear_max_iterations = 1;
        c.momentum_linear_solver.krylov = KrylovModel::GMRES;
        c.momentum_linear_solver.preconditioner = PreconditionerModel::Jacobi;
        c.momentum_linear_solver.allow_fallback = false;

        const Mesh mesh = make_channel_mesh(4, 4, 0.0);
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        U.fill(0.0);
        p.fill(17.0);
        EXPECT_THROW_WITH(
            solve_steady_incompressible(
                mesh, U, p, channel_velocity_bc(1.0), channel_pressure_bc(), c),
            std::runtime_error, "momentum solve did not converge");
        std::cout << "N9_7_ITEM9 linear failure propagated with pinned diagnostic\n";
    });

    // ------------------------------------------------------------------
    // Item 10 — nonlinear stagnation: a pathological momentum
    // under-relaxation (alpha_u = 1e-6) makes each nonlinear iteration move
    // the solution by ~1e-6 of the full correction. On this linear case
    // the outer-iteration scaling law with essentially exact inner Krylov
    // solves contracts the residual by ~(1 - alpha_u) per iteration, so
    // the per-iteration relative residual improvement (~1e-6) stays three
    // orders of magnitude below the production stagnation threshold
    // (1e-3 per iteration over the 10-iteration window). The production
    // monitor MUST terminate the solve with the STAGNATED verdict and a
    // non-empty reason, early — never a false converged and never an
    // unbounded iteration loop.
    // ------------------------------------------------------------------
    run_case("n9_7_item10_nonlinear_stagnation_detected", [] {
        auto c = controls_for(PressureVelocityAlgorithm::SIMPLEC);
        c.coupling.alpha_u = 1e-6;
        c.convergence.max_iterations = 400;

        const Mesh mesh = make_channel_mesh(12, 16, 0.0);
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        U.fill(0.0);
        p.fill(17.0);
        const auto result = solve_steady_incompressible(
            mesh, U, p, channel_velocity_bc(1.0), channel_pressure_bc(), c);
        EXPECT_FALSE(result.converged);
        EXPECT_TRUE(result.convergence_status == ConvergenceStatus::STAGNATED);
        EXPECT_FALSE(result.convergence_reason.empty());
        // Terminated by the stagnation detector, not by chance or by the
        // iteration cap: the 10-iteration window cannot fire before
        // iteration 11, and a detector termination must arrive long before
        // the 400-iteration maximum.
        EXPECT_TRUE(result.iterations >= 11);
        EXPECT_TRUE(result.iterations < 400);
        std::cout << "N9_7_ITEM10 status=" << to_string(result.convergence_status)
                  << " iterations=" << result.iterations
                  << " reason=" << result.convergence_reason << "\n";
    });

    // ------------------------------------------------------------------
    // Final gate — no false green: the same machinery that rejected the ten
    // corrupted cases still solves the reference Poiseuille case to the
    // unchanged physical gates.
    // ------------------------------------------------------------------
    run_case("n9_7_final_gate_reference_solution_still_accepted", [] {
        const auto& reference = reference_case();
        EXPECT_TRUE(reference.run.result.converged);
        EXPECT_FALSE(reference.run.result.history.empty());
        const auto& h = reference.run.result.history.back();
        EXPECT_TRUE(h.continuity_linf <= 1e-7);
        EXPECT_TRUE(h.momentum_residual <= 1e-7);
        const double l2 = poiseuille_l2(reference.run, 12, 16, 1.0, 0.1);
        EXPECT_TRUE(l2 <= 5e-3);
        std::cout << "N9_7_FINAL continuity_linf=" << h.continuity_linf
                  << " momentum_residual=" << h.momentum_residual
                  << " poiseuille_l2=" << l2 << "\n";
    });

    return run_all();
}
