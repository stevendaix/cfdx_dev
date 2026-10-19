// PCD Schur-action A/B benchmark on the REAL solved Ghia coupled systems.
//
// Purpose: measure, on the actual 32x32 Re=100 COUPLED+PCD cavity systems
// (not synthetic representative matrices), where the per-stage cost goes and
// what a stronger Schur action would buy. The measured question is the one
// left open by the pcd_schur.h contract note: the single-inverse PCD action
// -Fp^-1 is PROVISIONAL and was never qualified on a case with a non-zero
// mass flux; the Ghia Re=100 qualification is exactly such a case.
//
// Method:
//   1. Run the production solver on the 32x32 Re=100 cavity (identical mesh,
//      BCs and controls to test_ghia_cavity.cpp --pcd) for a few nonlinear
//      iterations, with the evidence-only coupled_assembly_callback capturing
//      the EXACT system handed to the Krylov solver. The hook fires at the
//      Krylov call point - after assembly, gauge treatment, finite checks
//      and preconditioner setup - so the captured matrix, RHS and initial
//      guess are the very objects production solves (no reconstruction, no
//      replication drift). The capture also carries the production PCD
//      pressure operators (gauge row already imposed) and the Rhie-Chow
//      Schur diagonal in both raw and gauge-pinned forms.
//   2. Replay two captured stages (assemblies #2 and #5: the early regime
//      where the measured production per-stage cost peaks) through the
//      documented PCD action variants, under identical conditions:
//      same captured matrix, same captured RHS, same captured production
//      initial guess, same outer FGMRES (restart 4096, budget 5000,
//      tolerance 1e-10, fixed KrylovControls).
//
// Validity guards (what makes this test PASS - they establish that the
// benchmark replays the production system; they do NOT establish any
// performance claim):
//   - at least five coupled assemblies were captured;
//   - the captured coupled system is gauge-coherent with the captured
//     state: the reference row is an identity row with reference RHS, and
//     the captured initial guess carries the reference value on the gauge
//     row (mirrors the production x0 seeding);
//   - the captured production PCD operators are reproduced by the
//     documented assembly path: rebuilding mass/Kp/Fp from the captured
//     linearization velocity via assemble_pcd_pressure_operators must
//     match the captured operators to 1e-12 (the replay itself uses the
//     captured production operators, so this check verifies the assembly
//     path rather than substituting for it);
//   - the baseline fp_ilu0 action converges every benchmarked stage to
//     the production 1e-10 true relative residual from the captured
//     production initial guess (known production behaviour: every stage
//     of the 32x32 qualification reaches ~9.8e-11);
//   - every converged variant agrees on the solution to 1e-6 relative
//     (a linear-solve consistency guard: no variant may "converge" to a
//     different solution).
// A non-baseline variant failing to converge is REPORTED, not failed: the
// measurement is the point. The test verdict below therefore separates
// benchmark VALIDITY (the guards above) from the PERFORMANCE report
// (printed measurements, never gating).
//
// Variants (all built from the existing PcdSchurApproximation class; NO
// production pcd_schur.h change is involved). The callback contract is
// audited at runtime: the laplacian_solve slot receives a callback that
// throws, so any unexpected invocation of the Kp slot inside the
// single-inverse action fails the solve loudly instead of silently
// measuring a different operator:
//   fp_ilu0   -Fp^-1, nested FGMRES(512)+ILU0, budget 200, tol 1e-8
//              (exactly the production --pcd configuration = baseline).
//   fp_amg    -Fp^-1, nested FGMRES(512)+NativeAMG, budget 200, tol 1e-8
//              (the nested-preconditioner axis of PR #965, measured on the
//              real matrix).
//   kp_amg    -Kp^-1, nested GMRES(40)+NativeAMG, budget 100, tol 1e-8
//              (the provisional alternative documented in pcd_schur.h; Kp
//              defaults used because that is how production would solve it).
//   fp_rc_scaled  -diag(RC)*Fp^-1 (scaling by the exact gauge-pinned
//              Rhie-Chow diagonal captured from the assembly). Measurement
//              variant: the algebraic derivation (left vs right scaling,
//              signs, gauge treatment) is NOT asserted here; the benchmark
//              measures whether the amplitude correction helps at all.
//   rc_diag   -diag(RC)^-1 with no nested solve (SIMPLE-style diagonal
//              Schur). Reference point; nested Krylov counters stay at
//              zero but Schur applications are counted and timed.
//   pcd_full  -Kp^-1*Fp*Mp^-1, the classic literature PCD action (Kay,
//              Loghin & Wathen; Elman, Silvester & Wathen 2001,
//              doi:10.1016/S0377-0427(00)00515-4). Composed inside the
//              convection_diffusion slot: one nested Kp solve
//              (GMRES(40)+NativeAMG, budget 100, tol 1e-8, mirroring
//              kp_amg because that is how production would solve Kp)
//              plus one forward Fp product plus the exact diagonal
//              Mp^-1 (Mp = diag(cell_volumes) by assembly, verified
//              fail-closed here instead of assumed). One inner
//              inverse, the same structural cost as production, NOT
//              the double-inverse compositions documented as
//              non-convergent in pcd_schur.h. Measurement variant:
//              the ordering, signs and gauge behaviour are NOT
//              asserted - the benchmark measures whether the classic
//              composition converges the real captured coupled system.
//
//   fp_ilu0_loose / fp_amg_loose  the NESTED TOLERANCE axis: the same
//              nested solvers as fp_ilu0 / fp_amg, but the nested
//              tolerance is 1e-4 instead of 1e-8 - the classic
//              preconditioner-level tolerance. Measurement question:
//              the outer coupled FGMRES is flexible precisely to
//              tolerate loose inner solves, so how many nested
//              iterations and wall ms does loosening save per
//              application, and does the coupled solve still converge?
//              The outer coupled tolerance (1e-10), the validity
//              guards and the baseline fp_ilu0 are unchanged; a loose
//              variant that fails to converge the coupled system is
//              reported, not sanctioned.
//
// Evidence printed (one line per solve, plus per-stage summaries):
//   PCD_ACTION_BENCHMARK stage=... action=... status=... converged=...
//   coupled_iterations=... true_rel_residual=... solve_wall_ms=...
//   schur_applications=... nested_iterations=... action_wall_ms=...
//   non_nested_wall_ms=... nested_setup_ms=... schur_setup_ms=...
//   PCD_ACTION_SUMMARY per-variant aggregates and baseline-relative
//   speedups follow the per-solve lines; the final verdict line is printed
//   only after all guards have been evaluated.
//
// Cost: bounded and reproducible - five production nonlinear iterations
// plus sixteen coupled solves plus two operator-assembly cross-checks, a few
// minutes in Release. Registered only under CFDX_ENABLE_LONG_VALIDATION.

#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

// Identical to test_ghia_cavity.cpp make_cavity_mesh (32x32 Re=100 cavity).
// Duplicated rather than shared so this benchmark cannot drift the production
// qualification test; the guards below establish replay fidelity instead.
Mesh make_cavity_mesh(std::size_t nx, std::size_t ny)
{
    if (nx < 4 || ny < 4) throw std::invalid_argument("cavity mesh is too small");
    Mesh mesh;
    mesh.points().resize((nx + 1) * (ny + 1) * 2);
    const auto point_id = [nx](std::size_t i, std::size_t j, std::size_t k) {
        return (j * (nx + 1) + i) * 2 + k;
    };
    for (std::size_t j=0;j<=ny;++j) for (std::size_t i=0;i<=nx;++i) {
        const double x=static_cast<double>(i)/static_cast<double>(nx);
        const double y=static_cast<double>(j)/static_cast<double>(ny);
        mesh.points().set(point_id(i,j,0),x,y,0.0);
        mesh.points().set(point_id(i,j,1),x,y,1.0);
    }

    std::map<std::vector<std::size_t>,std::size_t> face_map;
    std::vector<std::vector<std::size_t>> cell_faces(nx*ny);
    auto add_or_get_face=[&](std::initializer_list<std::size_t> vertices,std::size_t cell) {
        std::vector<std::size_t> sorted(vertices);
        std::sort(sorted.begin(),sorted.end());
        const auto found=face_map.find(sorted);
        if(found!=face_map.end()) {
            mesh.ownership().set_neighbour(found->second,static_cast<int>(cell));
            return found->second;
        }
        const std::size_t face=mesh.faces().n_faces();
        mesh.faces().push_face(vertices);
        mesh.ownership().resize(mesh.faces().n_faces());
        face_map.emplace(std::move(sorted),face);
        mesh.ownership().set_owner(face,cell);
        mesh.ownership().set_neighbour(face,FaceOwnership::BOUNDARY);
        return face;
    };

    for(std::size_t j=0;j<ny;++j) for(std::size_t i=0;i<nx;++i) {
        const std::size_t cell=j*nx+i;
        const std::size_t a=point_id(i,j,0), b=point_id(i+1,j,0);
        const std::size_t c=point_id(i+1,j+1,0), d=point_id(i,j+1,0);
        const std::size_t e=point_id(i,j,1), f=point_id(i+1,j,1);
        const std::size_t g=point_id(i+1,j+1,1), h=point_id(i,j+1,1);
        cell_faces[cell]={
            add_or_get_face({a,d,c,b},cell), add_or_get_face({e,f,g,h},cell),
            add_or_get_face({a,b,f,e},cell), add_or_get_face({d,h,g,c},cell),
            add_or_get_face({a,e,h,d},cell), add_or_get_face({b,c,g,f},cell)};
    }
    for(const auto& faces:cell_faces) mesh.cells().push_cell(faces);

    Patch bottom; bottom.name="bottom"; bottom.type=PatchType::WALL;
    Patch top; top.name="top"; top.type=PatchType::WALL;
    Patch left; left.name="left"; left.type=PatchType::WALL;
    Patch right; right.name="right"; right.type=PatchType::WALL;
    Patch front; front.name="front"; front.type=PatchType::WALL;
    Patch back; back.name="back"; back.type=PatchType::WALL;

    for(std::size_t f=0;f<mesh.n_faces();++f) {
        if(mesh.ownership().neighbour(f)>=0) continue;
        const auto& all_vertices=mesh.faces().vertices();
        const auto begin=all_vertices.begin()+static_cast<std::ptrdiff_t>(mesh.faces().face_offset(f));
        const auto end=begin+static_cast<std::ptrdiff_t>(mesh.faces().face_size(f));
        double x=0,y=0,z=0;
        for(auto it=begin;it!=end;++it) { x+=mesh.points().x(*it); y+=mesh.points().y(*it); z+=mesh.points().z(*it); }
        const double nverts=static_cast<double>(mesh.faces().face_size(f));
        x/=nverts; y/=nverts; z/=nverts;
        constexpr double tol=1e-12;
        if(std::abs(y)<tol) bottom.face_ids.push_back(f);
        else if(std::abs(y-1)<tol) top.face_ids.push_back(f);
        else if(std::abs(x)<tol) left.face_ids.push_back(f);
        else if(std::abs(x-1)<tol) right.face_ids.push_back(f);
        else if(std::abs(z)<tol) front.face_ids.push_back(f);
        else if(std::abs(z-1)<tol) back.face_ids.push_back(f);
        else throw std::runtime_error("cavity: unclassified boundary face");
    }
    mesh.boundary().add_patch(bottom); mesh.boundary().add_patch(top);
    mesh.boundary().add_patch(left); mesh.boundary().add_patch(right);
    mesh.boundary().add_patch(front); mesh.boundary().add_patch(back);
    return mesh;
}

double elapsed_ms(const std::chrono::steady_clock::time_point& start)
{
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();
}

// Independently recomputed ||b - A x|| / ||b||, mirroring the post-solve
// check in solve_coupled_momentum_continuity.
double true_relative_residual(const SparseMatrix& A, const Vector& b, const Vector& x)
{
    double residual_sq = 0.0;
    double rhs_sq = 0.0;
    for (std::size_t row = 0; row < A.n_rows(); ++row) {
        double ax = 0.0;
        const auto rb = A.row_offsets_data()[row];
        const auto re = A.row_offsets_data()[row + 1];
        for (std::size_t k = rb; k < re; ++k)
            ax += A.values_data()[k] * x(A.columns_data()[k]);
        const double r = b(row) - ax;
        residual_sq += r * r;
        rhs_sq += b(row) * b(row);
    }
    return std::sqrt(residual_sq) / (rhs_sq > 0.0 ? std::sqrt(rhs_sq) : 1.0);
}

// Largest elementwise |a-b| over the explicit pattern of two sparse
// matrices. Patterns are compared first: a pattern mismatch is a capture
// infidelity regardless of the values.
double sparse_max_abs_diff(const SparseMatrix& a, const SparseMatrix& b)
{
    if (a.n_rows() != b.n_rows() || a.nnz() != b.nnz())
        return std::numeric_limits<double>::infinity();
    double diff = 0.0;
    for (std::size_t k = 0; k < a.nnz(); ++k) {
        if (a.columns_data()[k] != b.columns_data()[k])
            return std::numeric_limits<double>::infinity();
        diff = std::max(diff, std::abs(a.values_data()[k] - b.values_data()[k]));
    }
    return diff;
}

// One captured coupled system: the EXACT objects handed to the production
// Krylov solver (captured at the Krylov call point), plus the production
// PCD pressure operators and both forms of the Rhie-Chow diagonal.
struct CapturedStage {
    SparseMatrix matrix;
    Vector rhs;
    Vector initial_guess;  // production x0: [U_old; p_old], gauge = reference value
    std::vector<double> rc_diagonal;       // gauge-pinned (reference cell = 1.0)
    std::vector<double> rc_diagonal_raw;    // as assembled, before the pin
    SparseMatrix pcd_mass;                 // production operators, gauge row
    SparseMatrix pcd_laplacian;            // already imposed by the assembler
    SparseMatrix pcd_convection_diffusion; //
    std::vector<double> velocity_old;  // component-major: [d*nc + c]
    std::size_t reference_cell = 0;
    bool has_fixed_pressure = false;
};

struct ActionStats {
    std::size_t applications = 0;  // Schur actions invoked by the outer Krylov
    std::size_t nested_iterations = 0;  // nested Krylov iterations (0 for rc_diag)
    double action_wall_ms = 0.0;    // total wall time inside the action
};

struct VariantSpec {
    const char* action;
    KrylovModel nested_krylov;
    PreconditionerModel nested_preconditioner;
    int nested_restart;
    std::size_t nested_budget;
    double nested_tolerance;
    bool use_kp_operator;   // solve Kp instead of Fp (documented provisional alternative)
    bool rc_scale_output;  // scale the nested solve output by diag(RC)
    bool diagonal_action;  // no nested solve: action = -diag(RC)^-1
    bool full_classic_action;  // classic literature composition -Kp^-1*Fp*Mp^-1
};

struct VariantResult {
    std::string action;
    int status = -1;
    bool converged = false;
    std::size_t coupled_iterations = 0;
    double true_rel_residual = std::numeric_limits<double>::infinity();
    double solve_wall_ms = 0.0;
    ActionStats action_stats;
    double nested_setup_ms = 0.0;
    double schur_setup_ms = 0.0;
    Vector x;
};

} // namespace

int main()
{
    constexpr std::size_t nx = 32;
    constexpr std::size_t ny = 32;
    constexpr double reynolds = 100.0;
    const std::size_t nc = nx * ny;

    Mesh mesh = make_cavity_mesh(nx, ny);
    Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
    Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
    U.fill(0.0); p.fill(0.0);

    VelocityBoundaryConditions ubc;
    ubc["bottom"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["left"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["right"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["front"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["back"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["top"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{1,0,0}};

    ScalarBoundaryConditions pbc;
    for(const char* name:{"bottom","top","left","right","front","back"})
        pbc[name]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

    // Mirror of the test_ghia_cavity.cpp --pcd controls, except the outer
    // nonlinear budget: five iterations are enough to capture the early
    // high-cost stages, and the benchmark never judges the physical result.
    IncompressibleSolverControls controls;
    controls.algorithm = PressureVelocityAlgorithm::COUPLED;
    controls.coupling.schur_model = CoupledSchurModel::PCD;
    controls.density = 1.0;
    controls.kinematic_viscosity = 1.0 / reynolds;
    controls.linear_max_iterations = 10000;
    controls.linear_tolerance = 1e-8;
    controls.pressure_reference_cell = (ny / 2) * nx + (nx / 2);
    controls.pressure_reference_value = 0.0;
    controls.coupling.alpha_u = 0.7;
    controls.coupling.alpha_p = 0.3;
    controls.coupling.coupled_max_iterations = 5000;
    controls.coupling.coupled_linear_tolerance = 1e-10;
    controls.coupled_linear_solver.krylov = KrylovModel::FGMRES;
    controls.coupled_linear_solver.preconditioner = PreconditionerModel::PCD;
    const std::size_t coupled_unknowns = 4 * nx * ny;
    controls.coupled_linear_solver.gmres_restart = static_cast<int>(
        coupled_unknowns <= 4096
            ? std::min<std::size_t>(coupled_unknowns,
                                    std::max<std::size_t>(1024,
                                                          controls.coupling.coupled_max_iterations))
            : 512);
    controls.pcd_convection_diffusion_max_iterations = 200;
    controls.pcd_convection_diffusion_linear_solver.gmres_restart =
        static_cast<int>(std::min<std::size_t>(
            coupled_unknowns / 4, std::max<std::size_t>(512, std::size_t{200})));
    controls.use_bounded_convection = true;
    controls.convection_scheme = ConvectionScheme::UPWIND;
    controls.convergence.max_iterations = 5;
    controls.convergence.relative_tolerance = 1e-8;
    controls.convergence.continuity_tolerance = 1e-8;

    std::vector<CapturedStage> stages;
    controls.coupled_assembly_callback =
        [&](const CoupledAssemblySnapshot& snapshot) {
            if (snapshot.matrix == nullptr || snapshot.rhs == nullptr ||
                snapshot.initial_guess == nullptr ||
                snapshot.pressure_schur_diagonal == nullptr ||
                snapshot.pressure_schur_diagonal_raw == nullptr ||
                snapshot.pcd_pressure_mass == nullptr ||
                snapshot.pcd_pressure_laplacian == nullptr ||
                snapshot.pcd_pressure_convection_diffusion == nullptr ||
                snapshot.velocity_old == nullptr)
                throw std::runtime_error(
                    "pcd_schur_action_benchmark: incomplete capture snapshot "
                    "(non-PCD coupled plan cannot be benchmarked)");
            CapturedStage stage;
            stage.matrix = *snapshot.matrix;
            stage.rhs = *snapshot.rhs;
            stage.initial_guess = *snapshot.initial_guess;
            stage.rc_diagonal = *snapshot.pressure_schur_diagonal;
            stage.rc_diagonal_raw = *snapshot.pressure_schur_diagonal_raw;
            stage.pcd_mass = *snapshot.pcd_pressure_mass;
            stage.pcd_laplacian = *snapshot.pcd_pressure_laplacian;
            stage.pcd_convection_diffusion =
                *snapshot.pcd_pressure_convection_diffusion;
            stage.velocity_old.resize(3 * nc);
            for (std::size_t d = 0; d < 3; ++d)
                for (std::size_t c = 0; c < nc; ++c)
                    stage.velocity_old[d * nc + c] =
                        snapshot.velocity_old->component_data(d)[c];
            stage.reference_cell = snapshot.reference_cell;
            stage.has_fixed_pressure = snapshot.has_fixed_pressure;
            stages.push_back(std::move(stage));
        };

    std::cerr << "PCD_ACTION_CAPTURE start grid=" << nx << "x" << ny
              << " Re=" << reynolds << " nonlinear_iterations=5\n";
    // Non-convergence after five iterations is expected and irrelevant: only
    // the captured solved systems are replayed below.
    (void)solve_steady_incompressible(mesh, U, p, ubc, pbc, controls);

    if (stages.size() < 5)
        throw std::runtime_error(
            "pcd_schur_action_benchmark: expected at least 5 captured coupled "
            "assemblies, got " + std::to_string(stages.size()));

    const FvGeometry geometry = build_fv_geometry(mesh);

    const std::array<VariantSpec, 8> variants = {{
        {"fp_ilu0", KrylovModel::FGMRES, PreconditionerModel::ILU0, 512, 200, 1e-8,
         false, false, false, false},
        {"fp_amg", KrylovModel::FGMRES, PreconditionerModel::NativeAMG, 512, 200, 1e-8,
         false, false, false, false},
        {"kp_amg", KrylovModel::GMRES, PreconditionerModel::NativeAMG, 40, 100, 1e-8,
         true, false, false, false},
        {"fp_rc_scaled", KrylovModel::FGMRES, PreconditionerModel::ILU0, 512, 200, 1e-8,
         false, true, false, false},
        {"rc_diag", KrylovModel::GMRES, PreconditionerModel::ILU0, 40, 100, 1e-8,
         false, false, true, false},
        {"pcd_full", KrylovModel::GMRES, PreconditionerModel::NativeAMG, 40, 100, 1e-8,
         false, false, false, true},
        {"fp_ilu0_loose", KrylovModel::FGMRES, PreconditionerModel::ILU0, 512, 200, 1e-4,
         false, false, false, false},
        {"fp_amg_loose", KrylovModel::FGMRES, PreconditionerModel::NativeAMG, 512, 200, 1e-4,
         false, false, false, false},
    }};

    const std::vector<std::size_t> benchmark_stages = {1, 4};  // assemblies #2 and #5
    std::size_t failures = 0;
    std::map<std::string, std::vector<VariantResult>> results_by_variant;

    for (const std::size_t stage_index : benchmark_stages) {
        const CapturedStage& stage = stages[stage_index];
        const std::size_t n = stage.matrix.n_rows();
        const std::size_t nv = 3 * nc;
        const std::size_t gauge_row = nv + stage.reference_cell;

        // ------------------------------------------------------------
        // Validity guard: gauge coherence of the captured system. With no
        // fixed pressure the reference row of the production coupled
        // matrix is an identity row with the reference RHS, and the
        // production initial guess carries the reference value on the
        // gauge row. A violation means the capture is NOT the solved
        // system and every measurement below is meaningless.
        // ------------------------------------------------------------
        if (!stage.has_fixed_pressure) {
            const auto row_begin = stage.matrix.row_offsets_data()[gauge_row];
            const auto row_end = stage.matrix.row_offsets_data()[gauge_row + 1];
            bool gauge_ok = (row_end - row_begin == 1) &&
                stage.matrix.columns_data()[row_begin] == gauge_row &&
                stage.matrix.values_data()[row_begin] == 1.0 &&
                stage.rhs(gauge_row) == controls.pressure_reference_value &&
                stage.initial_guess(gauge_row) == controls.pressure_reference_value;
            if (!gauge_ok) {
                std::cerr << "PCD_ACTION_GAUGE stage=" << (stage_index + 1)
                          << " gauge_row=" << gauge_row
                          << " row_nnz=" << (row_end - row_begin)
                          << " rhs=" << stage.rhs(gauge_row)
                          << " x0=" << stage.initial_guess(gauge_row)
                          << " - captured system NOT gauge-coherent\n";
                ++failures;
            } else {
                std::cerr << "PCD_ACTION_GAUGE stage=" << (stage_index + 1)
                          << " gauge_row=" << gauge_row
                          << " identity_row=1 rhs=reference x0=reference"
                          << " - gauge coherence OK\n";
            }
        }

        // ------------------------------------------------------------
        // Validity guard: the captured production PCD operators must be
        // reproducible by the documented assembly path from the captured
        // linearization velocity. The replay itself uses the CAPTURED
        // production operators (no reconstruction risk); this check
        // verifies that the assembly path and the captured production
        // operators describe the same mathematics. The rebuild wall time
        // is also the production per-assembly operator cost measurement.
        // ------------------------------------------------------------
        Field<double, Location::CELL> U_stage(mesh.n_cells(), "U_stage", "m/s", 3);
        for (std::size_t d = 0; d < 3; ++d)
            for (std::size_t c = 0; c < nc; ++c)
                U_stage.component_data(d)[c] = stage.velocity_old[d * nc + c];

        const auto rebuild_start = std::chrono::steady_clock::now();
        const auto rebuilt_ops = assemble_pcd_pressure_operators(
            mesh, geometry, U_stage, controls.density,
            controls.kinematic_viscosity, ubc, pbc, stage.reference_cell,
            ConvectionScheme::UPWIND);
        const double rebuild_ms = elapsed_ms(rebuild_start);

        const double mass_diff = sparse_max_abs_diff(
            stage.pcd_mass, rebuilt_ops.mass);
        const double kp_diff = sparse_max_abs_diff(
            stage.pcd_laplacian, rebuilt_ops.laplacian);
        const double fp_diff = sparse_max_abs_diff(
            stage.pcd_convection_diffusion, rebuilt_ops.convection_diffusion);
        const double ops_diff = std::max(mass_diff, std::max(kp_diff, fp_diff));
        std::cerr << "PCD_ACTION_OPERATORS stage=" << (stage_index + 1)
                  << std::scientific << std::setprecision(6)
                  << " mass_diff=" << mass_diff
                  << " kp_diff=" << kp_diff
                  << " fp_diff=" << fp_diff
                  << " rebuild_ms=" << rebuild_ms
                  << " (gate 1e-12) - " << (ops_diff <= 1e-12 ? "OK" : "MISMATCH")
                  << "\n";
        if (!(ops_diff <= 1e-12))
            ++failures;

        // Gauge-pinned Rhie-Chow diagonal for the rc variants (production
        // form: reference cell pinned to 1.0, gauge component preserved).
        const std::vector<double>& rc = stage.rc_diagonal;
        for (const VariantSpec& spec : variants) {
            if (spec.diagonal_action) {
                for (std::size_t c = 0; c < nc; ++c)
                    if (!(rc[c] > 0.0) || !std::isfinite(rc[c]))
                        throw std::runtime_error(
                            "pcd_schur_action_benchmark: invalid Rhie-Chow "
                            "diagonal entry at cell " + std::to_string(c));
            }
        }

        std::vector<VariantResult> stage_results;
        for (const VariantSpec& spec : variants) {
            const auto nested_setup_start = std::chrono::steady_clock::now();
            ActionStats action_stats;

            // Nested single-inverse solver for the PCD action. The production
            // policy (nested_schur_pressure_solve_usable) is reused verbatim:
            // a MAX_ITER nested solve with a finite result remains usable.
            PcdSchurApproximation::PressureSolve inverse_solve;
            if (spec.diagonal_action) {
                std::vector<double> rc_inv(nc);
                for (std::size_t c = 0; c < nc; ++c)
                    rc_inv[c] = 1.0 / rc[c];
                inverse_solve = [rc_inv, &action_stats](
                    const Vector& rhs, Vector& x) {
                    const auto apply_start = std::chrono::steady_clock::now();
                    if (x.size() != rhs.size())
                        x = Vector(rhs.size(), 0.0);
                    for (std::size_t i = 0; i < rhs.size(); ++i)
                        x(i) = rc_inv[i] * rhs(i);
                    ++action_stats.applications;
                    action_stats.action_wall_ms += elapsed_ms(apply_start);
                    return true;
                };
            } else if (spec.full_classic_action) {
                // Classic literature PCD action: S^-1 rhs ~= -Kp^-1 Fp
                // Mp^-1 rhs (Kay, Loghin & Wathen; Elman, Silvester &
                // Wathen 2001). apply() in pcd_schur.h provides the
                // overall negation, so the slot returns +Kp^-1 Fp Mp^-1
                // rhs - the same overall sign structure as the
                // production slot (+Fp^-1 rhs for the production action
                // -Fp^-1). One nested inverse (Kp) plus one forward
                // product (Fp) plus one diagonal solve (Mp^-1): the same
                // single-inner-inverse structure as production, NOT the
                // double-inverse compositions documented as non-convergent
                // in pcd_schur.h. The ordering, signs and gauge behaviour
                // are NOT asserted: this is a measurement variant, the
                // benchmark measures whether the classic composition
                // converges the real captured coupled system. Mp is
                // diagonal by construction (assemble_pcd_pressure_operators
                // builds diag(cell_volumes)); that structure is verified
                // fail-closed here instead of assumed.
                std::vector<double> mass_diag(nc);
                for (std::size_t r = 0; r < nc; ++r) {
                    bool diag_found = false;
                    for (std::size_t k = stage.pcd_mass.row_offsets_data()[r];
                         k < stage.pcd_mass.row_offsets_data()[r + 1]; ++k) {
                        const std::size_t col = stage.pcd_mass.columns_data()[k];
                        const double val = stage.pcd_mass.values_data()[k];
                        if (col != r)
                            throw std::runtime_error(
                                "pcd_schur_action_benchmark: the pressure "
                                "mass operator is not diagonal - the classic "
                                "PCD composition requires a diagonal Mp "
                                "(off-diagonal entry at cell " +
                                std::to_string(r) + ")");
                        if (!(val > 0.0) || !std::isfinite(val))
                            throw std::runtime_error(
                                "pcd_schur_action_benchmark: invalid pressure "
                                "mass diagonal at cell " + std::to_string(r));
                        mass_diag[r] = val;
                        diag_found = true;
                    }
                    if (!diag_found)
                        throw std::runtime_error(
                            "pcd_schur_action_benchmark: missing pressure "
                            "mass diagonal at cell " + std::to_string(r));
                }
                auto classic_preconditioner =
                    make_scalar_preconditioner(spec.nested_preconditioner);
                if (!classic_preconditioner ||
                    !classic_preconditioner->setup(stage.pcd_laplacian))
                    throw std::runtime_error(
                        "pcd_schur_action_benchmark: Kp preconditioner "
                        "setup failed for action " + std::string(spec.action));
                const auto nested_preconditioner =
                    std::shared_ptr<Preconditioner>(std::move(classic_preconditioner));
                const int nested_restart = std::min<int>(
                    spec.nested_restart,
                    static_cast<int>(stage.pcd_laplacian.n_rows()));
                const KrylovModel krylov = spec.nested_krylov;
                const std::size_t budget = spec.nested_budget;
                const double tolerance = spec.nested_tolerance;
                inverse_solve = [&stage, nested_preconditioner, nested_restart,
                                 budget, tolerance, krylov,
                                 mass_diag = std::move(mass_diag),
                                 &action_stats](const Vector& rhs, Vector& x) {
                    const auto apply_start = std::chrono::steady_clock::now();
                    const std::size_t np = rhs.size();
                    // z = Mp^-1 rhs (exact: Mp is diagonal, verified above)
                    Vector z(np, 0.0);
                    for (std::size_t i = 0; i < np; ++i)
                        z(i) = rhs(i) / mass_diag[i];
                    // y = Fp z (forward product with the CAPTURED production
                    // convection-diffusion operator)
                    Vector y(np, 0.0);
                    for (std::size_t row = 0; row < np; ++row) {
                        double ax = 0.0;
                        const auto rb = stage.pcd_convection_diffusion
                            .row_offsets_data()[row];
                        const auto re = stage.pcd_convection_diffusion
                            .row_offsets_data()[row + 1];
                        for (std::size_t k = rb; k < re; ++k)
                            ax += stage.pcd_convection_diffusion.values_data()[k]
                                * z(stage.pcd_convection_diffusion.columns_data()[k]);
                        y(row) = ax;
                    }
                    // x = Kp^-1 y (nested solve with the CAPTURED production
                    // laplacian; production nested-usable policy applies)
                    SolverResult result;
                    if (krylov == KrylovModel::GMRES)
                        result = solve_gmres(stage.pcd_laplacian, y, x,
                                             nested_restart, budget, tolerance,
                                             nested_preconditioner.get(), {},
                                             {}, false);
                    else
                        result = solve_fgmres(stage.pcd_laplacian, y, x,
                                              nested_restart, budget, tolerance,
                                              nested_preconditioner.get(), {},
                                              {}, false);
                    ++action_stats.applications;
                    action_stats.nested_iterations += result.iterations;
                    action_stats.action_wall_ms += elapsed_ms(apply_start);
                    return nested_schur_pressure_solve_usable(result, x, np);
                };
            } else {
                // The variant's nested operator: Fp (production baseline and
                // scaling variants) or Kp (documented provisional
                // alternative). Both come from the CAPTURED production
                // operators.
                const SparseMatrix& nested_matrix = spec.use_kp_operator
                    ? stage.pcd_laplacian
                    : stage.pcd_convection_diffusion;
                auto preconditioner =
                    make_scalar_preconditioner(spec.nested_preconditioner);
                if (!preconditioner || !preconditioner->setup(nested_matrix))
                    throw std::runtime_error(
                        "pcd_schur_action_benchmark: nested preconditioner "
                        "setup failed for action " + std::string(spec.action));
                const auto nested_preconditioner =
                    std::shared_ptr<Preconditioner>(std::move(preconditioner));
                const int nested_restart = std::min<int>(
                    spec.nested_restart, static_cast<int>(nested_matrix.n_rows()));
                const KrylovModel krylov = spec.nested_krylov;
                const std::size_t budget = spec.nested_budget;
                const double tolerance = spec.nested_tolerance;
                const bool scale_output = spec.rc_scale_output;
                inverse_solve = [nested_preconditioner, &nested_matrix, nested_restart,
                                 budget, tolerance, krylov, scale_output, &rc,
                                 &action_stats](const Vector& rhs, Vector& x) {
                    const auto solve_start = std::chrono::steady_clock::now();
                    SolverResult result;
                    if (krylov == KrylovModel::GMRES)
                        result = solve_gmres(nested_matrix, rhs, x, nested_restart,
                                             budget, tolerance,
                                             nested_preconditioner.get(), {}, {},
                                             false);
                    else
                        result = solve_fgmres(nested_matrix, rhs, x, nested_restart,
                                              budget, tolerance,
                                              nested_preconditioner.get(), {}, {},
                                              false);
                    if (scale_output)
                        for (std::size_t i = 0; i < x.size(); ++i)
                            x(i) *= rc[i];
                    ++action_stats.applications;
                    action_stats.nested_iterations += result.iterations;
                    action_stats.action_wall_ms += elapsed_ms(solve_start);
                    return nested_schur_pressure_solve_usable(
                        result, x, nested_matrix.n_rows());
                };
            }
            const double nested_setup_ms = elapsed_ms(nested_setup_start);

            // Callback-contract audit: the single-inverse PCD action must
            // never invoke the laplacian_solve slot (pcd_schur.h documents
            // apply() as -convection_diffusion_solve_ then a sign flip).
            // Passing a throwing callback there turns any unexpected
            // invocation into a loud failure instead of a silent
            // measurement of the wrong operator. The setup() contract
            // only requires the slot to be non-null, which a throwing
            // std::function satisfies.
            const PcdSchurApproximation::PressureSolve must_not_be_called =
                [](const Vector&, Vector&) {
                    throw std::runtime_error(
                        "pcd_schur_action_benchmark: the single-inverse PCD "
                        "action invoked the laplacian_solve slot - callback "
                        "contract violated");
                    return false;  // unreachable: gives the audit callback
                                   // the PressureSolve bool contract
                };

            // PCD action built from the CAPTURED production operators,
            // with the production gauge policy.
            auto pcd = std::make_unique<PcdSchurApproximation>(
                stage.pcd_mass, stage.pcd_laplacian,
                stage.pcd_convection_diffusion,
                must_not_be_called, inverse_solve, std::nullopt,
                stage.has_fixed_pressure
                    ? std::optional<std::size_t>()
                    : std::optional<std::size_t>(stage.reference_cell));

            CoupledBlockSchurOptions options;
            options.factorization = CoupledSchurFactorization::Full;
            options.velocity_approximation = CoupledSchurVelocityApproximation::Block;
            options.schur_approximation = CoupledSchurApproximationModel::PCD;

            const auto schur_setup_start = std::chrono::steady_clock::now();
            auto schur = std::make_unique<CoupledBlockSchurAMGPreconditioner>(
                nc, options,
                stage.has_fixed_pressure
                    ? std::numeric_limits<std::size_t>::max()
                    : stage.reference_cell);
            schur->set_pcd_schur(std::move(pcd));
            if (!schur->setup(stage.matrix))
                throw std::runtime_error(
                    "pcd_schur_action_benchmark: coupled Schur setup failed "
                    "for action " + std::string(spec.action));
            const double schur_setup_ms = elapsed_ms(schur_setup_start);

            // Identical outer coupled solve for every variant: same captured
            // matrix, same captured RHS, same captured production initial
            // guess, same fixed restart.
            const int coupled_restart = static_cast<int>(
                std::min<std::size_t>(4096, n));
            KrylovControls coupled_controls;
            coupled_controls.restart_min = coupled_restart;
            coupled_controls.restart_max = coupled_restart;
            coupled_controls.adaptive_restart = false;

            Vector x = stage.initial_guess;
            const auto solve_start = std::chrono::steady_clock::now();
            const SolverResult coupled = solve_fgmres(
                stage.matrix, stage.rhs, x, coupled_restart,
                controls.coupling.coupled_max_iterations,
                controls.coupling.coupled_linear_tolerance, schur.get(),
                coupled_controls);
            const double solve_wall_ms = elapsed_ms(solve_start);

            VariantResult vr;
            vr.action = spec.action;
            vr.status = static_cast<int>(coupled.status);
            vr.converged = coupled.status == SolverStatus::CONVERGED;
            vr.coupled_iterations = coupled.iterations;
            vr.true_rel_residual = true_relative_residual(stage.matrix, stage.rhs, x);
            vr.solve_wall_ms = solve_wall_ms;
            vr.action_stats = action_stats;
            vr.nested_setup_ms = nested_setup_ms;
            vr.schur_setup_ms = schur_setup_ms;
            vr.x = x;

            std::cout << std::scientific << std::setprecision(6)
                      << "PCD_ACTION_BENCHMARK"
                      << " stage=" << (stage_index + 1)
                      << " n=" << n
                      << " nnz=" << stage.matrix.nnz()
                      << " action=" << vr.action
                      << " status=" << vr.status
                      << " converged=" << (vr.converged ? 1 : 0)
                      << " coupled_iterations=" << vr.coupled_iterations
                      << " true_rel_residual=" << vr.true_rel_residual
                      << " solve_wall_ms=" << vr.solve_wall_ms
                      << " schur_applications=" << vr.action_stats.applications
                      << " nested_iterations=" << vr.action_stats.nested_iterations
                      << " action_wall_ms=" << vr.action_stats.action_wall_ms
                      << " non_nested_wall_ms="
                      << (vr.solve_wall_ms - vr.action_stats.action_wall_ms)
                      << " nested_setup_ms=" << vr.nested_setup_ms
                      << " schur_setup_ms=" << vr.schur_setup_ms
                      << "\n";
            std::cout.flush();

            // Baseline fidelity guard: the production -Fp^-1 + ILU0 action
            // must solve the captured system - from the captured production
            // initial guess - to the production gate. If it does not, the
            // capture diverged from production and every other number is
            // meaningless. (A wrong capture can still be solved to a tiny
            // residual; this guard is therefore paired with the gauge and
            // operator-coherence guards above, which constrain WHAT is
            // solved, while this guard constrains the replay's ability to
            // reproduce the production convergence on it.)
            if (vr.action == "fp_ilu0") {
                if (!vr.converged || !(vr.true_rel_residual <= 1e-10)) {
                    std::cerr << "PCD_ACTION_BASELINE stage=" << (stage_index + 1)
                              << " converged=" << (vr.converged ? 1 : 0)
                              << " true_rel_residual=" << vr.true_rel_residual
                              << " (gate 1e-10) - capture fidelity FAILED\n";
                    ++failures;
                } else {
                    std::cerr << "PCD_ACTION_BASELINE stage=" << (stage_index + 1)
                              << " converged=1"
                              << " true_rel_residual=" << vr.true_rel_residual
                              << " (gate 1e-10) - capture fidelity OK\n";
                }
            }

            // Retain the result for the stage-level consistency check and
            // the final performance aggregate.
            stage_results.push_back(vr);
            results_by_variant[vr.action].push_back(vr);
        }

        // Linear-solution consistency guard: every converged variant must
        // agree with every other converged variant on the solution of the
        // SAME captured system (identical matrix, RHS, production initial
        // guess and tolerances). A variant that "converges" to a different
        // solution is a broken action, not merely a slower one.
        for (std::size_t i = 0; i < stage_results.size(); ++i) {
            if (!stage_results[i].converged) continue;
            for (std::size_t j = i + 1; j < stage_results.size(); ++j) {
                if (!stage_results[j].converged) continue;
                double diff_sq = 0.0;
                double ref_sq = 0.0;
                for (std::size_t k = 0; k < stage_results[i].x.size(); ++k) {
                    const double d = stage_results[i].x(k) - stage_results[j].x(k);
                    diff_sq += d * d;
                    ref_sq += stage_results[i].x(k) * stage_results[i].x(k);
                }
                const double rel = std::sqrt(diff_sq) /
                    (ref_sq > 0.0 ? std::sqrt(ref_sq) : 1.0);
                if (!(rel <= 1e-6)) {
                    std::cerr << "PCD_ACTION_CONSISTENCY stage=" << (stage_index + 1)
                              << " a=" << stage_results[i].action
                              << " b=" << stage_results[j].action
                              << " rel_diff=" << rel
                              << " (gate 1e-6) - solutions DISAGREE\n";
                    ++failures;
                }
            }
        }
    }

    // ----------------------------------------------------------------
    // PERFORMANCE REPORT (never gating): per-variant aggregates over the
    // benchmarked stages and baseline-relative comparison. A variant that
    // did not converge on some stage is reported as such, not failed.
    // ----------------------------------------------------------------
    double baseline_total_ms = -1.0;
    const auto baseline_it = results_by_variant.find("fp_ilu0");
    if (baseline_it != results_by_variant.end()) {
        baseline_total_ms = 0.0;
        for (const VariantResult& vr : baseline_it->second)
            if (vr.converged)
                baseline_total_ms += vr.solve_wall_ms;
    }
    for (const auto& [name, results] : results_by_variant) {
        std::size_t converged_stages = 0;
        double total_ms = 0.0;
        std::size_t coupled_iterations = 0;
        std::size_t schur_applications = 0;
        std::size_t nested_iterations = 0;
        for (const VariantResult& vr : results) {
            if (vr.converged) {
                ++converged_stages;
                total_ms += vr.solve_wall_ms;
            }
            coupled_iterations += vr.coupled_iterations;
            schur_applications += vr.action_stats.applications;
            nested_iterations += vr.action_stats.nested_iterations;
        }
        std::cout << std::scientific << std::setprecision(6)
                  << "PCD_ACTION_SUMMARY"
                  << " action=" << name
                  << " stages_converged=" << converged_stages
                  << " of " << results.size()
                  << " total_solve_wall_ms=" << total_ms
                  << " coupled_iterations=" << coupled_iterations
                  << " schur_applications=" << schur_applications
                  << " nested_iterations=" << nested_iterations;
        if (name != "fp_ilu0" && baseline_total_ms > 0.0 && total_ms > 0.0)
            std::cout << " speedup_vs_baseline=" << baseline_total_ms / total_ms;
        std::cout << "\n";
    }
    std::cout.flush();

    // Final verdict, printed only after every guard has been evaluated.
    // VALIDITY is what this test gates: capture fidelity (gauge coherence,
    // operator coherence, baseline reproduction) and solution consistency.
    // PERFORMANCE is reported above and never gates the verdict.
    if (failures == 0) {
        std::cout << "PCD_ACTION_VALIDITY OK: capture is the solved production "
                     "system; baseline reproduces the production gate\n";
        std::cout << "PASS: pcd_schur_action_benchmark\n";
        return 0;
    }
    std::cout << "PCD_ACTION_VALIDITY FAILED: " << failures
              << " guard failure(s) - measurements must NOT be interpreted\n";
    std::cout << "FAIL: pcd_schur_action_benchmark\n";
    return 1;
}
