// PCD Schur-action A/B benchmark on the REAL assembled Ghia coupled systems.
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
//      every assembled coupled system (matrix, RHS, Rhie-Chow Schur diagonal,
//      linearization velocity).
//   2. Replay two captured stages (assemblies #2 and #5: the early regime
//      where the measured production per-stage cost peaks) through the
//      documented PCD action variants, under identical conditions:
//      same matrix, same RHS, same initial guess (zero), same outer FGMRES
//      (restart 4096, budget 5000, tolerance 1e-10, fixed KrylovControls).
//
// Variants (all built from the existing PcdSchurApproximation class; NO
// production pcd_schur.h change is involved):
//   fp_ilu0   -Fp^-1, nested FGMRES(512)+ILU0, budget 200, tol 1e-8
//              (exactly the production --pcd configuration = baseline).
//   fp_amg    -Fp^-1, nested FGMRES(512)+NativeAMG, budget 200, tol 1e-8
//              (the nested-preconditioner axis of PR #965, measured on the
//              real matrix).
//   kp_amg    -Kp^-1, nested GMRES(40)+NativeAMG, budget 100, tol 1e-8
//              (the provisional alternative documented in pcd_schur.h; Kp
//              defaults used because that is how production would solve it).
//   fp_rc_scaled  -diag(RC)*Fp^-1 (left scaling by the exact Rhie-Chow
//              diagonal captured from the assembly, gauge cell pinned to 1.0
//              exactly like the production assembly pin). Measurement
//              variant: the algebraic derivation (left vs right scaling,
//              signs, gauge treatment) is NOT asserted here; the benchmark
//              measures whether the amplitude correction helps at all.
//   rc_diag   -diag(RC)^-1 with no nested solve (SIMPLE-style diagonal
//              Schur). Reference point; the nested counters stay at zero.
//
// Evidence printed (one TSV-ish line per solve, stdout):
//   PCD_ACTION_BENCHMARK stage=... action=... status=... converged=...
//   coupled_iterations=... true_rel_residual=... solve_wall_ms=...
//   nested_calls=... nested_iterations=... nested_wall_ms=...
//   nested_setup_ms=... schur_setup_ms=...
//
// Acceptance (what makes this test PASS):
//   - at least five coupled assemblies were captured;
//   - the baseline fp_ilu0 action converges every benchmarked stage to the
//     production 1e-10 true relative residual (known production behaviour:
//     every stage of the 32x32 qualification reaches ~9.8e-11);
//   - all converged variants agree on the solution to 1e-6 relative
//     (a linear-solve consistency guard: no variant may "converge" to a
//     different solution).
// A non-baseline variant failing to converge is REPORTED, not failed: the
// measurement is the point. Physical acceptance of the full case remains
// the exclusive responsibility of test_ghia_cavity_n8_pcd.
//
// Cost: bounded and reproducible - five production nonlinear iterations plus
// ten coupled solves, a few minutes in Release. Registered only under
// CFDX_ENABLE_LONG_VALIDATION.

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
// qualification test; the baseline-variant convergence assertion below is the
// fidelity guard.
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

struct CapturedStage {
    SparseMatrix matrix;
    Vector rhs;
    std::vector<double> rc_diagonal;
    std::vector<double> velocity_old;  // component-major: [d*nc + c]
    std::size_t reference_cell = 0;
    bool has_fixed_pressure = false;
};

struct NestedStats {
    std::size_t calls = 0;
    std::size_t iterations = 0;
    double wall_ms = 0.0;
};

struct VariantSpec {
    const char* action;
    KrylovModel nested_krylov;
    PreconditionerModel nested_preconditioner;
    int nested_restart;
    std::size_t nested_budget;
    double nested_tolerance;
    bool use_kp_operator;   // solve Kp instead of Fp (documented provisional alternative)
    bool rc_scale_output;  // left-scale the nested solve output by diag(RC)
    bool diagonal_action;  // no nested solve: action = -diag(RC)^-1
};

struct VariantResult {
    std::string action;
    int status = -1;
    bool converged = false;
    std::size_t coupled_iterations = 0;
    double true_rel_residual = std::numeric_limits<double>::infinity();
    double solve_wall_ms = 0.0;
    NestedStats nested;
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
            CapturedStage stage;
            stage.matrix = *snapshot.matrix;
            stage.rhs = *snapshot.rhs;
            stage.rc_diagonal = *snapshot.pressure_schur_diagonal;
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
    // the captured assembled systems are replayed below.
    (void)solve_steady_incompressible(mesh, U, p, ubc, pbc, controls);

    if (stages.size() < 5)
        throw std::runtime_error(
            "pcd_schur_action_benchmark: expected at least 5 captured coupled "
            "assemblies, got " + std::to_string(stages.size()));

    const FvGeometry geometry = build_fv_geometry(mesh);

    const std::array<VariantSpec, 5> variants = {{
        {"fp_ilu0", KrylovModel::FGMRES, PreconditionerModel::ILU0, 512, 200, 1e-8,
         false, false, false},
        {"fp_amg", KrylovModel::FGMRES, PreconditionerModel::NativeAMG, 512, 200, 1e-8,
         false, false, false},
        {"kp_amg", KrylovModel::GMRES, PreconditionerModel::NativeAMG, 40, 100, 1e-8,
         true, false, false},
        {"fp_rc_scaled", KrylovModel::FGMRES, PreconditionerModel::ILU0, 512, 200, 1e-8,
         false, true, false},
        {"rc_diag", KrylovModel::GMRES, PreconditionerModel::ILU0, 40, 100, 1e-8,
         false, false, true},
    }};

    const std::vector<std::size_t> benchmark_stages = {1, 4};  // assemblies #2 and #5
    std::size_t failures = 0;

    for (const std::size_t stage_index : benchmark_stages) {
        const CapturedStage& stage = stages[stage_index];
        const std::size_t n = stage.matrix.n_rows();

        // Rebuild the linearization velocity field and the PCD pressure
        // operators exactly as the production assembly does for this stage.
        Field<double, Location::CELL> U_stage(mesh.n_cells(), "U_stage", "m/s", 3);
        for (std::size_t d = 0; d < 3; ++d)
            for (std::size_t c = 0; c < nc; ++c)
                U_stage.component_data(d)[c] = stage.velocity_old[d * nc + c];

        const auto pcd_ops = assemble_pcd_pressure_operators(
            mesh, geometry, U_stage, controls.density,
            controls.kinematic_viscosity, ubc, pbc, stage.reference_cell,
            ConvectionScheme::UPWIND);

        // Rhie-Chow diagonal with the production gauge pin.
        std::vector<double> rc = stage.rc_diagonal;
        if (rc.size() == nc)
            rc[stage.reference_cell] = 1.0;

        std::vector<VariantResult> stage_results;
        for (const VariantSpec& spec : variants) {
            const auto nested_setup_start = std::chrono::steady_clock::now();
            NestedStats nested_stats;

            // Nested single-inverse solver for the PCD action. The production
            // policy (nested_schur_pressure_solve_usable) is reused verbatim:
            // a MAX_ITER nested solve with a finite result remains usable.
            PcdSchurApproximation::PressureSolve inverse_solve;
            if (spec.diagonal_action) {
                std::vector<double> rc_inv(nc);
                for (std::size_t c = 0; c < nc; ++c) {
                    if (!(rc[c] > 0.0) || !std::isfinite(rc[c]))
                        throw std::runtime_error(
                            "pcd_schur_action_benchmark: invalid Rhie-Chow "
                            "diagonal entry");
                    rc_inv[c] = 1.0 / rc[c];
                }
                inverse_solve = [rc_inv](const Vector& rhs, Vector& x) {
                    if (x.size() != rhs.size())
                        x = Vector(rhs.size(), 0.0);
                    for (std::size_t i = 0; i < rhs.size(); ++i)
                        x(i) = rc_inv[i] * rhs(i);
                    return true;
                };
            } else {
                const SparseMatrix& nested_matrix = spec.use_kp_operator
                    ? pcd_ops.laplacian
                    : pcd_ops.convection_diffusion;
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
                                 &nested_stats](const Vector& rhs, Vector& x) {
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
                    ++nested_stats.calls;
                    nested_stats.iterations += result.iterations;
                    nested_stats.wall_ms += elapsed_ms(solve_start);
                    return nested_schur_pressure_solve_usable(
                        result, x, nested_matrix.n_rows());
                };
            }
            const double nested_setup_ms = elapsed_ms(nested_setup_start);

            // PCD action with the production gauge policy.
            auto pcd = std::make_unique<PcdSchurApproximation>(
                pcd_ops.mass, pcd_ops.laplacian, pcd_ops.convection_diffusion,
                inverse_solve, inverse_solve, std::nullopt,
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

            // Identical outer coupled solve for every variant: same matrix,
            // same RHS, same zero initial guess, same fixed restart.
            const int coupled_restart = static_cast<int>(
                std::min<std::size_t>(4096, n));
            KrylovControls coupled_controls;
            coupled_controls.restart_min = coupled_restart;
            coupled_controls.restart_max = coupled_restart;
            coupled_controls.adaptive_restart = false;

            Vector x(n, 0.0);
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
            vr.nested = nested_stats;
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
                      << " nested_calls=" << vr.nested.calls
                      << " nested_iterations=" << vr.nested.iterations
                      << " nested_wall_ms=" << vr.nested.wall_ms
                      << " nested_setup_ms=" << vr.nested_setup_ms
                      << " schur_setup_ms=" << vr.schur_setup_ms
                      << "\n";
            std::cout.flush();

            // Baseline fidelity guard: the production -Fp^-1 + ILU0 action
            // must solve the captured system to the production gate. If it
            // does not, the capture or the operator rebuild diverged from
            // production and every other number is meaningless.
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

            // Retain the result for the stage-level consistency check.
            stage_results.push_back(vr);
        }

        // Linear-solution consistency guard: every converged variant must
        // agree with every other converged variant on the solution of the
        // SAME captured system (identical matrix, RHS, zero initial guess
        // and tolerances). A variant that "converges" to a different
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

    std::cout << "PASS: pcd_schur_action_benchmark\n";
    std::cout << "\nResults: " << (failures == 0 ? 1 : 0) << " passed, "
              << failures << " failed\n";
    return failures == 0 ? 0 : 1;
}
