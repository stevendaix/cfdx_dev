#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "qualification_boundary_helpers.h"
#include "naca_momentum_trace.h"
#include "naca_relaxation_controls.h"
#include "cfdx/physics/wall_face_traction_reconstruction.h"
#include "cfdx/physics/wall_shear_quadratic_fit.h"
#include "naca_qualification_state_io.h"
#include <cstdlib>
#include <cstdio>
#include <cerrno>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

struct Forces { double cd,cl,cm,cdp,cdv,continuity; };

static Forces run(const std::string& path,std::size_t level,
                  bool growth_guard=false, double alpha_u=0.6,
                  double alpha_p=0.25, std::size_t max_iterations=8000)
{
    (void)level;  // the reported resolution is inferred from the actual mesh
    // Fail closed on every output that could be mistaken for physical
    // qualification after an early mesh/load/restart error. Validate *all*
    // destinations before deleting anything, and never truncate inputs.
    const char* output_names[] = {
        "CFDX_NACA_CONVERGED_FORCE_CSV",
        "CFDX_NACA_CONVERGED_FACE_FLUX_CSV",
        "CFDX_NACA_EXPORT_CONVERGED_STATE_CSV",
        "CFDX_NACA_FACE_PROFILES_CSV",
        "CFDX_NACA_FACE_TRACTION_RAW_CSV",
        "CFDX_NACA_PHYSICAL_HISTORY_CSV",
        "CFDX_NACA_MOMENTUM_BUDGET_CSV"};
    std::vector<std::string> outputs;
    for (const char* name : output_names) {
        const char* candidate = std::getenv(name);
        if (!candidate) continue;
        const std::string output(candidate);
        if (output.empty() || output == path ||
            (std::getenv("CFDX_NACA_IMPORT_INITIAL_STATE_CSV") &&
             output == std::getenv("CFDX_NACA_IMPORT_INITIAL_STATE_CSV")) ||
            (std::getenv("CFDX_NACA_IMPORT_FACE_FLUX_CSV") &&
             output == std::getenv("CFDX_NACA_IMPORT_FACE_FLUX_CSV")) ||
            std::find(outputs.begin(), outputs.end(), output) != outputs.end()) {
            throw std::invalid_argument(std::string("invalid or conflicting NACA output path: ") + name);
        }
        outputs.push_back(output);
    }
    for (const auto& output : outputs)
        if (std::remove(output.c_str()) != 0 && errno != ENOENT)
            throw std::runtime_error("cannot clear stale NACA diagnostic/certificate");
    Mesh mesh;
    if(!cfdx::io::read_mesh_hdf5(path,mesh))
        throw std::runtime_error("cannot read NACA0012 mesh: "+path);
    if(!mesh.topo_validate().ok)
        throw std::runtime_error("invalid NACA0012 topology: "+path);

    constexpr double rho=1.0, Uinf=1.0, Re=1000.0, chord=1.0;
    constexpr double nu=Uinf*chord/Re;
    Field<double,Location::CELL> U(mesh.n_cells(),"U","m/s",3);
    Field<double,Location::CELL> p(mesh.n_cells(),"p","Pa",1);
    U.fill(0.0); p.fill(0.0);
    const bool split_farfield = mesh.boundary().find("farfield_inlet") <
                                mesh.boundary().n_patches() &&
                                mesh.boundary().find("farfield_outlet") <
                                mesh.boundary().n_patches();
    VelocityBoundaryConditions ubc;
    if (split_farfield) {
        ubc["farfield_inlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,
                                  {Uinf,0,0}};
        ubc["farfield_outlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,
                                   {0,0,0}};
    } else {
        ubc["farfield"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{Uinf,0,0}};
    }
    ubc["airfoil"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["front"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    ubc["back"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    ScalarBoundaryConditions pbc;
    if (split_farfield) {
        pbc["farfield_inlet"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
        pbc["farfield_outlet"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    } else {
        pbc["farfield"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    }
    pbc["airfoil"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    pbc["front"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    pbc["back"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};

    IncompressibleSolverControls c;
    c.algorithm=PressureVelocityAlgorithm::SIMPLE;
    // Isolated A/B only. Both algorithms retain the exact same nonlinear
    // momentum, continuity and independent-face convergence criteria.
    if (const char* algorithm=std::getenv("CFDX_NACA_ALGORITHM")) {
        const std::string requested(algorithm);
        if (requested == "SIMPLEC")
            c.algorithm=PressureVelocityAlgorithm::SIMPLEC;
        else if (requested != "SIMPLE")
            throw std::invalid_argument("CFDX_NACA_ALGORITHM expects SIMPLE or SIMPLEC");
    }
    if (const char* pseudo=std::getenv("CFDX_NACA_PSEUDO_MOMENTUM")) {
        if (std::string(pseudo) != "1")
            throw std::invalid_argument("CFDX_NACA_PSEUDO_MOMENTUM must be 1");
        c.pseudo_momentum.enabled = true;
    }
    if (const char* cap=std::getenv("CFDX_NACA_PSEUDO_DIAGONAL_CAP")) {
        if (!c.pseudo_momentum.enabled)
            throw std::invalid_argument(
                "CFDX_NACA_PSEUDO_DIAGONAL_CAP requires pseudo momentum");
        std::size_t parsed = 0;
        const std::string raw(cap);
        const double fraction = std::stod(raw, &parsed);
        if (parsed != raw.size() || !std::isfinite(fraction) ||
            !(fraction > 0.0 && fraction <= 2.0))
            throw std::invalid_argument(
                "CFDX_NACA_PSEUDO_DIAGONAL_CAP expects a finite fraction in (0,2]");
        c.pseudo_momentum.maximum_diagonal_fraction = fraction;
    }
    c.density=rho; c.kinematic_viscosity=nu;
    c.linear_max_iterations=4000; c.linear_tolerance=1e-10;
    // Precision-only A/B knob. Tightening the algebraic solve can expose a
    // nonlinear residual floor without changing ANY physical stop criterion.
    // Forbid loosening the baseline linear tolerance in qualifications.
    if (const char* raw = std::getenv("CFDX_NACA_LINEAR_TOLERANCE")) {
        const std::string text(raw);
        std::size_t parsed = 0;
        double tolerance = 0.0;
        try { tolerance = std::stod(text, &parsed); }
        catch (const std::exception&) {
            throw std::invalid_argument("CFDX_NACA_LINEAR_TOLERANCE invalid number");
        }
        if (parsed != text.size() || !std::isfinite(tolerance) ||
            tolerance < 1e-14 || tolerance > 1e-10)
            throw std::invalid_argument(
                "CFDX_NACA_LINEAR_TOLERANCE requires [1e-14,1e-10]");
        c.linear_tolerance = tolerance;
    }
    c.convergence.max_iterations=max_iterations;
    c.convergence.relative_tolerance=1e-8;
    c.convergence.continuity_tolerance=1e-8;
    // Force/GCI output is forbidden until the independent face-transport
    // momentum residual reaches the same strict nonlinear tolerance.
    c.require_independent_face_residual_for_convergence = true;
    c.coupling.alpha_u=alpha_u; c.coupling.alpha_p=alpha_p;
    c.nonlinear_retry.reject_severe_physical_growth=growth_guard;
    cfdx::validation::apply_naca_relaxation_ab(c);
    // A/B only: the default NACA qualification is unchanged.
    if (const char* value = std::getenv("CFDX_NACA_SKEW_GRADIENT")) {
        if (std::string(value) != "1")
            throw std::invalid_argument("CFDX_NACA_SKEW_GRADIENT must be 1");
        c.pressure_skew_corrected_gradient = true;
    }
    // Experimental A/B only: the normal production and CI path is unchanged.
    if (const char* line_search = std::getenv("CFDX_NACA_PRESSURE_LINE_SEARCH")) {
        if (std::string(line_search) != "1")
            throw std::invalid_argument("CFDX_NACA_PRESSURE_LINE_SEARCH must be 1");
        c.pressure_continuity_line_search.enabled = true;
    }

    // A/B only: TVD is an existing bounded CFDX convection scheme. Keep the
    // unchanged UPWIND reference for all default NACA qualification gates.
    c.use_bounded_convection=true;
    c.convection_scheme=ConvectionScheme::UPWIND;
    if (const char* tvd = std::getenv("CFDX_NACA_CONVECTION_TVD")) {
        if (std::string(tvd) != "1")
            throw std::invalid_argument("CFDX_NACA_CONVECTION_TVD must be 1");
        c.convection_scheme=ConvectionScheme::TVD;
    }
    // Optional smooth TVD A/B variant. Never enable this on the default
    // NACA qualification path or jointly with the original MINMOD experiment.
    if (const char* vanleer = std::getenv("CFDX_NACA_CONVECTION_VANLEER")) {
        if (std::string(vanleer) != "1" || c.convection_scheme == ConvectionScheme::TVD)
            throw std::invalid_argument("NACA van Leer requires an exclusive value of 1");
        c.convection_scheme=ConvectionScheme::TVD_VANLEER;
    }
    // Alternative smooth limiter from the core numerics library, opt-in only.
    if (const char* albada = std::getenv("CFDX_NACA_CONVECTION_VANALBADA")) {
        if (std::string(albada) != "1" ||
            c.convection_scheme != ConvectionScheme::UPWIND)
            throw std::invalid_argument("NACA van Albada must be exclusive and equal 1");
        c.convection_scheme = ConvectionScheme::TVD_VANALBADA;
    }
    // Experimental multidimensional face-stencil MUSCL. Not the default TVD
    // qualification path. Mutual exclusion is intentional for A/B replay.
    if (const char* limited = std::getenv("CFDX_NACA_CONVECTION_CELL_LIMITED")) {
        if (std::string(limited) != "1" ||
            c.convection_scheme != ConvectionScheme::UPWIND)
            throw std::invalid_argument(
                "NACA cell-limited MUSCL must be exclusive and equal 1");
        c.convection_scheme = ConvectionScheme::MUSCL_CELL_LIMITED;
    }
    if (const char* continuation = std::getenv("CFDX_NACA_TVD_CONTINUATION")) {
        if (std::string(continuation) != "1" ||
            (c.convection_scheme != ConvectionScheme::TVD &&
             c.convection_scheme != ConvectionScheme::TVD_VANLEER &&
             c.convection_scheme != ConvectionScheme::TVD_VANALBADA &&
             c.convection_scheme != ConvectionScheme::MUSCL_CELL_LIMITED))
            throw std::invalid_argument(
                "CFDX_NACA_TVD_CONTINUATION requires TVD or van Leer and value 1");
        c.convection_continuation.enabled = true;
    }
    if (const char* affine = std::getenv("CFDX_NACA_AFFINE_FLUX_TRIALS")) {
        if (std::string(affine) != "1")
            throw std::invalid_argument("CFDX_NACA_AFFINE_FLUX_TRIALS must be 1");
        // The fast path requires the globally enabled line-search.
        if (!c.pressure_continuity_line_search.enabled)
            throw std::invalid_argument("affine flux trials require line-search");
        c.pressure_continuity_line_search.affine_flux_trials = true;
    }
    // Forensic only: select the exact difficult production pressure systems,
    // e.g. N256 iterations 20,21,240. Never alters a coefficient or tolerance.
    const char* replay_prefix = std::getenv("CFDX_NACA_PRESSURE_REPLAY_PREFIX");
    const char* replay_iters = std::getenv("CFDX_NACA_PRESSURE_REPLAY_ITERATIONS");
    if ((replay_prefix != nullptr) != (replay_iters != nullptr))
        throw std::invalid_argument("NACA pressure replay requires BOTH env variables");
    if (replay_prefix) {
        if (!*replay_prefix)
            throw std::invalid_argument("NACA pressure replay prefix is empty");
        c.pressure_replay_prefix = replay_prefix;
        c.pressure_replay_iterations = cfdx::core::parse_linear_capture_iterations(
            replay_iters);
    }
    c.pressure_reference_cell=0;
    c.pressure_reference_value=0.0;
    c.diagnostics.iteration_trace=true;
    c.diagnostics.iteration_trace_frequency=growth_guard ? 10 : 200;

    const auto geometry=build_fv_geometry(mesh);
    if (const char* checkpoint = std::getenv("CFDX_NACA_CHECKPOINT_STATE_CSV")) {
        const char* frequency = std::getenv("CFDX_NACA_CHECKPOINT_EVERY");
        if (!frequency || !*checkpoint)
            throw std::invalid_argument("NACA checkpoints require a nonempty path and frequency");
        const std::string raw(frequency);
        std::size_t n_read = 0;
        const auto interval = std::stoull(raw, &n_read);
        if (n_read != raw.size() || interval < 10 || interval > 10000)
            throw std::invalid_argument("NACA checkpoint frequency must be in [10,10000]");
        c.iteration_output_callback = [checkpoint_path=std::string(checkpoint),
                                       interval, &geometry](
            std::size_t iteration, double, const Mesh&, const auto& velocity,
            const auto& pressure) {
            if (iteration % interval != 0) return true;
            std::vector<cfdx::validation::NacaStateRow> values;
            values.reserve(geometry.cell_centres.size());
            for (std::size_t i=0;i<geometry.cell_centres.size();++i)
                values.push_back({geometry.cell_centres[i],
                                  velocity.component_data(0)[i],
                                  velocity.component_data(1)[i], pressure(i)});
            cfdx::validation::naca_write_state_csv_atomic(checkpoint_path, values);
            std::cout << "NACA_CHECKPOINT iteration=" << iteration
                      << " status=STATE_ONLY_NOT_CONVERGED\n";
            return true;
        };
    }
    // Cellwise signed momentum evidence: exported only on a requested
    // iteration, never used by the solve or to promote convergence.
    std::ofstream face_budget_output;
    if (const char* budget_path = std::getenv("CFDX_NACA_MOMENTUM_BUDGET_CSV")) {
        if (!*budget_path)
            throw std::invalid_argument("NACA momentum budget path is empty");
        std::size_t every = 10;
        if (const char* interval = std::getenv("CFDX_NACA_MOMENTUM_BUDGET_EVERY")) {
            const std::string raw(interval);
            std::size_t parsed = 0;
            every = std::stoull(raw, &parsed);
            if (parsed != raw.size() || every < 1 || every > 100000)
                throw std::invalid_argument("invalid NACA momentum budget frequency");
        }
        face_budget_output.open(budget_path);
        if (!face_budget_output)
            throw std::runtime_error("cannot open NACA momentum budget output");
        face_budget_output << std::setprecision(17);
        face_budget_output << "iteration,component,rank,cell,x,y,z,"
            "convection,viscous,transpose_viscous,pressure,body_force,"
            "bounded_divergence,mass_divergence,collocated_difference,"
            "physical_residual,bounded_operator_residual,worst_residual\n";
        std::size_t all_cells = 0;
        if (const char* selection = std::getenv("CFDX_NACA_MOMENTUM_BUDGET_ALL_CELLS")) {
            if (std::string(selection) != "1")
                throw std::invalid_argument("NACA full cell budget requires value 1");
            all_cells = 1;
        }
        c.momentum_face_budget_every = every;
        c.momentum_face_budget_callback =
            [&face_budget_output, all_cells](std::size_t iteration, std::size_t d,
                                  const cfdx::physics::MomentumFaceBudget& budget,
                                  const cfdx::physics::FvGeometry& geom) {
            std::vector<std::size_t> order(budget.cells.size());
            for (std::size_t cell = 0; cell < order.size(); ++cell)
                order[cell] = cell;
            const std::size_t top = all_cells
                ? order.size() : std::min<std::size_t>(12, order.size());
            std::partial_sort(order.begin(), order.begin()+static_cast<std::ptrdiff_t>(top), order.end(),
                [&](std::size_t a, std::size_t b) {
                    const double av = std::abs(budget.cells[a].physical_residual);
                    const double bv = std::abs(budget.cells[b].physical_residual);
                    return av == bv ? a < b : av > bv;
                });
            for (std::size_t rank = 0; rank < top; ++rank) {
                const auto cell = order[rank];
                const auto& p = geom.cell_centres[cell];
                const auto& b = budget.cells[cell];
                face_budget_output << iteration << ',' << d << ',' << rank
                    << ',' << cell << ',' << p.x << ',' << p.y << ',' << p.z
                    << ',' << b.convection << ',' << b.viscous
                    << ',' << b.transpose_viscous << ',' << b.pressure
                    << ',' << b.body_force << ',' << b.bounded_divergence
                    << ',' << b.mass_divergence
                    << ',' << b.collocated_difference
                    << ',' << b.physical_residual
                    << ',' << b.bounded_operator_residual
                    << ',' << budget.worst_absolute_residual << '\n';
            }
            face_budget_output.flush();
            if (!face_budget_output)
                throw std::runtime_error("cannot write NACA momentum budget");
        };
    }
    cfdx::validation::NacaMomentumTrace naca_trace(mesh);
    if (const char* forensic = std::getenv("CFDX_NACA_FORENSICS");
        forensic != nullptr && std::string(forensic) == "1") {
        c.momentum_assembly_forensic_callback =
            [&naca_trace](std::size_t iteration, const auto& pressure,
                          const auto& grad_pressure, const auto& velocity,
                          const auto& mass_flux, const auto& mu_eff) {
                naca_trace(iteration, pressure, grad_pressure, velocity,
                           mass_flux, mu_eff);
            };
        c.momentum_equation_forensic_callback =
            [&naca_trace](std::size_t iteration, std::size_t component,
                          const auto& equation, const auto& geometry,
                          const auto& grad_pressure, const auto& body_force) {
                naca_trace.equation(iteration, component, equation, geometry,
                                    grad_pressure, body_force);
            };
    }
    for(std::size_t i=0;i<mesh.n_cells();++i)
        U.component_data(0)[i]=Uinf;
    // Experimental coarse-to-fine continuation; strict mesh alignment.
    // The solver still needs to satisfy every ordinary nonlinear/force gate.
    if (const char* path = std::getenv("CFDX_NACA_IMPORT_INITIAL_STATE_CSV")) {
        const auto init = cfdx::validation::naca_read_state_csv(
            path, geometry.cell_centres);
        for (std::size_t i=0;i<init.size();++i) {
            U.component_data(0)[i] = init[i].ux;
            U.component_data(1)[i] = init[i].uy;
            U.component_data(2)[i] = 0.0;
            p(i) = init[i].pressure;
        }
        std::cout << "NACA_INITIAL_STATE: INTERPOLATED_GUESS_NOT_VALIDATED"
                  << " cells=" << init.size() << "\n";
    }

    // Optional SAME-MESH restart of the last authoritative Rhie-Chow mass
    // flux. An interpolated U,p state on a finer mesh cannot supply this file;
    // those runs must reproject the flux and reconverge independently.
    if (const char* path = std::getenv("CFDX_NACA_IMPORT_FACE_FLUX_CSV")) {
        if (!std::getenv("CFDX_NACA_IMPORT_INITIAL_STATE_CSV"))
            throw std::invalid_argument("NACA face-flux restart requires the matching U,p state");
        c.initial_face_mass_flux = cfdx::validation::naca_read_face_flux_csv(
            path, mesh, geometry);
        std::cout << "NACA_INITIAL_FACE_FLUX: SAME_MESH_GUESS_NOT_A_CERTIFICATE"
                  << " faces=" << c.initial_face_mass_flux.size() << "\n";
    }

    // Exercise the strict mathematical BC/FVM path before the production solve.
    // The legacy maps below are retained for the current steady-solver compatibility API.
    cfdx::physics::BoundaryConstraintMap strict_bc;
    if (split_farfield)
        cfdx::validation::add_velocity_dirichlet(
            strict_bc, "farfield_inlet", {Uinf, 0.0, 0.0});
    else
        cfdx::validation::add_velocity_dirichlet(
            strict_bc, "farfield", {Uinf, 0.0, 0.0});
    if (split_farfield) {
        cfdx::validation::add_velocity_flux_dependent(
            strict_bc, "farfield_outlet", {Uinf,0.0,0.0}, {0.0,0.0,0.0});
        cfdx::validation::add_pressure_neumann(strict_bc, "farfield_inlet", 0.0);
    }
    cfdx::validation::add_velocity_dirichlet(strict_bc, "airfoil", {0.0, 0.0, 0.0});
    cfdx::validation::add_pressure_dirichlet(
        strict_bc, split_farfield ? "farfield_outlet" : "farfield", 0.0);
    cfdx::validation::add_pressure_neumann(strict_bc, "airfoil", 0.0);
    const auto validation_flux = cfdx::validation::make_validation_face_flux(mesh, geometry, U);
    cfdx::validation::exercise_new_velocity_bc_contract(
        mesh, geometry, strict_bc, validation_flux, "NACA0012_RE1000_A0");

    const auto solve_started = std::chrono::steady_clock::now();
    const auto solve=solve_steady_incompressible(mesh,U,p,ubc,pbc,c);
    const double solve_elapsed_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      solve_started).count();
    std::cout << "NACA_SOLVE_WALL elapsed_seconds=" << solve_elapsed_seconds
              << " accepted_iterations=" << solve.history.size()
              << " converged=" << solve.converged
              << " algorithm=" << (c.algorithm == PressureVelocityAlgorithm::SIMPLEC
                     ? "SIMPLEC" : "SIMPLE")
              << " pseudo_momentum=" << c.pseudo_momentum.enabled << '\n';
    // Save physical nonlinear diagnostics even when the full CFD gate fails.
    // The official force/GCI checks still require solve.converged == true.
    // This opt-in creates no matrix dumps and can be used for long runs.
    const auto write_outer_history = [&](const std::string& path) {
        std::ofstream out(path, std::ios::trunc);
        if (!out) throw std::runtime_error("cannot create NACA physical history CSV");
        out << "iteration,continuity_normalized,momentum_rhs_relative,"
               "mass_normalized_imbalance,mass_local_linf,"
               "flux_velocity_mismatch_normalized,pressure_linear_iterations,"
               "momentum_linear_iterations,nonfinite_mass_faces,nonfinite_velocities,"
               "high_order_fraction,face_momentum_x_relative,"
               "face_momentum_y_relative,face_momentum_z_relative,"
               "face_momentum_x_linf,face_momentum_y_linf,face_momentum_z_linf,"
               "convection_scheme,momentum_reference_relative,"
               "momentum_x_rhs_relative,momentum_y_rhs_relative,"
               "momentum_z_rhs_relative,momentum_component_max,"
               "velocity_change_inf,pressure_change_inf,pressure_residual,"
               "alpha_u,alpha_p,pseudo_cfl,nonlinear_metric,"
               "algorithm,pseudo_capped_cells\n";
        out << std::setprecision(17);
        for (const auto& h : solve.history)
            out << h.iteration << ',' << h.continuity_normalized << ','
                << h.momentum_equation_residual_rhs_relative << ','
                << h.mass_normalized_imbalance << ',' << h.mass_local_linf << ','
                << h.flux_velocity_mismatch_normalized << ','
                << h.pressure_linear_iterations << ','
                << h.momentum_linear_iterations << ','
                << h.mass_nonfinite_faces << ','
                << h.boundedness_nonfinite_velocity << ','
                << h.high_order_fraction << ','
                << h.momentum_face_local_relative[0] << ','
                << h.momentum_face_local_relative[1] << ','
                << h.momentum_face_local_relative[2] << ','
                << h.momentum_face_local_linf[0] << ','
                << h.momentum_face_local_linf[1] << ','
                << h.momentum_face_local_linf[2] << ','
                << (c.convection_scheme == ConvectionScheme::TVD ? "TVD" :
                    c.convection_scheme == ConvectionScheme::TVD_VANLEER ? "TVD_VANLEER" :
                    c.convection_scheme == ConvectionScheme::TVD_VANALBADA ? "TVD_VANALBADA" :
                    c.convection_scheme == ConvectionScheme::MUSCL_CELL_LIMITED ? "MUSCL_CELL_LIMITED" :
                    c.convection_scheme == ConvectionScheme::UPWIND ? "UPWIND" :
                    "OTHER") << ','
                << h.momentum_equation_residual_relative << ','
                << h.momentum_equation_residual_rhs_relative_components[0] << ','
                << h.momentum_equation_residual_rhs_relative_components[1] << ','
                << h.momentum_equation_residual_rhs_relative_components[2] << ','
                << h.momentum_equation_residual_rhs_active_component_max << ','
                << h.velocity_change_inf << ',' << h.pressure_change_inf << ','
                << h.pressure_residual << ','
                << h.effective_alpha_u << ',' << h.effective_alpha_p << ','
                << (c.pseudo_momentum.enabled ? h.v17_pseudo_cfl : 0.0) << ','
                << h.nonlinear_convergence_metric << ','
                << (c.algorithm == PressureVelocityAlgorithm::SIMPLEC
                        ? "SIMPLEC" : "SIMPLE") << ','
                << h.v17_pseudo_capped_cells << '\n';
        out.flush();
        if (!out) throw std::runtime_error("failed writing NACA physical history CSV");
    };
    if (const char* path = std::getenv("CFDX_NACA_PHYSICAL_HISTORY_CSV")) {
        if (!*path) throw std::invalid_argument("empty NACA history path");
        write_outer_history(path);
        std::cout << "NACA_PHYSICAL_HISTORY path=" << path
                  << " accepted_iterations=" << solve.history.size()
                  << " converged=" << solve.converged
              << " algorithm=" << (c.algorithm == PressureVelocityAlgorithm::SIMPLEC
                     ? "SIMPLEC" : "SIMPLE")
              << " pseudo_momentum=" << c.pseudo_momentum.enabled << '\n';
    }

    if(!solve.converged) {
        std::ostringstream message;
        message << "NACA0012 solver did not converge: status="
                << cfdx::core::to_string(solve.convergence_status)
                << " reason=" << solve.convergence_reason
                << " iterations=" << solve.iterations;
        if(!solve.history.empty()) {
            const auto& last=solve.history.back();
            message << " momentum_rhs_relative="
                    << last.momentum_equation_residual_rhs_relative
                    << " momentum_reference_relative="
                    << last.momentum_equation_residual_relative
                    << " pressure_residual=" << last.pressure_residual
                    << " continuity=" << last.continuity_normalized
                    << " dU=" << last.velocity_change_inf
                    << " dP=" << last.pressure_change_inf
                    << " alpha_u=" << last.effective_alpha_u
                    << " alpha_p=" << last.effective_alpha_p;
        }
        throw std::runtime_error(message.str());
    }

    if (const char* path = std::getenv("CFDX_NACA_EXPORT_CONVERGED_STATE_CSV")) {
        std::vector<cfdx::validation::NacaStateRow> state;
        state.reserve(mesh.n_cells());
        for (std::size_t i=0;i<mesh.n_cells();++i)
            state.push_back({geometry.cell_centres[i], U.component_data(0)[i],
                             U.component_data(1)[i], p(i)});
        cfdx::validation::naca_write_state_csv(path, state);
        std::cout << "NACA_CONVERGED_STATE: EXPORTED_ONLY_AFTER_PASS cells="
                  << state.size() << "\n";
    }

    if (const char* path = std::getenv("CFDX_NACA_CONVERGED_FACE_FLUX_CSV")) {
        cfdx::validation::naca_write_face_flux_csv(
            path, mesh, geometry, solve.authoritative_mass_flux);
        std::cout << "NACA_CONVERGED_FACE_FLUX: EXPORTED_ONLY_AFTER_PASS"
                  << " faces=" << mesh.n_faces() << "\n";
    }

    ScalarBoundaryConditions gx,gy;
    if (split_farfield) {
        gx["farfield_inlet"]={ScalarBoundaryType::FIXED_VALUE,Uinf,0};
        gx["farfield_outlet"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    } else gx["farfield"]={ScalarBoundaryType::FIXED_VALUE,Uinf,0};
    gx["airfoil"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    gx["front"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    gx["back"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    if (split_farfield) {
        gy["farfield_inlet"]={ScalarBoundaryType::FIXED_VALUE,0,0};
        gy["farfield_outlet"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    } else gy["farfield"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    gy["airfoil"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    gy["front"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    gy["back"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    Field<double,Location::CELL> ux(mesh.n_cells(),"ux","m/s",1), uy(mesh.n_cells(),"uy","m/s",1);
    for(std::size_t i=0;i<mesh.n_cells();++i){ux(i)=U.component_data(0)[i];uy(i)=U.component_data(1)[i];}
    // Post-processing A/B: wall tractions need gradients at physical faces.
    // The rank-certified boundary-aware WLS gradient is affine-exact on a
    // skewed wall stencil; keep the historical Gauss path for comparison.
    const bool wall_wls = std::getenv("CFDX_NACA_WALL_GRADIENT_WLS") != nullptr;
    if (wall_wls && std::string(std::getenv("CFDX_NACA_WALL_GRADIENT_WLS")) != "1")
        throw std::invalid_argument("CFDX_NACA_WALL_GRADIENT_WLS must be 1");
    const auto gux=wall_wls
        ? compute_boundary_aware_wls_gradient(ux,mesh,geometry,gx)
        : gauss_gradient_with_boundary(ux,mesh,geometry,gx);
    const auto guy=wall_wls
        ? compute_boundary_aware_wls_gradient(uy,mesh,geometry,gy)
        : gauss_gradient_with_boundary(uy,mesh,geometry,gy);
    // Optional face-traction A/B; no modification to the historical force
    // gate.  This evaluation remains experimental until MMS + grid studies.
    bool face_traction = false;
    if (const char* mode = std::getenv("CFDX_NACA_WALL_TRACTION_FACE")) {
        if (std::string(mode) != "1")
            throw std::invalid_argument("CFDX_NACA_WALL_TRACTION_FACE must be 1");
        face_traction = true;
    }
    if (wall_wls && !face_traction)
        throw std::invalid_argument(
            "CFDX_NACA_WALL_GRADIENT_WLS requires face-traction post-processing");
    // A/B diagnostic: the production pressure BC uses dp/dn=0 on the foil.
    // Cell-gradient normal extrapolation is informative, not a trustworthy
    // replacement for a physically derived wall-normal pressure BC.
    bool use_cell_normal_dp = false;
    if (const char* mode = std::getenv("CFDX_NACA_WALL_PRESSURE_NORMAL_MODE")) {
        if (!face_traction || (std::string(mode) != "zero" &&
                              std::string(mode) != "cell_gradient"))
            throw std::invalid_argument(
                "CFDX_NACA_WALL_PRESSURE_NORMAL_MODE requires face traction and zero or cell_gradient");
        use_cell_normal_dp = std::string(mode) == "cell_gradient";
    }
    Field<double, Location::CELL> wall_pressure_gradient;
    if (face_traction)
        wall_pressure_gradient = wall_wls
            ? compute_boundary_aware_wls_gradient(p, mesh, geometry, pbc)
            : gauss_gradient_with_boundary(
                p, mesh, geometry, pbc, c.pressure_skew_corrected_gradient);

    const auto& patch=mesh.boundary().patch(mesh.boundary().find("airfoil"));
    std::ofstream face_profiles;
    if (const char* output = std::getenv("CFDX_NACA_FACE_PROFILES_CSV")) {
        if (!face_traction || std::string(output).empty())
            throw std::invalid_argument("NACA Cp/Cf profiles require face traction and a path");
        face_profiles.open(output);
        if (!face_profiles)
            throw std::runtime_error("cannot open NACA face profiles CSV");
        face_profiles << std::setprecision(17);
        face_profiles << "face,x_over_c,y_over_c,Cp,Cf_tangent,Cd_pressure_face,Cd_viscous_face\n";
    }
    std::ofstream face_traction_raw;
    if (const char* output = std::getenv("CFDX_NACA_FACE_TRACTION_RAW_CSV")) {
        if (!face_traction || std::string(output).empty())
            throw std::invalid_argument("NACA raw wall traction requires face traction and a path");
        face_traction_raw.open(output);
        if (!face_traction_raw)
            throw std::runtime_error("cannot create NACA raw wall traction CSV");
        face_traction_raw << std::setprecision(std::numeric_limits<double>::max_digits10);
        face_traction_raw << "face,x,y,z,nx,ny,nz,area,pressure,"
                             "tau_n_x,tau_n_y,tau_n_z,dynamic_pressure,"
                             "reference_area,reference_chord,x_ref,y_ref\n";
    }
    double Fx=0,Fy=0,Mz=0,Fxp=0,Fyp=0,Fxv=0,Fyv=0;
    const double q=0.5*rho*Uinf*Uinf;
    for(const auto f:patch.face_ids){
        const std::size_t o=mesh.ownership().owner(f);
        const auto Sf=geometry.face_area_vectors[f];
        const auto n=Sf.normalized();
        double pp=p(o);
        double txx=2.0*rho*nu*gux.component_data(0)[o];
        double tyy=2.0*rho*nu*guy.component_data(1)[o];
        double txy=rho*nu*(gux.component_data(1)[o]+guy.component_data(0)[o]);
        if (face_traction) {
            const auto wall = reconstruct_wall_face_traction(
                geometry.face_centres[f] - geometry.cell_centres[o], n,
                p(o),
                {wall_pressure_gradient.component_data(0)[o],
                 wall_pressure_gradient.component_data(1)[o],
                 wall_pressure_gradient.component_data(2)[o]},
                {ux(o), uy(o), 0.0}, {0.0,0.0,0.0},
                {{{gux.component_data(0)[o],gux.component_data(1)[o],
                   gux.component_data(2)[o]},
                  {guy.component_data(0)[o],guy.component_data(1)[o],
                   guy.component_data(2)[o]},
                  {0.0,0.0,0.0}}}, rho*nu,
                use_cell_normal_dp ?
                    (wall_pressure_gradient.component_data(0)[o]*n.x +
                     wall_pressure_gradient.component_data(1)[o]*n.y +
                     wall_pressure_gradient.component_data(2)[o]*n.z) : 0.0);
            pp = wall.pressure;
            // Existing total-force integration uses the traction components
            // of the symmetric stress tensor contracted with n.
            // Reconstruct stress components by evaluating the wall gradient;
            // the returned traction is the physically authoritative vector.
            const double face_txv = wall.viscous_on_fluid_normal.x;
            const double face_tyv = wall.viscous_on_fluid_normal.y;
            const double tx = (pp*n.x-face_txv)*Sf.mag();
            const double ty = (pp*n.y-face_tyv)*Sf.mag();
            Fx += tx; Fy += ty;
            Fxp += pp*n.x*Sf.mag(); Fyp += pp*n.y*Sf.mag();
            Fxv -= face_txv*Sf.mag(); Fyv -= face_tyv*Sf.mag();
            const auto r = geometry.face_centres[f]-Vec3{0.25,0,0};
            Mz += r.x*ty-r.y*tx;
            if (face_traction_raw) {
                const auto x_face = geometry.face_centres[f];
                // Current native force coefficient divides by q for the
                // unit chord and unit span case. Export this area explicitly
                // so the postprocessor detects any future normalization drift.
                face_traction_raw << f << ',' << x_face.x << ',' << x_face.y
                                  << ',' << x_face.z << ',' << n.x << ',' << n.y
                                  << ',' << n.z << ',' << Sf.mag() << ',' << pp
                                  << ',' << wall.viscous_on_fluid_normal.x
                                  << ',' << wall.viscous_on_fluid_normal.y
                                  << ',' << wall.viscous_on_fluid_normal.z
                                  << ',' << q << ",1," << chord << ",0.25,0\n";
            }
            if (face_profiles) {
                Vec3 tangent{-n.y,n.x,0.0};
                if (tangent.x < 0.0) tangent *= -1.0;
                const double cf_tangent =
                    -(face_txv*tangent.x + face_tyv*tangent.y) / q;
                face_profiles << f << ',' << geometry.face_centres[f].x/chord
                              << ',' << geometry.face_centres[f].y/chord
                              << ',' << pp/q << ',' << cf_tangent
                              << ',' << (pp*n.x*Sf.mag())/q
                              << ',' << (-face_txv*Sf.mag())/q << '\n';
            }
            continue;
        }
        // Face normal convention: n = Sf/|Sf| points out of the fluid domain
        // (into the body) on wall faces, because build_fv_geometry()
        // enforces Sf.(Cf - C_owner) > 0 on boundary faces. With this
        // fluid-outward n, the force exerted by the fluid on the body is
        //   F = int (p n - tau n) dS
        // i.e. pressure pushes along +p n and viscous traction acts along
        // -tau n (equivalent to the body->fluid normal contract of
        // src/cfdx/physics/forces.h).
        const double tx=(pp*n.x-txx*n.x-txy*n.y)*Sf.mag();
        const double ty=(pp*n.y-txy*n.x-tyy*n.y)*Sf.mag();
        Fx+=tx; Fy+=ty;
        Fxp+=pp*n.x*Sf.mag(); Fyp+=pp*n.y*Sf.mag();
        Fxv+=-(txx*n.x+txy*n.y)*Sf.mag();
        Fyv+=-(txy*n.x+tyy*n.y)*Sf.mag();
        const auto r=geometry.face_centres[f]-Vec3{0.25,0,0};
        Mz+=r.x*ty-r.y*tx;
    }
    // Independent, experimental two-depth wall-shear diagnostic. Does not
    // change the qualified force coefficient or hide missing near-wall cells.
    if (const char* fit = std::getenv("CFDX_NACA_WALL_SHEAR_QUADRATIC")) {
        if (std::string(fit) != "1")
            throw std::invalid_argument("CFDX_NACA_WALL_SHEAR_QUADRATIC must be 1");
        const auto* offsets = mesh.cells().offsets_data();
        const auto* faces = mesh.cells().faces_data();
        double integrated_cd_viscous = 0.0;
        std::ofstream quadratic_profiles;
        if (const char* path = std::getenv("CFDX_NACA_WALL_SHEAR_QUADRATIC_CSV")) {
            quadratic_profiles.open(path);
            if (!quadratic_profiles)
                throw std::runtime_error("cannot open NACA quadratic wall shear CSV");
            quadratic_profiles << std::setprecision(std::numeric_limits<double>::max_digits10)
                               << "face,x_over_c,y_over_c,owner,second_cell,"
                                  "h1,h2,Cf_tangent,Cd_viscous_face\n";
        }
        std::size_t accepted = 0;
        std::size_t missing = 0;
        for (const auto f : patch.face_ids) {
            const std::size_t owner = mesh.ownership().owner(f);
            const Vec3 sf = geometry.face_area_vectors[f];
            const Vec3 outward = sf.normalized();
            Vec3 tangent{-outward.y, outward.x, 0.0};
            if (tangent.x < 0.0) tangent *= -1.0;
            if (!(tangent.mag2() > 0.0)) {
                ++missing;
                continue;
            }
            tangent = tangent.normalized();
            const Vec3 face = geometry.face_centres[f];
            const Vec3 inward = outward * -1.0;
            const double h1 = (geometry.cell_centres[owner]-face).dot(inward);
            double best_score = std::numeric_limits<double>::infinity();
            std::size_t deeper = mesh.n_cells();
            if (!(h1 > 0.0)) { ++missing; continue; }
            for (Offset k = offsets[owner]; k < offsets[owner + 1]; ++k) {
                const auto other_face = faces[k];
                const auto neighbour = mesh.ownership().neighbour(other_face);
                if (neighbour < 0) continue;
                const auto candidate = mesh.ownership().owner(other_face) == owner
                    ? static_cast<std::size_t>(neighbour)
                    : mesh.ownership().owner(other_face);
                if (candidate == owner) continue;
                const Vec3 displacement = geometry.cell_centres[candidate] - face;
                const double h2 = displacement.dot(inward);
                if (!(h2 >= 1.5*h1 && h2 <= 6.0*h1)) continue;
                const Vec3 skew = displacement - inward * h2;
                const double score = std::abs(h2/h1 - 2.0) + 0.3*skew.mag()/h2;
                if (score < best_score) {best_score=score; deeper=candidate;}
            }
            if (deeper == mesh.n_cells()) { ++missing; continue; }
            const auto sample = [&](std::size_t cell) {
                return WallShearSample{
                    geometry.cell_centres[cell],
                    {ux(cell),uy(cell),0.0},
                    {{{gux.component_data(0)[cell],gux.component_data(1)[cell],
                       gux.component_data(2)[cell]},
                      {guy.component_data(0)[cell],guy.component_data(1)[cell],
                       guy.component_data(2)[cell]}, {0.0,0.0,0.0}}}};
            };
            try {
                const auto value = reconstruct_wall_shear_quadratic(
                    face, outward, tangent, {0,0,0},
                    sample(owner), sample(deeper), rho*nu);
                const double drag_face = value.shear_on_body_tangent *
                    tangent.x * sf.mag() / q;
                integrated_cd_viscous += drag_face;
                if (quadratic_profiles)
                    quadratic_profiles << f << ',' << face.x/chord << ','
                        << face.y/chord << ',' << owner << ',' << deeper << ','
                        << value.first_distance << ',' << value.second_distance << ','
                        << value.shear_on_body_tangent/q << ',' << drag_face << '\n';
                ++accepted;
            } catch (const std::invalid_argument&) {
                ++missing;
            }
        }
        if (quadratic_profiles && (!quadratic_profiles.good() || missing != 0))
            throw std::runtime_error("NACA quadratic wall shear CSV incomplete");
        std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
                  << "NACA_WALL_SHEAR_QUADRATIC"
                  << " covered_faces=" << accepted
                  << " missing_faces=" << missing
                  << " total_faces=" << patch.face_ids.size()
                  << " Cd_viscous_partial=" << integrated_cd_viscous
                  << " status=" << (missing == 0 ? "COVERED_NOT_VALIDATED" : "INSUFFICIENT_COVERAGE")
                  << "\n";
    }
    if (face_profiles && !face_profiles.good())
        throw std::runtime_error("failed to write NACA face profiles CSV");
    if (face_traction_raw && (!face_traction_raw.flush() || !face_traction_raw.good()))
        throw std::runtime_error("failed to write NACA raw wall traction CSV");
    const double cd=Fx/q, cl=Fy/q, cm=Mz/q;
    const auto& last=solve.history.back();
    // Optional, machine-readable single-grid certificate. The three-grid
    // Richardson/GCI gate must still validate mesh-family similarity and
    // cannot infer convergence from the existence of this file alone.
    if (const char* output = std::getenv("CFDX_NACA_CONVERGED_FORCE_CSV")) {
        std::ofstream out(output, std::ios::trunc);
        if (!out) throw std::runtime_error("cannot open converged NACA force certificate");
        // The NACA O-grid generator has 2*N+1 airfoil faces per mesh;
        // infer the surface resolution only for that verified topology.
        const auto airfoil_faces = patch.face_ids.size();
        if (airfoil_faces < 3 || airfoil_faces % 2 == 0)
            throw std::runtime_error("unexpected NACA O-grid wall-face topology for GCI metadata");
        const auto surface_resolution = (airfoil_faces - 1) / 2;
        out << "case,h,converged,res_mom,res_mass,cd,cdp,cdv,cl,cm\n";
        out << std::setprecision(std::numeric_limits<double>::max_digits10)
            << 'N' << surface_resolution << ','
            << chord / static_cast<double>(surface_resolution)
            << ",yes," << last.momentum_equation_residual_rhs_relative << ','
            << std::max({last.continuity_normalized, last.mass_normalized_imbalance,
                         last.mass_local_linf}) << ',' << cd << ',' << Fxp/q
            << ',' << Fxv/q << ',' << cl << ',' << cm << '\n';
        out.flush();
        if (!out) throw std::runtime_error("failed writing converged force certificate");
        std::cout << "NACA_CONVERGED_FORCE_CERTIFICATE path=" << output
                  << " qualification=ONE_GRID_ONLY\n";
    }
    std::cout<<"NACA0012_RE1000_A0 level="<<((patch.face_ids.size()-1)/2)
             <<" Cd="<<cd<<" Cl="<<cl<<" Cm="<<cm
             <<" Cd_pressure="<<Fxp/q<<" Cd_viscous="<<Fxv/q
             <<" continuity="<<last.continuity_linf
             <<" iterations="<<solve.iterations<<"\n";
    return {cd,cl,cm,Fxp/q,Fxv/q,last.continuity_linf};
}

static void run_quick_contract(const std::string& path)
{
    Mesh mesh;
    if(!cfdx::io::read_mesh_hdf5(path,mesh))
        throw std::runtime_error("cannot read NACA0012 mesh: "+path);
    const auto topo=mesh.topo_validate();
    if(!topo.ok)
        throw std::runtime_error("invalid NACA0012 topology: "+path);

    const auto geometry=build_fv_geometry(mesh);
    Field<double,Location::CELL> U(mesh.n_cells(),"U","m/s",3);
    U.fill(0.0);
    if (const char* path = std::getenv("CFDX_NACA_EXPORT_CENTERS_CSV"))
        cfdx::validation::naca_write_centers_csv(
            path, geometry.cell_centres);
    for(std::size_t i=0;i<mesh.n_cells();++i)
        U.component_data(0)[i]=1.0;

    cfdx::physics::BoundaryConstraintMap strict_bc;
    const auto split_in = mesh.boundary().find("farfield_inlet");
    const auto split_out = mesh.boundary().find("farfield_outlet");
    const bool split_farfield = split_in < mesh.boundary().n_patches() &&
                                split_out < mesh.boundary().n_patches();
    cfdx::validation::add_velocity_dirichlet(
        strict_bc, split_farfield ? "farfield_inlet" : "farfield", {1.0, 0.0, 0.0});
    if (split_farfield) {
        cfdx::validation::add_velocity_flux_dependent(
            strict_bc, "farfield_outlet", {1.0,0.0,0.0}, {0.0,0.0,0.0});
        cfdx::validation::add_pressure_neumann(strict_bc, "farfield_inlet", 0.0);
    }
    cfdx::validation::add_velocity_dirichlet(strict_bc, "airfoil", {0.0, 0.0, 0.0});
    cfdx::validation::add_pressure_dirichlet(
        strict_bc, split_farfield ? "farfield_outlet" : "farfield", 0.0);
    cfdx::validation::add_pressure_neumann(strict_bc, "airfoil", 0.0);
    const auto validation_flux =
        cfdx::validation::make_validation_face_flux(mesh, geometry, U);
    cfdx::validation::exercise_new_velocity_bc_contract(
        mesh, geometry, strict_bc, validation_flux, "NACA0012_RE1000_A0");

    const auto airfoil_id=mesh.boundary().find("airfoil");
    const auto farfield_id=mesh.boundary().find("farfield");
    if(airfoil_id>=mesh.boundary().n_patches() ||
       (!split_farfield && farfield_id>=mesh.boundary().n_patches()))
        throw std::runtime_error("NACA0012 required boundary patches are missing");

    std::cout<<"NACA0012_QUICK_CONTRACT: PASS cells="<<mesh.n_cells()
             <<" airfoil_faces="<<mesh.boundary().patch(airfoil_id).size()
             <<" farfield_faces="<<(split_farfield ?
                 mesh.boundary().patch(split_in).size() +
                 mesh.boundary().patch(split_out).size() :
                 mesh.boundary().patch(farfield_id).size())
             <<" solver_campaign_pending\n";
}

int main(int argc,char**argv)
{
    try {
    // Isolated quantitative r24/r36/r48/r72 experiment. The complete 3-grid
    // qualification remains separate and may not be promoted by this mode.
    // argc==6 dispatch preserves the 3-arg --single replay defined below.
    if (argc == 6 && std::string(argv[1]) == "--single") {
        if (argc != 6)
            throw std::invalid_argument(
                "usage: test_naca0012_qualification --single mesh.h5 alpha_u alpha_p max_iterations");
        const auto alpha_u = std::stod(argv[3]);
        const auto alpha_p = std::stod(argv[4]);
        const auto maxit = static_cast<std::size_t>(std::stoull(argv[5]));
        if (!(std::isfinite(alpha_u) && alpha_u > 0. && alpha_u <= 1. &&
              std::isfinite(alpha_p) && alpha_p > 0. && alpha_p <= 1. &&
              maxit > 0))
            throw std::invalid_argument("invalid NACA diagnostic continuation controls");
        (void)run(argv[2],64,true,alpha_u,alpha_p,maxit);
        return 0;
    }
    const bool quick=argc>=2 && std::string(argv[1])=="--quick";
    if(quick) {
        if(argc!=3)
            throw std::invalid_argument("usage: test_naca0012_qualification --quick mesh_n64.h5");
        run_quick_contract(argv[2]);
        return 0;
    }

    // Isolated replay for one refinement level. This does not qualify
    // mesh convergence; it preserves the per-level physical acceptance gates.
    if(argc==3 && std::string(argv[1])=="--single") {
        const Forces x=run(argv[2],0);
        if(!(x.continuity<1e-7))
            throw std::runtime_error("NACA0012 continuity gate failed");
        if(!std::isfinite(x.cd) || !std::isfinite(x.cl))
            throw std::runtime_error("NACA0012 force integration produced non-finite QoI");
        if(std::abs(x.cl)>2e-2)
            throw std::runtime_error("NACA0012 alpha=0 symmetry gate failed");
        if(!(x.cdp>0.0))
            throw std::runtime_error("NACA0012 pressure drag sign gate failed");
        std::cout<<"NACA0012_SINGLE: PASS (one level only; not a mesh qualification)\n";
        return 0;
    }

    const std::size_t need=3;
    if(argc!=4)
        throw std::invalid_argument("usage: test_naca0012_qualification [--single mesh.h5 | mesh_n64.h5 mesh_n128.h5 mesh_n256.h5]");
    std::vector<Forces> r;
    for(std::size_t i=0;i<need;++i) r.push_back(run(argv[i+1],64u<<i));
    for(const auto& x:r) {
        if(!(x.continuity<1e-7)) throw std::runtime_error("NACA0012 continuity gate failed");
        if(!std::isfinite(x.cd) || !std::isfinite(x.cl)) throw std::runtime_error("NACA0012 force integration produced non-finite QoI");
        if(std::abs(x.cl)>2e-2) throw std::runtime_error("NACA0012 alpha=0 symmetry gate failed");
        // Regression gate for the wall-normal orientation of the pressure
        // force integration: at alpha=0 the airfoil pressure drag is
        // strictly positive. This catches the historical sign inversion
        // where Fxp accumulated -p n with n pointing into the body.
        if(!(x.cdp>0.0)) throw std::runtime_error("NACA0012 pressure drag sign gate failed (Cd_pressure must be > 0 at alpha=0)");
    }
    if(!quick) {
        const double d01=std::abs(r[1].cd-r[0].cd), d12=std::abs(r[2].cd-r[1].cd);
        std::cout<<"NACA0012 Cd mesh differences="<<d01<<","<<d12<<" literature Re=1000 alpha=0 Cd~0.12\n";
        if(!(d12<=d01)) throw std::runtime_error("NACA0012 drag is not stabilizing under refinement");
    }
    std::cout<<"NACA0012_QUALIFICATION: DIAGNOSTIC solver_converged forces_conservation_symmetry_reference_gate_pending\n";
    } catch (const std::exception& error) {
        std::cerr << "NACA0012_QUALIFICATION: FAIL: " << error.what() << '\n';
        return 1;
    }
}
