// M1 — Navier-Stokes Incompressible — SIMPLE Algorithm
//
// Spécification CFDX v0.7 §36-38 :
//   - Navier-Stokes incompressible: ∂u/∂t + ∇·(u⊗u) = -∇p/ρ + ν∇²u + f
//   - Continuité: ∇·u = 0
//   - Couplage presseur-vitesse: SIMPLE
//
// Steps:
//   1. Assemble momentum equation (convection + diffusion + pressure gradient)
//   2. Solve momentum for u* (predicted velocity, using p^old)
//   3. Assemble pressure correction equation (from continuity)
//   4. Solve pressure correction p'
//   5. Update velocity: u^{n+1} = u* + (D/ν) * ∇p'
//   6. Update pressure: p^{n+1} = p^old + α_p * p'
//   7. Iterate until convergence

#pragma once

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/linalg/linear_system.h"
#include "cfdx/core/numerics/source_term.h"
#include "cfdx/core/numerics/temporal.h"
#include <vector>
#include <string>
#include <cstdio>

namespace cfdx {
namespace solvers {

struct NSParameters {
    double rho = 1.0;
    double nu = 1e-3;
    double dt = 1e-2;
    double tol_velocity = 1e-6;
    double tol_pressure = 1e-6;
    int max_iterations = 100;
    int max_inner = 50;
    double under_relaxation_u = 0.8;
    double under_relaxation_p = 0.2;
    bool use_local_time_stepping = false;
    TimeScheme time_scheme = TimeScheme::EULER_IMPLICIT;
    std::string convection_scheme = "upwind";
};

struct NSResult {
    bool converged = false;
    int iterations = 0;
    double residual_velocity = 0.0;
    double residual_pressure = 0.0;
    std::string message;
};

NSResult solve_navier_stokes_simple(const Mesh& mesh,
                                     Field<double, Location::CELL>& u,
                                     Field<double, Location::CELL>& v,
                                     Field<double, Location::CELL>& w,
                                     Field<double, Location::CELL>& p,
                                     const NSParameters& params);

double compute_velocity_divergence(const Mesh& mesh,
                                    const Field<double, Location::CELL>& u,
                                    const Field<double, Location::CELL>& v,
                                    const Field<double, Location::CELL>& w);

void compute_pressure_gradient(const Mesh& mesh,
                                const Field<double, Location::CELL>& p,
                                std::vector<double>& pg_u,
                                std::vector<double>& pg_v,
                                std::vector<double>& pg_w);

}  // namespace solvers
}  // namespace cfdx
