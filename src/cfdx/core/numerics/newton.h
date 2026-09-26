#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/preconditioner.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>

namespace cfdx::core {

// Generic Newton solve for a cell field. The callback assembles the residual
// field R(phi) and its sparse Jacobian J = dR/dphi. The field is flattened in
// component-major order: dof = component*n_cells + cell.
struct NewtonFieldControls {
    std::size_t max_iterations = 25;
    double tolerance = 1e-10;
    double linear_tolerance = 1e-12;
    std::size_t linear_max_iterations = 500;
    int gmres_restart = 30;
    double line_search_min = 1e-6;
};

struct NewtonFieldResult {
    SolverStatus status = SolverStatus::NOT_APPLICABLE;
    std::size_t iterations = 0;
    double residual = std::numeric_limits<double>::infinity();
    double residual_relative = std::numeric_limits<double>::infinity();
};

using NewtonFieldAssembler = std::function<void(
    const Field<double, Location::CELL>&,
    Field<double, Location::CELL>&,
    SparseMatrix&)>;

inline SolverResult solve_newton_field(
    Field<double, Location::CELL>& phi,
    const NewtonFieldAssembler& assemble,
    const NewtonFieldControls& controls = {})
{
    if (!assemble || phi.size() == 0 || phi.dimension() == 0 ||
        controls.max_iterations == 0 || !(controls.tolerance > 0.0) ||
        !(controls.linear_tolerance > 0.0) ||
        controls.linear_max_iterations == 0 ||
        controls.gmres_restart <= 0 ||
        !(controls.line_search_min > 0.0 && controls.line_search_min <= 1.0))
        throw std::invalid_argument("solve_newton_field: invalid controls or field");

    const std::size_t n = phi.size() * phi.dimension();
    Field<double, Location::CELL> residual(
        phi.size(), phi.name() + "_residual", phi.metadata().unit, phi.dimension());

    SolverResult result;
    double initial_norm = 0.0;

    auto flatten = [n, &phi](const Field<double, Location::CELL>& f, Vector& v) {
        if (f.size() != phi.size() || f.dimension() != phi.dimension())
            throw std::invalid_argument("solve_newton_field: field shape mismatch");
        for (std::size_t d = 0; d < f.dimension(); ++d)
            for (std::size_t c = 0; c < f.size(); ++c)
                v(d * f.size() + c) = f.component_data(d)[c];
        (void)n;
    };

    auto unflatten_add = [](Field<double, Location::CELL>& f, const Vector& delta, double alpha) {
        for (std::size_t d = 0; d < f.dimension(); ++d)
            for (std::size_t c = 0; c < f.size(); ++c)
                f.component_data(d)[c] += alpha * delta(d * f.size() + c);
    };

    Vector r(n, 0.0);
    for (std::size_t iteration = 0; iteration < controls.max_iterations; ++iteration) {
        SparseMatrix jacobian(n, n);
        assemble(phi, residual, jacobian);
        flatten(residual, r);
        const double norm = r.norm2();
        if (!std::isfinite(norm))
            throw std::runtime_error("solve_newton_field: non-finite nonlinear residual");

        if (iteration == 0) initial_norm = norm;
        const double scale = std::max(initial_norm, 1e-300);
        if (norm <= controls.tolerance * scale) {
            result.status = SolverStatus::CONVERGED;
            result.iterations = iteration;
            result.residual = norm;
            result.residual_relative = norm / scale;
            return result;
        }

        jacobian.finalize();
        Vector rhs(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) rhs(i) = -r(i);
        Vector delta(n, 0.0);
        JacobiPreconditioner preconditioner;
        if (!preconditioner.setup(jacobian))
            throw std::runtime_error("solve_newton_field: invalid Jacobian diagonal");

        LinearOperator op;
        op.size = n;
        op.apply = [&jacobian](const Vector& x, Vector& y) { y = jacobian.matvec(x); };
        const auto linear = solve_gmres(
            op, rhs, delta, controls.gmres_restart, controls.linear_max_iterations,
            controls.linear_tolerance, &preconditioner);
        if (linear.status != SolverStatus::CONVERGED)
            throw std::runtime_error("solve_newton_field: Newton linear solve failed");

        // Backtracking on the nonlinear residual prevents a full Newton step
        // from increasing the residual for strongly nonlinear field equations.
        const Field<double, Location::CELL> old = phi;
        double alpha = 1.0;
        bool accepted = false;
        while (alpha >= controls.line_search_min) {
            phi = old;
            unflatten_add(phi, delta, alpha);
            assemble(phi, residual, jacobian);
            flatten(residual, r);
            const double trial_norm = r.norm2();
            if (std::isfinite(trial_norm) && trial_norm < norm) {
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if (!accepted) {
            phi = old;
            result.status = SolverStatus::DIVERGED;
            result.iterations = iteration + 1;
            result.residual = norm;
            result.residual_relative = norm / scale;
            return result;
        }
    }

    result.status = SolverStatus::MAX_ITER;
    result.iterations = controls.max_iterations;
    result.residual = r.norm2();
    result.residual_relative = result.residual / std::max(initial_norm, 1e-300);
    return result;
}

} // namespace cfdx::core
