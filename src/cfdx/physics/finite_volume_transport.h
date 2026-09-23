        solution(0) += controls.relaxation * (candidate_value - solution(0));
        return {
            residual <= controls.tolerance
                ? cfdx::core::SolverStatus::CONVERGED
                : cfdx::core::SolverStatus::MAX_ITER_REACHED,
            1, residual, residual
        };
    }

    cfdx::core::Vector candidate = solution;
    auto result = cfdx::core::solve_bicgstab(
        equation.matrix, equation.rhs, candidate,
        controls.max_iterations, controls.tolerance);

    // Momentum matrices can move between nearly symmetric diffusion-dominated
    // states and mildly nonsymmetric convection-dominated states. BiCGStab may
    // either stagnate or break down on the former even though the linear
    // system is well posed. Retry from the original nonlinear iterate with
    // restarted GMRES, then with CG as a final robust path. Every retry starts
    // from the original iterate so no unconverged Krylov state is injected.
    if (result.status != cfdx::core::SolverStatus::CONVERGED) {
        candidate = solution;
        result = cfdx::core::solve_gmres(
            equation.matrix, equation.rhs, candidate,
            64, controls.max_iterations, controls.tolerance);
    }

    if (result.status != cfdx::core::SolverStatus::CONVERGED) {
        candidate = solution;
        result = cfdx::core::solve_cg(