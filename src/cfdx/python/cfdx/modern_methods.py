"""Modern numerical techniques, with explicit CFDX implementation status.

This is an educational/diagnostic catalog, not a solver configuration backend.
Only adapter-mapped selections are claimed as adapter mappings; they are NOT
promoted to tested/production-ready numerical implementations by this module.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class MethodInfo:
    key: str
    family: str
    label: str
    integration: str
    purpose: str
    qualification: str
    reference: str


METHODS: tuple[MethodInfo, ...] = (
    MethodInfo(
        "SIMPLE",
        "Pressure–velocity",
        "SIMPLE",
        "Adapter mapping exists",
        "Segregated pressure/velocity updates.",
        "Independently re-evaluate continuity, momentum, face-flux and solution-change norms.",
        "https://doc.cfd.direct/openfoam/user-guide-v13/fvsolution",
    ),
    MethodInfo(
        "SIMPLEC",
        "Pressure–velocity",
        "SIMPLEC",
        "Adapter mapping exists",
        "Pressure correction using a modified momentum coefficient treatment.",
        "Compare physical residual and conservation with SIMPLE on the same case.",
        "https://doc.cfd.direct/openfoam/user-guide-v13/fvsolution",
    ),
    MethodInfo(
        "COUPLED",
        "Pressure–velocity",
        "Block COUPLED",
        "Adapter mapping exists",
        "Joint pressure/velocity system; Schur selection must be explicit.",
        "Verify block operator, Rhie–Chow flux, PCD/Schur and nonlinear residual on the returned state.",
        "https://petsc.org/release/manual/ksp/",
    ),
    MethodInfo(
        "PCD",
        "Preconditioners",
        "PCD Schur",
        "Adapter mapping exists",
        "Pressure convection-diffusion approximation for coupled Schur solves.",
        "Check true coupled-system residual and preconditioner reuse under anisotropy.",
        "https://petsc.org/release/manual/ksp/",
    ),
    MethodInfo(
        "AMG",
        "Preconditioners",
        "Algebraic multigrid",
        "Adapter mapping exists",
        "Accelerate elliptic/diffusive linear systems.",
        "Test anisotropic meshes, rebuild conditions and true residual after preconditioning.",
        "https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html",
    ),
    MethodInfo(
        "FGMRES",
        "Linear algebra",
        "Flexible GMRES",
        "Adapter mapping exists",
        "Krylov method supporting a changing preconditioner.",
        "Compare against GMRES/CG for appropriate operator properties and equal true residual thresholds.",
        "https://petsc.org/release/manual/ksp/",
    ),
    MethodInfo(
        "WLS",
        "Discretization",
        "Weighted least squares",
        "Adapter mapping exists",
        "Reconstruct gradients on non-orthogonal/skewed meshes.",
        "Use affine/quadratic manufactured solutions and boundary face reconstruction tests.",
        "https://www.openfoam.com/documentation/guides/latest/doc/guide-schemes-gradient.html",
    ),
    MethodInfo(
        "TVD",
        "Discretization",
        "Bounded TVD convection",
        "Adapter mapping exists",
        "Limit numerical oscillations while retaining second-order accuracy where possible.",
        "Require converged full TVD physical residual, not merely a first-order continuation stage.",
        "https://www.openfoam.com/documentation/guides/latest/doc/guide-schemes-divergence.html",
    ),
    MethodInfo(
        "NEWTON_KRYLOV",
        "Nonlinear acceleration",
        "Newton–Krylov",
        "Research candidate — no native integration certified",
        "Use a Jacobian-vector action derived from the full nonlinear finite-volume residual.",
        "Prove Jv via directional finite differences and full conservation before enabling in production.",
        "https://petsc.org/release/manual/snes/",
    ),
    MethodInfo(
        "ANDERSON",
        "Nonlinear acceleration",
        "Anderson mixing",
        "Research candidate — no native integration certified",
        "Accelerate a nonlinear fixed-point iteration using a limited iterate history.",
        "Benchmark against unchanged tolerances, preserve safeguards, and reject nonphysical iterates.",
        "https://petsc.org/release/manual/snes/",
    ),
    MethodInfo(
        "NGMRES",
        "Nonlinear acceleration",
        "Nonlinear GMRES",
        "Research candidate — no native integration certified",
        "Combine fixed-point iterates to minimize an independently evaluated nonlinear residual.",
        "Require repeated-restart tests, residual decreases and conservation on final physical scheme.",
        "https://petsc.org/release/manual/snes/",
    ),
)


def methods_for_family(family: str = "") -> tuple[MethodInfo, ...]:
    if family:
        return tuple(method for method in METHODS if method.family == family)
    return METHODS


def method_by_key(key: str) -> MethodInfo:
    for method in METHODS:
        if method.key == key:
            return method
    raise KeyError(key)


def method_families() -> tuple[str, ...]:
    return tuple(dict.fromkeys(method.family for method in METHODS))


def search_methods(
    query: str = "", *, family: str = "", integration: str = "all"
) -> tuple[MethodInfo, ...]:
    """Search trustworthy catalog descriptions without changing solver state.

    `integration` filters a documented mapping or research status, never a
    production approval.  Unknown filters are rejected, not silently ignored.
    """
    if integration not in {"all", "mapped", "research"}:
        raise ValueError(f"invalid method integration filter: {integration!r}")
    needle = query.casefold().strip()
    methods = methods_for_family(family)
    return tuple(
        method
        for method in methods
        if (
            not needle
            or needle
            in " ".join(
                (
                    method.key,
                    method.label,
                    method.family,
                    method.integration,
                    method.purpose,
                    method.qualification,
                )
            ).casefold()
        )
        and (
            integration == "all"
            or (
                integration == "mapped"
                and method.integration == "Adapter mapping exists"
            )
            or (
                integration == "research" and "Research candidate" in method.integration
            )
        )
    )


def compare_methods(left: str, right: str) -> str:
    """Build a reversible, evidence-aware comparison, not a solver selection."""
    if left == right:
        raise ValueError("comparison needs two different numerical methods")
    a, b = method_by_key(left), method_by_key(right)
    return (
        "METHOD COMPARISON (READ ONLY)\n"
        "Adapter mapping is not certification; no current-HEAD benchmark is included.\n\n"
        f"A: {a.label} ({a.family})\n"
        f"Integration: {a.integration}\n"
        f"Purpose: {a.purpose}\n"
        f"Required evidence: {a.qualification}\n"
        f"Reference: {a.reference}\n\n"
        f"B: {b.label} ({b.family})\n"
        f"Integration: {b.integration}\n"
        f"Purpose: {b.purpose}\n"
        f"Required evidence: {b.qualification}\n"
        f"Reference: {b.reference}\n\n"
        "Compare on identical case, mesh, final scheme, tolerances and measured residuals."
    )
