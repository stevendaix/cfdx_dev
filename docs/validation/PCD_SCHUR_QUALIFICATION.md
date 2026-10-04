# PCD Schur approximation — N8 implementation contract

## Scope

This document defines the CFDX algebraic implementation of the pressure-convection-diffusion (PCD) Schur inverse approximation.

For the coupled pressure-velocity operator

    A = [ Auu  G ]
        [ D    C ],

the exact pressure Schur complement is

    S = C - D Auu^{-1} G.

For incompressible flow, the PCD approximation used here is

    S ~= -Kp Mp^{-1} Fp,

and therefore

    S^{-1} ~= -Fp^{-1} Mp Kp^{-1}.

The implementation is an algebraic Schur approximation. It does not claim that
a particular pressure discretisation is universally correct.

## Operator contract

The caller supplies three pressure-space operators:

- Mp: pressure mass operator;
- Kp: pressure diffusion/Laplacian operator;
- Fp: pressure convection-diffusion operator.

The caller also supplies independent solves for Kp and Fp.

The apply sequence is exactly:

1. z = Kp^{-1} rhs;
2. y = Mp z;
3. pressure = -Fp^{-1} y.

No dense inverse is formed and no pressure operator is inferred from Auu, G
or D.

This separation is intentional: PCD depends on the pressure convection and
diffusion discretisations, so constructing those operators implicitly from the
saddle-point blocks would hide a numerical assumption.

## Sign convention

CFDX defines the saddle-point Schur complement as

    S = C - D Auu^{-1} G.

For the usual incompressible convention with C = 0, the Schur complement is
negative pressure-Laplacian-like. The PCD inverse consequently carries the
explicit leading minus sign.

The unit test checks the complete action against an independent dense
calculation of -Fp^{-1} Mp Kp^{-1}.

## Lifecycle

setup() validates:

- saddle-point block dimensions;
- pressure operator dimensions;
- finite pressure-operator values;
- presence of both pressure solve callbacks;
- pressure null-space dimension when a policy is configured.

update_values(blocks) preserves the common SchurApproximation contract.

Because BlockOperator intentionally contains only the saddle-point blocks,
PCD exposes an explicit update_pressure_values(...) overload for its three
additional pressure operators. This path:

- accepts value-only updates when the CSR graphs are unchanged;
- rejects any pressure-operator graph change;
- updates the retained operator references only after all validation succeeds.

A symbolic graph change therefore requires a fresh setup()/new PCD object
rather than silently reusing stale symbolic assumptions.

## Null-space policy

When a pressure null-space projector is supplied:

- an incompatible input right-hand side is rejected;
- compatible input is projected before the pressure solves;
- intermediate and final pressure vectors are projected;
- no incompatible component is silently manufactured or accepted.

The policy is optional because not every pressure formulation has the same
gauge/null-space structure.

## Current verification

The implementation PR provides deterministic unit coverage for:

- exact PCD algebraic action;
- method name/registration surface;
- dimension rejection;
- pressure graph-change rejection;
- value-only pressure updates;
- null-space compatibility.

This is implementation and algebraic verification, not production qualification.

## Qualification still required

PCD must not be marked N8-qualified solely from these unit tests. The remaining
qualification work is:

- representative Couette, Poiseuille and Ghia matrices;
- controlled skew/non-orthogonal matrices;
- independent true-residual measurement;
- FGMRES iteration and convergence evidence;
- pressure Kp/Fp solve setup/solve cost;
- memory/NNZ accounting;
- null-space/gauge coverage on production pressure systems;
- numerical-update lifecycle in the production dispatcher;
- measured comparison against the exact-Schur oracle;
- evidence-based acceptance envelope before automatic selection.

No arbitrary approximation-quality threshold is introduced by this implementation.
