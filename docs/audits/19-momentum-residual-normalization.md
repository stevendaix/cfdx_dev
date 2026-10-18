# METHOD SHEET — Momentum residual normalization in steady incompressible convergence

**Status:** DRAFT — candidate correction, exact-HEAD CI and `validation-total` required before qualification.

**Scope:** acceptance gate and nonlinear monitor in `solve_steady_incompressible()`. This sheet does not change momentum/pressure discretization, pressure correction, boundary conditions, the configured tolerance, the outer-iteration budget, or physical QoI gates.

**Related work:** N12 coupled pressure-velocity MMS and issue #905; candidate fix PR #956.

## 1. Mathematical contract

For a momentum equation (A_u U=b_u), let (r_u=A_u U-b_u). CFDX records two distinct ratios:

- RHS-scaled residual: (\eta_b=\|r_u\|/\max(\|b_u\|,\epsilon)).
- First-iteration ratio: (\eta_0=\|r_u\|/\max(\|r_u^{(0)}\|,\epsilon)).

They answer different questions. The RHS-scaled residual is the acceptance metric for the momentum equation. The first-iteration ratio is useful to describe progress relative to a particular starting state, but it is not a robust absolute acceptance criterion when the initial residual is already at round-off.

## 2. Failure mode and correction

If (\|r_u^{(0)}\|\approx\epsilon_{machine}), dividing a later round-off residual by that near-zero initial value can produce a ratio substantially greater than one (the N12 diagnostic reported approximately 27.8) even while the raw and RHS-scaled residuals and continuity are at numerical noise. Requiring that ratio to satisfy the same tolerance as the RHS-scaled residual can therefore reject an algebraically converged state because of an ill-conditioned denominator.

The candidate change uses the RHS-scaled momentum residual in the nonlinear monitor and acceptance gate. The first-iteration ratio remains recorded in iteration history for diagnostics and adaptive-relaxation analysis. Continuity, pressure, velocity-change, pressure-change and QoI gates remain in force.

## 3. Discretization and implementation contract

No discrete operator is changed. The existing RHS scale and residual assembly are retained; only the interpretation of the two already-recorded residual metrics in the acceptance path is corrected. The RHS-scaled residual must remain finite, and all other configured gates must still pass. This is not permission to bypass a failed physical or conservation criterion.

## 4. Verification plan

- Run the existing N12 coupled pressure-velocity MMS, including the near-exact initial-state diagnostic.
- Run Release and DebugSanitizers builds, unit/fast regression tests, N9 qualification, and the full `validation-total` campaign on the exact candidate SHA.
- Inspect the reported raw residual, RHS scale, RHS-relative residual, first-iteration ratio, continuity, and velocity/pressure changes.
- Confirm that no tolerance, iteration budget, expected-failure classification, or physical oracle changed.
- Compare any remaining long-campaign failures against the same base SHA before attributing them to this change.

## 5. References

- Y. Saad and M. H. Schultz, “GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems,” *SIAM Journal on Scientific and Statistical Computing*, 7(3), 1986, pp. 856–869. DOI: 10.1137/0907058. (Residual minimization context; not a justification for the nonlinear acceptance gate itself.)
- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere, 1980. DOI: 10.1201/9781482234213. (Iterative solution and convergence monitoring context.)

This sheet documents the intended numerical contract; it does not itself constitute validation evidence.
