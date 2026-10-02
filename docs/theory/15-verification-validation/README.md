# 15 — Verification and Validation

**Status: REPOSITORY-GROUNDED educational chapter. Authoritative maturity remains in the V&V evidence layer.**

## Definitions

**Code verification:** correct implementation of the intended algorithm.

**Solution verification:** sufficient numerical resolution and convergence for the stated question.

**Validation:** comparison against an independent physical observation or trusted reference.

**Qualification:** evidence that a declared capability is demonstrated over a defined scope.

These questions are not interchangeable.

## Evidence chain

$$
requirement
\rightarrow mathematical\ contract
\rightarrow implementation
\rightarrow verification
\rightarrow solution\ verification
\rightarrow validation
\rightarrow qualification.
$$

## Manufactured solutions

Choose u_e, derive the forcing required by the continuous PDE and compare the computed solution against u_e. The exact solution and source must be independent of the production discretisation.

## Refinement

For

$$
E_h\sim Ch^p,
$$

$$
p_{obs}=\frac{\ln(E_h/E_{rh})}{\ln r}.
$$

The study must identify the population, ratio, norm and boundary treatment.

## Iterative convergence

A case must not PASS solely because an iteration limit was reached or a solver-reported residual decreased. Independent evidence can include true residual, mass imbalance, local conservation, boundedness and QoI stabilisation.

## Conservation

For internal interfaces, opposite-oriented fluxes must cancel. Global balances should be reconstructed independently where practical.

## Validation

A validation result identifies the independent reference, uncertainty, model discrepancy and comparison population.

## Reproducibility

Evidence retains exact commit, case/setup, mesh identity, numerical-method selection, compiler/build mode, hardware/backend, raw machine-readable metrics, report and acceptance rule.

The V&V governance and validation documents remain the authoritative status source.


## Scientific explanation standard

For every important equation, algorithm or data-model rule, document the motivation; definitions/units/sign conventions; derivation or formal rationale; discrete/FVM representation where applicable; actual CFDX execution path; example; error/limitation/sensitivity analysis; exact implementation files; executable verification; benchmark/V&V evidence; and stable bibliography/reference identifiers. Tables summarize explanations and do not replace them. Keep **Implemented / Verified / Validated / Qualified** distinct; documentation maturity never promotes numerical maturity.