# Turbulence model implementation hardening

This document records model-level equations implemented or exposed by CFDX. It is an implementation/verification document, not solver-level V&V.

## Standard k-epsilon

The transported model uses the standard production/destruction structure:
- nu_t = C_mu k^2 / epsilon
- P_k = rho nu_t S^2
- k destruction = rho epsilon
- epsilon production = C1 P_k epsilon/k
- epsilon destruction = C2 rho epsilon^2/k

Positivity floors are applied to k and epsilon by the transport-solver layer before algebraic evaluation. Public model helpers validate finite, physically admissible inputs. Explicit `_unchecked` kernels are available for prevalidated hot loops and retain debug assertions.

## RNG k-epsilon

The RNG epsilon production coefficient is evaluated from

C1* = C1 - eta (1 - eta/eta0) / (1 + beta eta^3)

with eta = S k / epsilon.

The transport solver uses the dedicated RNG diffusion coefficients and C_mu.

## Realizable k-epsilon

The model-level kernel exposes the variable C_mu formulation

C_mu = 1 / (A0 + As U* k/epsilon)

where U* combines strain and rotation invariants and As is obtained from the bounded third strain invariant. The checked helper rejects invalid inputs in every build. Its explicit `_unchecked` counterpart requires finite, non-negative strain/rotation/k and strictly positive epsilon/A0, with debug assertions guarding that contract.

The transported model is available through `solve_realizable_kepsilon_transport`.
It uses the variable C_mu closure, the realizable epsilon production coefficient
`max(0.43, eta/(eta+5))`, and the regularized epsilon destruction denominator
`k + sqrt(nu epsilon)`. The default A0, C2, sigma_k and sigma_epsilon values
match OpenFOAM's `realizableKE` implementation. Rotation magnitude and the
normalized third strain invariant can be supplied when the velocity-gradient
invariants are available; otherwise the solver uses the strain-only form.

## k-omega / SST

SST blending is explicit for sigma_k, sigma_omega, beta and gamma. The existing SST solver separately evaluates:
- F1/F2 blending;
- cross-diffusion;
- production limiting;
- blended diffusion coefficients;
- Menter eddy-viscosity limiter.

The eddy-viscosity denominator remains max(a1 omega, F2 S).

## Spalart-Allmaras

The existing transport implementation retains explicit SA constants and the modified-vorticity, production, destruction, trip-function and gradient-diffusion terms. Positivity of nu-tilde is enforced after each scalar solve.

## LES / hybrid models

Smagorinsky, WALE and the algebraic DES/DDES/IDDES length-scale helpers are model kernels. They are not claimed to be complete LES/DES solver implementations merely because the algebraic helpers exist.

## Wall treatment

Wall y+ is explicitly classified into:
- viscosity-affected region: y+ <= 5;
- buffer region: 5 < y+ < 30;
- logarithmic region: y+ >= 30.

This classification is diagnostic/model infrastructure. The checked classifier rejects non-finite or negative y+. Prevalidated hot loops may call `classify_wall_y_plus_unchecked`. It does not claim a solver-level wall-function V&V result.

## Validation boundary

Quantitative channel, flat-plate, cavity, DNS/experimental comparisons, grid convergence, observed order, and VMFL PASS promotion remain under the validation workstream (#118). This PR only adds equation-level model kernels and deterministic regression tests.

## Remaining model gaps

Comparison with the current SU2 and OpenFOAM model families leaves these major
items deliberately unadvertised as complete transport models:

- SA-negative and the SA rotation/compressibility/QCR variants;
- transition transport such as Langtry-Menter gamma-Re-theta;
- Reynolds-stress transport and non-linear eddy-viscosity RANS models;
- transported LES and hybrid DES/DDES/IDDES models with shielding and wall treatment;
- complete turbulence wall functions coupled to boundary production terms.

The existing WALE, DES, DDES and IDDES functions remain algebraic closures or
length-scale kernels until their transport, boundary conditions and validation
cases are implemented.

## Upstream references

- [OpenFOAM realizableKE source](https://cpp.openfoam.org/v13/realizableKE_8C_source.html)
- [OpenFOAM turbulence model overview](https://doc.openfoam.com/2212/tools/processing/models/turbulence/)
- [SU2 turbulence and transition options](https://su2code.github.io/docs_v7/Physical-Definition/)
