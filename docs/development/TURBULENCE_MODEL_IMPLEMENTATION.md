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

The model-selection layer now dispatches all registered LES/hybrid identities
(Smagorinsky, WALE, dynamic one-equation LES, DES, DDES and IDDES) through a
single canonical transport API. Smagorinsky, WALE and DES-family closures have
deterministic algebraic kernels, cell-scale requirements, and explicit
shielding/length-scale branches. Dynamic LES exposes its dynamic coefficient
contract with bounded coefficient evaluation.

DES/DDES/IDDES remain dependent on a transient 3-D scale-resolving driver and
base-RANS coupling for production use; the selector therefore does not falsely
promote these kernels to solver-ready status.

## Wall treatment

Wall y+ is explicitly classified into:
- viscosity-affected region: y+ <= 5;
- buffer region: 5 < y+ < 30;
- logarithmic region: y+ >= 30.

This classification is diagnostic/model infrastructure. The checked classifier rejects non-finite or negative y+. Prevalidated hot loops may call `classify_wall_y_plus_unchecked`. It does not claim a solver-level wall-function V&V result.

## Validation boundary

Quantitative channel, flat-plate, cavity, DNS/experimental comparisons, grid convergence, observed order, and VMFL PASS promotion remain under the validation workstream (#118). This PR only adds equation-level model kernels and deterministic regression tests.

## Corrections and model extensions

The canonical selection schema now exposes composable correction families rather
than creating a new model ID for every combination:

- rotation / streamline-curvature;
- compressibility;
- wall roughness;
- production limiting;
- Kato-Launder rotation-aware production;
- SA QCR.

The C++ layer validates the correction controls and provides deterministic
invariant-based kernels. This follows the separation used by mature CFD
interfaces: SU2, for example, exposes base turbulence models separately from
their corrections, while Fluent documents curvature/rotation corrections for
SA and multiple two-equation families. citeturn0search16turn0search18

The remaining implementation work is solver-level rather than selector-level:
SA-negative/Edwards and exact literature-specific SA/SST correction variants,
complete wall-function boundary production, full dynamic LES transport/filter
operations, and complete DES/DDES/IDDES transient coupling. These must not be
declared validated until the dedicated V&V campaign is prepared.


The existing WALE, DES, DDES and IDDES functions remain algebraic closures or
length-scale kernels until their transport, boundary conditions and validation
cases are implemented.

## Upstream references

- [OpenFOAM realizableKE source](https://cpp.openfoam.org/v13/realizableKE_8C_source.html)
- [OpenFOAM turbulence model overview](https://doc.openfoam.com/2212/tools/processing/models/turbulence/)
- [SU2 turbulence and transition options](https://su2code.github.io/docs_v7/Physical-Definition/)
