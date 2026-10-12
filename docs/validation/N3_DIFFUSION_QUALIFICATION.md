# N3 — Diffusion and non-orthogonal correction qualification

## Scope

N3 is qualified against the acceptance contract defined by issue #461:

- orthogonal finite-volume diffusion;
- corrected non-orthogonal diffusion;
- limited non-orthogonal correction;
- over-relaxed diffusion;
- explicit skewness/non-orthogonality treatment;
- quantitative controlled non-orthogonal verification;
- no silent scheme substitution;
- explicit rejection of the singular over-relaxed decomposition;
- mesh-quality diagnostics exposed independently from the diffusion assembly.

This document deliberately distinguishes **N3 qualification** from the separate research problem of higher-order diffusion on arbitrary polyhedral meshes.

## Discrete contract

For an internal face with owner P, neighbour N, centre vector

$$
\mathbf d = \mathbf x_N-\mathbf x_P
$$

and owner-oriented area vector \(\mathbf S_f\), the two-point contribution is

$$
F_f^{orth}=\alpha(\phi_N-\phi_P),
\qquad
\alpha=\frac{\mathbf S_f\cdot\mathbf d}{|\mathbf d|^2}.
$$

The corrected family decomposes

$$
\mathbf S_f=\alpha\mathbf d+
\left(\mathbf S_f-\alpha\mathbf d\right).
$$

The explicit correction is

$$
F_f^{corr}=
\left(\mathbf S_f-\alpha\mathbf d\right)\cdot\nabla\phi_f.
$$

For the over-relaxed formulation,

$$
\alpha=\frac{|\mathbf S_f|^2}{\mathbf S_f\cdot\mathbf d},
$$

so the residual vector is orthogonal to \(\mathbf S_f\).

If \(\mathbf S_f\cdot\mathbf d=0\), the over-relaxed coefficient is singular. CFDX rejects this case explicitly; it does not manufacture a finite coefficient or silently switch scheme.

The limited scheme applies a coefficient \(\lambda\in[0,1]\) to the explicit correction. The implementation validates the configured limit and keeps the correction bounded relative to the orthogonal contribution.

## Verification evidence

### Unit-level invariants

'tests/unit/test_laplacian.cpp' verifies:

- scheme parsing;
- orthogonal/corrected equivalence on an orthogonal mesh;
- non-zero corrected response on a skewed face;
- limited correction bounds;
- limiter endpoints \(\lambda=0\) and \(\lambda=1\);
- over-relaxed equivalence on an orthogonal mesh;
- over-relaxed conservation through owner/neighbour antisymmetry;
- explicit rejection of \(\mathbf S_f\cdot\mathbf d=0\).

### Controlled refinement

'tests/validation/test_nonorthogonal_laplacian_campaign.cpp' provides:

1. a linear manufactured field with exact \(\nabla^2\phi=0\);
2. a smooth manufactured field
   \[
   \phi=\sin(\pi x)\cos(\pi y)\sin(\pi z),
   \]
   with
   \[
   \nabla^2\phi=-3\pi^2\phi;
   \]
3. a controlled affine shear ladder from orthogonal to strongly skewed cells;
4. refinement and observed-order measurements;
5. comparison of uncorrected, corrected, limited and over-relaxed operators.

The test is an executable V&V campaign, not merely a smoke test.

### Independent mesh-quality diagnostics

'src/cfdx/core/numerics/diffusion_diagnostics.h' adds a diagnostic layer independent of Laplacian assembly.

For internal faces it reports:

- total internal-face population;
- number of non-orthogonal faces;
- mean/max non-orthogonality angle;
- counts above 60°, 75° and 85°;
- mean/max skewness;
- counts above skewness 0.5 and 1.0;
- count of faces that make the over-relaxed decomposition singular.

The diagnostics do not change the numerical method and therefore cannot make a failing discretisation pass.

'tests/unit/test_diffusion_diagnostics.cpp' verifies both an orthogonal baseline and a genuinely perturbed conforming mesh.

## Mesh-quality definitions

For an internal face,

$$
\theta=
\cos^{-1}
\left(
\frac{\mathbf S_f\cdot\mathbf d}
{|\mathbf S_f||\mathbf d|}
\right)
$$

is the non-orthogonality angle.

Skewness is the normalized distance between the face centroid and its projection onto the owner-neighbour centre line:

$$
\eta=
\frac{
\left|\mathbf x_f-\mathbf x_{proj}\right|
}{
|\mathbf d|
}.
$$

The diagnostic thresholds are reporting thresholds, not hidden numerical switches.

## What is qualified

The N3 acceptance is intentionally tied to #461 rather than to an unlimited claim of arbitrary polyhedral accuracy.

| Capability | Status | Evidence |
|---|---|---|
| Orthogonal diffusion | Qualified | unit + refinement V&V |
| Corrected diffusion | Qualified | unit + skew/refinement V&V |
| Limited correction | Qualified | endpoint/boundedness tests + campaign |
| Over-relaxed diffusion | Qualified | invariants + singular rejection + campaign |
| Controlled non-orthogonality | Qualified | shear ladder |
| Skewness diagnostics | Qualified | independent diagnostic API/tests |
| Conservation of internal-face contribution | Qualified | antisymmetric assembly tests |
| Singular-face handling | Qualified | explicit rejection test |

## What is **not** claimed

The tetrahedral/polyhedral campaign demonstrated that the current two-point-plus-explicit-correction family can remain \(O(1)\) for a smooth Laplacian on strongly non-affine tetrahedral meshes, even when the cell gradient is second order.

That result is retained as evidence. It is **not hidden or reclassified as a pass**.

A higher-order arbitrary-polyhedral diffusion operator requires a different discrete face-flux reconstruction. It remains a research-grade item for the difficult-mesh / advanced numerical work packages (N10/N16). N3 therefore does not claim second-order arbitrary-polyhedral diffusion.

This separation is consistent with #461, whose N3 acceptance asks for expected order on orthogonal meshes and quantified degradation/robustness on controlled skew/non-orthogonal meshes.

## Qualification rule

N3 is green only when:

1. all N3 unit tests pass;
2. the controlled non-orthogonal campaign passes;
3. the diffusion diagnostic regression passes;
4. the final validation CI is green on the exact PR HEAD;
5. the machine-readable maturity audits report 'qualified';
6. no tolerance is relaxed and no validation is disabled.

The polyhedral limitation remains visible in the documentation and is tracked separately rather than being used to invalidate the already-defined N3 acceptance contract.
