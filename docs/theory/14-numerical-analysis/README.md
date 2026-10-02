> **Executable source:** [chapter.py](chapter.py)  
> This README is navigation and chapter contract. Quantitative theory belongs to the Python percent source.

# 14 — Numerical Analysis

**Status: REPOSITORY-GROUNDED.**

## Consistency

A discrete operator L_h is consistent with L if

$$
L_hu-Lu\rightarrow0
$$

as h tends to zero for sufficiently smooth u.

## Observed order

If

$$
E_h\approx Ch^p,
$$

then

$$
p_{obs}=\frac{\log(E_h/E_{h/2})}{\log2}.
$$

Observed order belongs to a documented mesh family, solution regularity, boundary treatment and error metric. It is not an unconditional property of a scheme name.

## Stability

For

$$
\frac{du}{dt}=Lu,
$$

time integration generates an amplification operator. Stability depends on its spectrum and step size. Nonlinear CFD also depends on mesh quality, relaxation and physical coefficients.

## Boundedness and monotonicity

Accuracy, stability, boundedness and conservation are distinct properties. A bounded low-order method is not automatically high-order, and a high-order reconstruction is not automatically bounded.

## Error decomposition

A useful engineering decomposition is

$$
e_{total}=e_{discretisation}+e_{iterative}+e_{roundoff}+e_{model}+e_{input}.
$$

Grid convergence is meaningful only when iterative error is controlled or quantified.

## MMS

For an exact field u_e, compute the continuous forcing and solve the modified problem. For weighted L2,

$$
L_2=
\sqrt{\frac{\sum_i|u_i-u_{e,i}|^2w_i}
{\sum_i|u_{e,i}|^2w_i}}.
$$

The weights and population must be reported.

## Richardson and GCI

The refinement ratio, observed order, solution values and uncertainty assumptions must all come from the same controlled family.

CFDX N2 already demonstrates that linear exactness is not proof of second-order convergence on arbitrary polyhedra.

## Evidence

Every convergence claim records mesh family, refinement ratio, norm, order, solver convergence, conservation, boundary treatment, exact software revision and negative results.


## Scientific explanation standard

Every major equation or method in this chapter must explain: physical/mathematical motivation; variables, units and sign conventions; derivation; FVM/discrete formulation; actual CFDX algorithm/code path; analytical or numerical example; errors, limitations and sensitivities; exact source files; executable verification tests; benchmark/V&V evidence; and bibliography with stable identifiers.

Comparison tables are summaries after the mathematics, never substitutes for it. Status must remain **Implemented / Verified / Validated / Qualified** and documentation must never promote numerical maturity.