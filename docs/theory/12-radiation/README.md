> **Executable source:** [chapter.py](chapter.py)  
> This README is navigation and chapter contract. Quantitative theory belongs to the Python percent source.

# 12 — Radiation

**Status: REPOSITORY-GROUNDED; supported model scope must come from executable evidence.**

For a blackbody,

$$
E_b=\sigma T^4.
$$

Radiation is nonlinear in temperature and enters the thermal energy balance through surface or volumetric exchange.

## Surface-to-surface radiation

Diffuse-gray exchange uses view factors and emissivities. Two fundamental geometric checks are

$$
\sum_jF_{ij}=1,
\qquad
A_iF_{ij}=A_jF_{ji},
$$

when the corresponding closed-surface assumptions hold.

## CFDX implementation families

The source tree contains radiation.h, radiation_models.h, radiation_s2s.h, radiation_solver.h, radiation_advanced.h and radiation_streaming.h.

Each family must be mapped to its mathematical formulation before assigning maturity.

## Coupling

Radiative heat flux enters the energy balance as a boundary/source contribution. Radiation and energy modules must share the same sign convention.

## Verification

Use blackbody emission, view-factor closure/reciprocity where implemented, symmetric two-surface exchange, radiative energy conservation and coupled thermal-radiation regression. A radiation unit test is not physical validation.


## Scientific explanation standard

Every major equation or method in this chapter must explain: physical/mathematical motivation; variables, units and sign conventions; derivation; FVM/discrete formulation; actual CFDX algorithm/code path; analytical or numerical example; errors, limitations and sensitivities; exact source files; executable verification tests; benchmark/V&V evidence; and bibliography with stable identifiers.

Comparison tables are summaries after the mathematics, never substitutes for it. Status must remain **Implemented / Verified / Validated / Qualified** and documentation must never promote numerical maturity.