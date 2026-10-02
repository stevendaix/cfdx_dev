# 15 — Verification and Validation

## 1. Verification versus validation
Verification asks whether CFDX solves the intended mathematical equations. Validation asks whether those equations and their implementation represent the physical system for the intended use.

## 2. Code verification
Unit/algebraic tests establish local properties. For manufactured solution \(u_M\),
\[
S_M=L(u_M),
\]
and the discrete consistency residual is
\[
R_h=L_h(u_M)-S_M.
\]
Refinement of \(R_h\) is stronger evidence than simply obtaining a small value on one mesh.

## 3. Solution verification
With reference solution \(u_{ref}\),
\[
E_h=\|u_h-u_{ref}\|.
\]
For refinement ratio \(r\),
\[
p=\frac{\ln(E_h/E_{h/r})}{\ln r}.
\]
The reference must itself be sufficiently accurate for the claimed quantity.

## 4. Richardson extrapolation
If
\[
u_h=u_0+Ch^p,
\]
then
\[
u_{ext}=u_h+\frac{u_h-u_{h/r}}{r^p-1}.
\]
The asymptotic assumption must be checked rather than presumed.

## 5. GCI
A common form is
\[
GCI=F_s\frac{|u_1-u_2|}{|u_1|}\frac{1}{r^p-1}.
\]
The safety factor and indexing convention must be fixed in the V&V procedure.

## 6. Iterative error
If \(u_h^*\) is the exact discrete solution,
\[
e_i=u_h-u_h^*.
\]
A small residual bounds algebraic error only under assumptions involving the operator conditioning; it does not prove discretisation accuracy.

## 7. Validation
For experiment \(D_{exp}\),
\[
E_v=\frac{D_{CFDX}-D_{exp}}{D_{exp}}.
\]
This difference must be interpreted with experimental and numerical uncertainties and with identical operating conditions.

## 8. Qualification
Qualification is a defined evidence decision for a specified capability. It is not synonymous with a test passing.

## 9. Evidence
Each result must retain equations, mesh, BCs, solver settings, tolerances, residuals, reference data, uncertainty, commit and machine-readable artifact.

**CFDX governance:** `docs/vv/`, `docs/validation/CFDX_VV_GOVERNANCE.md`. **References:** `docs/references/bibliography.bib` (Roache, ASME V&V 20, Ghia et al.).