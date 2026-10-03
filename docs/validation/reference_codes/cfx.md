# Ansys CFX reference dossier

## Evidence scope

CFX public documentation provides unusually detailed mathematical descriptions of its finite-volume/control-volume discretisation and advection schemes.

## Relevant documented capabilities

CFX documents:
- first-order Upwind advection;
- specified blend between first and second-order advection;
- Central Difference;
- bounded High Resolution advection;
- first- and second-order backward Euler transient schemes;
- control-volume gradients evaluated through a Gauss-divergence-theorem form at integration points.

The CFX High Resolution scheme uses a nonlinear boundedness construction and is therefore not assumed equivalent to any particular CFDX TVD limiter.

## Sources

- https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/cfx_mod/cfx_mod.html
- https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/cfx_thry/i1311648.html
- https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/cfx_mod/i1313685.html
- https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/cfx_thry/CIHEGJAA.html
- https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/pdf/Ansys_CFX-Solver_Theory_Guide.pdf

## N14 interpretation

CFX is a reference for conservative finite-volume/control-volume formulation, high-resolution bounded advection and transient discretisation. The matrix must preserve the distinction between common purpose and identical algebra.
