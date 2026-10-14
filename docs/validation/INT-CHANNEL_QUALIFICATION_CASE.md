# INT-CHANNEL — plane channel qualification

Reference: fully developed laminar plane-channel analytical solution.

Required evidence: at least three systematically refined meshes; no-slip walls; pressure-gradient or periodic driving; velocity profile, bulk flow and pressure-gradient/friction QoI; independent mass/flux conservation; nonlinear and linear convergence; observed mesh behaviour.

Mesh policy: all 2-D qualification meshes are generated directly by Python and written to CFDX-HDF5. Gmsh is deliberately not used for this case family.

**Boundary formulation (reformulated after diagnostic #894):** body-force driven plane channel (G=1, nu=0.1) with no-slip walls, zero-gradient inlet velocity and fixed outlet pressure p=0. The former all-Neumann pressure setup was ill posed: the pressure slope was undetermined and the pure-Neumann path eliminated the continuity row of the reference (inlet) cell, so the global mass imbalance was invisible to the pressure correction (mirror of the #889 root cause; the diagnostic measured umax=0.0849 against the analytical umax=1.25). The fixed outlet pressure removes the null mode of the pressure equation (fiche methode 15, section 4): the pressure operator is SPD, no reference row is eliminated, and the outlet boundary flux is driven by (p_bc - p_P). The developed body-force solution u = G/(2 nu) y(1-y) with uniform p = 0 satisfies every boundary condition exactly; no gate or tolerance was modified.

**Current status:** DIAGNOSTIC in the registry. The executable case consumes the Python-generated HDF5 meshes for the N=16/32/64 refinement series (PR #443) and is now well posed (fixed outlet pressure, gauge determined). Quantitative PASS still requires the validation-total campaign to execute successfully at a merged HEAD and retain the full convergence/QoI evidence.