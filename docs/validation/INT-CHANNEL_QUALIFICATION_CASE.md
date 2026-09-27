# INT-CHANNEL — plane channel qualification

Reference: fully developed laminar plane-channel analytical solution.

Required evidence: at least three systematically refined meshes; no-slip walls; pressure-gradient or periodic driving; velocity profile, bulk flow and pressure-gradient/friction QoI; independent mass/flux conservation; nonlinear and linear convergence; observed mesh behaviour.

Mesh policy: all 2-D qualification meshes are generated directly by Python and written to CFDX-HDF5. Gmsh is deliberately not used for this case family.

**Current status:** READY in the registry. PR #443 now makes the channel solver consume the Python-generated HDF5 meshes for the N=16/32/64 refinement series. Quantitative PASS still requires the campaign to execute successfully in CI and retain the full convergence/QoI evidence.