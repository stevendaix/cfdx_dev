# SEP-BFS — backward-facing-step qualification

Reference: Armaly et al. (1983), JFM 127, 473–496.

Required evidence: canonical 2-D geometry, frozen Reynolds number, inlet profile, no-slip walls, pressure outlet, controlled wall/shear-layer refinement, reattachment length, pressure recovery and wall shear, conservation and at least three mesh levels.

**Current status:** READY for implementation. The conformal 2-D mesh is generated directly in Python and written to CFDX-HDF5; Gmsh is not used. A Re=200 Armaly/Fluent-aligned target is frozen for the first solver campaign. The current PR still needs the real separated-flow solve, wall-shear extraction and reattachment-length gate.