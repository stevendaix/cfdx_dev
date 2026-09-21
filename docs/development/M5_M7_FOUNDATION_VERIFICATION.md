# M5–M7 verification status

## M5 — VOF

The first CFDX reference kernel is a conservative periodic 1-D volume-fraction transport update.

Acceptance criteria:
- 0 <= alpha <= 1;
- CFL <= 1 is enforced;
- periodic integral of alpha is conserved to floating-point tolerance;
- surface-tension force follows the CSF scalar relation sigma*kappa*|grad(alpha)|;
- contact-angle wall-normal component is checked against cos(theta).

This is a verification kernel, not yet the full polyhedral PLIC reconstruction. Full multi-dimensional interface reconstruction remains a separate solver-level work package.

## M6 — Dynamic mesh

The foundation applies point displacements to the native mesh and recomputes cell volumes from the moved topology.

Acceptance criteria:
- displacement vector has one entry per mesh point;
- non-finite motion is rejected;
- all resulting cell volumes remain positive;
- minimum volume ratio is reported and gated;
- topology remapping has explicit old-to-new/new-to-old validation.

A full remeshing algorithm and conservative field transfer are still separate implementation work.

## M7 — FSI

The foundation provides backend-neutral partitioned coupling utilities:
- interface work;
- displacement residual;
- bounded Aitken relaxation;
- displacement/force convergence gate.

The utilities deliberately do not embed a structural discretization. Fluid and structure solvers remain independent, with coupling at the interface as required by the architecture.

## Evidence level

These tests are Level-A kernel verification. They must not be interpreted as full application validation of VOF, dynamic-mesh or FSI capabilities.
