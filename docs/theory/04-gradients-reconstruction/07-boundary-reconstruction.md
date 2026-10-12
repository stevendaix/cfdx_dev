# 7. Boundary Reconstruction

A boundary face normally has no second cell-centre neighbour. Interior stencils therefore cannot simply be copied to boundaries.

A boundary policy must specify:

- available geometric information;
- prescribed versus derived boundary values;
- normal and tangential reconstruction;
- any ghost information;
- corner and curved-boundary treatment;
- behaviour for incomplete local stencils.

Boundary handling is part of the numerical method. For CFDX it must remain explicit so that verification can reproduce the stencil and boundary-condition implementations cannot silently alter gradient behaviour.
