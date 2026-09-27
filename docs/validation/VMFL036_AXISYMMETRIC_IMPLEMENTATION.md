# VMFL036 native axisymmetric implementation

This directory contains the native meridional kernel required to implement the
official VMFL036 definition as a 2-D axisymmetric calculation.

Implemented foundations:
1. explicit (x,r) geometry and revolution metrics;
2. deterministic meridional mesh representation;
3. conservative 2*pi*r face/volume metrics;
4. cylindrical no-swirl momentum terms;
5. axis regularity guards;
6. pressure-correction coefficient using the axisymmetric face measure;
7. explicit axis/wall/inlet/outlet/farfield BC vocabulary;
8. independent pressure/viscous sphere traction;
9. kernel verification including SU2-inspired axisymmetric operator checks;
10. explicit VMFL036 Re=100 and Ansys-property variants;
11. CI/CTest registration and roadmap traceability.

The existing generic 3-D solver is not silently reused as an axisymmetric solver.
The next integration step is to connect this metric/kernel to the production
finite-volume assembly and its pressure-velocity loop. A planar 2-D solve,
thin 3-D extrusion, or post-hoc 2*pi*r scaling is not accepted as VMFL036.
