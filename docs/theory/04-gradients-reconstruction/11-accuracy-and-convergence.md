# 11. Accuracy and Convergence

For mesh sizes `h_i` and errors `E_i`, the observed order between two levels is

[
p=
rac{log(E_1/E_2)}
{log(h_1/h_2)}.
]

A meaningful refinement study defines the reference solution, mesh family, characteristic size, error norm, refinement ratio, boundary treatment, convergence criteria and numerical options.

Linear-field exactness is useful, but it does not establish asymptotic second-order accuracy for general smooth fields on distorted polyhedral meshes.

Quadratic fields and manufactured solutions are therefore important. The campaign should separate discretisation, geometry, iterative, reconstruction and reference-solution errors where possible.
