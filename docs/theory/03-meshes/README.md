# 03 — Meshes

A CFD mesh is a discrete geometric model, not merely a list of coordinates. CFDX needs topology, orientation, geometry, boundary ownership and quality information before numerical operators are meaningful.

## 1. Topology

The basic entities are points \(P_i\), oriented faces \(f\), and cells \(P\). A face stores its vertex cycle and owner/neighbour relation. Boundary faces have an owner and no internal neighbour.

**CFDX paths:** `src/cfdx/core/mesh/mesh.h`, `src/cfdx/core/mesh/mesh.cpp`, `src/cfdx/core/mesh/boundary.cpp`.

## 2. Face orientation

For a face with area vector
\[
\mathbf S_f=A_f\mathbf n_f,
\]
the normal must be outward from the owner cell. For an internal neighbour,
\[
\mathbf S_{f,N}=-\mathbf S_{f,P}.
\]
A sign error here changes convection, diffusion and pressure forces simultaneously.

## 3. Cell volume

For a closed polyhedron, the volume can be evaluated through oriented tetrahedral decomposition or a divergence identity. A convenient identity is
\[
V=\frac13\oint_{\partial V}\mathbf x\cdot\mathbf n\,dA
\approx\frac13\sum_f\mathbf C_f\cdot\mathbf S_f.
\]
The implementation must use a geometrically consistent face-centre convention.

## 4. Geometric invariants

For every valid closed cell:
\[
V>0,\qquad\sum_f\mathbf S_f=0.
\]
For an internal face, owner and neighbour centres must be distinct unless a degenerate mesh is deliberately supported.

## 5. Quality

Useful measures include
\[
\theta_{nonorth}=\cos^{-1}
\left(
\frac{\mathbf d_{PN}\cdot\mathbf S_f}
{|\mathbf d_{PN}||\mathbf S_f|}
\right),
\]
skewness, aspect ratio, volume ratio and face-area imbalance. A quality metric is diagnostic, not a universal acceptance threshold: each numerical method has its own sensitivity.

## 6. Non-planar faces

A polyhedral face may not be exactly planar. Area-vector construction and face-centre calculation must use a documented convention. Numerical operators must not silently assume a planar polygon if the mesh importer permits general polyhedra.

## 7. Boundary patches

A patch groups boundary faces sharing a physical condition. Patch identity is therefore part of the physics model.

**CFDX paths:** `src/cfdx/core/mesh/boundary.cpp`, `src/cfdx/core/boundary/boundary_condition.h`.

## 8. Mesh import

The import chain is
\[
\text{external mesh}\rightarrow\text{topology mapping}
\rightarrow\text{orientation/geometry reconstruction}
\rightarrow\text{validation}\rightarrow\text{solver mesh}.
\]
Relevant importers include `src/cfdx/io/gmsh/gmsh_importer.h` and `src/cfdx/io/mesh/mesh_importer.h`.

## 9. Mesh rejection

A mesh should be rejected when mandatory invariants fail: invalid indices, open cells where closed cells are required, inconsistent owner/neighbour references, non-positive volume or irreconcilable face orientation. Silent repair is dangerous because it can change the physical problem.

## 10. Verification

A mesh verification report should contain entity counts, boundary counts, volume range, area range, closure residuals, orientation diagnostics, quality distributions and import provenance.

**Tests / tools:** `tests/unit/`, `apps/cfdx_convert/main.cpp`, `src/cfdx/runtime/dynamic_mesh.h`.
