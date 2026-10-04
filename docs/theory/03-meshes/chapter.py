# %% [markdown]
# 03 — Meshes and Polyhedral Geometry
#
# A finite-volume solver does not operate directly on CAD geometry. It operates on a discrete geometric/topological graph: cells, faces, vertices, owner/neighbour connectivity, oriented face area vectors, cell volumes and centres. Every numerical operator in later chapters consumes this representation.
#
# ## 3.1 Topology versus geometry
#
# Topology answers which points, faces and cells are connected. Geometry supplies coordinates, centres, areas, normals and volumes. A valid CFD mesh requires both layers to be mutually consistent.
#
# ## 3.2 Mesh entities and orientation
#
# For a face $f$,
# $
# \mathbf S_f=A_f\mathbf n_f,\qquad |\mathbf n_f|=1.
# $
# In CFDX, the cached face vector is oriented owner (\rightarrow) exterior. For an internal face shared by owner $P$ and neighbour $N$,
# $
# \mathbf S_{f,N}=-\mathbf S_{f,P}.
# $
# The face therefore has one geometric object and two algebraic orientations.
#
# The minimum connectivity relation is
# $
# f=(P,N,\{v_1,\ldots,v_m\}),
# $
# with $N$ absent for a boundary face. The owner/neighbour relation must not be reconstructed from face ordering heuristics at solver level.
#
# ## 3.3 Closure theorem
#
# For a closed cell,
# $
# \boxed{\sum_{f\in\partial V_P}\mathbf S_f=\mathbf0}.
# $
# More generally, applying Gauss' theorem to a constant vector field gives
# $
# \int_{V_P}\nabla\mathbf a\,dV
# =\oint_{\partial V_P}\mathbf a\otimes\mathbf n\,dA
# =\mathbf0,
# $
# hence the area-vector closure. This is simultaneously a geometric identity and a prerequisite for constant-field preservation.
#
# ## 3.4 Volume and centroid
#
# For a planar face and reference point (\mathbf x_0), an oriented polyhedral volume can be assembled from pyramidal contributions:
# $
# V_P=\frac13\sum_f
# A_f(\mathbf C_f-\mathbf x_0)\cdot\mathbf n_f.
# $
# The result is independent of (\mathbf x_0) for a closed, consistently oriented polyhedron.
#
# The first geometric moment can be written in terms of a consistent face decomposition. In implementation, the exact convention matters more than a generic centroid formula: face triangulation, orientation and signed sub-volumes must use the same convention as the volume routine.
#
# ## 3.5 Face area and polygon decomposition
#
# For a planar polygon triangulated around a reference vertex (\mathbf x_0),
# $
# \mathbf S_f=\frac12\sum_i
# (\mathbf x_i-\mathbf x_0)\times
# (\mathbf x_{i+1}-\mathbf x_0).
# $
# The polygon area is
# $
# A_f=|\mathbf S_f|.
# $
# This expression makes the orientation dependence explicit: reversing the vertex order changes the sign of (\mathbf S_f), not the physical area.
#
# ## 3.6 Centres and centre-to-centre vectors
#
# Define
# $
# \mathbf d_{PN}=\mathbf C_N-\mathbf C_P,\qquad
# d_{PN}=|\mathbf d_{PN}|.
# $
# A centre-to-face interpolation parameter is
# $
# w=\frac{|\mathbf C_f-\mathbf C_P|}{|\mathbf C_N-\mathbf C_P|}.
# $
# This is a geometric interpolation coordinate only when the face-centre location and reconstruction policy make that interpretation valid.
#
# ## 3.7 Non-orthogonality, skewness and aspect ratio
#
# A representative non-orthogonality angle is
# $
# \theta_f=
# \cos^{-1}\!\left(
# \frac{\mathbf d_{PN}\cdot\mathbf n_f}{|\mathbf d_{PN}|}\right).
# $
# A more direct decomposition used by diffusion schemes is
# $
# \mathbf S_f=\alpha_f\mathbf d_{PN}
# +\mathbf S_f^\perp,\qquad
# \alpha_f=\frac{\mathbf S_f\cdot\mathbf d_{PN}}
# {|\mathbf d_{PN}|^2}.
# $
# Then (\mathbf S_f^\perp\cdot\mathbf d_{PN}=0). This separates the two-point orthogonal contribution from the non-orthogonal correction.
#
# Skewness measures displacement between the interpolation location implied by the stencil and the actual face integration point. Aspect ratio measures directional stretching. They are distinct error mechanisms and should be reported separately.
#
# ## 3.8 Non-planar faces
#
# A non-coplanar polygon has no unique planar normal without a convention. CFDX must therefore define how area, centre and orientation are computed. A single physical internal face must have the same geometry for both adjacent cells, with only the orientation sign changing.
#
# ## 3.9 Degeneracy and quality gates
#
# Geometry validation should reject, or explicitly classify,
# $
# V_P\le V_{\min},\quad
# A_f\le A_{\min},\quad
# d_{PN}\le d_{\min},
# $
# non-finite coordinates, invalid indices, duplicate topology and inconsistent owner/neighbour orientation.
#
# Quality thresholds are acceptance criteria, not correction mechanisms. A threshold must never silently turn a negative or degenerate volume into a positive one.
#
# ## 3.10 Parallel partitioning
#
# Domain decomposition duplicates interface information as required by the local stencil. A global conservation audit must account for the partition interface exactly once:
# $
# \sum_P R_P=0
# $
# for a closed conserved system, including all owned-cell contributions and correctly paired inter-partition fluxes.
#
# ## 3.11 Geometry-to-discretisation chain
#
# $
# \boxed{
# \text{vertices}
# \rightarrow\text{faces/topology}
# \rightarrow
# (\mathbf S_f,A_f,\mathbf C_f)
# \rightarrow
# (V_P,\mathbf C_P)
# \rightarrow
# \text{gradients}
# \rightarrow
# \text{fluxes}
# \rightarrow
# A
# }
# $
# A geometry error can therefore appear later as a gradient, diffusion or conservation error.
#
# ## 3.12 CFDX implementation
#
# The relevant source family is [mesh](../../../src/cfdx/core/mesh/), [face_geometry.h](../../../src/cfdx/core/geometry/face_geometry.h), [cell_geometry.h](../../../src/cfdx/core/geometry/cell_geometry.h), [geometry_cache.h](../../../src/cfdx/core/geometry/geometry_cache.h), [mesh_quality.h](../../../src/cfdx/core/geometry/mesh_quality.h) and [mesh_validator.h](../../../src/cfdx/core/geometry/mesh_validator.h).
#
# ## 3.13 Verification ladder
#
# $
# \text{indices}
# \rightarrow\text{orientation}
# \rightarrow\text{closure}
# \rightarrow\text{positive volume}
# \rightarrow\text{known geometry}
# \rightarrow\text{quality metrics}
# \rightarrow\text{operator invariants}.
# $
#
# A mesh-quality number is not itself evidence that a numerical operator is accurate.
#
# ## 3.14 Executable closure test
# %%
import numpy as np
S = np.array([[1.,0.,0.],[-1.,0.,0.],[0.,1.,0.],[0.,-1.,0.],[0.,0.,1.],[0.,0.,-1.]])
assert np.allclose(S.sum(axis=0), 0.0)
