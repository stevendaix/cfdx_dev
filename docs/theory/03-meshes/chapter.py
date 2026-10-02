# %% [markdown]
# Meshes
#
# ## 3.1 Topology and orientation
#
# A finite-volume mesh is a topological graph of points, faces and cells. Every face has an orientation. For cell P,
# \[
# \mathbf S_f=A_f\mathbf n_f
# \]
# is the outward area vector.
#
# ## 3.2 Geometric closure
#
# For every closed cell,
# \[
# \sum_{f\in\partial V_P}\mathbf S_f=\mathbf0.
# \]
# This identity is required by the divergence theorem and is an immediate diagnostic for corrupted geometry.
#
# ## 3.3 Cell volume
#
# A polyhedral volume can be obtained by oriented tetrahedral decomposition. For a planar face and reference point \(\mathbf x_0\),
# \[
# V_P=\frac13\sum_fA_f(\mathbf x_f-\mathbf x_0)\cdot\mathbf n_f,
# \]
# with the exact sign convention defined by the implementation.
#
# ## 3.4 Neighbour geometry
#
# \[
# \mathbf d_{PN}=\mathbf C_N-\mathbf C_P,\qquad
# d_{PN}=|\mathbf d_{PN}|.
# \]
# A representative non-orthogonality angle is
# \[
# \theta_f=\cos^{-1}
# \left(\frac{\mathbf d_{PN}\cdot\mathbf n_f}{|\mathbf d_{PN}|}\right).
# \]
# Skewness measures the displacement of the face centre from the centre-to-centre interpolation line. Aspect ratio measures anisotropy. These metrics describe different failure modes.
#
# ## 3.5 Polyhedral cells
#
# Polyhedra require consistent vertex ordering, face orientation, face-centre and area-vector definitions. A non-planar polygon requires a deterministic geometric convention; both adjacent cells must use the same physical face with opposite orientation.
#
# ## 3.6 Mesh validity
#
# Before assembly, validate indices, ownership, face orientation, positive volume, finite coordinates and forbidden degeneracies. A mesh failure must be reported as a geometry error rather than allowed to become a later linear-solver failure.
#
# ## 3.7 Numerical consequences
#
# Mesh distortion changes truncation error, gradient conditioning, diffusion accuracy and solver convergence. Consequently, an accuracy claim must state the mesh family used.
#
# ## 3.8 CFDX traceability
#
# Mesh topology: src/cfdx/core/mesh/  
# Face geometry: src/cfdx/core/geometry/face_geometry.h  
# Cell geometry: src/cfdx/core/geometry/cell_geometry.h  
# Quality metrics: src/cfdx/core/geometry/mesh_quality.h  
# Validation: src/cfdx/core/geometry/mesh_validator.h
#
# %%
from __future__ import annotations
import numpy as np
S=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],float)
assert np.allclose(S.sum(axis=0),0.0)
