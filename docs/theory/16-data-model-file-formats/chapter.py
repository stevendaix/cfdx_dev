# %% [markdown]
# 16 — Data Model and File Formats
#
# ## 16.1 Three semantic products
#
# CFDX distinguishes authoritative setup, restart state and visualization:
# \[
# \boxed{\text{case.cfdx.h5}=\text{setup}},\qquad
# \boxed{\text{case.dat.h5}=\text{state/checkpoint}},
# \qquad
# \boxed{\text{case}_{t}.vtu=\text{visualization}}.
# \]
# A visualization file must not silently become the authoritative case definition.
#
# ## 16.2 Topology representation
#
# The v1 HDF5 case schema contains root topology datasets:
# \[
# points,\quad face\_vertices,\quad face\_offsets,\quad owner,\quad
# neighbour,\quad cell\_faces,\quad cell\_offsets.
# \]
# The face/cell arrays form CSR-like ragged representations:
# \[
# face\_vertices[face\_offsets_f:face\_offsets_{f+1}]
# \]
# lists the vertices of face \(f\), while
# \[
# cell\_faces[cell\_offsets_c:cell\_offsets_{c+1}]
# \]
# lists the faces attached to cell \(c\).
#
# ## 16.3 Topology invariants
#
# For every face:
# \[
# 0\le owner_f<N_{cell},
# \]
# and for internal faces
# \[
# 0\le neighbour_f<N_{cell},\qquad owner_f\ne neighbour_f.
# \]
# Boundary faces use the documented sentinel/optional boundary-patch representation. Offset arrays must be monotone and terminate at the corresponding flattened-array length.
#
# ## 16.4 Fields
#
# Cell scalar and vector fields are stored under
# \[
# /fields/scalar/<name>,\qquad
# /fields/vector/<name>.
# \]
# A field contract must specify location, component count, units, association to a mesh and semantic name.
#
# ## 16.5 Metadata and compatibility
#
# Root attributes encode the implemented schema/version and case metadata. Readers must distinguish unsupported schema versions from absent optional data. Compatibility must never be implemented by silently dropping physics/numerics data.
#
# ## 16.6 Checkpoint state
#
# A restart state is not merely the latest field snapshot. If the selected temporal scheme requires history:
# \[
# U^n,\ U^{n-1},\ldots
# \]
# and if adaptive controllers require state, their controller variables are part of restart completeness.
#
# ## 16.7 Restart invariant
#
# Let \(R\) be a restart read and \(W\) a checkpoint write. For supported data:
# \[
# R(W(x))\equiv x
# \]
# up to explicitly documented floating-point/serialization differences.
#
# More strongly, a restarted calculation should reproduce the same next-step discrete state under deterministic execution:
# \[
# \mathcal S(R(W(x)))\approx\mathcal S(x).
# \]
#
# ## 16.8 VTU
#
# VTU is a presentation/output format. Cell/point topology and fields must map unambiguously from CFDX's internal representation. Derived visualisation quantities should be marked as derived rather than mistaken for authoritative solver state.
#
# ## 16.9 Provenance
#
# A result should be traceable to
# \[
# (case\ revision,\ mesh\ revision,\ physics,\ numerics,\ solver,\ software\ revision).
# \]
# This is essential for V&V reproducibility.
#
# ## 16.10 Implementation traceability
#
# HDF5 implementation is under [src/cfdx/io/hdf5](../../../src/cfdx/io/hdf5/); case schema under [src/cfdx/io/cfdx_io](../../../src/cfdx/io/cfdx_io/); VTU output under [src/cfdx/io/vtu](../../../src/cfdx/io/vtu/). Exact dataset definitions must remain synchronized with the reader/writer and schema tests.
#
# ## 16.11 Verification
#
# Required checks: schema version handling, dimensions/dtypes, offset consistency, topology round-trip, field round-trip, checkpoint round-trip, restart continuation and VTU field/topology integrity.
#
# ## 16.12 Executable round-trip abstraction
# %%
import numpy as np
x=np.arange(12,dtype=float).reshape(3,4)
serialized=x.copy()
restored=serialized.copy()
assert np.array_equal(restored,x)
