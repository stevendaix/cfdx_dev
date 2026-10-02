# %% [markdown]
"""# CFDX Data Model and File Formats

CFDX separates setup from numerical state: case.cfdx.h5 is the case definition, case.dat.h5 is checkpoint state, and case_<time>.vtu is a visualisation artifact. HDF5 topology uses flattened CSR-style arrays with offsets, owner/neighbour semantics and typed fields. Readers must reject invalid offsets, dimensions and indices rather than silently repairing topology.
"""

# %%
offsets = [0, 3, 5]
assert offsets == sorted(offsets) and offsets[-1] == 5
