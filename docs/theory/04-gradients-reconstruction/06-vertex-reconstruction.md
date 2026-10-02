# 6. Vertex-Based Reconstruction

Vertex-based methods reconstruct or approximate field information at mesh vertices and use that information to recover a cell gradient.

They can exploit richer polyhedral connectivity, but introduce additional choices:

- vertex-to-cell weighting;
- vertex field reconstruction;
- boundary vertices;
- unusual topology;
- hanging or irregular connectivity.

A reproducible method definition must specify these choices explicitly. The label "vertex Green–Gauss" alone is not sufficient to reproduce an algorithm.

Verification must also distinguish the properties of vertex reconstruction from properties inherited from the test field or mesh geometry.
