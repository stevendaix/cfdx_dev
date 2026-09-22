# Mesh browser model

The GUI mesh browser must consume the existing CFDX HDF5 mesh interchange data rather than a hard-coded patch list. `read_mesh_catalog()` reads the real topology dimensions and boundary patch metadata (`boundary_patches`, `patch_face_ids`, `patch_face_offsets`) emitted by the existing mesh import pipeline.

The browser model exposes stable `MeshSelection("patch", index, name)` handles. It is metadata-only: numerical topology remains owned by the C++ Mesh/Patch structures.
