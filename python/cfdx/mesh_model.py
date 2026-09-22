"""Headless mesh catalog for the GUI/TUI application layer.

The reader consumes the existing CFDX HDF5 mesh interchange layout. It does
not duplicate Mesh topology; it exposes only metadata required by a browser.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py

from .setup_model import MeshSelection

_PATCH_TYPES = {0:"wall",1:"inlet",2:"outlet",3:"symmetry",4:"periodic",5:"interface",6:"empty",7:"unknown"}

@dataclass(frozen=True)
class MeshPatch:
    index: int
    name: str
    type: str
    face_count: int

    @property
    def selection(self) -> MeshSelection:
        return MeshSelection("patch", self.index, self.name)

@dataclass(frozen=True)
class MeshCatalog:
    n_points: int
    n_faces: int
    n_cells: int
    patches: tuple[MeshPatch, ...]

    @property
    def patch_names(self) -> tuple[str, ...]:
        return tuple(p.name for p in self.patches)

    def patch(self, index: int) -> MeshPatch:
        try:
            return self.patches[index]
        except IndexError as exc:
            raise IndexError(f"mesh patch index out of range: {index}") from exc

def _attr_int(h5: h5py.File, name: str, default: int = 0) -> int:
    value = h5.attrs.get(name, default)
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    return int(value)

def _parse_patch_metadata(raw: object) -> list[tuple[str, int]]:
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8")
    if not raw:
        return []
    entries: list[tuple[str, int]] = []
    for item in str(raw).split(";"):
        if not item:
            continue
        fields = item.split(":")
        if len(fields) < 4:
            raise ValueError(f"invalid boundary patch metadata: {item!r}")
        name = fields[0].strip()
        if not name:
            raise ValueError("boundary patch name is empty")
        entries.append((name, int(fields[3])))
    return entries

def read_mesh_catalog(path: Path) -> MeshCatalog:
    """Read mesh dimensions and real patch metadata from a CFDX HDF5 mesh."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    with h5py.File(path, "r") as h5:
        if "points" not in h5 or "owner" not in h5 or "neighbour" not in h5:
            raise ValueError("CFDX mesh is missing mandatory topology datasets")
        n_points = len(h5["points"])
        n_faces = len(h5["owner"])
        n_cells = _attr_int(h5, "n_cells", len(h5["cell_offsets"]) - 1 if "cell_offsets" in h5 else 0)
        if "patch_face_offsets" not in h5 or "patch_face_ids" not in h5:
            raise ValueError("CFDX mesh is missing boundary patch datasets")
        offsets = h5["patch_face_offsets"][:]
        face_ids = h5["patch_face_ids"][:]
        metadata = _parse_patch_metadata(h5.attrs.get("boundary_patches", ""))
        if len(offsets) != len(metadata) + 1:
            raise ValueError("boundary patch metadata/offset count mismatch")
        if len(offsets) == 0 or int(offsets[0]) != 0:
            raise ValueError("boundary patch offsets must start at zero")
        if any(int(offsets[i + 1]) < int(offsets[i]) for i in range(len(offsets) - 1)):
            raise ValueError("boundary patch offsets must be monotonic")
        if any(int(fid) < 0 or int(fid) >= n_faces for fid in face_ids):
            raise ValueError("boundary patch contains an out-of-range face id")
        patches=[]
        for i,(name,type_code) in enumerate(metadata):
            start,end=int(offsets[i]),int(offsets[i+1])
            if start < 0 or end < start or end > len(face_ids):
                raise ValueError(f"invalid face range for patch {name!r}")
            patches.append(MeshPatch(i,name,_PATCH_TYPES.get(type_code,"unknown"),end-start))
        return MeshCatalog(n_points,n_faces,n_cells,tuple(patches))
