"""Headless mesh catalog for the GUI/TUI application layer.

The reader consumes the CFDX HDF5 mesh interchange layout.  It exposes only
metadata required by the browser and keeps numerical topology owned by C++.
Both legacy single-region files and region-grouped HDF5 cases are supported.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py

from .setup_model import MeshSelection

_PATCH_TYPES = {
    0: "wall", 1: "inlet", 2: "outlet", 3: "symmetry",
    4: "periodic", 5: "interface", 6: "empty", 7: "unknown",
}


@dataclass(frozen=True)
class MeshCell:
    index: int
    stable_id: str

    @property
    def selection(self) -> MeshSelection:
        return MeshSelection("cell", self.index, None, self.stable_id)


@dataclass(frozen=True)
class MeshPatch:
    index: int
    name: str
    type: str
    face_count: int
    region: str = "default"

    @property
    def stable_id(self) -> str:
        prefix = "" if self.region == "default" else f"region:{self.region}/"
        return f"{prefix}patch:{self.name}"

    @property
    def selection(self) -> MeshSelection:
        return MeshSelection("patch", self.index, self.name, self.stable_id)


@dataclass(frozen=True)
class MeshRegion:
    name: str
    n_points: int
    n_faces: int
    n_cells: int
    patches: tuple[MeshPatch, ...]
    cell_ids: tuple[str, ...]

    @property
    def stable_id(self) -> str:
        return f"region:{self.name}"

    def cell(self, index: int) -> MeshCell:
        return MeshCell(index, self.cell_ids[index])


@dataclass(frozen=True)
class MeshCatalog:
    n_points: int
    n_faces: int
    n_cells: int
    patches: tuple[MeshPatch, ...]
    regions: tuple[MeshRegion, ...] = ()

    @property
    def patch_names(self) -> tuple[str, ...]:
        return tuple(p.name for p in self.patches)

    @property
    def cell_ids(self) -> tuple[str, ...]:
        if self.regions:
            if len(self.regions) == 1:
                return self.regions[0].cell_ids
            return tuple(cell_id for region in self.regions for cell_id in region.cell_ids)
        return tuple(f"cell:{i}" for i in range(self.n_cells))

    def cell(self, index: int) -> MeshCell:
        return MeshCell(index, self.cell_ids[index])

    def patch(self, index: int) -> MeshPatch:
        try:
            return self.patches[index]
        except IndexError as exc:
            raise IndexError(f"mesh patch index out of range: {index}") from exc


def _attr_int(h5: h5py.Group, name: str, default: int = 0) -> int:
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


def _region_groups(h5: h5py.File) -> list[tuple[str, h5py.Group]]:
    if "regions" in h5 and isinstance(h5["regions"], h5py.Group):
        return [
            (name, group)
            for name, group in h5["regions"].items()
            if isinstance(group, h5py.Group)
        ]
    if "points" in h5:
        return [("default", h5)]
    return []


def _read_region(name: str, h5: h5py.Group) -> MeshRegion:
    if "points" not in h5 or "owner" not in h5 or "neighbour" not in h5:
        raise ValueError(f"region {name!r} is missing mandatory topology datasets")
    n_points = len(h5["points"])
    n_faces = len(h5["owner"])
    n_cells = _attr_int(
        h5, "n_cells", len(h5["cell_offsets"]) - 1 if "cell_offsets" in h5 else 0
    )
    if n_cells < 0:
        raise ValueError(f"region {name!r} has invalid cell count")

    if "cell_ids" in h5:
        raw_ids = h5["cell_ids"][:]
        if len(raw_ids) != n_cells:
            raise ValueError(f"region {name!r} cell_ids length mismatch")
        ids: list[str] = []
        for value in raw_ids:
            if isinstance(value, bytes):
                value = value.decode("utf-8")
            value = str(value)
            if not value:
                raise ValueError(f"region {name!r} contains an empty cell ID")
            ids.append(f"region:{name}/cell:{value}" if name != "default" else f"cell:{value}")
        if len(ids) != len(set(ids)):
            raise ValueError(f"region {name!r} cell IDs must be unique")
    else:
        ids = [
            f"region:{name}/cell:{i}" if name != "default" else f"cell:{i}"
            for i in range(n_cells)
        ]

    if "patch_face_offsets" not in h5 or "patch_face_ids" not in h5:
        raise ValueError(f"region {name!r} is missing boundary patch datasets")
    offsets = h5["patch_face_offsets"][:]
    face_ids = h5["patch_face_ids"][:]
    metadata = _parse_patch_metadata(h5.attrs.get("boundary_patches", ""))
    if len(offsets) != len(metadata) + 1:
        raise ValueError(f"region {name!r} boundary patch metadata/offset count mismatch")
    if len(offsets) == 0 or int(offsets[0]) != 0:
        raise ValueError(f"region {name!r} patch offsets must start at zero")
    if any(int(offsets[i + 1]) < int(offsets[i]) for i in range(len(offsets) - 1)):
        raise ValueError(f"region {name!r} patch offsets must be monotonic")
    if any(int(fid) < 0 or int(fid) >= n_faces for fid in face_ids):
        raise ValueError(f"region {name!r} contains an out-of-range patch face ID")

    names = [item[0] for item in metadata]
    if len(names) != len(set(names)):
        raise ValueError(f"region {name!r} boundary patch names must be unique")

    patches = []
    for i, (patch_name, type_code) in enumerate(metadata):
        start, end = int(offsets[i]), int(offsets[i + 1])
        if start < 0 or end < start or end > len(face_ids):
            raise ValueError(f"invalid face range for patch {patch_name!r}")
        patches.append(
            MeshPatch(i, patch_name, _PATCH_TYPES.get(type_code, "unknown"),
                      end - start, name)
        )
    return MeshRegion(name, n_points, n_faces, n_cells, tuple(patches), tuple(ids))


def read_mesh_catalog(path: Path) -> MeshCatalog:
    """Read a real CFDX HDF5 mesh, including a multi-region catalog."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    with h5py.File(path, "r") as h5:
        groups = _region_groups(h5)
        if not groups:
            raise ValueError("CFDX mesh contains no readable region")
        regions = tuple(_read_region(name, group) for name, group in groups)
        if len(regions) == 1 and regions[0].name == "default":
            region = regions[0]
            return MeshCatalog(region.n_points, region.n_faces, region.n_cells,
                               region.patches, regions)
        patches = tuple(patch for region in regions for patch in region.patches)
        return MeshCatalog(
            sum(region.n_points for region in regions),
            sum(region.n_faces for region in regions),
            sum(region.n_cells for region in regions),
            patches,
            regions,
        )
