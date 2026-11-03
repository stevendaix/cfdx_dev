"""Read-only, metadata-first inspection of CFDX HDF5 artifacts.

The inspector deliberately reads HDF5 metadata only. It never materializes mesh
coordinates, connectivity, cell IDs, or field arrays, so opening a large case
does not imply loading its scientific payload into RAM.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import prod
from pathlib import Path
from typing import Any
import json

import h5py

_FORMATS = {"CFDX", "CFDX-HDF5-mesh-v1"}


def _text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    shape = getattr(value, "shape", ())
    if shape:
        return str(value.tolist() if hasattr(value, "tolist") else value)
    if hasattr(value, "item"):
        value = value.item()
    return str(value)


def _shape_text(shape: tuple[int, ...]) -> str:
    return " × ".join(str(int(size)) for size in shape) or "scalar"


@dataclass(frozen=True)
class Hdf5DatasetInfo:
    path: str
    shape: tuple[int, ...]
    dtype: str
    size_bytes: int
    chunks: tuple[int, ...] | None
    compression: str | None

    @property
    def display(self) -> str:
        return f"{self.path}  [{_shape_text(self.shape)} · {self.dtype}]"


@dataclass(frozen=True)
class CaseFileSummary:
    path: Path
    format: str | None
    schema_version: str | None
    is_cfdx: bool
    has_case_config: bool
    has_project_config: bool
    has_mesh: bool
    preview_supported: bool
    n_points: int | None
    n_faces: int | None
    n_cells: int | None
    regions: tuple[str, ...]
    scalar_fields: tuple[str, ...]
    vector_fields: tuple[str, ...]
    case_name: str | None
    physics: tuple[str, ...]
    numerics: tuple[str, ...]
    boundaries: tuple[str, ...]
    attributes: tuple[tuple[str, str], ...]
    datasets: tuple[Hdf5DatasetInfo, ...]
    datasets_truncated: bool = False

    @property
    def total_field_count(self) -> int:
        return len(self.scalar_fields) + len(self.vector_fields)


def _first_dataset(h5: h5py.File, candidates: tuple[str, ...]):
    for candidate in candidates:
        if candidate in h5 and isinstance(h5[candidate], h5py.Dataset):
            return h5[candidate]
    return None


def _region_names(h5: h5py.File) -> tuple[str, ...]:
    if "regions" not in h5 or not isinstance(h5["regions"], h5py.Group):
        return ()
    return tuple(sorted(
        name for name, value in h5["regions"].items()
        if isinstance(value, h5py.Group)
    ))


def _mesh_counts(h5: h5py.File, regions: tuple[str, ...]):
    if regions:
        counts = []
        for name in regions:
            group = h5["regions"][name]
            points = group.get("points")
            owner = group.get("owner")
            offsets = group.get("cell_offsets")
            if (
                not isinstance(points, h5py.Dataset)
                or not isinstance(owner, h5py.Dataset)
            ):
                return None, None, None
            raw_cells = group.attrs.get(
                "n_cells",
                offsets.shape[0] - 1 if isinstance(offsets, h5py.Dataset) else 0,
            )
            counts.append((int(points.shape[0]), int(owner.shape[0]), int(raw_cells)))
        return tuple(sum(item[i] for item in counts) for i in range(3))

    points = _first_dataset(h5, ("points", "mesh/points", "mesh/topology/points"))
    owner = _first_dataset(h5, ("owner", "mesh/owner", "mesh/topology/owner"))
    offsets = _first_dataset(
        h5,
        ("cell_offsets", "mesh/cell_offsets", "mesh/topology/cell_offsets"),
    )
    n_points = int(points.shape[0]) if points is not None and points.ndim else None
    n_faces = int(owner.shape[0]) if owner is not None and owner.ndim else None
    raw_cells = h5.attrs.get("n_cells")
    if raw_cells is None and offsets is not None and offsets.shape[0]:
        raw_cells = offsets.shape[0] - 1
    n_cells = int(raw_cells) if raw_cells is not None else None
    return n_points, n_faces, n_cells


def _field_names(h5: h5py.File, kind: str) -> tuple[str, ...]:
    path = f"fields/{kind}"
    if path not in h5 or not isinstance(h5[path], h5py.Group):
        return ()
    return tuple(sorted(
        name for name, value in h5[path].items()
        if isinstance(value, h5py.Dataset)
    ))


def _case_config(h5: h5py.File) -> dict[str, Any] | None:
    raw = None
    if "case/config" in h5:
        dataset = h5["case/config"]
        if not isinstance(dataset, h5py.Dataset) or dataset.size != 1:
            return None
        if dataset.dtype.kind not in "SUO":
            return None
        if dataset.id.get_storage_size() > 1_048_576:
            return None
        raw = dataset[()]
    elif "case_setup_json" in h5.attrs:
        raw = h5.attrs["case_setup_json"]
    if raw is None:
        return None
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8")
    try:
        value = json.loads(str(raw))
    except (TypeError, ValueError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def inspect_cfdx_file(path: str | Path, *, max_datasets: int = 300) -> CaseFileSummary:
    """Inspect CFDX HDF5 metadata without reading any scientific array values.

    max_datasets bounds the returned browser tree; dataset discovery stops
    once the limit is reached and reports how many entries were omitted.
    """
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    if max_datasets < 1:
        raise ValueError("max_datasets must be positive")

    with h5py.File(source, "r") as h5:
        attrs = tuple(
            sorted((str(key), _text(value)) for key, value in h5.attrs.items())
        )
        attr_map = dict(attrs)
        file_format = attr_map.get("format")
        schema = attr_map.get("schema_version")
        regions = _region_names(h5)
        n_points, n_faces, n_cells = _mesh_counts(h5, regions)
        has_mesh = n_points is not None and n_faces is not None
        config = _case_config(h5)
        has_project_config = "case/config" in h5 and config is not None
        scalar_fields = _field_names(h5, "scalar")
        vector_fields = _field_names(h5, "vector")

        datasets: list[Hdf5DatasetInfo] = []
        datasets_truncated = False

        def collect(name: str, obj: Any):
            nonlocal datasets_truncated
            if not isinstance(obj, h5py.Dataset):
                return None
            if len(datasets) >= max_datasets:
                datasets_truncated = True
                return 1
            shape = tuple(int(size) for size in obj.shape)
            size_bytes = int(prod(shape) * obj.dtype.itemsize)
            datasets.append(Hdf5DatasetInfo(
                path="/" + name,
                shape=shape,
                dtype=str(obj.dtype),
                size_bytes=size_bytes,
                chunks=tuple(int(size) for size in obj.chunks) if obj.chunks else None,
                compression=str(obj.compression) if obj.compression else None,
            ))

        h5.visititems(collect)
        is_cfdx = file_format in _FORMATS or config is not None or has_mesh
        preview_supported = (
            has_mesh and not regions
            and n_points is not None and n_faces is not None
            and n_points <= 100_000 and n_faces <= 100_000
            and n_cells is not None
            and n_cells <= 100_000
            and all(
                name in h5
                for name in (
                    "points",
                    "face_vertices",
                    "face_offsets",
                    "owner",
                    "neighbour",
                    "patch_face_ids",
                    "patch_face_offsets",
                )
            )
        )
        return CaseFileSummary(
            path=source.resolve(),
            format=file_format,
            schema_version=schema,
            is_cfdx=is_cfdx,
            has_case_config=config is not None,
            has_project_config=has_project_config,
            has_mesh=has_mesh,
            preview_supported=preview_supported,
            n_points=n_points,
            n_faces=n_faces,
            n_cells=n_cells,
            regions=regions,
            scalar_fields=scalar_fields,
            vector_fields=vector_fields,
            case_name=(
                str(config["name"])
                if config and config.get("name")
                else attr_map.get("source_case_name")
            ),
            physics=tuple(sorted(config["physics"]))
            if config and isinstance(config.get("physics"), dict)
            else (),
            numerics=tuple(sorted(config["numerics"]))
            if config and isinstance(config.get("numerics"), dict)
            else (),
            boundaries=tuple(sorted(config["boundaries"]))
            if config and isinstance(config.get("boundaries"), dict)
            else (),
            attributes=attrs,
            datasets=tuple(datasets),
            datasets_truncated=datasets_truncated,
        )
