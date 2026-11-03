"""Read-only, metadata-first inspection of CFDX HDF5 artifacts.

The inspector deliberately reads HDF5 metadata only. It never materializes mesh
coordinates, connectivity, cell IDs, or field arrays, so opening a large case
does not imply loading its scientific payload into RAM.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from math import prod
from pathlib import Path
from typing import Any

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
    preview_block_reason: str | None
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
    attributes_truncated: bool = False

    @property
    def total_field_count(self) -> int:
        return len(self.scalar_fields) + len(self.vector_fields)


def _safe_attribute(
    h5: h5py.Group,
    name: str,
    *,
    max_bytes: int = 4096,
    allow_variable_length: bool = False,
):
    """Read one attribute only when its stored payload is explicitly bounded."""
    if name not in h5.attrs:
        return None
    attribute_id = h5.attrs.get_id(name)
    if attribute_id.dtype.kind == "O" and not allow_variable_length:
        return None
    if attribute_id.get_storage_size() > max_bytes:
        return None
    value = h5.attrs[name]
    if isinstance(value, str) and len(value.encode("utf-8")) > max_bytes:
        return None
    return value


def _root_attributes(
    h5: h5py.File, *, max_count: int = 100, max_value_bytes: int = 4096
) -> tuple[tuple[tuple[str, str], ...], bool]:
    """Return a bounded, display-safe snapshot of root attributes."""
    attributes = []
    truncated = False
    for index, name in enumerate(h5.attrs):
        if index >= max_count:
            truncated = True
            break
        attribute_id = h5.attrs.get_id(name)
        size = attribute_id.get_storage_size()
        if attribute_id.dtype.kind == "O" and name not in {
            "format",
            "schema_version",
            "source_case_name",
        }:
            value = "<omitted: variable-length attribute>"
        elif size > max_value_bytes:
            value = f"<omitted: {size:,} bytes>"
        else:
            raw = h5.attrs[name]
            if isinstance(raw, str) and len(raw.encode("utf-8")) > max_value_bytes:
                value = "<omitted: oversized variable-length attribute>"
            else:
                value = _text(raw)
        attributes.append((str(name), value))
    return tuple(sorted(attributes)), truncated


def _first_dataset(h5: h5py.File, candidates: tuple[str, ...]):
    for candidate in candidates:
        if candidate in h5 and isinstance(h5[candidate], h5py.Dataset):
            return h5[candidate]
    return None


def _region_names(h5: h5py.File) -> tuple[str, ...]:
    if "regions" not in h5 or not isinstance(h5["regions"], h5py.Group):
        return ()
    return tuple(
        sorted(
            name
            for name, value in h5["regions"].items()
            if isinstance(value, h5py.Group)
        )
    )


def _mesh_counts(h5: h5py.File, regions: tuple[str, ...]):
    if regions:
        counts = []
        for name in regions:
            group = h5["regions"][name]
            points = group.get("points")
            owner = group.get("owner")
            offsets = group.get("cell_offsets")
            if not isinstance(points, h5py.Dataset) or not isinstance(
                owner, h5py.Dataset
            ):
                return None, None, None
            raw_cells = _safe_attribute(group, "n_cells")
            if raw_cells is None:
                raw_cells = (
                    offsets.shape[0] - 1 if isinstance(offsets, h5py.Dataset) else 0
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
    raw_cells = _safe_attribute(h5, "n_cells")
    if raw_cells is None and offsets is not None and offsets.shape[0]:
        raw_cells = offsets.shape[0] - 1
    n_cells = int(raw_cells) if raw_cells is not None else None
    return n_points, n_faces, n_cells


def _field_names(h5: h5py.File, kind: str) -> tuple[str, ...]:
    path = f"fields/{kind}"
    if path not in h5 or not isinstance(h5[path], h5py.Group):
        return ()
    return tuple(
        sorted(
            name for name, value in h5[path].items() if isinstance(value, h5py.Dataset)
        )
    )


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
        raw = _safe_attribute(
            h5,
            "case_setup_json",
            max_bytes=1_048_576,
            allow_variable_length=True,
        )
    if raw is None:
        return None
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8")
    if isinstance(raw, str) and len(raw.encode("utf-8")) > 1_048_576:
        return None
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
        attrs, attributes_truncated = _root_attributes(h5)
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
            datasets.append(
                Hdf5DatasetInfo(
                    path="/" + name,
                    shape=shape,
                    dtype=str(obj.dtype),
                    size_bytes=size_bytes,
                    chunks=tuple(int(size) for size in obj.chunks)
                    if obj.chunks
                    else None,
                    compression=str(obj.compression) if obj.compression else None,
                )
            )

        h5.visititems(collect)
        is_cfdx = file_format in _FORMATS or config is not None or has_mesh
        preview_datasets = (
            "points",
            "face_vertices",
            "face_offsets",
            "owner",
            "neighbour",
            "patch_face_ids",
            "patch_face_offsets",
        )
        preview_block_reason = None
        if not has_mesh:
            preview_block_reason = "mesh topology metadata is incomplete"
        elif regions:
            preview_block_reason = (
                "region-grouped meshes require the streaming renderer"
            )
        elif n_points is None or n_faces is None or n_cells is None:
            preview_block_reason = "point, face and cell counts must all be declared"
        elif max(n_points, n_faces, n_cells) > 100_000:
            preview_block_reason = "mesh exceeds the 100k entity safety limit"
        elif not all(
            name in h5 and isinstance(h5[name], h5py.Dataset)
            for name in preview_datasets
        ):
            preview_block_reason = "required surface topology datasets are missing"
        else:
            points_dataset = h5["points"]
            face_vertices_dataset = h5["face_vertices"]
            face_offsets_dataset = h5["face_offsets"]
            owner_dataset = h5["owner"]
            neighbour_dataset = h5["neighbour"]
            patch_ids_dataset = h5["patch_face_ids"]
            patch_offsets_dataset = h5["patch_face_offsets"]
            if points_dataset.ndim != 2 or points_dataset.shape[1] != 3:
                preview_block_reason = "points must have shape (n_points, 3)"
            elif face_vertices_dataset.ndim != 1:
                preview_block_reason = (
                    "face connectivity must be a one-dimensional array"
                )
            elif face_offsets_dataset.shape != (n_faces + 1,):
                preview_block_reason = "face offsets do not match the face count"
            elif owner_dataset.shape != (n_faces,) or neighbour_dataset.shape != (
                n_faces,
            ):
                preview_block_reason = (
                    "owner/neighbour arrays do not match the face count"
                )
            elif patch_offsets_dataset.ndim != 1 or patch_offsets_dataset.shape[0] < 1:
                preview_block_reason = "patch offsets are malformed"
            elif patch_ids_dataset.ndim != 1 or patch_offsets_dataset.shape[0] < 2:
                preview_block_reason = "patch face metadata is missing"
            else:
                source_bytes = sum(
                    int(prod(h5[name].shape) * h5[name].dtype.itemsize)
                    for name in preview_datasets
                )
                connectivity_entries = int(face_vertices_dataset.shape[0])
                if source_bytes > 32 * 1024 * 1024:
                    preview_block_reason = (
                        "surface topology exceeds the 32 MiB eager-load budget"
                    )
                elif connectivity_entries > 2_000_000:
                    preview_block_reason = (
                        "surface connectivity exceeds the 2M-entry eager-load budget"
                    )
        preview_supported = preview_block_reason is None
        return CaseFileSummary(
            path=source.resolve(),
            format=file_format,
            schema_version=schema,
            is_cfdx=is_cfdx,
            has_case_config=config is not None,
            has_project_config=has_project_config,
            has_mesh=has_mesh,
            preview_supported=preview_supported,
            preview_block_reason=preview_block_reason,
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
            attributes_truncated=attributes_truncated,
        )
