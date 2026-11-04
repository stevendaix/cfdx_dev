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
    size_bytes: int | None
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
    mesh_counts_complete: bool
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
    tree_truncated: bool = False
    attributes_truncated: bool = False
    regions_truncated: bool = False
    scalar_fields_truncated: bool = False
    vector_fields_truncated: bool = False
    file_signature: tuple[int, int, int, int] | None = None

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


def _internal_object(h5: h5py.Group, path: str):
    """Resolve only hard links inside the inspected HDF5 file."""
    parts = path.strip("/").split("/")
    if not parts or any(part in {"", ".", ".."} for part in parts):
        return None
    current = h5
    for index, part in enumerate(parts):
        try:
            link = current.get(part, getlink=True)
        except KeyError:
            return None
        if not isinstance(link, h5py.HardLink):
            return None
        try:
            current = current[part]
        except (KeyError, OSError, RuntimeError):
            return None
        if index < len(parts) - 1 and not isinstance(current, h5py.Group):
            return None
    return current


def _first_dataset(h5: h5py.File, candidates: tuple[str, ...]):
    for candidate in candidates:
        dataset = _internal_object(h5, candidate)
        if isinstance(dataset, h5py.Dataset):
            return dataset
    return None


def _region_names(
    h5: h5py.File, *, max_items: int = 300
) -> tuple[tuple[str, ...], bool]:
    region_group = _internal_object(h5, "regions")
    if not isinstance(region_group, h5py.Group):
        return (), False
    names = []
    truncated = False
    scanned = 0
    scan_limit = max_items * 4
    for name in region_group:
        scanned += 1
        if scanned > scan_limit:
            truncated = True
            break
        if not isinstance(region_group.get(name, getlink=True), h5py.HardLink):
            truncated = True
            continue
        value = region_group[name]
        if not isinstance(value, h5py.Group):
            truncated = True
            continue
        if len(names) >= max_items:
            truncated = True
            break
        names.append(name)
    return tuple(sorted(names)), truncated


def _nonnegative_count(value: Any) -> int | None:
    """Parse a stored entity count without accepting malformed or negative data."""
    if isinstance(value, bool):
        return None
    try:
        count = int(value)
    except (TypeError, ValueError, OverflowError):
        return None
    if count < 0 or (not isinstance(value, (str, bytes)) and count != value):
        return None
    return count


def _mesh_counts(h5: h5py.File, regions: tuple[str, ...]):
    if regions:
        counts = []
        for name in regions:
            group = _internal_object(h5, f"regions/{name}")
            if not isinstance(group, h5py.Group):
                return None, None, None
            points = _internal_object(group, "points")
            owner = _internal_object(group, "owner")
            offsets = _internal_object(group, "cell_offsets")
            if not isinstance(points, h5py.Dataset) or not isinstance(
                owner, h5py.Dataset
            ):
                return None, None, None
            raw_cells = _nonnegative_count(_safe_attribute(group, "n_cells"))
            if raw_cells is None and isinstance(offsets, h5py.Dataset):
                raw_cells = max(0, int(offsets.shape[0]) - 1)
            if raw_cells is None:
                return None, None, None
            counts.append((int(points.shape[0]), int(owner.shape[0]), raw_cells))
        return tuple(sum(item[i] for item in counts) for i in range(3))

    points = _first_dataset(h5, ("points", "mesh/points", "mesh/topology/points"))
    owner = _first_dataset(h5, ("owner", "mesh/owner", "mesh/topology/owner"))
    offsets = _first_dataset(
        h5,
        ("cell_offsets", "mesh/cell_offsets", "mesh/topology/cell_offsets"),
    )
    n_points = int(points.shape[0]) if points is not None and points.ndim else None
    n_faces = int(owner.shape[0]) if owner is not None and owner.ndim else None
    raw_cells = _nonnegative_count(_safe_attribute(h5, "n_cells"))
    if raw_cells is None and offsets is not None and offsets.shape[0]:
        raw_cells = max(0, int(offsets.shape[0]) - 1)
    n_cells = raw_cells
    return n_points, n_faces, n_cells


def _field_names(
    h5: h5py.File, kind: str, *, max_items: int = 300
) -> tuple[tuple[str, ...], bool]:
    path = f"fields/{kind}"
    field_group = _internal_object(h5, path)
    if not isinstance(field_group, h5py.Group):
        return (), False
    names = []
    truncated = False
    scanned = 0
    scan_limit = max_items * 4
    for name in field_group:
        scanned += 1
        if scanned > scan_limit:
            truncated = True
            break
        if not isinstance(field_group.get(name, getlink=True), h5py.HardLink):
            truncated = True
            continue
        value = field_group[name]
        if not isinstance(value, h5py.Dataset):
            truncated = True
            continue
        if len(names) >= max_items:
            truncated = True
            break
        names.append(name)
    return tuple(sorted(names)), truncated


def _case_config(h5: h5py.File) -> dict[str, Any] | None:
    raw = None
    project_config = _internal_object(h5, "case/config")
    if project_config is not None:
        dataset = project_config
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


@dataclass(frozen=True)
class Hdf5DatasetSample:
    path: str
    shape: tuple[int, ...]
    dtype: str
    total_values: int
    values: tuple[str, ...]


def read_hdf5_dataset_sample(
    path: str | Path, dataset_path: str, *, max_values: int = 32
) -> Hdf5DatasetSample:
    """Read a bounded numeric/fixed-string sample, never a full HDF5 array."""
    if not dataset_path.startswith("/") or dataset_path == "/":
        raise ValueError("dataset_path must be an absolute HDF5 dataset path")
    parts = dataset_path.strip("/").split("/")
    if any(part in {"", ".", ".."} for part in parts):
        raise ValueError("dataset_path is not a safe HDF5 path")
    if not 1 <= max_values <= 256:
        raise ValueError("max_values must be between 1 and 256")

    with h5py.File(Path(path), "r") as h5:
        dataset = _internal_object(h5, dataset_path)
        if dataset is None:
            raise ValueError(
                "dataset path must resolve through internal HDF5 hard links"
            )
        if not isinstance(dataset, h5py.Dataset):
            raise TypeError(f"{dataset_path} is not an HDF5 dataset")
        shape = tuple(int(size) for size in dataset.shape)
        dtype = dataset.dtype
        if dtype.kind not in "biufcS":
            raise ValueError(
                f"sampling dtype {dtype} is disabled; only numeric and fixed "
                "strings are supported"
            )
        if dtype.itemsize * max_values > 1_048_576:
            raise ValueError("dataset element width exceeds the 1 MiB sample budget")
        if dataset.is_virtual:
            raise ValueError("virtual datasets are not sampled")
        if dataset.external:
            raise ValueError("externally stored datasets are not sampled")
        if dataset.chunks is not None:
            chunk_bytes = prod(dataset.chunks) * dtype.itemsize
            if chunk_bytes > 1_048_576:
                raise ValueError(
                    "dataset chunk exceeds the 1 MiB decompression sample budget"
                )

        if not shape:
            values = (_text(dataset[()]),)
        elif any(size == 0 for size in shape):
            values = ()
        else:
            counts = [1] * len(shape)
            remaining = max_values
            for axis in reversed(range(len(shape))):
                counts[axis] = min(shape[axis], remaining)
                remaining = max(1, remaining // counts[axis])
            selection = tuple(slice(0, count) for count in counts)
            import numpy as np

            raw = np.asarray(dataset[selection]).reshape(-1)
            values = tuple(_text(value) for value in raw[:max_values])

        return Hdf5DatasetSample(
            path=dataset_path,
            shape=shape,
            dtype=str(dtype),
            total_values=int(prod(shape)),
            values=values,
        )


def inspect_cfdx_file(path: str | Path, *, max_datasets: int = 300) -> CaseFileSummary:
    """Inspect CFDX HDF5 metadata without reading any scientific array values.

    max_datasets bounds the returned browser tree; dataset discovery stops
    once the limit is reached and reports how many entries were omitted.
    """
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    if not 1 <= max_datasets <= 5000:
        raise ValueError("max_datasets must be between 1 and 5000")

    with h5py.File(source, "r") as h5:
        attrs, attributes_truncated = _root_attributes(h5)
        attr_map = dict(attrs)
        file_format = attr_map.get("format")
        schema = attr_map.get("schema_version")
        regions, regions_truncated = _region_names(h5)
        n_points, n_faces, n_cells = _mesh_counts(h5, regions)
        has_mesh = n_points is not None and n_faces is not None
        config = _case_config(h5)
        project_config = _internal_object(h5, "case/config")
        has_project_config = (
            isinstance(project_config, h5py.Dataset) and config is not None
        )
        scalar_fields, scalar_fields_truncated = _field_names(h5, "scalar")
        vector_fields, vector_fields_truncated = _field_names(h5, "vector")

        datasets: list[Hdf5DatasetInfo] = []
        datasets_truncated = False
        tree_truncated = False
        visited_objects = 0
        object_limit = min(20_000, max(1_000, max_datasets * 20))

        def collect(name: str, obj: Any):
            nonlocal datasets_truncated, tree_truncated, visited_objects
            visited_objects += 1
            if visited_objects > object_limit:
                tree_truncated = True
                return 1
            if not isinstance(obj, h5py.Dataset):
                return None
            if len(datasets) >= max_datasets:
                datasets_truncated = True
                return 1
            shape = tuple(int(size) for size in obj.shape)
            size_bytes = (
                None if obj.dtype.kind == "O" else int(prod(shape) * obj.dtype.itemsize)
            )
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
        region_link = h5.get("regions", getlink=True)
        region_object = _internal_object(h5, "regions")
        if not has_mesh:
            preview_block_reason = "mesh topology metadata is incomplete"
        elif region_link is not None and not isinstance(region_object, h5py.Group):
            preview_block_reason = "region catalog must be an internal HDF5 group"
        elif regions_truncated:
            preview_block_reason = (
                "region catalog is incomplete or contains non-internal links"
            )
        elif regions:
            preview_block_reason = (
                "region-grouped meshes require the streaming renderer"
            )
        elif n_points is None or n_faces is None or n_cells is None:
            preview_block_reason = "point, face and cell counts must all be declared"
        elif max(n_points, n_faces, n_cells) > 100_000:
            preview_block_reason = "mesh exceeds the 100k entity safety limit"
        elif not all(
            isinstance(_internal_object(h5, name), h5py.Dataset)
            for name in preview_datasets
        ):
            preview_block_reason = "required surface topology datasets are missing"
        else:
            points_dataset = _internal_object(h5, "points")
            face_vertices_dataset = _internal_object(h5, "face_vertices")
            face_offsets_dataset = _internal_object(h5, "face_offsets")
            owner_dataset = _internal_object(h5, "owner")
            neighbour_dataset = _internal_object(h5, "neighbour")
            patch_ids_dataset = _internal_object(h5, "patch_face_ids")
            patch_offsets_dataset = _internal_object(h5, "patch_face_offsets")
            integer_topology = (
                face_vertices_dataset,
                face_offsets_dataset,
                owner_dataset,
                neighbour_dataset,
                patch_ids_dataset,
                patch_offsets_dataset,
            )
            if points_dataset.ndim != 2 or points_dataset.shape[1] != 3:
                preview_block_reason = "points must have shape (n_points, 3)"
            elif points_dataset.dtype.kind != "f":
                preview_block_reason = (
                    "point coordinates must use a floating-point dtype"
                )
            elif any(dataset.dtype.kind not in "iu" for dataset in integer_topology):
                preview_block_reason = (
                    "surface connectivity datasets must use integer dtypes"
                )
            elif neighbour_dataset.dtype.kind != "i":
                preview_block_reason = (
                    "neighbour indices must use a signed integer dtype"
                )
            elif any(
                dataset.is_virtual or dataset.external
                for dataset in integer_topology + (points_dataset,)
            ):
                preview_block_reason = (
                    "virtual or externally stored topology datasets are not previewed"
                )
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
            elif patch_offsets_dataset.ndim != 1:
                preview_block_reason = "patch offsets are malformed"
            elif patch_ids_dataset.ndim != 1 or patch_offsets_dataset.shape[0] < 2:
                preview_block_reason = "patch face metadata is missing"
            else:
                source_bytes = sum(
                    int(
                        prod(_internal_object(h5, name).shape)
                        * _internal_object(h5, name).dtype.itemsize
                    )
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
            mesh_counts_complete=(
                not regions_truncated
                and n_points is not None
                and n_faces is not None
                and n_cells is not None
            ),
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
            tree_truncated=tree_truncated,
            attributes_truncated=attributes_truncated,
            regions_truncated=regions_truncated,
            scalar_fields_truncated=scalar_fields_truncated,
            vector_fields_truncated=vector_fields_truncated,
        )
