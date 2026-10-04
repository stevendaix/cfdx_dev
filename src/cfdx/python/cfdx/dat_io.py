"""Reader/writer for CFDX DAT numerical checkpoints.

DAT checkpoints support both the legacy text representation and an HDF5
container representation. The HDF5 representation keeps the same logical
checkpoint model while allowing large field arrays to be loaded efficiently.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math

import h5py
import numpy as np


@dataclass(frozen=True)
class DatField:
    name: str
    dimension: int
    values: list[float]

    @property
    def is_vector(self) -> bool:
        return self.dimension > 1

    def component(self, index: int) -> list[float]:
        if index < 0 or index >= self.dimension:
            raise IndexError(index)
        return self.values[index::self.dimension]


@dataclass(frozen=True)
class DatRestart:
    version: int
    cells: int
    iteration: int
    time: float
    fields: dict[str, DatField]
    # Persistent global cell identities, independent of MPI rank/local ordering.
    cell_ids: tuple[int, ...] | None = None


def _validate(restart: DatRestart) -> DatRestart:
    if restart.version not in (1, 2):
        raise ValueError(f"unsupported DAT version {restart.version}")
    if restart.cells < 0 or restart.iteration < 0:
        raise ValueError("invalid DAT checkpoint metadata")
    if not math.isfinite(restart.time):
        raise ValueError("time must be finite")
    if not restart.fields:
        raise ValueError("DAT restart contains no fields")
    if restart.cell_ids is not None:
        if len(restart.cell_ids) != restart.cells:
            raise ValueError("DAT checkpoint cell_ids count differs from cells")
        if len(set(restart.cell_ids)) != len(restart.cell_ids):
            raise ValueError("DAT checkpoint contains duplicate cell ids")
        if any(cell_id < 0 for cell_id in restart.cell_ids):
            raise ValueError("DAT checkpoint contains negative cell ids")
    for name, field in restart.fields.items():
        if not name or field.dimension <= 0:
            raise ValueError(f"invalid field {name!r}")
        if len(field.values) != restart.cells * field.dimension:
            raise ValueError(f"field {name!r} has an invalid value count")
        if not all(math.isfinite(v) for v in field.values):
            raise ValueError(f"field {name!r} contains non-finite values")
    return restart


def _read_text(path: Path) -> DatRestart:
    tokens = path.read_text(encoding="utf-8").split()
    if not tokens:
        raise ValueError("empty DAT restart")
    pos = 0

    def take(label: str) -> str:
        nonlocal pos
        if pos >= len(tokens):
            raise ValueError(f"unexpected end of DAT while reading {label}")
        value = tokens[pos]
        pos += 1
        return value

    if take("DAT magic") != "CFDX-DAT":
        raise ValueError("unsupported DAT format")
    try:
        version = int(take("DAT version"))
    except ValueError as exc:
        raise ValueError("invalid DAT version") from exc
    if version not in (1, 2):
        raise ValueError(f"unsupported DAT version {version}")

    if take("cells key") != "cells":
        raise ValueError("expected 'cells'")
    try:
        cells = int(take("cell count"))
    except ValueError as exc:
        raise ValueError("invalid cell count") from exc

    if take("iteration key") != "iteration":
        raise ValueError("expected 'iteration'")
    try:
        iteration = int(take("iteration"))
    except ValueError as exc:
        raise ValueError("invalid iteration") from exc

    if take("time key") != "time":
        raise ValueError("expected 'time'")
    try:
        time = float(take("time"))
    except ValueError as exc:
        raise ValueError("invalid time") from exc

    fields: dict[str, DatField] = {}
    while pos < len(tokens):
        if take("field key") != "field":
            raise ValueError("expected 'field'")
        name = take("field name")
        if name in fields:
            raise ValueError(f"duplicate field name: {name!r}")
        try:
            dimension = int(take(f"{name} dimension"))
        except ValueError as exc:
            raise ValueError(f"invalid dimension for field {name!r}") from exc
        values: list[float] = []
        for _ in range(cells * dimension):
            try:
                values.append(float(take(f"{name} values")))
            except ValueError as exc:
                raise ValueError(f"invalid value in field {name!r}") from exc
        fields[name] = DatField(name, dimension, values)

    return _validate(DatRestart(version, cells, iteration, time, fields))


def _read_hdf5(path: Path) -> DatRestart:
    with h5py.File(path, "r") as h5:
        if h5.attrs.get("format", "") not in ("CFDX-DAT", b"CFDX-DAT"):
            raise ValueError("not a CFDX DAT HDF5 checkpoint")
        version = int(h5.attrs.get("version", 2))
        cells = int(h5.attrs["cells"])
        iteration = int(h5.attrs["iteration"])
        time = float(h5.attrs["time"])
        cell_ids: tuple[int, ...] | None = None
        if "cell_ids" in h5:
            ids = np.asarray(h5["cell_ids"][()], dtype=np.uint64)
            if ids.ndim != 1:
                raise ValueError("DAT HDF5 cell_ids must be one-dimensional")
            cell_ids = tuple(int(value) for value in ids)
        if "fields" not in h5:
            raise ValueError("DAT HDF5 checkpoint has no fields group")
        fields: dict[str, DatField] = {}
        for name, dataset in h5["fields"].items():
            values = np.asarray(dataset[()], dtype=float)
            if values.ndim == 1:
                dimension = 1
                flat = values
            elif values.ndim == 2:
                dimension = int(values.shape[1])
                flat = values.reshape(-1)
            else:
                raise ValueError(f"invalid HDF5 shape for field {name!r}")
            if values.shape[0] != cells:
                raise ValueError(f"field {name!r} cell count mismatch")
            fields[name] = DatField(name, dimension, flat.tolist())
    return _validate(DatRestart(version, cells, iteration, time, fields, cell_ids))


def read_dat_restart(path: str | Path) -> DatRestart:
    """Read a CFDX DAT checkpoint, detecting text or HDF5 representation."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    if h5py.is_hdf5(path):
        return _read_hdf5(path)
    return _read_text(path)


def write_dat_hdf5(path: str | Path, restart: DatRestart) -> Path:
    """Write a CFDX DAT checkpoint as an HDF5 container."""
    path = Path(path)
    _validate(restart)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-DAT"
        h5.attrs["version"] = restart.version
        h5.attrs["cells"] = restart.cells
        h5.attrs["iteration"] = restart.iteration
        h5.attrs["time"] = restart.time
        if restart.cell_ids is not None:
            h5.create_dataset("cell_ids", data=np.asarray(restart.cell_ids, dtype=np.uint64))
        group = h5.create_group("fields")
        for name, field in restart.fields.items():
            values = np.asarray(field.values, dtype=np.float64).reshape(
                restart.cells, field.dimension
            )
            if field.dimension == 1:
                values = values[:, 0]
            group.create_dataset(name, data=values)
    return path


def remap_dat_restart(
    restart: DatRestart, target_cell_ids: list[int] | tuple[int, ...]
) -> DatRestart:
    """Reorder checkpoint fields from persistent source IDs to target cell IDs.

    This is rank-independent: the source and target local orderings may differ
    and may represent different MPI decompositions, provided the same global
    cell IDs are present.
    """
    _validate(restart)
    if restart.cell_ids is None:
        raise ValueError("checkpoint has no persistent cell ids")
    target = tuple(int(value) for value in target_cell_ids)
    if len(set(target)) != len(target):
        raise ValueError("target cell ids contain duplicates")
    if any(value < 0 for value in target):
        raise ValueError("target cell ids contain negative values")
    source_index = {cell_id: i for i, cell_id in enumerate(restart.cell_ids)}
    missing = [cell_id for cell_id in target if cell_id not in source_index]
    if missing:
        raise ValueError(f"target cell id missing from checkpoint: {missing[0]}")
    fields: dict[str, DatField] = {}
    for name, field in restart.fields.items():
        values: list[float] = []
        for cell_id in target:
            src = source_index[cell_id]
            begin = src * field.dimension
            end = begin + field.dimension
            values.extend(field.values[begin:end])
        fields[name] = DatField(field.name, field.dimension, values)
    return _validate(
        DatRestart(restart.version, len(target), restart.iteration, restart.time, fields, target)
    )
