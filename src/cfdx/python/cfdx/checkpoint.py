"""Versioned, atomic persistence for orchestration checkpoint metadata."""
from __future__ import annotations

import json
import math
import os
import tempfile
from pathlib import Path

from .session import Checkpoint

CHECKPOINT_SCHEMA_VERSION = 1
_CHECKPOINT_FIELDS = {
    "schema_version",
    "case_revision",
    "mesh_revision",
    "physics_revision",
    "numerics_revision",
    "iteration",
    "time",
}


def _validate_integer(value: object, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"checkpoint field {field!r} must be an integer")
    if value < 0:
        raise ValueError(f"checkpoint field {field!r} must be non-negative")
    return value


def _validate_time(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("checkpoint field 'time' must be a number")
    time = float(value)
    if not math.isfinite(time):
        raise ValueError("checkpoint field 'time' must be finite")
    return time


def save_checkpoint(checkpoint: Checkpoint, path: Path) -> None:
    """Persist checkpoint metadata atomically.

    The temporary file is created in the destination directory so os.replace()
    is atomic on the same filesystem. A partially written checkpoint is never
    exposed at the destination path.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "case_revision": checkpoint.case_revision,
        "mesh_revision": checkpoint.mesh_revision,
        "physics_revision": checkpoint.physics_revision,
        "numerics_revision": checkpoint.numerics_revision,
        "iteration": checkpoint.iteration,
        "time": checkpoint.time,
    }
    # Validate the serialized representation before publishing it.
    _checkpoint_from_data(payload)

    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
        text=True,
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def _checkpoint_from_data(data: object) -> Checkpoint:
    if not isinstance(data, dict) or set(data) != _CHECKPOINT_FIELDS:
        raise ValueError("invalid CFDX checkpoint fields")

    schema_version = _validate_integer(data["schema_version"], "schema_version")
    if schema_version != CHECKPOINT_SCHEMA_VERSION:
        raise ValueError(
            f"unsupported CFDX checkpoint schema version {schema_version}"
        )

    return Checkpoint(
        case_revision=_validate_integer(data["case_revision"], "case_revision"),
        mesh_revision=_validate_integer(data["mesh_revision"], "mesh_revision"),
        physics_revision=_validate_integer(data["physics_revision"], "physics_revision"),
        numerics_revision=_validate_integer(data["numerics_revision"], "numerics_revision"),
        iteration=_validate_integer(data["iteration"], "iteration"),
        time=_validate_time(data["time"]),
    )


def load_checkpoint(path: Path) -> Checkpoint:
    """Load a checkpoint and reject unknown schema versions or invalid metadata."""
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    return _checkpoint_from_data(data)
