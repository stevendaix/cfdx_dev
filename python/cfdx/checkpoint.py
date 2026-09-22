"""Portable checkpoint metadata persistence for orchestration state."""
from __future__ import annotations

import json
from pathlib import Path

from .session import Checkpoint


def save_checkpoint(checkpoint: Checkpoint, path: Path) -> None:
    path.write_text(json.dumps(checkpoint.__dict__, indent=2, sort_keys=True) + "
", encoding="utf-8")


def load_checkpoint(path: Path) -> Checkpoint:
    data = json.loads(path.read_text(encoding="utf-8"))
    expected = {"case_revision", "mesh_revision", "physics_revision", "numerics_revision", "iteration", "time"}
    if set(data) != expected:
        raise ValueError("invalid CFDX checkpoint fields")
    return Checkpoint(
        case_revision=int(data["case_revision"]),
        mesh_revision=int(data["mesh_revision"]),
        physics_revision=int(data["physics_revision"]),
        numerics_revision=int(data["numerics_revision"]),
        iteration=int(data["iteration"]),
        time=float(data["time"]),
    )
