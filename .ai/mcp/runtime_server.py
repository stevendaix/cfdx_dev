"""Read-only CFDX Runtime MCP server.

This first Runtime MCP surface inspects existing CFDX case/checkpoint artifacts.
It never creates, modifies, executes, or deletes cases or solver state.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import h5py
import numpy as np
from mcp.server import MCPServer
from mcp.types import ToolAnnotations


def _root(configured: str | None) -> Path:
    return Path(configured or os.environ.get("CFDX_RUNTIME_ROOT", ".")).resolve()


def _safe_path(root: Path, value: str) -> Path:
    if not value or Path(value).is_absolute():
        raise ValueError("path must be a non-empty root-relative path")
    candidate = (root / value).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise ValueError("path must remain inside the configured runtime root") from exc
    return candidate


def _scalar(value: Any) -> Any:
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _attributes(h5: h5py.File) -> dict[str, Any]:
    return {str(key): _scalar(value) for key, value in h5.attrs.items()}


def _datasets(group: h5py.Group, prefix: str = "") -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name, item in group.items():
        path = f"{prefix}/{name}" if prefix else name
        if isinstance(item, h5py.Dataset):
            rows.append({"path": path, "shape": list(item.shape), "dtype": str(item.dtype), "size": int(item.size)})
        elif isinstance(item, h5py.Group):
            rows.extend(_datasets(item, path))
    return rows


def _inspect_case(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"ok": False, "errors": [f"case file does not exist: {path.name}"]}
    if not path.name.lower().endswith(".cfdx.h5"):
        return {"ok": False, "errors": ["case must use the canonical .cfdx.h5 extension"]}
    try:
        with h5py.File(path, "r") as h5:
            if _scalar(h5.attrs.get("format")) != "CFDX":
                return {"ok": False, "errors": ["not a CFDX case artifact"]}
            raw_config = h5["case/config"][()]
            if isinstance(raw_config, bytes):
                raw_config = raw_config.decode("utf-8")
            config = json.loads(raw_config)
            return {"ok": True, "artifact": "case", "path": path.name, "attributes": _attributes(h5), "configuration": config, "datasets": _datasets(h5), "runtime_state_present": "runtime" in h5}
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        return {"ok": False, "errors": [f"invalid CFDX case artifact: {exc}"]}


def _inspect_checkpoint(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"ok": False, "errors": [f"checkpoint file does not exist: {path.name}"]}
    if not path.name.lower().endswith(".dat.h5"):
        return {"ok": False, "errors": ["checkpoint must use the canonical .dat.h5 extension"]}
    try:
        with h5py.File(path, "r") as h5:
            if _scalar(h5.attrs.get("format")) != "CFDX-DAT":
                return {"ok": False, "errors": ["not a CFDX DAT HDF5 checkpoint"]}
            if "fields" not in h5:
                return {"ok": False, "errors": ["checkpoint has no fields group"]}
            fields = [{"name": name, "shape": list(dataset.shape), "dtype": str(dataset.dtype), "size": int(dataset.size)} for name, dataset in h5["fields"].items()]
            attributes = _attributes(h5)
            required_metadata = ("version", "cells", "iteration", "time")
            missing_metadata = [name for name in required_metadata if name not in attributes]
            if missing_metadata:
                return {"ok": False, "errors": [f"checkpoint is missing required metadata: {missing_metadata}"]}
            try:
                cells = int(attributes["cells"])
                iteration = int(attributes["iteration"])
                time_value = float(attributes["time"])
            except (TypeError, ValueError) as exc:
                return {"ok": False, "errors": [f"checkpoint metadata is invalid: {exc}"]}
            if cells < 0 or iteration < 0 or not np.isfinite(time_value):
                return {"ok": False, "errors": ["checkpoint metadata contains invalid cells, iteration, or time"]}
            cell_ids_present = "cell_ids" in h5
            if cell_ids_present and h5["cell_ids"].ndim != 1:
                return {"ok": False, "errors": ["checkpoint cell_ids must be one-dimensional"]}
            if cell_ids_present and h5["cell_ids"].shape[0] != cells:
                return {"ok": False, "errors": ["checkpoint cell_ids count differs from cells"]}
            if cell_ids_present:
                cell_ids = np.asarray(h5["cell_ids"][()], dtype=np.int64)
                if np.any(cell_ids < 0):
                    return {"ok": False, "errors": ["checkpoint contains negative cell ids"]}
                if np.unique(cell_ids).size != cell_ids.size:
                    return {"ok": False, "errors": ["checkpoint contains duplicate cell ids"]}
            for name, dataset in h5["fields"].items():
                if dataset.ndim not in (1, 2):
                    return {"ok": False, "errors": [f"field {name!r} has unsupported rank: {dataset.ndim}"]}
                if dataset.shape[0] != cells:
                    return {"ok": False, "errors": [f"field {name!r} cell count mismatch"]}
                if dataset.ndim == 2 and dataset.shape[1] <= 0:
                    return {"ok": False, "errors": [f"field {name!r} has invalid dimension"]}
                if not np.all(np.isfinite(dataset[()])):
                    return {"ok": False, "errors": [f"field {name!r} contains non-finite values"]}
            return {"ok": True, "artifact": "checkpoint", "path": path.name, "attributes": attributes, "metadata": {"version": int(attributes["version"]), "cells": cells, "iteration": iteration, "time": time_value}, "cell_ids_present": cell_ids_present, "fields": fields}
    except (OSError, KeyError, TypeError, ValueError) as exc:
        return {"ok": False, "errors": [f"invalid CFDX DAT checkpoint: {exc}"]}


def _field_summary(path: Path, field_name: str, component: int | None) -> dict[str, Any]:
    if not field_name:
        return {"ok": False, "errors": ["field_name must be non-empty"]}
    try:
        with h5py.File(path, "r") as h5:
            if _scalar(h5.attrs.get("format")) != "CFDX-DAT":
                return {"ok": False, "errors": ["not a CFDX DAT HDF5 checkpoint"]}
            fields = h5.get("fields")
            if fields is None or field_name not in fields:
                return {"ok": False, "errors": [f"field does not exist: {field_name}"]}
            dataset = fields[field_name]
            if dataset.ndim not in (1, 2):
                return {"ok": False, "errors": [f"field has unsupported rank: {dataset.ndim}"]}
            dimension = 1 if dataset.ndim == 1 else int(dataset.shape[1])
            if component is not None and (component < 0 or component >= dimension):
                return {"ok": False, "errors": [f"component must be in [0, {dimension - 1}]"]}
            values = dataset[()]
            if dataset.ndim == 2:
                values = values[:, component] if component is not None else values.reshape(-1)
            values = values.astype("float64", copy=False)
            finite = values[np.isfinite(values)] if hasattr(values, "dtype") else values
            nonfinite_count = int(values.size - finite.size)
            result: dict[str, Any] = {
                "ok": True,
                "artifact": "checkpoint-field",
                "checkpoint": path.name,
                "field": field_name,
                "dimension": dimension,
                "component": component,
                "cell_count": int(dataset.shape[0]),
                "value_count": int(values.size),
                "finite": nonfinite_count == 0,
                "nonfinite_count": nonfinite_count,
            }
            if finite.size:
                result.update({
                    "min": float(finite.min()),
                    "max": float(finite.max()),
                    "mean": float(finite.mean()),
                    "l2_norm": float((finite * finite).sum() ** 0.5),
                })
            return result
    except (OSError, KeyError, TypeError, ValueError) as exc:
        return {"ok": False, "errors": [f"invalid CFDX DAT field: {exc}"]}


def _checkpoint_compare(
    left_path: Path, right_path: Path, field_name: str, component: int | None
) -> dict[str, Any]:
    if not field_name:
        return {"ok": False, "errors": ["field_name must be non-empty"]}
    try:
        with h5py.File(left_path, "r") as left, h5py.File(right_path, "r") as right:
            if _scalar(left.attrs.get("format")) != "CFDX-DAT" or _scalar(right.attrs.get("format")) != "CFDX-DAT":
                return {"ok": False, "errors": ["both artifacts must be CFDX DAT HDF5 checkpoints"]}
            left_fields, right_fields = left.get("fields"), right.get("fields")
            if left_fields is None or right_fields is None:
                return {"ok": False, "errors": ["both checkpoints must contain a fields group"]}
            if field_name not in left_fields or field_name not in right_fields:
                return {"ok": False, "errors": [f"field does not exist in both checkpoints: {field_name}"]}
            a, b = left_fields[field_name], right_fields[field_name]
            if a.shape != b.shape or a.ndim not in (1, 2):
                return {"ok": False, "errors": ["field shapes are incompatible"]}
            dimension = 1 if a.ndim == 1 else int(a.shape[1])
            if component is not None and not 0 <= component < dimension:
                return {"ok": False, "errors": [f"component must be in [0, {dimension - 1}]"]}
            av, bv = a[()], b[()]
            if a.ndim == 2:
                av = av[:, component] if component is not None else av.reshape(-1)
                bv = bv[:, component] if component is not None else bv.reshape(-1)
            av, bv = av.astype("float64", copy=False), bv.astype("float64", copy=False)
            if not (np.all(np.isfinite(av)) and np.all(np.isfinite(bv))):
                return {"ok": False, "errors": ["field comparison requires finite values in both checkpoints"]}
            delta = bv - av
            diff_l2 = float(np.linalg.norm(delta.ravel()))
            reference_l2 = float(np.linalg.norm(av.ravel()))
            return {
                "ok": True,
                "artifact": "checkpoint-field-comparison",
                "left_checkpoint": left_path.name,
                "right_checkpoint": right_path.name,
                "field": field_name,
                "dimension": dimension,
                "component": component,
                "value_count": int(av.size),
                "max_abs_difference": float(np.max(np.abs(delta))) if delta.size else 0.0,
                "l2_difference": diff_l2,
                "reference_l2": reference_l2,
                "relative_l2_difference": diff_l2 / reference_l2 if reference_l2 else (0.0 if diff_l2 == 0.0 else None),
            }
    except (OSError, KeyError, TypeError, ValueError) as exc:
        return {"ok": False, "errors": [f"invalid CFDX DAT comparison: {exc}"]}


def create_server(root: str | None = None) -> MCPServer:
    runtime_root = _root(root)
    server = MCPServer("CFDX Runtime MCP", instructions="Read-only CFDX runtime artifact inspection. These tools inspect case/checkpoint files only; they never execute CFDX, mutate artifacts, or claim solver success.")
    annotations = ToolAnnotations(read_only_hint=True, open_world_hint=False)

    @server.tool(name="case.inspect", title="Inspect CFDX case", annotations=annotations)
    def case_inspect(case_path: str) -> dict[str, Any]:
        """Inspect a root-relative .cfdx.h5 case without loading solver state."""
        try:
            return _inspect_case(_safe_path(runtime_root, case_path))
        except ValueError as exc:
            return {"ok": False, "errors": [str(exc)]}

    @server.tool(name="checkpoint.compare", title="Compare CFDX checkpoints", annotations=annotations)
    def checkpoint_compare(left_checkpoint_path: str, right_checkpoint_path: str, field_name: str, component: int | None = None) -> dict[str, Any]:
        """Compare one field between two checkpoints without modifying runtime state."""
        try:
            left = _safe_path(runtime_root, left_checkpoint_path)
            right = _safe_path(runtime_root, right_checkpoint_path)
            for path in (left, right):
                if not path.is_file() or not path.name.lower().endswith(".dat.h5"):
                    return {"ok": False, "errors": ["both checkpoints must be existing canonical .dat.h5 artifacts"]}
            return _checkpoint_compare(left, right, field_name, component)
        except ValueError as exc:
            return {"ok": False, "errors": [str(exc)]}

    @server.tool(name="checkpoint.field.inspect", title="Inspect CFDX checkpoint field", annotations=annotations)
    def checkpoint_field_inspect(checkpoint_path: str, field_name: str, component: int | None = None) -> dict[str, Any]:
        """Inspect scalar/vector field statistics without returning field arrays."""
        try:
            path = _safe_path(runtime_root, checkpoint_path)
            if not path.is_file() or not path.name.lower().endswith(".dat.h5"):
                return {"ok": False, "errors": ["checkpoint must be an existing canonical .dat.h5 artifact"]}
            return _field_summary(path, field_name, component)
        except ValueError as exc:
            return {"ok": False, "errors": [str(exc)]}

    @server.tool(name="checkpoint.inspect", title="Inspect CFDX checkpoint", annotations=annotations)
    def checkpoint_inspect(checkpoint_path: str) -> dict[str, Any]:
        """Inspect a root-relative .dat.h5 checkpoint without executing CFDX."""
        try:
            return _inspect_checkpoint(_safe_path(runtime_root, checkpoint_path))
        except ValueError as exc:
            return {"ok": False, "errors": [str(exc)]}

    return server


mcp = create_server()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=None)
    parser.add_argument("--transport", choices=("stdio", "streamable-http"), default="stdio")
    args = parser.parse_args()
    create_server(str(args.root) if args.root is not None else None).run(transport=args.transport)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
