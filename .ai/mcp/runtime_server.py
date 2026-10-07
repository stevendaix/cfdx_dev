"""Read-only CFDX Runtime MCP server.

This first Runtime MCP surface inspects existing CFDX case/checkpoint artifacts.
It never creates, modifies, executes, or deletes cases or solver state.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

import h5py
import numpy as np

_SRC_ROOT = Path(__file__).resolve().parents[2] / "src" / "cfdx" / "python"
if str(_SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(_SRC_ROOT))

from mcp.server import MCPServer
from mcp.types import ToolAnnotations
from cfdx.case_io import validate_case_bundle


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


def _validate_case(path: Path) -> dict[str, Any]:
    try:
        result = validate_case_bundle(path)
        return {"ok": True, "artifact": "case-validation", "path": path.name, "checks": result}
    except (OSError, KeyError, TypeError, ValueError) as exc:
        return {"ok": False, "errors": [f"invalid CFDX case artifact: {exc}"]}


_CONVERGENCE_REQUIRED = (
    "format",
    "schema_version",
    "converged",
    "iterations",
    "reference_momentum_residual",
    "convergence_status",
    "convergence_reason",
    "history",
)

_CONVERGENCE_HISTORY_REQUIRED = (
    "iteration",
    "momentum_residual",
    "pressure_residual",
    "continuity_l1",
    "continuity_linf",
    "continuity_normalized",
    "momentum_equation_residual",
    "momentum_equation_residual_relative",
    "velocity_change_inf",
    "pressure_change_inf",
    "effective_alpha_u",
    "effective_alpha_p",
    "nonlinear_convergence_metric",
    "momentum_linear_iterations",
    "pressure_linear_iterations",
    "pressure_correctors_used",
    "linear_tolerance_used",
    "corrected_flux_continuity_linf",
    "reconstructed_velocity_continuity_linf",
    "flux_velocity_mismatch_linf",
    "momentum_equation_residual_components",
    "momentum_equation_residual_internal",
    "momentum_equation_residual_boundary",
    "pressure_gradient_linf",
    "pressure_gradient_l2",
    "momentum_residual_cell",
    "momentum_residual_no_pressure",
    "momentum_pressure_contribution",
    "momentum_residual_patch",
    "mass_boundary_flux",
    "mass_global_cell_balance",
    "mass_local_l1",
    "mass_local_linf",
    "mass_local_l2",
    "mass_normalized_imbalance",
    "mass_worst_cell",
    "mass_nonfinite_faces",
    "boundedness_nonfinite_velocity",
    "momentum_conservation_residual",
    "momentum_conservation_normalized",
    "momentum_conservation_worst_cell",
)

_CONVERGENCE_NUMERIC_FIELDS = {
    "momentum_residual",
    "pressure_residual",
    "continuity_l1",
    "continuity_linf",
    "continuity_normalized",
    "momentum_equation_residual",
    "momentum_equation_residual_relative",
    "velocity_change_inf",
    "pressure_change_inf",
    "effective_alpha_u",
    "effective_alpha_p",
    "nonlinear_convergence_metric",
    "linear_tolerance_used",
    "corrected_flux_continuity_linf",
    "reconstructed_velocity_continuity_linf",
    "flux_velocity_mismatch_linf",
    "momentum_equation_residual_internal",
    "momentum_equation_residual_boundary",
    "pressure_gradient_linf",
    "pressure_gradient_l2",
    "momentum_residual_no_pressure",
    "momentum_pressure_contribution",
    "mass_boundary_flux",
    "mass_global_cell_balance",
    "mass_local_l1",
    "mass_local_linf",
    "mass_local_l2",
    "mass_normalized_imbalance",
}

_CONVERGENCE_INTEGER_FIELDS = {
    "iteration",
    "momentum_linear_iterations",
    "pressure_linear_iterations",
    "pressure_correctors_used",
    "momentum_residual_cell",
    "mass_worst_cell",
    "mass_nonfinite_faces",
    "boundedness_nonfinite_velocity",
}


def _reject_nonfinite(value: str) -> None:
    raise ValueError(f"convergence artifact contains non-finite JSON constant: {value}")


def _validate_convergence(path: Path, history_limit: int) -> dict[str, Any]:
    if not path.is_file():
        return {"ok": False, "errors": [f"convergence artifact does not exist: {path.name}"]}
    if not path.name.lower().endswith("convergence.json"):
        return {"ok": False, "errors": ["convergence artifact must use the canonical convergence.json filename"]}
    if history_limit < 0 or history_limit > 10000:
        return {"ok": False, "errors": ["history_limit must be between 0 and 10000"]}
    try:
        with path.open("r", encoding="utf-8") as stream:
            document = json.load(stream, parse_constant=_reject_nonfinite)
        if not isinstance(document, dict):
            return {"ok": False, "errors": ["convergence artifact root must be an object"]}
        missing = [key for key in _CONVERGENCE_REQUIRED if key not in document]
        if missing:
            return {"ok": False, "errors": [f"convergence artifact is missing required keys: {missing}"]}
        if document["format"] != "CFDX-CONVERGENCE":
            return {"ok": False, "errors": ["not a CFDX-CONVERGENCE artifact"]}
        if document["schema_version"] != 1:
            return {"ok": False, "errors": [f"unsupported convergence schema version: {document['schema_version']}"]}
        if not isinstance(document["converged"], bool):
            return {"ok": False, "errors": ["converged must be a boolean"]}
        if not isinstance(document["iterations"], int) or isinstance(document["iterations"], bool) or document["iterations"] < 0:
            return {"ok": False, "errors": ["iterations must be a non-negative integer"]}
        if not isinstance(document["reference_momentum_residual"], (int, float)) or isinstance(document["reference_momentum_residual"], bool):
            return {"ok": False, "errors": ["reference_momentum_residual must be numeric"]}
        if not np.isfinite(float(document["reference_momentum_residual"])):
            return {"ok": False, "errors": ["reference_momentum_residual must be finite"]}
        if not isinstance(document["convergence_status"], int) or isinstance(document["convergence_status"], bool):
            return {"ok": False, "errors": ["convergence_status must be an integer"]}
        if not isinstance(document["convergence_reason"], str):
            return {"ok": False, "errors": ["convergence_reason must be a string"]}
        history = document["history"]
        if not isinstance(history, list):
            return {"ok": False, "errors": ["history must be an array"]}
        for index, record in enumerate(history):
            if not isinstance(record, dict):
                return {"ok": False, "errors": [f"history record {index} must be an object"]}
            missing = [key for key in _CONVERGENCE_HISTORY_REQUIRED if key not in record]
            if missing:
                return {"ok": False, "errors": [f"history record {index} is missing required keys: {missing}"]}
            for key in _CONVERGENCE_NUMERIC_FIELDS:
                value = record[key]
                if not isinstance(value, (int, float)) or isinstance(value, bool) or not np.isfinite(float(value)):
                    return {"ok": False, "errors": [f"history record {index} field {key!r} must be finite numeric"]}
            for key in _CONVERGENCE_INTEGER_FIELDS:
                value = record[key]
                if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                    return {"ok": False, "errors": [f"history record {index} field {key!r} must be a non-negative integer"]}
            for key in ("momentum_residual_patch",):
                if not isinstance(record[key], str):
                    return {"ok": False, "errors": [f"history record {index} field {key!r} must be a string"]}
            for key in ("momentum_equation_residual_components", "momentum_conservation_residual", "momentum_conservation_normalized"):
                value = record[key]
                if not isinstance(value, list) or len(value) != 3 or any(
                    not isinstance(item, (int, float)) or isinstance(item, bool) or not np.isfinite(float(item))
                    for item in value
                ):
                    return {"ok": False, "errors": [f"history record {index} field {key!r} must be a finite numeric array of length 3"]}
            worst_cell = record["momentum_conservation_worst_cell"]
            if not isinstance(worst_cell, list) or len(worst_cell) != 3 or any(
                not isinstance(item, int) or isinstance(item, bool) or item < 0 for item in worst_cell
            ):
                return {"ok": False, "errors": [f"history record {index} field 'momentum_conservation_worst_cell' must be a non-negative integer array of length 3"]}
        return {
            "ok": True,
            "artifact": "convergence-history",
            "path": path.name,
            "format": document["format"],
            "schema_version": document["schema_version"],
            "converged": document["converged"],
            "iterations": document["iterations"],
            "reference_momentum_residual": float(document["reference_momentum_residual"]),
            "convergence_status": document["convergence_status"],
            "convergence_reason": document["convergence_reason"],
            "history_count": len(history),
            "history": history[:history_limit],
            "history_truncated": len(history) > history_limit,
        }
    except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
        return {"ok": False, "errors": [f"invalid CFDX convergence artifact: {exc}"]}



_EXECUTION_REQUIRED = (
    "format",
    "schema_version",
    "process_exit_code",
    "converged",
    "iterations",
    "convergence_status",
    "convergence_reason",
    "artifacts",
)


def _validate_execution(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"ok": False, "errors": [f"execution artifact does not exist: {path.name}"]}
    if path.name.lower() != "execution.json":
        return {"ok": False, "errors": ["execution artifact must use the canonical execution.json filename"]}
    try:
        with path.open("r", encoding="utf-8") as stream:
            document = json.load(stream, parse_constant=_reject_nonfinite)
        if not isinstance(document, dict):
            return {"ok": False, "errors": ["execution artifact root must be an object"]}
        missing = [key for key in _EXECUTION_REQUIRED if key not in document]
        if missing:
            return {"ok": False, "errors": [f"execution artifact is missing required keys: {missing}"]}
        if document["format"] != "CFDX-EXECUTION":
            return {"ok": False, "errors": ["not a CFDX-EXECUTION artifact"]}
        if document["schema_version"] != 1:
            return {"ok": False, "errors": [f"unsupported execution schema version: {document['schema_version']}"]}
        exit_code = document["process_exit_code"]
        if not isinstance(exit_code, int) or isinstance(exit_code, bool) or exit_code < 0:
            return {"ok": False, "errors": ["process_exit_code must be a non-negative integer"]}
        if not isinstance(document["converged"], bool):
            return {"ok": False, "errors": ["converged must be a boolean"]}
        if not isinstance(document["iterations"], int) or isinstance(document["iterations"], bool) or document["iterations"] < 0:
            return {"ok": False, "errors": ["iterations must be a non-negative integer"]}
        if not isinstance(document["convergence_status"], int) or isinstance(document["convergence_status"], bool):
            return {"ok": False, "errors": ["convergence_status must be an integer"]}
        if not isinstance(document["convergence_reason"], str):
            return {"ok": False, "errors": ["convergence_reason must be a string"]}
        artifacts = document["artifacts"]
        if not isinstance(artifacts, dict):
            return {"ok": False, "errors": ["artifacts must be an object"]}
        for key in ("restart_dat", "convergence_json"):
            if key not in artifacts or not isinstance(artifacts[key], bool):
                return {"ok": False, "errors": [f"artifacts.{key} must be a boolean"]}
        return {
            "ok": True,
            "artifact": "execution-summary",
            "path": path.name,
            "format": document["format"],
            "schema_version": document["schema_version"],
            "process_exit_code": exit_code,
            "converged": document["converged"],
            "iterations": document["iterations"],
            "convergence_status": document["convergence_status"],
            "convergence_reason": document["convergence_reason"],
            "artifacts": {
                "restart_dat": artifacts["restart_dat"],
                "convergence_json": artifacts["convergence_json"],
            },
        }
    except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
        return {"ok": False, "errors": [f"invalid CFDX execution artifact: {exc}"]}



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

    @server.tool(name="case.validate", title="Validate CFDX case", annotations=annotations)
    def case_validate(case_path: str) -> dict[str, Any]:
        """Validate a root-relative .cfdx.h5 case using the canonical CFDX bundle contract."""
        try:
            path = _safe_path(runtime_root, case_path)
            if not path.is_file() or not path.name.lower().endswith(".cfdx.h5"):
                return {"ok": False, "errors": ["case must be an existing canonical .cfdx.h5 artifact"]}
            return _validate_case(path)
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

    @server.tool(name="convergence.inspect", title="Inspect CFDX convergence history", annotations=annotations)
    def convergence_inspect(convergence_path: str, history_limit: int = 100) -> dict[str, Any]:
        """Inspect canonical convergence.json execution evidence without modifying runtime state."""
        try:
            path = _safe_path(runtime_root, convergence_path)
            return _validate_convergence(path, history_limit)
        except ValueError as exc:
            return {"ok": False, "errors": [str(exc)]}

    @server.tool(name="execution.inspect", title="Inspect CFDX execution summary", annotations=annotations)
    def execution_inspect(execution_path: str) -> dict[str, Any]:
        """Inspect canonical execution.json evidence without modifying runtime state."""
        try:
            path = _safe_path(runtime_root, execution_path)
            return _validate_execution(path)
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
