"""Fail-closed, run-scoped reading of CFDX's persisted solver evidence.

`execution.json` is the overall solver verdict (including energy/turbulence
stages), whereas `convergence.json` holds the incompressible flow history.
Neither a clean process exit, a shrinking residual, nor an old file qualifies
an execution. This module never invents numerical tolerances.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

MAX_EVIDENCE_BYTES = 64 * 1024 * 1024
MetricValue = tuple[str, float]
Signature = tuple[int, int, int, int]


class EvidenceError(ValueError):
    """The emitted evidence is absent, stale, malformed or contradictory."""


def file_signature(path: Path) -> Signature | None:
    """Read file identity without following symlinks or opening its contents."""
    if path.is_symlink():
        raise EvidenceError(f"symbolic-link evidence is not supported: {path.name}")
    try:
        stat = path.stat()
    except FileNotFoundError:
        return None
    if not path.is_file():
        raise EvidenceError(f"evidence path is not a regular file: {path.name}")
    return (
        int(stat.st_dev),
        int(stat.st_ino),
        int(stat.st_size),
        int(stat.st_mtime_ns),
    )


def baseline_signatures(directory: Path | None) -> dict[str, Signature | None]:
    """Snapshot pre-run artifacts. A subsequent artifact must differ."""
    if directory is None:
        return {}
    root = Path(directory)
    return {
        name: file_signature(root / name)
        for name in ("execution.json", "convergence.json")
    }


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise EvidenceError(f"duplicate evidence JSON key: {key}")
        result[key] = value
    return result


def _reject_constant(value: str) -> object:
    raise EvidenceError(f"non-finite JSON numeric literal: {value}")


def _read_json(
    path: Path, baseline: Mapping[str, Signature | None], max_bytes: int
) -> tuple[dict, str]:
    before = file_signature(path)
    if before is None:
        raise EvidenceError(f"{path.name} is missing")
    if before == baseline.get(path.name):
        raise EvidenceError(f"{path.name} belongs to a previous run (unchanged)")
    if before[2] <= 0 or before[2] > max_bytes:
        raise EvidenceError(
            f"{path.name} exceeds the bounded evidence read ({max_bytes} bytes)"
        )
    raw = path.read_bytes()
    if len(raw) != before[2] or file_signature(path) != before:
        raise EvidenceError(f"{path.name} changed while it was inspected")
    try:
        payload = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=_unique_object,
            parse_constant=_reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise EvidenceError(f"{path.name} is not valid UTF-8 JSON: {exc}") from exc
    if not isinstance(payload, dict):
        raise EvidenceError(f"{path.name} must contain a JSON object")
    return payload, hashlib.sha256(raw).hexdigest()


def _nonnegative_int(value: object, key: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise EvidenceError(f"{key} must be a non-negative integer")
    return value


def _format_check(data: dict, fmt: str, filename: str) -> None:
    if data.get("format") != fmt or data.get("schema_version") != 1:
        raise EvidenceError(f"unsupported {filename} format or schema_version")
    if not isinstance(data.get("converged"), bool):
        raise EvidenceError(f"{filename} must have a boolean converged field")
    _nonnegative_int(data.get("iterations"), f"{filename}.iterations")


@dataclass(frozen=True, slots=True)
class ConvergenceEvidence:
    """Observed solver evidence, not a grid/physics qualification certificate."""

    overall_converged: bool
    flow_converged: bool
    iterations: int
    reason: str
    convergence_status: int
    final_metrics: tuple[MetricValue, ...]
    flow_history_count: int
    execution_sha256: str
    convergence_sha256: str
    output_directory: str
    source: str = "CFDX-EXECUTION v1 + CFDX-CONVERGENCE v1"

    def to_dict(self) -> dict[str, object]:
        return {
            "overall_converged": self.overall_converged,
            "flow_converged": self.flow_converged,
            "iterations": self.iterations,
            "reason": self.reason,
            "convergence_status": self.convergence_status,
            "final_metrics": dict(self.final_metrics),
            "flow_history_count": self.flow_history_count,
            "execution_sha256": self.execution_sha256,
            "convergence_sha256": self.convergence_sha256,
            "output_directory": self.output_directory,
            "source": self.source,
            "independent_physical_qualification": "not_evaluated",
        }


_V8_CERTIFICATE_KEYS = {
    "evaluated",
    "algorithm",
    "schur",
    "momentum_global_relative",
    "momentum_max_active_component_relative",
    "momentum_component_relative",
    "continuity",
    "k",
    "omega",
    "eddy_viscosity_closure",
    "merit",
    "sst",
    "sst_transport_committed",
    "bounded",
    "all_finite",
    "converged",
}
_V8_NUMERIC_KEYS = (
    "momentum_global_relative",
    "momentum_max_active_component_relative",
    "continuity",
    "k",
    "omega",
    "eddy_viscosity_closure",
    "merit",
)
_V8_BOOLEAN_KEYS = (
    "evaluated",
    "sst",
    "sst_transport_committed",
    "bounded",
    "all_finite",
    "converged",
)
_V8_ALGORITHMS = {"simple", "simplec", "piso", "pimple", "fractional_step", "coupled"}


def _optional_nonnegative_finite(value: object, key: str) -> None:
    if value is None:
        return
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value < 0
    ):
        raise EvidenceError(f"{key} must be null or a finite non-negative number")


def _validate_global_nonlinear_v8(value: object, history_index: int) -> None:
    """Validate the one nested certificate object emitted by the V8 writer."""
    key = f"history[{history_index}].global_nonlinear_v8"
    if not isinstance(value, dict):
        raise EvidenceError(f"{key} must be an object")
    if set(value) != _V8_CERTIFICATE_KEYS:
        missing = sorted(_V8_CERTIFICATE_KEYS - set(value))
        extra = sorted(set(value) - _V8_CERTIFICATE_KEYS)
        raise EvidenceError(f"{key} has invalid keys (missing={missing}, extra={extra})")
    for name in _V8_BOOLEAN_KEYS:
        if not isinstance(value[name], bool):
            raise EvidenceError(f"{key}.{name} must be a boolean")
    algorithm = value["algorithm"]
    if not isinstance(algorithm, str) or algorithm not in _V8_ALGORITHMS:
        raise EvidenceError(f"{key}.algorithm is unsupported")
    schur = value["schur"]
    if algorithm == "coupled":
        if not isinstance(schur, str) or not schur:
            raise EvidenceError(f"{key}.schur must name the resolved coupled method")
    elif schur is not None:
        raise EvidenceError(f"{key}.schur must be null for a non-coupled algorithm")
    for name in _V8_NUMERIC_KEYS:
        _optional_nonnegative_finite(value[name], f"{key}.{name}")
    components = value["momentum_component_relative"]
    if not isinstance(components, list) or len(components) != 3:
        raise EvidenceError(f"{key}.momentum_component_relative must have 3 entries")
    for component in components:
        _optional_nonnegative_finite(
            component, f"{key}.momentum_component_relative"
        )
    if value["sst_transport_committed"] and not value["sst"]:
        raise EvidenceError(f"{key}.sst_transport_committed requires sst")
    if value["converged"] and not (
        value["evaluated"] and value["all_finite"] and value["bounded"]
    ):
        raise EvidenceError(f"{key}.converged requires evaluated, finite, bounded state")


def load_run_evidence(
    output_dir: str | Path,
    *,
    process_exit_code: int,
    baseline: Mapping[str, Signature | None] | None = None,
    max_bytes: int = MAX_EVIDENCE_BYTES,
) -> ConvergenceEvidence:
    """Validate both native solver artifacts, including the overall stage gate.

    `baseline` is mandatory when verifying an execution; callers must capture
    it *before* starting the child process. This rejects stale work directories.
    """
    if baseline is None or set(baseline) != {"execution.json", "convergence.json"}:
        raise EvidenceError("pre-run artifact baseline is required")
    if (
        isinstance(max_bytes, bool)
        or not isinstance(max_bytes, int)
        or max_bytes < 1024
    ):
        raise ValueError("max_bytes must be at least 1024")
    if process_exit_code not in (0, 1, 2):
        raise EvidenceError("not a recognized CFDX production-solver exit code")
    root = Path(output_dir)
    execution, execution_hash = _read_json(root / "execution.json", baseline, max_bytes)
    convergence, convergence_hash = _read_json(
        root / "convergence.json", baseline, max_bytes
    )
    _format_check(execution, "CFDX-EXECUTION", "execution.json")
    _format_check(convergence, "CFDX-CONVERGENCE", "convergence.json")
    if execution.get("process_exit_code") != process_exit_code:
        raise EvidenceError(
            "execution.json process_exit_code disagrees with the process"
        )
    if execution["iterations"] != convergence["iterations"]:
        raise EvidenceError("execution and convergence iteration counts disagree")
    declared = execution.get("artifacts")
    if not isinstance(declared, dict) or declared.get("convergence_json") is not True:
        raise EvidenceError("execution.json does not certify its convergence artifact")
    overall, flow = execution["converged"], convergence["converged"]
    if overall and not flow:
        raise EvidenceError("overall solver cannot converge if flow did not converge")
    if overall != (process_exit_code == 0):
        raise EvidenceError(
            "overall convergence flag disagrees with production exit semantics"
        )
    reason = execution.get("convergence_reason")
    status = execution.get("convergence_status")
    if not isinstance(reason, str) or len(reason) > 4096:
        raise EvidenceError("invalid convergence reason")
    if not isinstance(status, int) or isinstance(status, bool):
        raise EvidenceError("invalid convergence status")
    rows = convergence.get("history")
    if not isinstance(rows, list):
        raise EvidenceError("flow history is missing")
    if overall and not rows:
        raise EvidenceError("converged result has no final flow metrics")
    last_iteration = -1
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            raise EvidenceError(f"history[{index}] is not an object")
        iteration = _nonnegative_int(
            row.get("iteration"), f"history[{index}].iteration"
        )
        if iteration <= last_iteration or iteration > convergence["iterations"]:
            raise EvidenceError(
                "history iterations are not strictly increasing or exceed run count"
            )
        for metric, value in row.items():
            if metric == "iteration" or value is None or isinstance(value, (bool, str)):
                continue
            if metric == "global_nonlinear_v8":
                _validate_global_nonlinear_v8(value, index)
            elif isinstance(value, (int, float)):
                if not math.isfinite(value):
                    raise EvidenceError(f"non-finite history[{index}].{metric}")
            elif isinstance(value, list):
                if any(
                    not isinstance(item, (int, float))
                    or isinstance(item, bool)
                    or not math.isfinite(item)
                    for item in value
                    if item is not None
                ):
                    raise EvidenceError(f"invalid numeric history[{index}].{metric}")
            else:
                raise EvidenceError(f"unsupported history[{index}].{metric} value")
        last_iteration = iteration
    if overall and last_iteration != convergence["iterations"]:
        raise EvidenceError(
            "converged history does not reach its declared final iteration"
        )
    latest = rows[-1] if rows else {}
    metrics: list[MetricValue] = []
    for name, value in latest.items():
        if name == "iteration" or isinstance(value, bool):
            continue
        if isinstance(value, (float, int)):
            number = float(value)
            if not math.isfinite(number):
                raise EvidenceError(f"non-finite final flow metric {name}")
            metrics.append((name, number))
    return ConvergenceEvidence(
        overall_converged=overall,
        flow_converged=flow,
        iterations=execution["iterations"],
        reason=reason,
        convergence_status=status,
        final_metrics=tuple(sorted(metrics)),
        flow_history_count=len(rows),
        execution_sha256=execution_hash,
        convergence_sha256=convergence_hash,
        output_directory=str(root.resolve()),
    )
