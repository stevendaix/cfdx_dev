"""Compare TH-03 HDF5 geometry with the Driver-Seegmiller reference domain.

This is a geometry-only diagnostic. It never authorizes physical qualification.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import h5py
import numpy as np

H_REF = 0.0127


def gate_reference_geometry(data: dict) -> dict:
    required = ("h_m", "h_upstream", "h_downstream", "x_min_h", "x_max_h")
    if any(key not in data for key in required):
        raise ValueError("incomplete TH-03 geometry record")
    values = np.asarray([data[key] for key in required], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("nonfinite TH-03 geometry")
    checks = {
        "step_height_12p7_mm": abs(data["h_m"] - H_REF) <= 1e-10,
        "upstream_height_8H": abs(data["h_upstream"] - 8.0) <= 1e-7,
        "downstream_height_9H": abs(data["h_downstream"] - 9.0) <= 1e-7,
        "long_inlet_to_minus_130H": abs(data["x_min_h"] + 130.0) <= 1e-7,
        "outlet_at_plus_50H": abs(data["x_max_h"] - 50.0) <= 1e-7,
    }
    return {
        "checks": checks,
        "pass_reference_domain": all(checks.values()),
        "physical_validation_authorized": False,
    }


def read_case(path: Path) -> dict:
    with h5py.File(path, "r") as h:
        points = np.asarray(h["points"][:], dtype=float)
        if not np.isfinite(points).all():
            raise ValueError("nonfinite HDF5 mesh points")
        offsets = np.asarray(h["face_offsets"][:], dtype=np.int64)
        vertices = np.asarray(h["face_vertices"][:], dtype=np.int64)
        patch_offsets = np.asarray(h["patch_face_offsets"][:], dtype=np.int64)
        patch_ids = np.asarray(h["patch_face_ids"][:], dtype=np.int64)
        patch_names = [
            x.split(":")[0] for x in str(h.attrs["boundary_patches"]).split(";")
        ]
        patches = {
            name: patch_ids[patch_offsets[i] : patch_offsets[i + 1]]
            for i, name in enumerate(patch_names)
        }
        if "lowerWall" not in patches or len(patches["lowerWall"]) == 0:
            raise ValueError("TH-03 HDF5 case lacks lowerWall")
        lower_centres = []
        for face_id in patches["lowerWall"]:
            ids = vertices[offsets[face_id] : offsets[face_id + 1]]
            lower_centres.append(points[ids].mean(axis=0))
        lower = np.asarray(lower_centres)
        upstream = lower[(lower[:, 0] < 0.0) & (np.abs(lower[:, 1] - H_REF) < 1e-7)]
        downstream = lower[(lower[:, 0] > 0.0) & (np.abs(lower[:, 1]) < 1e-7)]
        if not len(upstream) or not len(downstream):
            raise ValueError("cannot locate upstream/downstream step walls")
        y_up = float(np.median(upstream[:, 1]))
        y_down = float(np.median(downstream[:, 1]))
        height = y_up - y_down
        if not height > 0.0:
            raise ValueError("invalid step height")
        top = float(points[:, 1].max())
        return {
            "case_file": str(path),
            "n_cells": int(h.attrs["n_cells"]),
            "n_faces": int(h.attrs["n_faces"]),
            "mesh_hash": str(h.attrs.get("mesh_hash", "")),
            "source_case_path": str(h.attrs.get("source_case_path", "")),
            "h_m": height,
            "h_upstream": (top - y_up) / height,
            "h_downstream": (top - y_down) / height,
            "x_min_h": float(points[:, 0].min() / height),
            "x_max_h": float(points[:, 0].max() / height),
        }


def pressure_failure_evidence(
    logfile: Path, execution: Path, relative_target: float = 1e-11
) -> dict:
    """Extract actual native pressure-correction rejection attempts; never certify."""
    source = Path(logfile).read_text(errors="replace")
    rejected = re.findall(
        r"^CFDX nonlinear iteration \\d+ rejected;[^\\n]*$", source, re.MULTILINE
    )
    pattern = re.compile(
        r"pressure-correction solve did not converge "
        r"\\(status=(\\d+), iterations=(\\d+), "
        r"residual=([0-9eE+\\-.]+), relative=([0-9eE+\\-.]+)\\)"
    )
    raw = [match.groups() for line in rejected if (match := pattern.search(line))]
    if not raw:
        raise ValueError("no native pressure-rejection evidence")
    execution_data = json.loads(Path(execution).read_text())
    if (
        execution_data.get("converged") is not False
        or execution_data.get("process_exit_code") != 2
    ):
        raise ValueError("execution record does not prove a rejected, unconverged run")
    trials = []
    for index, (status, iterations, residual, relative) in enumerate(raw, 1):
        row = {
            "trial": index,
            "status": int(status),
            "iterations": int(iterations),
            "residual_absolute": float(residual),
            "residual_relative": float(relative),
        }
        if not (
            math.isfinite(row["residual_absolute"])
            and math.isfinite(row["residual_relative"])
            and row["residual_absolute"] >= 0.0
            and row["residual_relative"] > 0.0
        ):
            raise ValueError("nonfinite or invalid native pressure residual")
        trials.append(row)
    return {
        "schema": "CFDX-PRESSURE-FAIL-REPLAY/V78",
        "target_relative": relative_target,
        "process_exit_code": 2,
        "outer_iterations_accepted": 0,
        "converged": False,
        "trials": trials,
        "min_relative": min(row["residual_relative"] for row in trials),
        "physical_validation_authorized": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_h5", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pressure-log", type=Path)
    parser.add_argument("--pressure-execution", type=Path)
    args = parser.parse_args()
    rows = []
    for path in args.case_h5:
        row = read_case(path)
        row["geometry_gate"] = gate_reference_geometry(row)
        row["geometry_status"] = (
            "REFERENCE_DOMAIN_MATCH"
            if row["geometry_gate"]["pass_reference_domain"]
            else "NOT_REFERENCE_DOMAIN"
        )
        row["physical_validation_authorized"] = False
        rows.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if bool(args.pressure_log) != bool(args.pressure_execution):
        raise ValueError("provide both --pressure-log and --pressure-execution")
    report = {
        "schema": "CFDX-TH03-GEOMETRY-COMPARABILITY/V78",
        "cases": rows,
        "geometry_only_not_physics": True,
    }
    if args.pressure_log:
        report["pressure_failure_evidence"] = pressure_failure_evidence(
            args.pressure_log, args.pressure_execution
        )
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    for row in rows:
        print(
            row["case_file"],
            row["geometry_status"],
            "upstream=",
            row["h_upstream"],
            "downstream=",
            row["h_downstream"],
            "physical_validation_authorized=0",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
