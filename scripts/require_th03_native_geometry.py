"""Fail-closed TH-03 native OpenFOAM geometry provenance gate.

This gate proves that the HDF5 mesh is the exact topology/geometry imported
from the referenced OpenFOAM constant/polyMesh. Matching bounding boxes alone
are insufficient. A pass is geometry provenance only, never physical validation.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "src/cfdx/python") not in sys.path:
    sys.path.insert(0, str(ROOT / "src/cfdx/python"))

from cfdx.io.hdf5_writer import _integrity_hashes

from scripts.generate_th03_openfoam_blockmesh import H, load_openfoam_polymesh


def _text(value: object) -> str:
    return value.decode() if isinstance(value, bytes) else str(value)


def native_case_path(source_case_path: object) -> Path:
    source = _text(source_case_path)
    prefix = "OpenFOAM polyMesh, "
    if not source.startswith(prefix):
        raise ValueError(
            "TH-03 physical initialisation requires native OpenFOAM polyMesh; "
            f"case provenance is {source!r}. Reconstructed six-block meshes "
            "are diagnostic-only."
        )
    raw_path = source[len(prefix):].strip()
    if not raw_path:
        raise ValueError("native OpenFOAM case path is empty")
    case = Path(raw_path).expanduser().resolve(strict=True)
    polymesh = case / "constant" / "polyMesh"
    required = ("points", "faces", "owner", "neighbour", "boundary")
    missing = [name for name in required if not (polymesh / name).is_file()]
    if missing:
        raise ValueError(f"native OpenFOAM polyMesh is incomplete: {missing}")
    return case


def reference_domain(mesh: dict) -> dict:
    points = np.asarray(mesh["points"], dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError("native polyMesh has invalid point coordinates")
    lower = next(
        (p for p in mesh["boundary_patch_assignment"] if p["name"] == "lowerWall"),
        None,
    )
    if lower is None or not lower["face_ids"]:
        raise ValueError("native polyMesh has no lowerWall patch faces")
    offsets = np.asarray(mesh["face_offsets"], dtype=np.int64)
    vertices = np.asarray(mesh["face_vertices"], dtype=np.int64)
    centres = []
    for face_id in lower["face_ids"]:
        ids = vertices[offsets[int(face_id)]:offsets[int(face_id) + 1]]
        centres.append(points[ids].mean(axis=0))
    centres = np.asarray(centres)
    upstream = centres[(centres[:, 0] < 0.0) & (np.abs(centres[:, 1] - H) <= 1e-7 * H)]
    downstream = centres[(centres[:, 0] > 0.0) & (np.abs(centres[:, 1]) <= 1e-7 * H)]
    if not len(upstream) or not len(downstream):
        raise ValueError("cannot identify upstream and downstream lower-wall segments")
    y_up = float(np.median(upstream[:, 1]))
    y_down = float(np.median(downstream[:, 1]))
    step = y_up - y_down
    if not step > 0.0:
        raise ValueError("invalid native step height")
    x_min, x_max = float(points[:, 0].min() / step), float(points[:, 0].max() / step)
    y_top = float(points[:, 1].max())
    checks = {
        "step_height_12p7_mm": abs(step - 0.0127) <= 1e-10,
        "upstream_height_8H": abs((y_top - y_up) / step - 8.0) <= 1e-7,
        "downstream_height_9H": abs((y_top - y_down) / step - 9.0) <= 1e-7,
        "inlet_at_minus_130H": abs(x_min + 130.0) <= 1e-7,
        "outlet_at_plus_50H": abs(x_max - 50.0) <= 1e-7,
    }
    if not all(checks.values()):
        raise ValueError(f"native OpenFOAM mesh does not match TH-03 reference domain: {checks}")
    return {
        "step_height_m": step,
        "upstream_height_over_H": (y_top - y_up) / step,
        "downstream_height_over_H": (y_top - y_down) / step,
        "x_min_over_H": x_min,
        "x_max_over_H": x_max,
        "checks": checks,
    }


def require_native_geometry(case_h5: Path) -> dict:
    case_h5 = Path(case_h5).resolve(strict=True)
    with h5py.File(case_h5, "r") as h:
        source = native_case_path(h.attrs.get("source_case_path", ""))
        if _text(h.attrs.get("source_format", "")) != "openfoam":
            raise ValueError("HDF5 source_format is not OpenFOAM")
        if _text(h.attrs.get("source_version", "")) != "v2606":
            raise ValueError("HDF5 source_version is not the pinned OpenFOAM v2606")
        stored = {
            key: _text(h.attrs.get(key, ""))
            for key in ("topology_hash", "geometry_hash", "mesh_hash")
        }
        if any(not value for value in stored.values()):
            raise ValueError("HDF5 case is missing native mesh integrity hashes")
    native = load_openfoam_polymesh(source)
    topology, geometry, mesh_hash = _integrity_hashes(native)
    actual = {
        "topology_hash": topology,
        "geometry_hash": geometry,
        "mesh_hash": mesh_hash,
    }
    if actual != stored:
        raise ValueError(
            "HDF5 mesh differs from the referenced native OpenFOAM polyMesh: "
            f"stored={stored}, reimported={actual}"
        )
    domain = reference_domain(native)
    return {
        "schema": "CFDX-TH03-NATIVE-GEOMETRY/1",
        "case_h5": str(case_h5),
        "openfoam_case": str(source),
        "geometry_provenance": "NATIVE_OPENFOAM_POLYMESH_REIMPORTED_AND_HASH_MATCHED",
        "hashes": actual,
        "reference_domain": domain,
        "physical_qualification": "NOT_AUTHORIZED_BY_GEOMETRY_GATE",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-h5", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = require_native_geometry(args.case_h5)
    rendered = json.dumps(result, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
