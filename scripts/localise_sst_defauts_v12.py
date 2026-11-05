"""Localise les pires residus SST retournes, sans qualification physique.

Lit la sortie HDF5 d'un cas CFDX et le journal du solveur V12.
Le point de cellule est une moyenne de sommets (diagnostic seulement).
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import h5py
import numpy as np

TRACE = re.compile(r"^SST_TRANSPORT_DIAGNOSTIC\s+(.*)$", re.MULTILINE)
FIELD = re.compile(r"([a-zA-Z_]+)=([^\s]+)")
TH03_STEP_HEIGHT_M = 0.0127  # Cas Driver-Seegmiller/OpenFOAM; ne pas inferer depuis la profondeur extrudee.


def localisation(mesh: Path, log: Path) -> dict:
    raw = log.read_text(encoding="utf-8", errors="replace")
    matches = TRACE.findall(raw)
    if not matches:
        raise ValueError("missing SST_TRANSPORT_DIAGNOSTIC: no localization")
    records = [dict(FIELD.findall(match)) for match in matches]
    qualified_records = [
        row
        for row in records
        if "k_local_worst_cell" in row
        and "omega_local_worst_cell" in row
        and "iteration" in row
    ]
    # A late retry may restart at iteration=1. Prefer the most advanced
    # inspected SST Picard state, never mislabel a fresh retry as iter=100.
    # Multiple outer retries can all reach the same Picard iteration budget.
    # On a tie, use the LAST attempt, not a stale earlier retry.
    selected = max(
        enumerate(qualified_records),
        key=lambda pair: (int(pair[1]["iteration"]), pair[0]),
        default=None,
    )
    if selected is None:
        raise ValueError("SST trace lacks local worst-cell indices")
    record = selected[1]
    with h5py.File(mesh, "r") as file:
        xyz = np.asarray(file["points"], dtype=float)
        face_offsets = np.asarray(file["face_offsets"], dtype=np.int64)
        face_vertices = np.asarray(file["face_vertices"], dtype=np.int64)
        cell_offsets = np.asarray(file["cell_offsets"], dtype=np.int64)
        cell_faces = np.asarray(file["cell_faces"], dtype=np.int64)
        patch_offsets = np.asarray(file["patch_face_offsets"], dtype=np.int64)
        patch_fids = np.asarray(file["patch_face_ids"], dtype=np.int64)
        headers = str(file.attrs["boundary_patches"]).split(";")
        patch_names = [v.split(":")[0] for v in headers]
        patch_codes = [int(v.split(":")[3]) for v in headers]
        h = TH03_STEP_HEIGHT_M
        if not math.isfinite(h) or h <= 0:
            raise ValueError("non-positive span used to normalize position")

        def face_point(face_id: int) -> np.ndarray:
            a, b = face_offsets[face_id : face_id + 2]
            return np.mean(xyz[face_vertices[a:b]], axis=0)

        wall_faces = [
            (name, face_point(int(fid)))
            for p, (name, code) in enumerate(zip(patch_names, patch_codes))
            if code == 0
            for fid in patch_fids[patch_offsets[p] : patch_offsets[p + 1]]
        ]
        if not wall_faces:
            raise ValueError("no physical WALL faces in HDF5")
        results = {}
        for variable in ("k", "omega"):
            index = int(record[f"{variable}_local_worst_cell"])
            if not 0 <= index < len(cell_offsets) - 1:
                raise ValueError(f"{variable}: invalid cell index {index}")
            faces = cell_faces[cell_offsets[index] : cell_offsets[index + 1]]
            vertices = sorted(
                {
                    int(v)
                    for fid in faces
                    for v in face_vertices[face_offsets[fid] : face_offsets[fid + 1]]
                }
            )
            point = np.mean(xyz[vertices], axis=0)
            name, distance = min(
                (
                    (name, float(np.linalg.norm(point - center)))
                    for name, center in wall_faces
                ),
                key=lambda pair: pair[1],
            )
            defect = float(record[f"{variable}_local_relative"])
            if not math.isfinite(defect) or defect < 0:
                raise ValueError("nonfinite local defect")
            results[variable] = {
                "worst_cell": index,
                "defect_local_relative": defect,
                "point_over_h": [float(v / h) for v in point],
                "nearest_wall_patch": name,
                "nearest_wall_face_center_distance_over_h": distance / h,
            }
    return {
        "status": "DIAGNOSTIC_ONLY_NOT_QUALIFIED",
        "iteration": int(record["iteration"]),
        "h_m": h,
        "fields": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = localisation(args.mesh, args.log)
    payload = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
