#!/usr/bin/env bash
set -euo pipefail

# End-to-end OpenFOAM icoFoam/cavity import qualification.
# Requires sourced OpenFOAM environment, cfdx-io, Python, h5py and numpy.
# Does not claim OpenFOAM/CFDX numerical equivalence.

if ! command -v foamVersion >/dev/null 2>&1; then
  echo "ERROR: source the OpenFOAM environment first (foamVersion not found)." >&2
  exit 2
fi
if ! command -v cfdx-io >/dev/null 2>&1; then
  echo "ERROR: cfdx-io is not on PATH; install/build CFDX Python I/O first." >&2
  exit 2
fi

: "${FOAM_TUTORIALS:?FOAM_TUTORIALS is unset; source the OpenFOAM environment}"
tutorial="$(find "$FOAM_TUTORIALS" -type d -path '*/icoFoam/cavity/cavity' -print -quit)"
if [[ -z "$tutorial" ]]; then
  echo "ERROR: icoFoam/cavity/cavity tutorial not found under $FOAM_TUTORIALS" >&2
  exit 2
fi

root="${1:-$PWD/openfoam-icofoam-qualification}"
mkdir -p "$root"
root="$(cd "$root" && pwd)"
case_dir="$root/openfoam-cavity"
converted="$root/converted"
evidence="$root/evidence"
if [[ -e "$case_dir" || -e "$converted" || -e "$evidence" ]]; then
  echo "ERROR: output paths already exist under $root; choose a fresh directory." >&2
  exit 2
fi

printf '%s\n' "OpenFOAM version: $(foamVersion)" | tee "$root/provenance.txt"
printf '%s\n' "Tutorial source: $tutorial" | tee -a "$root/provenance.txt"
printf '%s\n' "Tutorial root: $FOAM_TUTORIALS" | tee -a "$root/provenance.txt"
cp -a "$tutorial" "$case_dir"

cd "$case_dir"
if [[ -x ./Allrun ]]; then
  echo "Reviewing and executing tutorial Allrun script."
  ./Allrun
else
  echo "Running standard icoFoam cavity tutorial commands."
  blockMesh
  icoFoam
fi

latest_time="$(find . -mindepth 1 -maxdepth 1 -type d -printf '%f\n' \
  | awk '/^[0-9]+([.][0-9]+)?$/ {print}' | sort -g | tail -n 1)"
if [[ -z "$latest_time" || ! -f "$latest_time/U" || ! -f "$latest_time/p" ]]; then
  echo "ERROR: no numeric time directory containing both U and p was found." >&2
  exit 3
fi
printf '%s\n' "Selected time: $latest_time" | tee -a "$root/provenance.txt"
cd "$root"

cfdx-io convert "$case_dir" --results --output "$converted"

python - "$converted/case.cfdx.h5" "$evidence" "$latest_time" <<'PY'
import json
import math
import sys
from pathlib import Path

import h5py
import numpy as np

case_path = Path(sys.argv[1])
evidence_dir = Path(sys.argv[2])
selected_time = sys.argv[3]
evidence_dir.mkdir(parents=True, exist_ok=False)

with h5py.File(case_path, "r") as h5:
    counts = {
        "n_points": int(h5.attrs["n_points"]),
        "n_faces": int(h5.attrs["n_faces"]),
        "n_cells": int(h5.attrs["n_cells"]),
    }
    if min(counts.values()) <= 0:
        raise SystemExit(f"Invalid non-positive mesh counts: {counts}")
    for name in ("owner", "neighbour"):
        if name not in h5:
            raise SystemExit(f"Missing required topology dataset {name}")
    owner = np.asarray(h5["owner"][:])
    neighbour = np.asarray(h5["neighbour"][:])
    if owner.size != counts["n_faces"] or neighbour.size != counts["n_faces"]:
        raise SystemExit(
            f"Topology length mismatch: faces={counts['n_faces']}, "
            f"owner={owner.size}, neighbour={neighbour.size}"
        )
    if owner.size and (owner.min() < 0 or owner.max() >= counts["n_cells"]):
        raise SystemExit("Owner array contains out-of-range cell indices")
    internal = neighbour[neighbour >= 0]
    if internal.size and (internal.min() < 0 or internal.max() >= counts["n_cells"]):
        raise SystemExit("Neighbour array contains out-of-range cell indices")

    fields = {}
    for group_name, field_name in (("vector", "U"), ("scalar", "p")):
        path = f"fields/{group_name}/{field_name}"
        if path not in h5:
            raise SystemExit(f"Missing imported field {path}")
        values = np.asarray(h5[path][:])
        expected_shape = (counts["n_cells"], 3) if group_name == "vector" else (counts["n_cells"],)
        if values.shape != expected_shape:
            raise SystemExit(f"{path}: expected shape {expected_shape}, got {values.shape}")
        if not np.isfinite(values).all():
            raise SystemExit(f"{path} contains non-finite values")
        fields[field_name] = {
            "dataset": path,
            "shape": list(values.shape),
            "finite": True,
            "min": float(values.min()) if values.size else None,
            "max": float(values.max()) if values.size else None,
        }

    patches = h5.attrs.get("boundary_patches", "")
    if isinstance(patches, bytes):
        patches = patches.decode("utf-8", errors="replace")
    evidence = {
        "qualification": "openfoam_icofoam_cavity_reader_import",
        "status": "import_structure_and_finiteness_passed",
        "selected_time": selected_time,
        "mesh": counts,
        "topology": {
            "owner_entries": int(owner.size),
            "neighbour_entries": int(neighbour.size),
            "owner_indices_in_range": True,
            "neighbour_internal_indices_in_range": True,
        },
        "fields": fields,
        "boundary_patches_metadata": str(patches),
        "limitations": [
            "This report checks imported structure and finite values, not exact source-to-import equality.",
            "CFDX production solver initialization is not run by this script.",
            "OpenFOAM kinematic p is not mapped to CFDX pressure in Pa.",
            "This is single-phase icoFoam qualification, not VOF qualification.",
        ],
    }
(evidence_dir / "reader_import_evidence.json").write_text(
    json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
print(json.dumps(evidence, indent=2, sort_keys=True))
PY

echo "Reader import checks passed. Evidence: $evidence/reader_import_evidence.json"
echo "This is not a solver-backed qualification; run the dedicated CFDX solver test separately."
