#!/usr/bin/env bash
set -euo pipefail

# Reproduce yushan-lab/openfoam-cfd-minilab's 40x40 icoFoam cavity,
# import native U/p fields into CFDX, and compare imported values to source.
# Requires a sourced OpenFOAM environment, git, Python, h5py and numpy.

repo_url="https://github.com/yushan-lab/openfoam-cfd-minilab.git"
root="${1:-$PWD/openfoam-cfd-minilab-qualification}"
mkdir -p "$root"
root="$(cd "$root" && pwd)"
source_dir="$root/source"
run_output="$root/run-output"
converted="$root/converted"
evidence="$root/evidence"

for path in "$source_dir" "$run_output" "$converted" "$evidence"; do
  if [[ -e "$path" ]]; then
    echo "ERROR: output path already exists: $path (choose a fresh output directory)" >&2
    exit 2
  fi
done

if ! command -v foamVersion >/dev/null 2>&1; then
  echo "ERROR: source an OpenFOAM environment first (foamVersion not found)." >&2
  exit 2
fi
if ! command -v git >/dev/null 2>&1; then
  echo "ERROR: git is required." >&2
  exit 2
fi

repo_commit="$(git ls-remote "$repo_url" refs/heads/main | awk 'NR==1 {print $1}')"
if [[ -z "$repo_commit" ]]; then
  echo "ERROR: unable to resolve the public repository main commit." >&2
  exit 2
fi
git clone --no-checkout "$repo_url" "$source_dir"
git -C "$source_dir" checkout --detach "$repo_commit"

case_dir="$source_dir/cases/lid_driven_cavity"
if [[ ! -f "$source_dir/scripts/run_cavity.sh" || ! -f "$case_dir/system/controlDict" ]]; then
  echo "ERROR: expected cavity case or run script is missing from repository." >&2
  exit 2
fi

mkdir -p "$run_output"
{
  echo "Repository: $repo_url"
  echo "Repository commit: $repo_commit"
  echo "OpenFOAM version: $(foamVersion)"
  echo "Case: cases/lid_driven_cavity"
  echo "Mesh resolution: 40x40x1"
  echo "Run command: MESH_RESOLUTION=40 OUTPUT_ROOT=$run_output SAVE_FINAL_FIELDS=1 bash scripts/run_cavity.sh"
} | tee "$root/provenance.txt"

(
  cd "$source_dir"
  MESH_RESOLUTION=40 \
    OUTPUT_ROOT="$run_output" \
    SAVE_FINAL_FIELDS=1 \
    bash scripts/run_cavity.sh
)

latest_time="$(
  find "$case_dir" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' \
    | awk '/^[0-9]+([.][0-9]+)?$/ {print}' | sort -g | tail -n 1
)"
if [[ -z "$latest_time" || ! -s "$case_dir/$latest_time/U" || ! -s "$case_dir/$latest_time/p" ]]; then
  echo "ERROR: run did not produce a numeric final time containing both U and p." >&2
  exit 3
fi
echo "Selected time: $latest_time" | tee -a "$root/provenance.txt"

package_root="$PWD/src/cfdx/python"
if [[ ! -d "$package_root/cfdx/io" ]]; then
  echo "ERROR: run this script from the CFDX repository root; package path not found: $package_root" >&2
  exit 2
fi
export PYTHONPATH="$package_root${PYTHONPATH:+:$PYTHONPATH}"
python -m cfdx.io.cli convert "$case_dir" --results --output "$converted"

mkdir -p "$evidence"
python - "$case_dir" "$latest_time" "$converted/case.cfdx.h5" "$evidence" "$repo_commit" <<'PY'
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np

case_dir = Path(sys.argv[1])
time_name = sys.argv[2]
h5_path = Path(sys.argv[3])
evidence_dir = Path(sys.argv[4])
repo_commit = sys.argv[5]
time_dir = case_dir / time_name


def read_field(path: Path, kind: str, components: int) -> tuple[np.ndarray, str]:
    text = path.read_text(encoding="utf-8")
    text = re.sub(r"/\*.*?\*/|//[^\n]*", "", text, flags=re.S)
    dimensions_match = re.search(r"\bdimensions\s*\[([^\]]+)\]\s*;", text)
    if not dimensions_match:
        raise SystemExit(f"{path}: missing dimensions declaration")
    dimensions = " ".join(dimensions_match.group(1).split())
    if kind == "vector":
        pattern = r"\binternalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;"
        matches = re.search(pattern, text, flags=re.S)
        if not matches:
            raise SystemExit(f"{path}: expected ASCII nonuniform List<vector> internalField")
        count = int(matches.group(1))
        values = np.asarray(
            [tuple(float(v) for v in row) for row in re.findall(
                r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)", matches.group(2)
            )],
            dtype=np.float64,
        )
        if values.shape != (count, components):
            raise SystemExit(f"{path}: declared {count} vectors but parsed {values.shape}")
    else:
        pattern = r"\binternalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;"
        matches = re.search(pattern, text, flags=re.S)
        if not matches:
            raise SystemExit(f"{path}: expected ASCII nonuniform List<scalar> internalField")
        count = int(matches.group(1))
        values = np.asarray(
            [float(v) for v in re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", matches.group(2))],
            dtype=np.float64,
        )
        if values.shape != (count,):
            raise SystemExit(f"{path}: declared {count} scalars but parsed {values.shape}")
    if not np.isfinite(values).all():
        raise SystemExit(f"{path}: source field contains non-finite values")
    return values, dimensions


source_u, dims_u = read_field(time_dir / "U", "vector", 3)
source_p, dims_p = read_field(time_dir / "p", "scalar", 1)

with h5py.File(h5_path, "r") as h5:
    n_cells = int(h5.attrs["n_cells"])
    if source_u.shape != (n_cells, 3) or source_p.shape != (n_cells,):
        raise SystemExit(
            f"Source field sizes do not match imported mesh: cells={n_cells}, "
            f"U={source_u.shape}, p={source_p.shape}"
        )
    imported_u = np.asarray(h5["fields/vector/U"][:])
    imported_p = np.asarray(h5["fields/scalar/p"][:])
    owner = np.asarray(h5["owner"][:])
    neighbour = np.asarray(h5["neighbour"][:])
    n_faces = int(h5.attrs["n_faces"])
    if owner.shape != (n_faces,) or neighbour.shape != (n_faces,):
        raise SystemExit("Imported owner/neighbour sizes do not match n_faces")
    if owner.size and (owner.min() < 0 or owner.max() >= n_cells):
        raise SystemExit("Imported owner indices are out of range")
    internal_neighbour = neighbour[neighbour >= 0]
    if internal_neighbour.size and internal_neighbour.max() >= n_cells:
        raise SystemExit("Imported neighbour indices are out of range")
    np.testing.assert_array_equal(imported_u, source_u, err_msg="Imported U differs from source U")
    np.testing.assert_array_equal(imported_p, source_p, err_msg="Imported p differs from source p")

    report = {
        "qualification": "yushan_lab_openfoam_cfd_minilab_icofoam_cavity",
        "status": "exact_source_to_import_field_values_passed",
        "repository": "https://github.com/yushan-lab/openfoam-cfd-minilab",
        "repository_commit": repo_commit,
        "selected_time": time_name,
        "mesh": {
            "n_points": int(h5.attrs["n_points"]),
            "n_faces": n_faces,
            "n_cells": n_cells,
        },
        "topology": {
            "owner_entries": int(owner.size),
            "neighbour_entries": int(neighbour.size),
            "owner_indices_in_range": True,
            "internal_neighbour_indices_in_range": True,
        },
        "fields": {
            "U": {
                "dimensions": dims_u,
                "expected_dimensions": "0 1 -1 0 0 0 0",
                "shape": list(imported_u.shape),
                "exact_array_equality": True,
                "max_absolute_difference": float(np.max(np.abs(imported_u - source_u))),
            },
            "p": {
                "dimensions": dims_p,
                "expected_dimensions": "0 2 -2 0 0 0 0",
                "shape": list(imported_p.shape),
                "exact_array_equality": True,
                "max_absolute_difference": float(np.max(np.abs(imported_p - source_p))),
                "note": "OpenFOAM p is kinematic pressure; this does not authorize mapping to CFDX pressure in Pa.",
            },
        },
        "limitations": [
            "This compares cell-centred ASCII internalField values and basic topology arrays.",
            "It does not qualify boundary-condition mapping or numerical equivalence of solvers.",
            "It does not execute the CFDX production solver.",
            "Only the 40x40x1 base cavity case is imported; the four-grid study remains separate validation evidence.",
        ],
    }
    if dims_u != "0 1 -1 0 0 0 0":
        raise SystemExit(f"Unexpected OpenFOAM U dimensions: {dims_u}")
    if dims_p != "0 2 -2 0 0 0 0":
        raise SystemExit(f"Unexpected OpenFOAM p dimensions: {dims_p}")

(evidence_dir / "source_import_comparison.json").write_text(
    json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
print(json.dumps(report, indent=2, sort_keys=True))
PY

echo "Exact source-to-import comparison passed."
echo "Evidence: $evidence/source_import_comparison.json"
echo "This does not qualify production-solver initialization or OpenFOAM/CFDX numerical equivalence."
