"""Check physical symmetry and force integration of zero-incidence NACA Cp/Cf.

This is a consistency gate, NOT an experimental aerodynamic reference.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

FIELDS = (
    "face", "x_over_c", "y_over_c", "Cp", "Cf_tangent",
    "Cd_pressure_face", "Cd_viscous_face",
)


def analyze(path: Path, expected_pressure: float | None = None,
            expected_viscous: float | None = None,
            symmetry_tolerance: float = 1e-5) -> dict:
    if not math.isfinite(symmetry_tolerance) or symmetry_tolerance <= 0:
        raise ValueError("invalid symmetry tolerance")
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or any(f not in reader.fieldnames for f in FIELDS):
            raise ValueError("missing profile columns")
        raw = list(reader)
    if len(raw) < 32:
        raise ValueError("not enough wall-face samples")
    rows = []
    seen = set()
    for row in raw:
        face = int(row["face"])
        vals = {key: float(row[key]) for key in FIELDS[1:]}
        if face in seen or not all(math.isfinite(v) for v in vals.values()):
            raise ValueError("duplicated face or nonfinite profile")
        seen.add(face)
        rows.append(vals)
    upper = sorted((r for r in rows if r["y_over_c"] > 0),
                   key=lambda r: r["x_over_c"])
    lower = sorted((r for r in rows if r["y_over_c"] < 0),
                   key=lambda r: r["x_over_c"])
    if len(upper) != len(lower) or len(upper) < 12:
        raise ValueError("upper and lower surfaces are not paired")
    cp_err = 0.0
    cf_err = 0.0
    for a, b in zip(upper, lower):
        if abs(a["x_over_c"] - b["x_over_c"]) > 1e-8:
            raise ValueError("wall face x/c does not pair")
        if abs(a["y_over_c"] + b["y_over_c"]) > 1e-8:
            raise ValueError("wall face y/c not symmetric")
        cp_err = max(cp_err, abs(a["Cp"] - b["Cp"]))
        cf_err = max(cf_err, abs(a["Cf_tangent"] - b["Cf_tangent"]))
    if cp_err > symmetry_tolerance or cf_err > symmetry_tolerance:
        raise ValueError("Cp or Cf symmetry failed")
    cdp = sum(row["Cd_pressure_face"] for row in rows)
    cdv = sum(row["Cd_viscous_face"] for row in rows)
    for computed, expected in ((cdp, expected_pressure), (cdv, expected_viscous)):
        if expected is not None and (
            not math.isfinite(expected) or not math.isclose(
                computed, expected, rel_tol=2e-7, abs_tol=2e-9)
        ):
            raise ValueError("wall-integrated and reported force disagree")
    return {
        "status": "SYMMETRY_AND_FORCE_CONSISTENT",
        "physical_reference_validation": "NOT_PERFORMED",
        "faces": len(rows),
        "paired_faces_per_side": len(upper),
        "cp_max_pair_difference": cp_err,
        "cf_max_pair_difference": cf_err,
        "cd_pressure_face_integrated": cdp,
        "cd_viscous_face_integrated": cdv,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path)
    parser.add_argument("--cd-pressure", type=float)
    parser.add_argument("--cd-viscous", type=float)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        report = analyze(args.csv, args.cd_pressure, args.cd_viscous)
    except (ValueError, TypeError, KeyError) as exc:
        print("NACA_WALL_PROFILE_CONTRACT: FAIL:", exc)
        return 1
    content = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(content, encoding="utf-8")
    print(content, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
