"""Independent surface-integral audit from raw NACA face geometry and stresses.

C++ exports raw p, tau*n, face area/normals; Python integrates F=(p*n-tau*n)dA
and moment without using the C++ accumulated forces. This detects postprocess
normalization/sign/omitted-face mistakes. It does not independently re-solve the
fluid equations or prove the wall-gradient reconstruction is accurate.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path

COLUMNS = (
    "face", "x", "y", "z", "nx", "ny", "nz", "area", "pressure",
    "tau_n_x", "tau_n_y", "tau_n_z", "dynamic_pressure", "reference_area",
    "reference_chord", "x_ref", "y_ref",
)


def read_rows(path: Path) -> list[dict[str, float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != COLUMNS:
            raise ValueError("wrong NACA raw face schema")
        rows = []
        seen = set()
        for line, raw in enumerate(reader, 2):
            if None in raw or any(v is None or not v.strip() for v in raw.values()):
                raise ValueError(f"truncated or empty raw face row {line}")
            try:
                idx = int(raw["face"])
                values = {k: float(raw[k]) for k in COLUMNS[1:]}
            except ValueError as error:
                raise ValueError(f"invalid face row {line}: {error}") from error
            if idx < 0 or idx in seen:
                raise ValueError("negative or duplicate wall face ID")
            seen.add(idx)
            if not all(math.isfinite(v) for v in values.values()):
                raise ValueError("nonfinite face data")
            if values["area"] <= 0:
                raise ValueError("nonpositive face area")
            n2 = sum(values[f"n{axis}"] ** 2 for axis in "xyz")
            if not math.isclose(n2, 1.0, abs_tol=1e-10, rel_tol=1e-10):
                raise ValueError("normal is not a unit vector")
            if any(values[key] <= 0 for key in (
                "dynamic_pressure", "reference_area", "reference_chord"
            )):
                raise ValueError("invalid force normalization")
            rows.append(values)
    if len(rows) < 32:
        raise ValueError("incomplete wall surface (fewer than 32 faces)")
    return rows


def integrate(rows: list[dict[str, float]],
              certificate: dict[str, str] | None = None,
              *, force_tolerance: float = 1e-8,
              closure_tolerance: float = 1e-8) -> dict:
    if (not math.isfinite(force_tolerance) or force_tolerance <= 0 or
        not math.isfinite(closure_tolerance) or closure_tolerance <= 0):
        raise ValueError("invalid force or geometry tolerance")
    first = rows[0]
    fixed = ("dynamic_pressure", "reference_area", "reference_chord", "x_ref", "y_ref")
    for row in rows:
        for key in fixed:
            if not math.isclose(row[key], first[key], abs_tol=1e-12, rel_tol=1e-12):
                raise ValueError("inconsistent reference normalization or force origin")
        # The NACA solver is the 2-D unit-span extrusion; normal on airfoil
        # wall must lie in x-y, unlike a cap/farfield normal.
        if abs(row["nz"]) > 1e-8:
            raise ValueError("NACA wall contains a non-planar surface normal")
    area = math.fsum(row["area"] for row in rows)
    closure = {axis: math.fsum(row[f"n{axis}"] * row["area"] for row in rows)
               for axis in "xyz"}
    normalized_closure = math.sqrt(sum(v*v for v in closure.values())) / area
    if normalized_closure > closure_tolerance:
        raise ValueError(f"incomplete or incorrectly oriented closed surface: {normalized_closure:.5g}")
    scale = first["dynamic_pressure"] * first["reference_area"]
    origin_x, origin_y = first["x_ref"], first["y_ref"]
    def p_force(row, axis):
        return row["pressure"] * row[f"n{axis}"] * row["area"]
    def v_force(row, axis):
        return -row[f"tau_n_{axis}"] * row["area"]
    fp = {axis: math.fsum(p_force(r, axis) for r in rows) for axis in "xy"}
    fv = {axis: math.fsum(v_force(r, axis) for r in rows) for axis in "xy"}
    moment = math.fsum(
        (r["x"] - origin_x) * (p_force(r, "y") + v_force(r, "y")) -
        (r["y"] - origin_y) * (p_force(r, "x") + v_force(r, "x"))
        for r in rows
    )
    results = {
        "cd": (fp["x"] + fv["x"]) / scale,
        "cdp": fp["x"] / scale,
        "cdv": fv["x"] / scale,
        "cl": (fp["y"] + fv["y"]) / scale,
        "cm": moment / (scale * first["reference_chord"]),
    }
    if not math.isclose(results["cd"], results["cdp"] + results["cdv"],
                        rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("internal force decomposition inconsistency")
    if certificate is not None:
        for key, actual in results.items():
            try:
                claimed = float(certificate[key])
            except (ValueError, KeyError, TypeError) as error:
                raise ValueError(f"invalid certificate {key}") from error
            if not math.isfinite(claimed) or abs(actual - claimed) > force_tolerance * max(1.0, abs(actual)):
                raise ValueError(f"{key}: independent surface integral disagrees with C++ certificate")
    return {
        "status": "INDEPENDENT_SURFACE_FORCE_CONSISTENT" if certificate else
                  "INDEPENDENT_SURFACE_INTEGRAL_ONLY",
        "physical_reference_validation": "NOT_PERFORMED",
        "faces": len(rows),
        "total_area": area,
        "closed_surface_relative_defect": normalized_closure,
        "pressure_resultant": fp,
        "viscous_resultant": fv,
        "coefficients": results,
        "reference_area": first["reference_area"],
        "reference_chord": first["reference_chord"],
        "pressure_normal_BC": "CHECK_SOURCE_CONFIGURATION",
    }


def read_certificate(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("invalid force certificate header")
        rows = list(reader)
    if len(rows) != 1 or rows[0].get("converged", "").lower() not in ("yes", "true", "1"):
        raise ValueError("a single converged force certificate is required")
    return rows[0]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-faces", type=Path, required=True)
    parser.add_argument("--certificate", type=Path,
                        help="Compare against the converged C++ coefficients")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        report = integrate(read_rows(args.raw_faces),
                           read_certificate(args.certificate) if args.certificate else None)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        staging = args.output.with_name(args.output.name + ".tmp")
        staging.write_text(json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        os.replace(staging, args.output)
    except (OSError, ValueError, KeyError, OverflowError) as error:
        args.output.unlink(missing_ok=True)
        print(f"NACA_SURFACE_FORCE_ORACLE FAIL: {error}")
        return 2
    print(json.dumps(report, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
