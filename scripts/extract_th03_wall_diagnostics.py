"""Extract diagnostic BFS Cf and reattachment from CFDX *native* face CSV.

This uses CFDX's explicitly first-order one-sided wall shear. It is a
mesh-sensitivity diagnostic and cannot by itself validate Driver--Seegmiller.
The reference-comparison oracle still requires full three-mesh evidence.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path


class WallDiagnosticError(ValueError):
    pass


def calculate_wall_diagnostics(
    face_csv: Path,
    *,
    patch: str,
    density: float,
    bulk_speed: float,
    step_height: float,
    step_x: float = 0.0,
) -> tuple[list[dict[str, float]], float]:
    if (
        any(not math.isfinite(x) for x in (density, bulk_speed, step_height, step_x))
        or min(density, abs(bulk_speed), step_height) <= 0
    ):
        raise WallDiagnosticError(
            "invalid density, inlet speed, step height or step origin"
        )
    grouped: dict[float, list[tuple[float, float]]] = defaultdict(list)
    with face_csv.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        required = {"patch", "x_m", "area_m2", "wall_shear_x_Pa", "wall_shear_method"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise WallDiagnosticError("CFDX face CSV columns missing")
        for row in reader:
            if row["patch"] != patch:
                continue
            if row["wall_shear_method"] != "first_order_one_sided_diagnostic":
                raise WallDiagnosticError(
                    "wall-shear reconstruction missing or not supported"
                )
            try:
                x = float(row["x_m"])
                area = float(row["area_m2"])
                shear = float(row["wall_shear_x_Pa"])
            except (TypeError, ValueError) as error:
                raise WallDiagnosticError("non-numeric wall face output") from error
            if any(not math.isfinite(a) for a in (x, area, shear)) or area <= 0:
                raise WallDiagnosticError("non-finite face position, area or shear")
            x_over_h = (x - step_x) / step_height
            if x_over_h >= 0:
                # Group only faces at the same streamwise station. For an
                # unstructured wall this remains a raw-face diagnostic; a
                # production comparison needs explicit binning/error bounds.
                grouped[round(x_over_h, 12)].append((shear, area))
    if len(grouped) < 3:
        raise WallDiagnosticError("at least 3 downstream wall stations are mandatory")
    profiles = []
    for x, samples in sorted(grouped.items()):
        mean_shear = sum(s * a for s, a in samples) / sum(a for _, a in samples)
        cf = 2.0 * mean_shear / (density * bulk_speed * bulk_speed)
        profiles.append({"x": x, "value": cf})
    crossings = []
    for previous, current in itertools.pairwise(profiles):
        # The separated zone has Cf<0; the first negative->positive
        # crossing is the definition. Zero-valued isolated wall samples
        # are not assumed to be positive (avoid accidental false root).
        if previous["value"] < 0 < current["value"]:
            xa, xb = previous["x"], current["x"]
            ca, cb = previous["value"], current["value"]
            crossings.append(xa - ca * (xb - xa) / (cb - ca))
    if len(crossings) != 1:
        raise WallDiagnosticError(
            f"expected one negative-to-positive Cf crossing, found {len(crossings)}"
        )
    return profiles, crossings[0]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--face-csv", required=True, type=Path)
    parser.add_argument("--patch", required=True)
    parser.add_argument("--density", required=True, type=float)
    parser.add_argument("--bulk-speed", required=True, type=float)
    parser.add_argument("--step-height", required=True, type=float)
    parser.add_argument("--step-x", type=float, default=0.0)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        profiles, root = calculate_wall_diagnostics(
            args.face_csv,
            patch=args.patch,
            density=args.density,
            bulk_speed=args.bulk_speed,
            step_height=args.step_height,
            step_x=args.step_x,
        )
    except (OSError, WallDiagnosticError) as error:
        parser.exit(2, f"TH-03 wall diagnostic FAIL: {error}\n")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "cf_diagnostic.csv").open(
        "w", newline="", encoding="utf-8"
    ) as target:
        writer = csv.DictWriter(target, fieldnames=("x", "value"))
        writer.writeheader()
        writer.writerows(profiles)
    (args.output_dir / "reattachment_diagnostic.json").write_text(
        json.dumps(
            {
                "case_id": "TH-03",
                "reattachment_x_over_h": root,
                "physical_qualification": "NOT_QUALIFIED",
                "wall_shear_method": "first_order_one_sided_diagnostic",
                "wall_face_stations": len(profiles),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"TH-03 first-order diagnostic x_r/h = {root:.6g} (NOT VALIDATED)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
