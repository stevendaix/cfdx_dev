"""Independent three-grid wall-shear method-qualification gate.

A converged CFD solve is necessary but not sufficient. This intentionally
requires three converged levels and agreement between two independent wall
shear reconstructions; it never manufactures a literature reference value.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


LEVELS = (64, 128, 256)


def verify(payload: dict, method_agreement_limit: float = 0.05) -> dict:
    if not (math.isfinite(method_agreement_limit)
            and 0 < method_agreement_limit < 1):
        raise ValueError("invalid method-agreement limit")
    if payload.get("case") != "NACA0012_RE1000_A0":
        raise ValueError("wrong physical case")
    if not isinstance(payload.get("results"), list):
        raise ValueError("missing results list")
    rows = payload["results"]
    if len(rows) != 3 or tuple(sorted(int(row["N"]) for row in rows)) != LEVELS:
        raise ValueError("three distinct NACA mesh levels are required")
    rows = sorted(rows, key=lambda row: int(row["N"]))
    methods = ("cd_viscous_face", "cd_viscous_quadratic")
    series = {m: [] for m in methods}
    for row in rows:
        if row.get("status") != "PASS" or row.get("converged") is not True:
            raise ValueError(f"unconverged level N={row['N']}")
        continuity = float(row["continuity"])
        if not math.isfinite(continuity) or continuity >= 1e-7 or continuity < 0:
            raise ValueError(f"continuity gate failed for N={row['N']}")
        if int(row.get("covered_wall_faces", 0)) != int(row.get("wall_faces", -1)):
            raise ValueError(f"wall sample coverage incomplete for N={row['N']}")
        if int(row["wall_faces"]) <= 0:
            raise ValueError("empty wall sample")
        for m in methods:
            v = float(row[m])
            if not math.isfinite(v) or v <= 0:
                raise ValueError(f"invalid {m} on N={row['N']}")
            series[m].append(v)
    studies = {}
    for m, values in series.items():
        d1 = values[1] - values[0]
        d2 = values[2] - values[1]
        if d1*d2 <= 0 or abs(d2) >= abs(d1):
            raise ValueError(f"{m} is not in a monotone convergent sequence")
        # At fixed two-dimensional topology the nominal h refinement is 2.
        p = math.log(abs(d1/d2)) / math.log(2)
        if not math.isfinite(p) or p <= 0:
            raise ValueError(f"{m} has no positive observed order")
        extrapolated = values[2] + d2 / (2**p - 1)
        studies[m] = {
            "observed_order": p,
            "richardson_extrapolated": extrapolated,
            "fine_change": d2,
            "fine_value": values[2],
        }
    discrepancy = abs(series[methods[0]][-1] - series[methods[1]][-1]) / max(
        abs(series[methods[0]][-1]), abs(series[methods[1]][-1]))
    if discrepancy > method_agreement_limit:
        raise ValueError(f"independent wall-shear reconstructions disagree: {discrepancy}")
    return {
        "status": "GRID_AND_METHOD_CONSISTENT",
        "reference_validation": "NOT_PERFORMED",
        "method_discrepancy_on_n256": discrepancy,
        "studies": studies,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        result = verify(json.loads(args.input.read_text(encoding="utf-8")))
    except (ValueError, TypeError, KeyError, OverflowError) as error:
        print(f"NACA_WALL_SHEAR_INDEPENDENT: FAIL: {error}")
        return 1
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("NACA_WALL_SHEAR_INDEPENDENT: " + result["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
