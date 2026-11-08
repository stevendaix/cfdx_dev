"""Conservative bin-wise NACA drag decomposition across DIFFERENT wall samplings.

Each face belongs to one x/c bin based on its centre (diagnostic only).
The face-integrated force channels must exactly reproduce the published
converged force certificate; do not compare raw face indices across meshes.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

BINS = (0.0, 0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 1.000000001)


def read_profile(path: Path, certificate: Path) -> dict:
    with path.open(encoding="utf-8", newline="") as file:
        rows = list(csv.DictReader(file))
    with certificate.open(encoding="utf-8", newline="") as file:
        frows = list(csv.DictReader(file))
    if (len(frows) != 1 or not rows
            or frows[0].get("converged") not in ("yes", "true", "1")):
        raise ValueError("missing full face-profile data or unconverged certificate")
    if len({row["face"] for row in rows}) != len(rows):
        raise ValueError("duplicate face index")
    channels = ("x_over_c", "y_over_c", "Cp", "Cf_tangent",
                "Cd_pressure_face", "Cd_viscous_face")
    sums = {"cdp": 0.0, "cdv": 0.0}
    region = {side: [] for side in ("upper", "lower", "te_closure")}
    for row in rows:
        fields = {key: float(row[key]) for key in channels}
        if not all(math.isfinite(value) for value in fields.values()):
            raise ValueError("nonfinite wall profile")
        x = fields["x_over_c"]
        if x < -1e-6 or x > 1.000001:
            raise ValueError("wall position outside NACA chord")
        # TE thickness closure lies at exactly x=1 and is not split with N.
        if x >= 1.0 - 1e-12:
            side = "te_closure"
        else:
            side = "upper" if fields["y_over_c"] >= 0 else "lower"
        region[side].append(fields)
        sums["cdp"] += fields["Cd_pressure_face"]
        sums["cdv"] += fields["Cd_viscous_face"]
    expected = {key: float(frows[0][key]) for key in ("cd", "cdp", "cdv")}
    if not all(math.isfinite(value) for value in expected.values()):
        raise ValueError("nonfinite force certificate")
    if (not math.isclose(sums["cdp"], expected["cdp"], rel_tol=2e-7, abs_tol=2e-9)
            or not math.isclose(sums["cdv"], expected["cdv"], rel_tol=2e-7, abs_tol=2e-9)
            or not math.isclose(expected["cd"], sums["cdp"] + sums["cdv"],
                                rel_tol=2e-7, abs_tol=2e-9)):
        raise ValueError("face forces do not reconstruct drag certificate")
    result = {"force": expected, "bins": [], "n_wall_faces": len(rows)}
    for side, entries in region.items():
        for start, end in zip(BINS[:-1], BINS[1:]):
            values = [r for r in entries if start <= r["x_over_c"] < end]
            result["bins"].append({
                "side": side, "range": [start, end], "face_count": len(values),
                "cdp": sum(v["Cd_pressure_face"] for v in values),
                "cdv": sum(v["Cd_viscous_face"] for v in values),
            })
    return result


def compare(a: dict, b: dict) -> dict:
    if len(a["bins"]) != len(b["bins"]):
        raise ValueError("bin layout mismatch")
    delta = {comp: float(b["force"][comp]) - float(a["force"][comp])
             for comp in ("cd", "cdp", "cdv")}
    rows = []
    for x, y in zip(a["bins"], b["bins"]):
        if x["side"] != y["side"] or x["range"] != y["range"]:
            raise ValueError("mismatched bin definitions")
        rows.append({
            "side": x["side"], "range": x["range"],
            "delta_cdp": y["cdp"] - x["cdp"],
            "delta_cdv": y["cdv"] - x["cdv"],
            "n_faces_a": x["face_count"], "n_faces_b": y["face_count"],
        })
    for comp in ("cdp", "cdv"):
        total = sum(r["delta_" + comp] for r in rows)
        if not math.isclose(total, delta[comp], rel_tol=1e-8, abs_tol=1e-9):
            raise ValueError("integrated bin difference does not close")
    if not math.isclose(delta["cd"], delta["cdp"] + delta["cdv"], abs_tol=1e-9):
        raise ValueError("force decomposition does not close")
    return {
        "status": "TWO_MESH_PROFILE_DIAGNOSTIC_ONLY",
        "delta_force": delta, "bins": rows,
        "pressure_viscous_cancellation_fraction": (
            1.0 - abs(delta["cd"])
            / max(abs(delta["cdp"]) + abs(delta["cdv"]), 1e-30)
        ),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("a", "b"):
        p.add_argument(f"--{key}-profiles", type=Path, required=True)
        p.add_argument(f"--{key}-force", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    try:
        r = compare(read_profile(args.a_profiles, args.a_force),
                    read_profile(args.b_profiles, args.b_force))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(r, indent=2) + "\n")
    except (ValueError, TypeError, KeyError, OSError) as err:
        args.output.unlink(missing_ok=True)
        print("CROSS_SURFACE_PROFILE: BLOCKED", err)
        return 2
    print(json.dumps({"status": r["status"], "delta_force": r["delta_force"],
                      "cancellation": r["pressure_viscous_cancellation_fraction"]},
                     indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
