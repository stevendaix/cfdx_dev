"""Fail-closed NACA 3-grid verification. Never manufactures missing CFD data.

CSV columns: case,h,converged,res_mom,res_mass,cd,cdp,cdv,cl,cm
h must represent the same characteristic mesh spacing along a geometrically
similar family (not just the number of radial layers). Rows are sorted by h.
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

REQUIRED = (
    "case",
    "h",
    "converged",
    "res_mom",
    "res_mass",
    "cd",
    "cdp",
    "cdv",
    "cl",
    "cm",
)


def qualify(rows, residual_gate=1e-8):
    if not math.isfinite(residual_gate) or residual_gate <= 0.0:
        raise ValueError("physical residual gate must be finite and strictly positive")
    if len(rows) != 3:
        raise ValueError("GCI requires exactly 3 converged geometrically similar grids")
    for row in rows:
        if row["converged"].lower() not in ("yes", "true", "1"):
            raise ValueError(
                f"{row['case']}: physical nonlinear convergence NOT attained"
            )
        for key in REQUIRED[1:]:
            if key == "converged":
                continue
            if not math.isfinite(float(row[key])):
                raise ValueError(f"{row['case']}: nonfinite {key}")
        if float(row["res_mom"]) < 0 or float(row["res_mass"]) < 0:
            raise ValueError(f"{row['case']}: negative residual norm")
        if max(float(row["res_mom"]), float(row["res_mass"])) > residual_gate:
            raise ValueError(
                f"{row['case']}: physical residual exceeds {residual_gate}"
                continue
            if not math.isfinite(float(row[key])):
                raise ValueError(f"{row['case']}: nonfinite {key}")
        if max(float(row["res_mom"]), float(row["res_mass"])) > residual_gate:
            raise ValueError(
                f"{row['case']}: physical residual exceeds {residual_gate}"
            )
        if abs(float(row["cd"]) - float(row["cdp"]) - float(row["cdv"])) > 1e-8:
            raise ValueError(f"{row['case']}: Cd != Cd_pressure+Cd_viscous")
    ordered = sorted(rows, key=lambda row: float(row["h"]), reverse=True)
    h = [float(row["h"]) for row in ordered]
    if not 0 < h[2] < h[1] < h[0]:
        raise ValueError(
            "mesh characteristic sizes must be positive and strictly decreasing"
        )
    r21, r32 = h[0] / h[1], h[1] / h[2]
    if abs(r21 - r32) / max(r21, r32) > 0.05:
        raise ValueError(
            "unequal refinement ratios: this estimator requires near-uniform r"
        )
    # Do not replace two DIFFERENT mesh ratios by a geometric mean: the
    # Richardson order is then biased, even when their mismatch is below 5%.
    # With f(h) = f_exact + C*h**p, the difference ratio is
    # (h0**p-h1**p)/(h1**p-h2**p), not a power of mean(r21,r32).
    result = {
        "refinement_ratios": [r21, r32],
        "meshes": [row["case"] for row in ordered],
        "quantities": {},
        "estimator": "three_grid_unequal_ratio_power_law_candidate",
    }
    log21, log32 = math.log(r21), math.log(r32)

    def difference_ratio(order):
        # expm1 avoids cancellation at small orders.
        return math.expm1(order * log21) / (-math.expm1(-order * log32))

    def observed_order(magnitude_ratio):
        if magnitude_ratio <= log21 / log32:
            return None  # would require non-positive observed order
        lo, hi = 1.e-9, 8.0
        if magnitude_ratio >= difference_ratio(hi):
            return None
        for _ in range(80):
            mid = (lo + hi) / 2.0
            if difference_ratio(mid) < magnitude_ratio:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2.0

    for key in ("cd", "cdp", "cdv", "cl", "cm"):
        f0, f1, f2 = (float(row[key]) for row in ordered)
        d01, d12 = f0 - f1, f1 - f2
        if d01 == 0.0 or d12 == 0.0 or d01 * d12 <= 0:
            result["quantities"][key] = {
                "status": "nonmonotone_or_plateau", "GCI": None,
            }
            continue
        order = observed_order(abs(d01 / d12))
        if order is None or not math.isfinite(order):
            result["quantities"][key] = {
                "status": "invalid_observed_order", "GCI": None,
            }
            continue
        denominator = math.expm1(order * log32)
        extrapolated = f2 + (f2 - f1) / denominator
        absolute_gci = 1.25 * abs(f2 - f1) / denominator
        if not (math.isfinite(extrapolated) and math.isfinite(absolute_gci)):
            result["quantities"][key] = {
                "status": "invalid_refinement_extrapolation", "GCI": None,
            }
            continue
        result["quantities"][key] = {
            "status": "asymptotic_candidate",
            "order": order,
            "extrapolated": extrapolated,
            "GCI_absolute": absolute_gci,
            "GCI_relative": (absolute_gci / abs(f2) if abs(f2) > 1e-14 else None),
            "fine_pair_refinement_ratio": r32,
        }
    return result

            continue
        p = math.log(abs(e21 / e32)) / math.log(r)
        if not (0 < p < 8):
            result["quantities"][key] = {
                "status": "invalid_observed_order",
                "order": p,
                "GCI": None,
            }
            continue
        denom = r**p - 1.0
        extrapolated = f3 + (f3 - f2) / denom
        absolute_gci = 1.25 * abs(f3 - f2) / denom
        result["quantities"][key] = {
            "status": "asymptotic_candidate",
            "order": p,
            "extrapolated": extrapolated,
            "GCI_absolute": absolute_gci,
            "GCI_relative": (absolute_gci / abs(f3) if abs(f3) > 1e-14 else None),
        }
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv_path", type=Path)
    ap.add_argument("--residual-gate", type=float, default=1e-8)
    args = ap.parse_args()
    with args.csv_path.open(newline="") as f:
        reader = csv.DictReader(f)
        missing = set(REQUIRED) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"missing fields: {sorted(missing)}")
        rows = list(reader)
    report = qualify(rows, args.residual_gate)
    print(json.dumps(report, indent=2, sort_keys=True))
    if report["quantities"]["cd"]["status"] != "asymptotic_candidate":
        raise ValueError(
            "drag does not demonstrate monotone convergence; Cd GCI refused"
        )


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError) as exc:
        print(f"QUALIFICATION FAILED: {exc}", file=sys.stderr)
        raise SystemExit(2)
