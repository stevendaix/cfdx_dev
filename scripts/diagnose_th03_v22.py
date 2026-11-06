#!/usr/bin/env python3
"""Independent V22 TH-03 failure analysis. Never certifies an unconverged state.

This script measures spatial resolution, momentum hot spots, a sampled
limit-cycle signature, and two independent near-wall velocity reconstructions.
The reattachment estimates are NOT a substitute for mesh-converged wall Cf.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.qualify_th03_v21 import NU, H, read_ascii_vtu, wall_sign_change


def quadratic_wall_derivative(y1: float, u1: float, y2: float, u2: float) -> float:
    """Fit U_t(y)=a*y+b*y^2 through no-slip wall and two layer samples."""
    if (
        not np.isfinite([y1, u1, y2, u2]).all()
        or y1 <= 0
        or not (1.5 * y1 <= y2 <= 6 * y1)
    ):
        raise ValueError("ill-conditioned or non-finite quadratic wall stencil")
    return float(((u1 / y1) * y2 - (u2 / y2) * y1) / (y2 - y1))


def periodicity(
    history: list[dict], *, field: str = "momentum_equation_residual_rhs_relative"
) -> dict:
    values = np.asarray([r.get(field, float("nan")) for r in history], dtype=float)[
        -60:
    ]
    if len(values) < 18 or not np.isfinite(values).all():
        return {"status": "insufficient_data", "points": len(values)}
    denominator = max(float(np.mean(np.abs(values))), 1e-300)
    scores = {
        lag: float(np.sqrt(np.mean((values[lag:] - values[:-lag]) ** 2)) / denominator)
        for lag in range(1, min(12, len(values) // 3) + 1)
    }
    best = min(scores, key=scores.get)
    return {
        "status": "repeating_cycle_suspected"
        if scores[best] < 1e-3
        else "no_certified_cycle",
        "best_lag": best,
        "relative_rms_periodicity_error": scores[best],
        "scores_by_lag": {str(k): v for k, v in scores.items()},
        "points": len(values),
    }


def momentum_hotspots(
    history: list[dict], centres: np.ndarray, bounds: np.ndarray
) -> dict:
    counts: Counter[int] = Counter()
    for r in history:
        if "momentum_residual_cell" not in r:
            continue
        cell = int(r["momentum_residual_cell"])
        if cell < 0 or cell >= len(centres):
            raise ValueError("momentum residual cell index outside VTU mesh")
        counts[cell] += 1
    if not counts:
        return {"status": "missing_momentum_cell_history"}
    top = []
    for cell, hits in counts.most_common(10):
        dx = float((bounds[cell, 1] - bounds[cell, 0]) / H)
        dy = float((bounds[cell, 3] - bounds[cell, 2]) / H)
        if dx <= 0 or dy <= 0:
            raise ValueError("invalid cell x/y geometry")
        top.append(
            {
                "cell": cell,
                "count": hits,
                "x_over_h": float(centres[cell, 0] / H),
                "y_over_h": float(centres[cell, 1] / H),
                "dx_over_h": dx,
                "dy_over_h": dy,
                "bounding_box_aspect_ratio": max(dx / dy, dy / dx),
            }
        )
    up = sum(hits for cell, hits in counts.items() if centres[cell, 0] < 0)
    n = sum(counts.values())
    return {"total_observations": n, "upstream_fraction": up / n, "top": top}


def compare_wall_reconstructions(
    centres: np.ndarray, bounds: np.ndarray, ux: np.ndarray
) -> dict:
    """Use only columns in the horizontal lower block: x/H in (0,8)."""
    if ux.shape != (len(centres),):
        raise ValueError("Ux vector has wrong size")
    candidates = np.where(
        (centres[:, 0] > 0)
        & (centres[:, 0] < 8 * H)
        & np.isclose(bounds[:, 2], 0, atol=1e-10, rtol=0)
    )[0]
    rows = []
    for first in candidates:
        col = np.where(
            (np.abs(bounds[:, 0] - bounds[first, 0]) < 1e-10)
            & (np.abs(bounds[:, 1] - bounds[first, 1]) < 1e-10)
            & (bounds[:, 2] >= -1e-10)
        )[0]
        col = col[np.argsort(centres[col, 1])]
        if len(col) < 2 or col[0] != first:
            raise ValueError("missing second wall-normal cell; no silent fallback")
        second = int(col[1])
        y1, y2 = float(centres[first, 1]), float(centres[second, 1])
        if not (np.isfinite([y1, y2, ux[first], ux[second]]).all()):
            raise ValueError("nonfinite near-wall state")
        linear = float(ux[first] / y1)
        quadratic = quadratic_wall_derivative(
            y1, float(ux[first]), y2, float(ux[second])
        )
        rows.append(
            {
                "x_over_h": float(centres[first, 0] / H),
                "first_cell": int(first),
                "second_cell": second,
                "distance_ratio": y2 / y1,
                "linear_cf_proxy": 2 * NU * linear / 44.2**2,
                "quadratic_cf_proxy": 2 * NU * quadratic / 44.2**2,
                "ux_first": float(ux[first]),
            }
        )
    rows.sort(key=lambda r: r["x_over_h"])
    x = np.array([r["x_over_h"] for r in rows])
    a = np.array([r["linear_cf_proxy"] for r in rows])
    b = np.array([r["quadratic_cf_proxy"] for r in rows])
    if len(rows) < 3:
        raise ValueError("too few downstream wall stations")
    ra = wall_sign_change(x, a)
    rb = wall_sign_change(x, b)
    step = float(np.max(np.diff(x)))
    return {
        "method": "two near-wall velocity extrapolations; NOT exact wall traction",
        "n_columns_0_to_8h": len(rows),
        "max_dx_between_centres_over_h": step,
        "linear_reattachment_x_over_h_proxy": ra,
        "quadratic_reattachment_x_over_h_proxy": rb,
        "delta_last_root_over_h": rb[-1] - ra[-1] if ra and rb else None,
        "experimental_half_width_over_h": 0.10,
        "resolution_smaller_than_experimental_half_width": step <= 0.10,
        "rows": rows,
    }


def upstream_reference_probe(centres: np.ndarray, ux: np.ndarray) -> dict:
    upstream = centres[:, 0] < 0
    if not np.any(upstream):
        return {"status": "no_upstream_cells"}
    xs = np.unique(np.round(centres[upstream, 0] / H, 7))
    nearest_x = float(xs[np.argmin(np.abs(xs + 4))])
    ids = np.where(np.abs(centres[:, 0] / H - nearest_x) < 1e-5)[0]
    ids = ids[np.argsort(centres[ids, 1])]
    # The lower upstream wall is at y/H=1 and the upper one at y/H=9.
    # These are simple samples; experimental upstream Uref and delta_99
    # must be compared using the same x/H=-4 location and normalization.
    mid = int(ids[np.argmin(abs(centres[ids, 1] / H - 5))])
    return {
        "target_x_over_h": -4.0,
        "nearest_sample_x_over_h": nearest_x,
        "station_offset_over_h": nearest_x + 4,
        "center_channel_u_x_at_nearest_sample_m_s": float(ux[mid]),
        "n_cross_stream_cells": len(ids),
        "must_compare_full_upstream_profile": True,
    }


def analyse(vtu: Path, history_json: Path) -> dict:
    centres, bounds, _, fields = read_ascii_vtu(vtu)
    j = json.loads(history_json.read_text())
    history = j.get("history", [])
    if not isinstance(history, list):
        raise TypeError("history must be array")
    wall = compare_wall_reconstructions(centres, bounds, fields["u_x"])
    hotspots = momentum_hotspots(history, centres, bounds)
    period = periodicity(history)
    return {
        "source_vtu": str(vtu),
        "history_json": str(history_json),
        "reference": {
            "experiment_x_over_h": 6.26,
            "experimental_uncertainty_over_h": 0.10,
            "NASA_SSTm_quasisteady_x_over_h": 6.50,
            "experiment_uses_skin_friction_zero": True,
        },
        "convergence": {
            "solver_reported_converged": bool(j.get("converged", False)),
            "last_momentum_relative": history[-1].get(
                "momentum_equation_residual_rhs_relative"
            )
            if history
            else None,
            "last_continuity_normalized": history[-1].get("continuity_normalized")
            if history
            else None,
            "periodicity": period,
            "momentum_hotspots": hotspots,
        },
        "near_wall_reconstruction": wall,
        "upstream": upstream_reference_probe(centres, fields["u_x"]),
        "physical_qualification": "NOT_QUALIFIED",
        "limitations": [
            "state is not converged",
            "two layer wall shear is a diagnostic, not a verified face traction",
            "V21 uses first-order upwind momentum and surrogate mesh",
            "grid independence, y+, inlet profiles and Cp/Cf reference agreement not demonstrated",
        ],
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vtu", type=Path, required=True)
    p.add_argument("--history", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    report = analyse(args.vtu, args.history)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        "momentum hotspot upstream fraction:",
        report["convergence"]["momentum_hotspots"].get("upstream_fraction"),
    )
    print(
        "cycle:",
        report["convergence"]["periodicity"]["status"],
        report["convergence"]["periodicity"].get("best_lag"),
    )
    print(
        "wall Cf-proxy roots:",
        report["near_wall_reconstruction"]["linear_reattachment_x_over_h_proxy"],
        report["near_wall_reconstruction"]["quadratic_reattachment_x_over_h_proxy"],
    )
    print("physical qualification:", report["physical_qualification"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
