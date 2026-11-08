#!/usr/bin/env python3
"""Fail-closed cellwise NACA momentum budget audit (never a CFD certificate)."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path

COLUMNS = (
    "iteration",
    "component",
    "rank",
    "cell",
    "x",
    "y",
    "z",
    "convection",
    "viscous",
    "transpose_viscous",
    "pressure",
    "body_force",
    "bounded_divergence",
    "mass_divergence",
    "collocated_difference",
    "physical_residual",
    "bounded_operator_residual",
    "worst_residual",
)
# Extended records certify per-cell arithmetic closure of the accepted-state
# signed face budget against the independently evaluated production CSR row.
# Legacy files are read for diagnosis, but cannot claim matrix verification.
EXTENDED_COLUMNS = (
    "temporal",
    "algebraic_residual",
    "matrix_face_difference",
    "matrix_face_scale",
)

TERMS = (
    "convection",
    "viscous",
    "transpose_viscous",
    "pressure",
    "body_force",
    "bounded_divergence",
    "mass_divergence",
    "collocated_difference",
    "physical_residual",
    "bounded_operator_residual",
    "worst_residual",
)


def _close(left: float, right: float) -> bool:
    # Machine-level arithmetic consistency, NOT a physics tolerance gate.
    return abs(left - right) <= 2e-10 * max(1.0, abs(left), abs(right))


def read_budget(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or len(set(reader.fieldnames)) != len(
            reader.fieldnames
        ):
            raise ValueError("missing or duplicated cell-budget CSV headers")
        if not set(COLUMNS).issubset(reader.fieldnames):
            raise ValueError("missing required cell-budget columns")
        extended_present = set(EXTENDED_COLUMNS) & set(reader.fieldnames)
        if extended_present and extended_present != set(EXTENDED_COLUMNS):
            raise ValueError("incomplete extended matrix/face audit columns")
        has_extended = bool(extended_present)
        rows: list[dict] = []
        for line, data in enumerate(reader, 2):
            if None in data or any(v is None or not v for v in data.values()):
                raise ValueError(f"malformed cell-budget CSV row {line}")
            row: dict = {}
            for key in ("iteration", "component", "rank", "cell"):
                if not data[key].isdecimal():
                    raise ValueError(f"invalid integer {key} at line {line}")
                row[key] = int(data[key])
            if row["iteration"] < 1 or row["component"] not in (0, 1, 2):
                raise ValueError(f"invalid iteration or component at line {line}")
            for key in (
                "x",
                "y",
                "z",
                *TERMS,
                *(EXTENDED_COLUMNS if has_extended else ()),
            ):
                val = float(data[key])
                if not math.isfinite(val):
                    raise ValueError(f"nonfinite {key} at line {line}")
                row[key] = val
            if row["worst_residual"] < 0:
                raise ValueError(f"negative worst residual at line {line}")
            conv = row["convection"] + row["viscous"] + row["transpose_viscous"]
            phys = conv + row["pressure"] + row["body_force"] + row.get("temporal", 0.0)
            if not _close(phys, row["physical_residual"]):
                raise ValueError(
                    f"cell transport/source sum inconsistent at line {line}"
                )
            if not _close(
                phys + row["bounded_divergence"], row["bounded_operator_residual"]
            ):
                raise ValueError(f"bounded operator sum inconsistent at line {line}")
            if has_extended:
                scale = row["matrix_face_scale"]
                error = row["matrix_face_difference"]
                if scale < 0:
                    raise ValueError(f"negative matrix/face scale at line {line}")
                if not _close(
                    row["algebraic_residual"] - row["bounded_operator_residual"], error
                ):
                    raise ValueError(
                        f"matrix/face residual difference inconsistent at line {line}"
                    )
                if abs(error) > 2e-10 * max(scale, 1e-15):
                    raise ValueError(f"matrix/face operator mismatch at line {line}")
            rows.append(row)
    if not rows:
        raise ValueError("empty cell-budget file")
    groups: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for row in rows:
        groups[row["iteration"], row["component"]].append(row)
    for key, group in groups.items():
        if {r["rank"] for r in group} != set(range(len(group))):
            raise ValueError(f"noncontiguous ranks for iteration/component {key}")
        if len({r["cell"] for r in group}) != len(group):
            raise ValueError(f"duplicate cells for iteration/component {key}")
        ordered = sorted(group, key=lambda x: x["rank"])
        if any(
            abs(a["physical_residual"]) + 1e-12 < abs(b["physical_residual"])
            for a, b in itertools.pairwise(ordered)
        ):
            raise ValueError(f"ranks not ordered by local residual for {key}")
        if not _close(
            abs(ordered[0]["physical_residual"]), ordered[0]["worst_residual"]
        ):
            raise ValueError(f"top residual disagrees with reported maximum for {key}")
        if any(
            not _close(r["worst_residual"], ordered[0]["worst_residual"])
            for r in ordered
        ):
            raise ValueError(f"inconsistent worst residual within {key}")
    return rows


def sampled_sign_alternation(
    rows: list[dict], component: int, cell: int, window: int = 8
) -> dict:
    """Detect a stable alternating SIGNED residual at a fixed cell.

    A unit sampling stride is required before suggesting a period-2 iterate.
    Sparse odd-stride samples can support alternating samples only.
    This never proves nonlinear stability or a solver physical qualification.
    """
    series = sorted(
        (r for r in rows if r["component"] == component and r["cell"] == cell),
        key=lambda r: r["iteration"],
    )[-window:]
    result = {
        "component": "xyz"[component],
        "cell": cell,
        "samples": len(series),
        "status": "INSUFFICIENT_EVIDENCE",
    }
    if len(series) < 6:
        return result
    steps = [b["iteration"] - a["iteration"] for a, b in itertools.pairwise(series)]
    if not (steps[0] > 0 and all(step == steps[0] for step in steps)):
        result["status"] = "IRREGULAR_SAMPLING"
        return result
    signed = [r["physical_residual"] for r in series]
    amplitude = sum(abs(v) for v in signed) / len(signed)
    if not (math.isfinite(amplitude) and amplitude > 1e-14):
        result["status"] = "TOO_SMALL_TO_IDENTIFY"
        return result
    alternating = all(a * b < 0 for a, b in itertools.pairwise(signed))
    stable_amplitude = max(abs(abs(v) - amplitude) for v in signed) <= 0.15 * amplitude
    result.update(
        {
            "stride": steps[0],
            "first_iteration": series[0]["iteration"],
            "last_iteration": series[-1]["iteration"],
            "mean_absolute_residual": amplitude,
            "signed_residuals": signed,
        }
    )
    if alternating and stable_amplitude:
        result["status"] = (
            "PERIOD_TWO_SIGNATURE"
            if steps[0] == 1
            else "ALTERNATING_SAMPLES_PERIOD_UNPROVEN"
        )
    else:
        result["status"] = "NO_STABLE_ALTERNATION"
    return result


def period_two_term_decomposition(
    rows: list[dict], component: int, cell: int, window: int = 8
) -> dict:
    """Attribute a measured two-step signed residual to face/source contributions.

    The split is arithmetic, not causal: changing one term may change every
    other term through nonlinear pressure/velocity coupling. In particular the
    collocated-flux delta contains more than the Rhie-Chow correction.
    """
    signature = sampled_sign_alternation(rows, component, cell, window)
    if signature["status"] != "PERIOD_TWO_SIGNATURE":
        return {"signature": signature, "term_decomposition_available": False}
    selected = sorted(
        (r for r in rows if r["component"] == component and r["cell"] == cell),
        key=lambda r: r["iteration"],
    )[-window:]
    positive = [r for r in selected if r["physical_residual"] > 0.0]
    negative = [r for r in selected if r["physical_residual"] < 0.0]
    if not positive or not negative:
        return {"signature": signature, "term_decomposition_available": False}

    def mean(entries: list[dict], key: str) -> float:
        return sum(r.get(key, 0.0) for r in entries) / len(entries)

    keys = (
        "convection",
        "viscous",
        "transpose_viscous",
        "pressure",
        "body_force",
        "temporal",
        "bounded_divergence",
        "collocated_difference",
        "physical_residual",
        "bounded_operator_residual",
    )
    contributions = {}
    for key in keys:
        pos, neg = mean(positive, key), mean(negative, key)
        contributions[key] = {
            "positive_state_mean": pos,
            "negative_state_mean": neg,
            "positive_minus_negative": pos - neg,
            "half_amplitude": (pos - neg) / 2.0,
        }
    main_terms = (
        "convection",
        "viscous",
        "transpose_viscous",
        "pressure",
        "body_force",
        "temporal",
    )
    full_delta = contributions["physical_residual"]["positive_minus_negative"]
    sum_delta = sum(contributions[key]["positive_minus_negative"] for key in main_terms)
    # Signed closure of BOTH halves of the cycle (not of two global maxima).
    closure_error = full_delta - sum_delta
    scale = max(
        abs(full_delta),
        sum(abs(contributions[key]["positive_minus_negative"]) for key in main_terms),
        1e-15,
    )
    if abs(closure_error) > 2e-10 * scale:
        raise ValueError("period-two contribution sum does not close")
    return {
        "signature": signature,
        "term_decomposition_available": True,
        "positive_samples": len(positive),
        "negative_samples": len(negative),
        "contributions": contributions,
        "physical_cycle_delta": full_delta,
        "physical_cycle_half_amplitude": full_delta / 2.0,
        "term_delta_closure_error": closure_error,
        "terms_by_absolute_cycle_amplitude": sorted(
            main_terms,
            key=lambda key: abs(contributions[key]["half_amplitude"]),
            reverse=True,
        ),
        "causality_warning": "Term difference is a signed arithmetic attribution, not a unique physical cause. Flux difference against interpolated velocity is not pure Rhie-Chow.",
    }


def diagnose(rows: list[dict]) -> dict:
    groups: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for row in rows:
        groups[row["iteration"], row["component"]].append(row)
    latest = max(key[0] for key in groups)
    summary = []
    for iteration, component in sorted(groups):
        top = min(groups[iteration, component], key=lambda row: row["rank"])
        terms = {
            key: top[key]
            for key in (
                "convection",
                "viscous",
                "transpose_viscous",
                "pressure",
                "body_force",
                "temporal",
                "bounded_divergence",
                "collocated_difference",
            )
            if key in top
        }
        dominant = sorted(terms, key=lambda key: abs(terms[key]), reverse=True)
        summary.append(
            {
                "iteration": iteration,
                "component": "xyz"[component],
                "worst_cell": top["cell"],
                "x": top["x"],
                "y": top["y"],
                "convection": top["convection"],
                "viscous": top["viscous"],
                "pressure": top["pressure"],
                "collocated_flux_difference": top["collocated_difference"],
                "physical_residual": top["physical_residual"],
                "bounded_operator_residual": top["bounded_operator_residual"],
                "temporal": top.get("temporal", 0.0),
                "algebraic_residual": top.get("algebraic_residual"),
                "matrix_face_difference": top.get("matrix_face_difference"),
                "dominant_terms_by_absolute_magnitude": dominant,
                "signed_term_budget": terms,
                "absolute_contribution_sum": sum(abs(v) for v in terms.values()),
            }
        )
    # Only sample at the same fixed cell; rank-0 rows can switch cells and
    # generate a false oscillation. Prefer the latest top-ranked y hotspot.
    latest_y = groups.get((latest, 1), [])
    hotspot_cell = next((r["cell"] for r in latest_y if r["rank"] == 0), None)
    sampled_oscillation = (
        sampled_sign_alternation(rows, 1, hotspot_cell)
        if hotspot_cell is not None
        else None
    )
    return {
        "status": "DIAGNOSTIC_ONLY_NOT_A_SOLVER_CERTIFICATE",
        "transverse_hotspot_sign_alternation": sampled_oscillation,
        "transverse_hotspot_cycle_terms": (
            period_two_term_decomposition(rows, 1, hotspot_cell)
            if hotspot_cell is not None
            else None
        ),
        "latest_sampled_iteration": latest,
        "sampled_cell_rows": len(rows),
        "matrix_face_closure_checked": all(
            all(name in row for name in EXTENDED_COLUMNS) for row in rows
        ),
        "max_matrix_face_error": (
            max(abs(row["matrix_face_difference"]) for row in rows)
            if all("matrix_face_difference" in row for row in rows)
            else None
        ),
        "budget_series": summary,
        "solver_convergence_proven": False,
        "three_grid_gci_proven": False,
        "sst_thermal_radiation_qualification_proven": False,
        "flux_note": "collocated difference = authoritative corrected flux minus directly interpolated velocity flux; not pure Rhie-Chow",
        "closure_note": "A(U)U-b(U) versus independently accumulated signed face fluxes at the same accepted state. Does not establish an independent pressure-gradient or limiter discretization oracle.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("budget", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = diagnose(read_budget(args.budget))
    serialized = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        if args.output.resolve() == args.budget.resolve():
            parser.error("cannot overwrite source CSV")
        args.output.write_text(serialized, encoding="utf-8")
    else:
        print(serialized, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
