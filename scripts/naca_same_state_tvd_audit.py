"""Fail-closed diagnostic of NACA accepted-state TVD / H/A / Rhie-Chow traces.

This validates identities and compares *signed* changes in selected cells.
It does NOT validate physical convergence or mistake the accepted-state
reconstruction for the pre-correction Picard assembly at its earlier iterate.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

FACE_COLUMNS = (
    "iteration",
    "face",
    "owner",
    "neighbour",
    "mass_flux",
    "upwind_convection",
    "tvd_convection",
    "deferred_tvd_flux",
    "predictor_hbya_flux",
    "rc_relaxed_flux",
    "rc_full_flux",
    "conservative_flux",
    "rc_pressure_response",
    "pressure_relaxation_split",
    "simplec_extra_flux",
    "decomposition_error",
    "internal_pressure_operator",
    "pressure_operator_error",
)
PICARD_COLUMNS = (
    "iteration",
    "face",
    "owner",
    "neighbour",
    "high_order_fraction",
    "assembly_mass_flux",
    "upwind_flux",
    "selected_flux",
    "deferred_tvd_flux",
)
LINEAR_COLUMNS = (
    "iteration",
    "component",
    "alpha_u",
    "linear_tolerance",
    "relaxed_true_l2",
    "relaxed_true_relative",
    "unrelaxed_original_l2",
    "unrelaxed_original_relative",
    "relaxation_displacement_l2",
    "velocity_change_linf",
    "relaxed_true_linf",
    "original_operator_linf",
    "worst_original_cell",
)


def close(a: float, b: float, tolerance: float = 2e-10) -> bool:
    return abs(a - b) <= tolerance * max(1.0, abs(a), abs(b))


def parse_int(text: str, label: str, allow_negative: bool = False) -> int:
    try:
        v = int(text)
    except (ValueError, TypeError) as err:
        raise ValueError(f"{label}: invalid integer {text!r}") from err
    if str(v) != text or (v < 0 and not allow_negative):
        raise ValueError(f"{label}: invalid index {text!r}")
    return v


def parse_float(text: str, label: str, allow_nan: bool = False) -> float:
    try:
        value = float(text)
    except (ValueError, TypeError) as err:
        raise ValueError(f"{label}: invalid float {text!r}") from err
    if not math.isfinite(value) and not (allow_nan and math.isnan(value)):
        raise ValueError(f"{label}: nonfinite {text!r}")
    return value


def csv_records(path: Path, required: tuple[str, ...]):
    if not path.is_file():
        raise ValueError(f"missing evidence: {path}")
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or any(
            name not in reader.fieldnames for name in required
        ):
            raise ValueError(f"{path}: missing required columns")
        for n, row in enumerate(reader, 2):
            if None in row or any(row.get(key) is None for key in required):
                raise ValueError(f"{path}:{n}: incomplete row")
            yield n, row


def load_faces(path: Path, independent: bool = True):
    faces = defaultdict(dict)
    max_err = {
        "deferred": 0.0,
        "predictor_pressure": 0.0,
        "same_state_flux": 0.0,
        "independent_internal_pressure": 0.0,
    }
    first_topology = {}
    for lineno, r in csv_records(path, FACE_COLUMNS):
        pre = f"{path}:{lineno}"
        it = parse_int(r["iteration"], f"{pre} iteration")
        face = parse_int(r["face"], f"{pre} face")
        owner = parse_int(r["owner"], f"{pre} owner")
        neighbour = parse_int(r["neighbour"], f"{pre} neighbour", True)
        if it < 1 or neighbour < -1 or face in faces[it]:
            raise ValueError(f"{pre}: duplicate/invalid iteration or face")
        if face in first_topology and first_topology[face] != (owner, neighbour):
            raise ValueError(f"{pre}: face topology changed")
        first_topology[face] = (owner, neighbour)
        row = {"iteration": it, "face": face, "owner": owner, "neighbour": neighbour}
        for key in FACE_COLUMNS[4:]:
            row[key] = parse_float(
                r[key],
                f"{pre} {key}",
                allow_nan=(
                    neighbour == -1
                    and key in ("internal_pressure_operator", "pressure_operator_error")
                ),
            )

        def match(a, b, name, pre=pre):
            err = abs(a - b)
            max_err[name] = max(max_err[name], err)
            if not close(a, b):
                raise ValueError(f"{pre}: {name} mismatch ({a} != {b})")

        match(
            row["tvd_convection"] - row["upwind_convection"],
            row["deferred_tvd_flux"],
            "deferred",
        )
        match(row["mass_flux"], row["conservative_flux"], "same_state_flux")
        match(
            row["predictor_hbya_flux"] - row["rc_relaxed_flux"],
            row["rc_pressure_response"],
            "predictor_pressure",
        )
        match(
            row["rc_full_flux"] - row["rc_relaxed_flux"],
            row["pressure_relaxation_split"],
            "predictor_pressure",
        )
        match(
            row["mass_flux"] - row["rc_full_flux"],
            row["simplec_extra_flux"],
            "predictor_pressure",
        )
        match(
            row["predictor_hbya_flux"]
            - row["rc_pressure_response"]
            + row["pressure_relaxation_split"]
            + row["simplec_extra_flux"],
            row["mass_flux"],
            "predictor_pressure",
        )
        match(row["decomposition_error"], 0.0, "predictor_pressure")
        if neighbour >= 0 and independent:
            match(
                row["predictor_hbya_flux"] - row["rc_full_flux"],
                row["internal_pressure_operator"],
                "independent_internal_pressure",
            )
            match(
                row["pressure_operator_error"], 0.0, "independent_internal_pressure"
            )
        faces[it][face] = row
    if not faces:
        raise ValueError(f"{path}: empty face evidence")
    if sorted(faces) != list(range(min(faces), max(faces) + 1)):
        raise ValueError(f"{path}: missing iterations")
    patterns = {frozenset(rows) for rows in faces.values()}
    if len(patterns) != 1:
        raise ValueError(
            f"{path}: face set changed; incomplete accepted-state evidence"
        )
    return faces, max_err


def load_picard(path: Path, accepted: dict[int, dict]):
    picard = defaultdict(dict)
    max_errors = {
        "picard_deferred": 0.0,
        "previous_accepted_flux_lag": 0.0,
        "previous_accepted_tvd_lag": 0.0,
    }
    for lineno, r in csv_records(path, PICARD_COLUMNS):
        pre = f"{path}:{lineno}"
        it = parse_int(r["iteration"], f"{pre} iteration")
        face = parse_int(r["face"], f"{pre} face")
        owner = parse_int(r["owner"], f"{pre} owner")
        neighbour = parse_int(r["neighbour"], f"{pre} neighbour", True)
        if it < 1 or neighbour < -1 or face in picard[it]:
            raise ValueError(f"{pre}: duplicate/invalid Picard face")
        row = {"owner": owner, "neighbour": neighbour}
        for key in PICARD_COLUMNS[4:]:
            row[key] = parse_float(r[key], f"{pre} {key}")
        if not (0.0 <= row["high_order_fraction"] <= 1.0):
            raise ValueError(f"{pre}: invalid continuation fraction")
        err = abs(row["selected_flux"] - row["upwind_flux"] - row["deferred_tvd_flux"])
        max_errors["picard_deferred"] = max(max_errors["picard_deferred"], err)
        if not close(err, 0.0):
            raise ValueError(f"{pre}: Picard deferred correction mismatch")
        if it not in accepted or face not in accepted[it]:
            raise ValueError(f"{pre}: no matching accepted-state face")
        a = accepted[it][face]
        if (owner, neighbour) != (a["owner"], a["neighbour"]):
            raise ValueError(f"{pre}: Picard face topology mismatch")
        picard[it][face] = row
    if set(picard) != set(accepted) or any(
        set(picard[it]) != set(accepted[it]) for it in accepted
    ):
        raise ValueError(f"{path}: incomplete Picard trace for accepted iteration")
    # Phi_old must be EXACTLY the previously accepted authoritative RC mass
    # flux, but initial iteration starts with a separately initialized field.
    # The deferred TVD flux also agrees for constant convection continuation.
    for it in sorted(picard):
        if it - 1 not in accepted:
            continue
        for face, r in picard[it].items():
            a = accepted[it - 1][face]
            err = abs(r["assembly_mass_flux"] - a["mass_flux"])
            max_errors["previous_accepted_flux_lag"] = max(
                max_errors["previous_accepted_flux_lag"], err
            )
            if not close(err, 0.0):
                raise ValueError(
                    f"iteration {it} face {face}: stale/wrong Picard mass flux"
                )
            # At identical full-order scheme the reconstruction must be the
            # previous accepted-state reconstruction of the same U_y and phi.
            if r["high_order_fraction"] == 1.0:
                err = abs(r["deferred_tvd_flux"] - a["deferred_tvd_flux"])
                max_errors["previous_accepted_tvd_lag"] = max(
                    max_errors["previous_accepted_tvd_lag"], err
                )
                if not close(err, 0.0):
                    raise ValueError(
                        f"iteration {it} face {face}: inconsistent deferred Picard reconstruction"
                    )
    return picard, max_errors


def load_linear(path: Path, iterations: set[int]):
    seen = set()
    final = {}
    for lineno, r in csv_records(path, LINEAR_COLUMNS):
        pre = f"{path}:{lineno}"
        it = parse_int(r["iteration"], f"{pre} iteration")
        c = parse_int(r["component"], f"{pre} component")
        if c > 2 or (it, c) in seen:
            raise ValueError(f"{pre}: duplicate or invalid component")
        seen.add((it, c))
        floats = {
            key: parse_float(r[key], f"{pre} {key}") for key in LINEAR_COLUMNS[2:-1]
        }
        parse_int(r["worst_original_cell"], f"{pre} worst_original_cell")
        alpha = floats["alpha_u"]
        tol = floats["linear_tolerance"]
        if not (0 < alpha <= 1 and 1e-14 <= tol <= 1e-6):
            raise ValueError(f"{pre}: invalid alpha/linear tolerance")
        if floats["relaxed_true_relative"] > max(1.001 * tol, 1e-14):
            raise ValueError(
                f"{pre}: relaxed TRUE residual fails tolerance certificate"
            )
        if it == max(iterations):
            final[c] = floats
    expected = {(it, c) for it in iterations for c in range(3)}
    if seen != expected:
        raise ValueError(
            f"{path}: missing/extra linear certificates for accepted iterations"
        )
    return final


def signed_sum(rows: dict[int, dict], cell: int, field: str) -> float:
    return math.fsum(
        (1 if r["owner"] == cell else -1) * r[field]
        for r in rows.values()
        if r["owner"] == cell or r["neighbour"] == cell
    )


def analyze(
    face_csv: Path,
    linear_csv: Path,
    cells: tuple[int, ...] = (184, 201),
    period_start: int = 93,
    period_end: int = 100,
    picard_csv: Path | None = None,
) -> dict:
    faces, errors = load_faces(face_csv)
    present = set(faces)
    if period_start >= period_end or any(
        i not in present for i in range(period_start, period_end + 1)
    ):
        raise ValueError("missing iterations for requested cycle window")
    linear = load_linear(linear_csv, present)
    picard = None
    if picard_csv is not None:
        picard, picard_errors = load_picard(picard_csv, faces)
        errors.update(picard_errors)
    reported = []
    for cell in cells:
        for i in (period_start, period_end):
            if not any(
                r["owner"] == cell or r["neighbour"] == cell for r in faces[i].values()
            ):
                raise ValueError(
                    f"cell {cell}: missing complete incidence, iteration {i}"
                )
        cycle = []
        for it in range(period_start, period_end + 1):
            cycle.append(
                {
                    "iteration": it,
                    "upwind": signed_sum(faces[it], cell, "upwind_convection"),
                    "tvd": signed_sum(faces[it], cell, "tvd_convection"),
                    "deferred": signed_sum(faces[it], cell, "deferred_tvd_flux"),
                }
            )
        # Paired consecutive states isolate the alteration in accepted-state
        # face flux. The signs follow the face owner-neighbour orientation.
        a, b = faces[period_end - 1], faces[period_end]
        changes = {
            key: signed_sum(b, cell, key) - signed_sum(a, cell, key)
            for key in (
                "upwind_convection",
                "tvd_convection",
                "deferred_tvd_flux",
                "mass_flux",
                "predictor_hbya_flux",
                "rc_pressure_response",
            )
        }
        closure = (
            changes["tvd_convection"]
            - changes["upwind_convection"]
            - changes["deferred_tvd_flux"]
        )
        if not close(closure, 0.0):
            raise ValueError("cellwise flux delta lacks closure")
        deltas = []
        for face, row in b.items():
            if row["owner"] == cell or row["neighbour"] == cell:
                direction = 1 if row["owner"] == cell else -1
                deltas.append(
                    {
                        "face": face,
                        "deferred_tvd_flux_change": direction
                        * (row["deferred_tvd_flux"] - a[face]["deferred_tvd_flux"]),
                        "upwind_flux_change": direction
                        * (row["upwind_convection"] - a[face]["upwind_convection"]),
                        "pressure_response_change": direction
                        * (
                            row["rc_pressure_response"]
                            - a[face]["rc_pressure_response"]
                        ),
                    }
                )
        deltas.sort(
            key=lambda item: abs(item["deferred_tvd_flux_change"]), reverse=True
        )
        # A signed share is meaningful only if the total delta is not tiny.
        share = (
            changes["deferred_tvd_flux"] / changes["tvd_convection"]
            if abs(changes["tvd_convection"]) > 1e-15
            else None
        )
        even = [x["deferred"] for x in cycle if x["iteration"] % 2 == 0]
        odd = [x["deferred"] for x in cycle if x["iteration"] % 2 == 1]
        mean_even = math.fsum(even) / len(even)
        mean_odd = math.fsum(odd) / len(odd)
        period_two_amplitude = abs(mean_even - mean_odd)
        odd_drift = max(odd) - min(odd)
        even_drift = max(even) - min(even)
        picard_window = []
        if picard is not None:
            for it in range(period_start, period_end + 1):
                picard_window.append(
                    {
                        "iteration": it,
                        "assembly_deferred": signed_sum(
                            picard[it], cell, "deferred_tvd_flux"
                        ),
                        "accepted_deferred": signed_sum(
                            faces[it], cell, "deferred_tvd_flux"
                        ),
                    }
                )
        reported.append(
            {
                "cell": cell,
                "picard_vs_accepted_window": picard_window,
                "last_two_iteration_delta": changes,
                "signed_deferred_share": share,
                "period_two_even_odd_difference": period_two_amplitude,
                "within_parity_drift": max(odd_drift, even_drift),
                "period_two_local_diagnostic": period_two_amplitude
                > 10 * max(odd_drift, even_drift, 1e-14),
                "most_influential_faces": deltas[:6],
            }
        )
    return {
        "trace_consistency": "PASS",
        "physical_qualification": "NOT_EVALUATED",
        "face_csv": str(face_csv),
        "linear_csv": str(linear_csv),
        "picard_csv": str(picard_csv) if picard_csv else None,
        "iterations": len(faces),
        "faces_per_iteration": len(next(iter(faces.values()))),
        "max_absolute_closure_errors": errors,
        "linear_last_iteration": {str(k): v for k, v in linear.items()},
        "cells": reported,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--faces", required=True, type=Path)
    parser.add_argument("--linear", required=True, type=Path)
    parser.add_argument("--picard", type=Path, default=None)
    parser.add_argument("--cell", type=int, action="append", default=[])
    parser.add_argument("--start", type=int, default=93)
    parser.add_argument("--end", type=int, default=100)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        report = analyze(
            args.faces,
            args.linear,
            tuple(args.cell or (184, 201)),
            args.start,
            args.end,
            args.picard,
        )
    except (ValueError, OSError) as err:
        parser.exit(1, f"SAME_STATE_AUDIT_FAIL: {err}\n")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print("SAME_STATE_AUDIT_PASS: operator identities only, NOT physical qualification")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
