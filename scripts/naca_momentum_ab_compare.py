"""Fail-closed matched-run NACA TVD/UPWIND momentum diagnostic (not qualification)."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from naca_momentum_budget_diagnosis import (
    diagnose,
    period_two_term_decomposition,
    read_budget,
)

METRICS = (
    "momentum_rhs_relative",
    "momentum_y_rhs_relative",
    "face_momentum_y_linf",
    "continuity_normalized",
    "pressure_residual",
    "velocity_change_inf",
)


def read_history(path: Path, expected_scheme: str) -> list[dict]:
    rows: list[dict] = []
    needed = {
        "iteration",
        "convection_scheme",
        "alpha_u",
        "alpha_p",
        "high_order_fraction",
        *METRICS,
    }
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or not needed.issubset(reader.fieldnames):
            raise ValueError(f"missing required NACA history columns: {path}")
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"duplicated NACA history columns: {path}")
        for line, item in enumerate(reader, 2):
            if None in item or any(x is None or x == "" for x in item.values()):
                raise ValueError(f"malformed history row {line}")
            if item["convection_scheme"] != expected_scheme:
                raise ValueError(f"history row {line} is not {expected_scheme}")
            if not item["iteration"].isdecimal():
                raise ValueError(f"invalid history iteration row {line}")
            record = {"iteration": int(item["iteration"])}
            if record["iteration"] != len(rows) + 1:
                raise ValueError(
                    f"missing, repeated, or reordered history iteration {line}"
                )
            for key in (*METRICS, "alpha_u", "alpha_p", "high_order_fraction"):
                record[key] = float(item[key])
                if not math.isfinite(record[key]):
                    raise ValueError(f"nonfinite {key} at row {line}")
            for key in METRICS:
                if record[key] < 0.0:
                    raise ValueError(f"negative {key} at row {line}")
            if not (
                0.0 < record["alpha_u"] <= 1.0
                and 0.0 < record["alpha_p"] <= 1.0
                and 0.0 <= record["high_order_fraction"] <= 1.0
            ):
                raise ValueError(f"invalid relaxation/order fraction at row {line}")
            rows.append(record)
    if len(rows) < 8:
        raise ValueError("need at least eight accepted iterations in each campaign")
    return rows


def compare(
    tvd_history: list[dict],
    upwind_history: list[dict],
    tvd_budget: list[dict],
    upwind_budget: list[dict],
) -> dict:
    if len(tvd_history) != len(upwind_history):
        raise ValueError(
            "TVD/UPWIND runs must have identical accepted iteration counts"
        )
    for tvd, up in zip(tvd_history, upwind_history):
        if any(
            not math.isclose(tvd[key], up[key], rel_tol=1e-13, abs_tol=1e-14)
            for key in ("alpha_u", "alpha_p")
        ):
            raise ValueError(
                "TVD/UPWIND runs used different pressure/velocity relaxation"
            )
        if tvd["high_order_fraction"] < 1.0 - 1e-12:
            raise ValueError("TVD requires full high-order fraction for this A/B")
    if tvd_budget[-1]["iteration"] != len(tvd_history) or upwind_budget[-1][
        "iteration"
    ] != len(upwind_history):
        raise ValueError("budget evidence does not reach the final accepted iteration")
    tvd_report = diagnose(tvd_budget)
    up_report = diagnose(upwind_budget)
    last_tvd, last_up = tvd_history[-1], upwind_history[-1]
    latest_y = [
        r
        for r in tvd_budget
        if r["iteration"] == len(tvd_history) and r["component"] == 1 and r["rank"] == 0
    ]
    if len(latest_y) != 1:
        raise ValueError(
            "missing or duplicated TVD worst transverse cell at final iteration"
        )
    cell = latest_y[0]["cell"]
    return {
        "status": "DIAGNOSTIC_ONLY_NOT_PHYSICAL_QUALIFICATION",
        "accepted_iterations": len(tvd_history),
        "same_relaxation": True,
        "final_tvd_high_order_fraction": last_tvd["high_order_fraction"],
        "alpha_u": last_tvd["alpha_u"],
        "alpha_p": last_tvd["alpha_p"],
        "final": {
            metric: {
                "TVD": last_tvd[metric],
                "UPWIND": last_up[metric],
                "TVD_over_UPWIND": (
                    last_tvd[metric] / last_up[metric]
                    if last_up[metric] != 0.0
                    else None
                ),
            }
            for metric in METRICS
        },
        "tvd_transverse_hotspot": {
            "cell": cell,
            "x": latest_y[0]["x"],
            "y": latest_y[0]["y"],
            "signed_term_cycle": period_two_term_decomposition(tvd_budget, 1, cell),
        },
        "matrix_face_max_error": {
            "TVD": tvd_report["max_matrix_face_error"],
            "UPWIND": up_report["max_matrix_face_error"],
        },
        "tvd_convergence_proven": False,
        "upwind_convergence_proven": False,
        "explanation": (
            "Matched history compares two nonlinear iterations, not solutions. "
            "The pressure-linear residual does not establish nonlinear convergence. "
            "The cycle-term split is arithmetic attribution, not proof of the "
            "limiter or Rhie-Chow as the sole cause."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mesh",
        type=Path,
        required=True,
        help="Common mesh for both runs; its SHA-256 is recorded",
    )
    parser.add_argument("--tvd-history", type=Path, required=True)
    parser.add_argument("--upwind-history", type=Path, required=True)
    parser.add_argument("--tvd-budget", type=Path, required=True)
    parser.add_argument("--upwind-budget", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = compare(
        read_history(args.tvd_history, "TVD"),
        read_history(args.upwind_history, "UPWIND"),
        read_budget(args.tvd_budget),
        read_budget(args.upwind_budget),
    )
    if not args.mesh.is_file():
        parser.error("common mesh file does not exist")
    sha = hashlib.sha256()
    with args.mesh.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(block)
    result["common_mesh_path"] = str(args.mesh.resolve())
    result["common_mesh_sha256"] = sha.hexdigest()
    result["provenance_limit"] = (
        "Common mesh hash is recorded from this CLI argument; histories do not "
        "embed a mesh hash. Caller must execute both runs on this exact file. "
        "Matched histories alone cannot cryptographically prove mesh identity."
    )
    data = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        output = args.output.resolve()
        sources = (
            args.tvd_history,
            args.upwind_history,
            args.tvd_budget,
            args.upwind_budget,
        )
        if output in {source.resolve() for source in sources}:
            parser.error("output would overwrite a source file")
        args.output.write_text(data, encoding="utf-8")
    else:
        print(data, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
