#!/usr/bin/env python3
"""V15 SST non-convergence audit: matrix residual arithmetic, wall and step merit.

Observational only. Never infers TH-03 physical qualification from inner solver
acceptance, a y+ proxy, or a reduction in one nonlinear merit function.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.diagnose_th03_v14 import _number, _parsed_rows
from scripts.diagnose_th03_v14 import analyse as analyse_v14


def analyse(
    log_path: Path,
    pressure_prefix: Path | None = None,
    pressure_target: float = 1e-11,
    sst_tolerance: float = 1e-8,
) -> dict:
    if not math.isfinite(sst_tolerance) or sst_tolerance <= 0:
        raise ValueError("SST tolerance must be positive and finite")
    result = analyse_v14(log_path, pressure_prefix, pressure_target)
    text = log_path.read_text(encoding="utf-8", errors="replace")
    rows = _parsed_rows(text, "SST_TRANSPORT_DIAGNOSTIC")
    walls = _parsed_rows(text, "SST_WALL_DISCRETE_AUDIT")
    line = _parsed_rows(text, "SST_NONLINEAR_LINESEARCH")
    samples = []
    for row in rows:
        values = {
            key: _number(row, key)
            for key in (
                "iteration",
                "k_local_relative",
                "omega_local_relative",
                "k_local_compensated",
                "omega_local_compensated",
                "k_local_roundoff_gap",
                "omega_local_roundoff_gap",
                "compensated_gate",
            )
        }
        if any(values[k] is None for k in values):
            continue  # The V14 log may not carry V15's precision data.
        gap = values["k_local_roundoff_gap"]
        native = values["k_local_relative"] if not values["compensated_gate"] else None
        robust = values["k_local_compensated"]
        # For native-gate traces, the maximum per-row gap is an upper bound
        # on the difference between the two maximum residuals.
        roundoff_cannot_explain = (
            native is not None
            and native > sst_tolerance
            and gap < native - sst_tolerance
            and robust > sst_tolerance
        )
        samples.append(
            {**values, "roundoff_alone_cannot_explain_k_gate": roundoff_cannot_explain}
        )
    invalid_lines = []
    accepted = 0
    rejected = 0
    for item in line:
        alpha = _number(item, "alpha")
        old = _number(item, "old_merit")
        new = _number(item, "new_merit")
        flag = _number(item, "accepted")
        if alpha is None or old is None or new is None or flag not in (0.0, 1.0):
            invalid_lines.append("missing-or-nonfinite")
            continue
        if flag:
            accepted += 1
            if alpha <= 0 or alpha > 1 or new >= old:
                invalid_lines.append("accepted-without-strict-reduction")
        else:
            rejected += 1
            if alpha != 0 or new != old:
                invalid_lines.append("rejected-but-state-changed")
    result["v15"] = {
        "status": "DIAGNOSTIC_ONLY_NOT_QUALIFIED",
        "precision_samples": len(samples),
        "last_precision_sample": samples[-1] if samples else None,
        "wall_audit_samples": len(walls),
        "last_wall_audit": walls[-1] if walls else None,
        "line_search_samples": len(line),
        "line_search_accepted_samples": accepted,
        "line_search_rejected_samples": rejected,
        "line_search_contract_failures": invalid_lines,
        "sst_tolerance_unchanged": sst_tolerance,
        "physical_qualification": "NOT_QUALIFIED",
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--sst-tolerance", type=float, default=1e-8)
    args = parser.parse_args()
    result = analyse(args.log, sst_tolerance=args.sst_tolerance)
    payload = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
