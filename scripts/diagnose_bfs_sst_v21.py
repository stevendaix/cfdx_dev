"""V21 BFS owner-cell synchronization audit, without physical certification.

Read-only parser for V19/V20/V21 logs. A vanishing owner-row gap proves only
that the 69 owner constraints in this captured SST iterate are satisfied;
it cannot certify interior k/omega transport nor global U/p convergence.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

from diagnose_bfs_sst_v20 import analyze as analyze_v20

_OWNER = ("iteration", "owners", "initial_gap", "trial_gap", "final_gap",
          "final_relative", "worst_cell", "line_accepted")
_BOUNDARY = ("iteration", "owners", "max_absolute_mismatch", "worst_cell")
_INNER = ("iteration", "k_returned_nonlinear_raw", "omega_returned_nonlinear_raw",
          "k_returned_row_scaled", "omega_returned_row_scaled", "k_local_relative",
          "omega_local_relative", "k_local_worst_cell", "omega_local_worst_cell")


def _kv(line: str) -> dict[str, str]:
    return dict(re.findall(r"([A-Za-z_]+)=([^\s]+)", line))


def _get_record(line: str, fields: tuple[str, ...], integer_fields: tuple[str, ...]) -> dict:
    raw = _kv(line)
    if not set(fields).issubset(raw):
        raise ValueError("incomplete V21 SST diagnostic: " + line.split(" ", 1)[0])
    out = {}
    for field in fields:
        try:
            out[field] = int(raw[field]) if field in integer_fields else float(raw[field])
        except ValueError as exc:
            raise ValueError("invalid numeric diagnostic: " + field) from exc
        if not math.isfinite(out[field]) or out[field] < 0:
            raise ValueError("negative/non-finite diagnostic: " + field)
    return out


def analyze(text: str) -> dict:
    parent = analyze_v20(text)
    owner = []
    boundary = []
    transport = []
    for line in text.splitlines():
        if line.startswith("SST_OWNER_SYNC_V21 "):
            row = _get_record(line, _OWNER,
                              ("iteration", "owners", "worst_cell", "line_accepted"))
            if row["iteration"] < 1 or row["line_accepted"] not in (0, 1):
                raise ValueError("invalid owner iteration or acceptance flag")
            owner.append(row)
        elif line.startswith("SST_WALL_OWNER_CONSTRAINT_V18 "):
            boundary.append(_get_record(line, _BOUNDARY,
                                        ("iteration", "owners", "worst_cell")))
        elif line.startswith("SST_TRANSPORT_DIAGNOSTIC "):
            transport.append(_get_record(line, _INNER,
                                         ("iteration", "k_local_worst_cell",
                                          "omega_local_worst_cell")))
    stable = bool(owner) and all(
        item["owners"] > 0 and item["line_accepted"] == 1
        and item["final_relative"] <= 64 * sys.float_info.epsilon
        for item in owner)
    # A successful wall projection MUST NOT turn a failed SST / BFS into a
    # passing physical qualification. Even an empty log remains unqualified.
    return {"schema": "CFDX-BFS-SST-V21", "owner_projection_stable": stable,
            "owner_sync": owner, "owner_residual": boundary,
            "transport": transport, "row_budgets": parent["row_budgets"],
            "committed_states": parent["committed_states"],
            "three_state_contraction_proved": False,
            "driver_seegmiller_qualified": False,
            "reason": "Wall constraint consistency is not global SST convergence or NASA qualification"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = analyze(args.log.read_text(encoding="utf-8", errors="replace"))
    except (OSError, ValueError) as exc:
        result = {"schema": "CFDX-BFS-SST-V21", "error": str(exc),
                  "owner_projection_stable": False,
                  "three_state_contraction_proved": False,
                  "driver_seegmiller_qualified": False}
        status = 3
    else:
        status = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n",
                           encoding="utf-8")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
