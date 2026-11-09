#!/usr/bin/env python3
"""Fail-closed evidence of multiple *committed* physical state contractions.

Not a BFS experimental qualification: true reattachment/Cf/profiles and ALL
physical convergence gates must be assessed separately. Ignore candidates.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


class QualificationError(ValueError):
    pass


def certify(log: str, *, minimum_commits: int = 3, max_ratio: float = 0.995) -> dict:
    if minimum_commits < 3 or not 0 < max_ratio < 1 or not math.isfinite(max_ratio):
        raise QualificationError("invalid proof controls")
    records = []
    for line in log.splitlines():
        if not line.startswith("GLOBAL_CONTRACTION_COMMITTED_V14 "):
            continue
        fields = dict(re.findall(r"(\w+)=([^\s]+)", line))
        if not {"iteration", "committed", "contractions", "merit", "converged", "qualified"} <= fields.keys():
            raise QualificationError("incomplete committed contraction record")
        try:
            r = {
                "iteration": int(fields["iteration"]),
                "committed": int(fields["committed"]),
                "contractions": int(fields["contractions"]),
                "merit": float(fields["merit"]),
                "converged": int(fields["converged"]),
                "qualified": int(fields["qualified"]),
            }
        except (ValueError, OverflowError) as exc:
            raise QualificationError("malformed contraction record") from exc
        if r["iteration"] <= 0 or not math.isfinite(r["merit"]) or r["merit"] < 0:
            raise QualificationError("invalid physical merit")
        if r["converged"] not in (0, 1) or r["qualified"] != 0:
            raise QualificationError("unproven physical qualification asserted")
        if records:
            old = records[-1]
            if r["iteration"] <= old["iteration"] or r["committed"] != old["committed"] + 1:
                raise QualificationError("duplicate/noncommitted/skip iteration")
            # A reset of the streak is allowed only on an explicitly observed
            # operator change; this proof requires one identical operator.
            if r["contractions"] != old["contractions"] + 1:
                raise QualificationError("operator changed or contraction ledger reset")
            if not (r["merit"] <= max_ratio * old["merit"] * (1 + 1e-12)):
                raise QualificationError("merit did not contract")
        else:
            if r["committed"] != 1 or r["contractions"] != 0:
                raise QualificationError("missing first accepted baseline")
        records.append(r)
    if len(records) < minimum_commits:
        raise QualificationError(f"only {len(records)} committed states; need {minimum_commits}")
    return {
        "schema": "CFDX-BFS-GLOBAL-CONTRACTION-V14",
        "committed": len(records),
        "contractions": len(records) - 1,
        "initial_merit": records[0]["merit"],
        "final_merit": records[-1]["merit"],
        "strict_contraction_proved": True,
        "physical_converged": bool(records[-1]["converged"]),
        "bfs_driver_seegmiller_qualified": False,
        "states": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-commits", type=int, default=3)
    args = parser.parse_args()
    try:
        result = certify(args.log.read_text(encoding="utf-8", errors="replace"),
                         minimum_commits=args.min_commits)
    except (OSError, QualificationError) as exc:
        result = {"schema": "CFDX-BFS-GLOBAL-CONTRACTION-V14", "strict_contraction_proved": False,
                  "bfs_driver_seegmiller_qualified": False, "reason": str(exc)}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        return 3
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
