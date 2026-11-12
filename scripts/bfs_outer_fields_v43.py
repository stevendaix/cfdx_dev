"""Fail-closed CFDX V43 audit of a periodic *field* trajectory.

An oscillating norm is not evidence for an oscillating state vector. This
reader requires physically committed states with the matching V43 field trace,
checks its own n-4 and n-1 hypothesis, and records whether a putative period
four is supported for U, gauge-invariant p, conservative phi, k, omega, nut.
No trace is ever counted as a converged BFS solution.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re

RECORD = re.compile(r"^OUTER_FIELD_PERIOD_V43\s+(.*)$")
COMMITTED = re.compile(r"^GLOBAL_NONLINEAR_COMMITTED_V13\s+(.*)$")
FIELDS = ("Ux_four", "Uy_four", "Uz_four", "p_gauge_four", "phi_four",
          "k_four", "omega_four", "nut_four")


def parse_kv(blob: str) -> dict[str, str]:
    result = {}
    for item in blob.split():
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        if key in result:
            raise ValueError(f"duplicate diagnostic key {key}")
        result[key] = value
    return result


def audit(text: str, require_field_cycle: bool = False) -> dict:
    records = []
    committed: dict[int, float] = {}
    for line in text.splitlines():
        if match := COMMITTED.match(line):
            values = parse_kv(match.group(1))
            index = int(values["iteration"])
            merit = float(values["merit"])
            if not math.isfinite(merit) or merit < 0 or index in committed:
                raise ValueError("invalid or repeated committed physical merit")
            if values.get("qualified") != "0":
                raise ValueError("unsubstantiated physical qualification claim")
            committed[index] = merit
        if match := RECORD.match(line):
            row = parse_kv(match.group(1))
            index = int(row["iteration"])
            measurements = [float(row[name]) for name in FIELDS]
            four = float(row["max_four"])
            previous = float(row["max_previous"])
            if any(not math.isfinite(value) or value < 0
                   for value in [*measurements, four, previous]):
                raise ValueError("nonfinite / negative committed field comparison")
            if abs(four - max(measurements)) > 1e-10 * max(1.0, four):
                raise ValueError("V43 field comparisons are internally inconsistent")
            if row.get("physical_qualified") != "0":
                raise ValueError("field cycle observation cannot certify physical convergence")
            if int(row["same_phase_count"]) < 0:
                raise ValueError("negative phase count")
            records.append({"iteration": index,
                            "relative_four_back": dict(zip(FIELDS, measurements)),
                            "maximum_four_back": four,
                            "maximum_previous": previous,
                            "same_phase_count": int(row["same_phase_count"]),
                            "period_four_candidate": row["period_four_candidate"] == "1",
                            "worst_channel": int(row["worst_channel"]),
                            "worst_index": int(row["worst_index"])})
    if len(committed) < 5 or not records:
        raise ValueError("missing committed physical states or V43 field period evidence")
    indices = sorted(committed)
    if indices != list(range(indices[0], indices[-1] + 1)):
        raise ValueError("committed physical iteration history has gaps")
    if len(records) != len(committed) - 4:
        raise ValueError("missing n-4 comparison or uncommitted field sample")
    if [r["iteration"] for r in records] != indices[4:]:
        raise ValueError("field trace not paired to accepted iteration history")
    chain = 0
    for record in records:
        near = (record["maximum_four_back"] <= 1e-3 and
                record["maximum_previous"] >=
                    max(1e-4, 10.0 * record["maximum_four_back"]))
        chain = chain + 1 if near else 0
        if record["same_phase_count"] != chain or record["period_four_candidate"] != (chain >= 8):
            raise ValueError("inconsistent V43 four-phase confirmation counter")
    physical_cycle = chain >= 8
    if require_field_cycle and not physical_cycle:
        raise ValueError("field-period-four proof not established")
    return {
        "contract": "FIELD_PERIOD_DIAGNOSTIC_ONLY",
        "physical_converged": False,
        "field_period_four_candidate": physical_cycle,
        "consecutive_four_step_returns": chain,
        "committed_states": len(committed),
        "first_iteration": indices[0],
        "last_iteration": indices[-1],
        "merit_first": committed[indices[0]],
        "merit_last": committed[indices[-1]],
        "last_report": records[-1],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--require-field-cycle", action="store_true")
    args = parser.parse_args()
    data = audit(args.log.read_text(), args.require_field_cycle)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    print(f"FIELD_PERIOD_V43 {data['contract']} "
          f"period_four_candidate={int(data['field_period_four_candidate'])} "
          f"committed={data['committed_states']}")


if __name__ == "__main__":
    main()
