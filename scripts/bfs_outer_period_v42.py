"""Detect repeating *global residual merits* in committed BFS CFDX states.

This is an evidence-only observer. It does not prove periodicity of the full
fields and can never certify global CFD convergence. No solver flags change.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

_EVENT = re.compile(r"([A-Za-z_][A-Za-z_0-9]*)=([^ ]+)")
_TAG = "GLOBAL_NONLINEAR_COMMITTED_V13 "


def audit(path: Path, *, min_cycles: int = 4, max_period: int = 8,
          stability: float = 0.01, min_amplitude: float = 0.1) -> dict:
    path = Path(path)
    if not path.is_file() or not path.stat().st_size:
        raise ValueError("missing or empty CFDX production log")
    if min_cycles < 3 or max_period < 2 or not (0 < stability < 1) or not (0 < min_amplitude < 1):
        raise ValueError("invalid periodicity detection controls")
    records: list[dict] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(_TAG):
            continue
        fields = dict(_EVENT.findall(line))
        required = ("iteration", "U", "p_continuity", "k", "omega", "merit", "converged")
        if any(field not in fields for field in required):
            raise ValueError("incomplete committed-state evidence")
        index = int(fields["iteration"])
        measurements = {key: float(fields[key]) for key in required[1:-1]}
        if index <= 0 or any(not math.isfinite(v) or v < 0 for v in measurements.values()):
            raise ValueError("invalid global iteration or physical residual")
        if fields["converged"] not in ("0", "1"):
            raise ValueError("invalid global qualification marker")
        if records and index != records[-1]["iteration"] + 1:
            raise ValueError("non-contiguous committed global history")
        records.append({"iteration": index, **measurements,
                        "converged": fields["converged"] == "1"})
    if not records:
        raise ValueError("no committed global nonlinear states")
    # Verify that this evidence never falsely equates periodic residuals to
    # complete field-state cycles or Newton line-search accepted steps.
    result = {
        "schema": "CFDX_BFS_GLOBAL_MERIT_PERIOD_V42/1",
        "log_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "committed_states": len(records),
        "last_iteration": records[-1]["iteration"],
        "last_global_merit": records[-1]["merit"],
        "full_field_cycle_proven": False,
        "global_physical_convergence_proven": False,
        "qualification": "MERIT_TRACE_ONLY",
        "period": None,
    }
    if records[-1]["converged"]:
        result["status"] = "LAST_STATE_QUALIFIED_NOT_A_STALL"
        return result
    for period in range(2, max_period + 1):
        window = min_cycles * period
        if len(records) < window:
            continue
        tail = records[-window:]
        merit = [r["merit"] for r in tail]
        phases = [merit[j::period] for j in range(period)]
        means = [sum(v)/len(v) for v in phases]
        if not all(v > 0 for v in means):
            continue
        spread = max((max(v) - min(v))/m for v, m in zip(phases, means))
        amplitude = (max(means)-min(means))/max(means)
        if spread <= stability and amplitude >= min_amplitude:
            result.update({
                "status": "REPEATING_UNCONVERGED_GLOBAL_MERIT",
                "period": period,
                "phase_means": means,
                "relative_phase_spread_max": spread,
                "relative_cycle_amplitude": amplitude,
                "window_first_iteration": tail[0]["iteration"],
                "window_last_iteration": tail[-1]["iteration"],
            })
            return result
    result["status"] = "NO_PERIOD_DETECTED"
    return result


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--log", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--require-period", type=int, choices=tuple(range(2,9)))
    args = p.parse_args()
    try:
        outcome = audit(args.log)
        if args.require_period is not None and outcome["period"] != args.require_period:
            raise ValueError(f"period {args.require_period} not supported by input")
    except (ValueError, KeyError, OverflowError) as exc:
        p.error(str(exc))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(outcome, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("BFS_GLOBAL_MERIT_PERIOD_V42 status=" + outcome["status"] +
          " physical_qualified=0")


if __name__ == "__main__":
    main()
