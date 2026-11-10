"""BFS SST iteration forensic comparison, with fail-closed qualification.

The log may contain V19 convergence records or optional V20 row budgets.
This tool reports numerical evidence; it NEVER certifies NASA observables.
"""
from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path

_KEYS = ("iteration", "field", "cell", "convection", "diffusion", "production",
         "destruction", "cross", "balance", "closure_error", "source_identity")


def parse_kv(line: str) -> dict[str, str]:
    return dict(re.findall(r"([A-Za-z_]+)=([^\s]+)", line))


def analyze(text: str) -> dict:
    records: list[dict] = []
    committed: list[dict] = []
    for line in text.splitlines():
        if line.startswith("SST_ASSEMBLY_V20 "):
            raw = parse_kv(line)
            if not set(_KEYS).issubset(raw):
                raise ValueError("incomplete SST assembly record")
            sample = {"iteration": int(raw["iteration"]), "field": raw["field"],
                      "cell": int(raw["cell"])}
            if sample["iteration"] < 1 or sample["cell"] < 0 or sample["field"] not in ("k", "omega"):
                raise ValueError("invalid SST row identifier")
            for key in _KEYS[3:]:
                sample[key] = float(raw[key])
                if not math.isfinite(sample[key]):
                    raise ValueError("non-finite SST row budget")
            records.append(sample)
        elif line.startswith("GLOBAL_NONLINEAR_COMMITTED_V13 "):
            raw = parse_kv(line)
            required = ("iteration", "U", "p_continuity", "k", "omega", "merit", "converged")
            if not set(required).issubset(raw):
                raise ValueError("incomplete global committed state")
            state = {"iteration": int(raw["iteration"])}
            for key in required[1:-1]:
                state[key] = float(raw[key])
                if not math.isfinite(state[key]):
                    raise ValueError("nonfinite committed residual")
            state["converged"] = int(raw["converged"])
            if state["converged"] not in (0, 1):
                raise ValueError("invalid converged marker")
            if committed and state["iteration"] <= committed[-1]["iteration"]:
                raise ValueError("nonmonotone committed iterations")
            committed.append(state)
    return {"schema": "CFDX-BFS-SST-V20", "row_budgets": records,
            "committed_states": committed,
            "three_state_contraction_proved": False,
            "driver_seegmiller_qualified": False,
            "reason": "Independent global physical and NASA qualification required"}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--log", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    try:
        record = analyze(args.log.read_text(encoding="utf-8", errors="replace"))
    except (OSError, ValueError) as exc:
        record = {"schema": "CFDX-BFS-SST-V20", "error": str(exc),
                  "three_state_contraction_proved": False,
                  "driver_seegmiller_qualified": False}
        status = 3
    else:
        status = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
