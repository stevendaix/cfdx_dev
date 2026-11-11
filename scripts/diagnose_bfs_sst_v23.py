"""Compare V22 / V23 SST Krylov candidates, never certify BFS or NASA physics."""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

PAIRS = re.compile(r"([a-zA-Z_][a-zA-Z_0-9]*)=([^\s]+)")


def first_retry_samples(text: str) -> dict:
    steps: list[dict] = []
    diag: list[dict] = []
    gates: list[dict] = []
    for line in text.splitlines():
        if "CFDX nonlinear iteration 1 rejected" in line and steps:
            break  # Compare identical FIRST outer attempt only.
        if line.startswith("SST_JFNK_V22 "):
            steps.append(dict(PAIRS.findall(line)))
        if line.startswith("SST_TRANSPORT_DIAGNOSTIC "):
            diag.append(dict(PAIRS.findall(line)))
        if line.startswith("SST_ACCEPTANCE_GATES "):
            gates.append(dict(PAIRS.findall(line)))
    if not steps or not diag or not gates:
        raise ValueError("missing independent physical SST diagnostics")
    def sample10(rows: list[dict]) -> dict:
        return next((r for r in rows if r.get("iteration") == "10"), {})
    s, d, g = sample10(steps), sample10(diag), sample10(gates)
    if not (s and d and g):
        raise ValueError("no complete physical SST iteration 10")
    if s.get("selected") == "jfnk_block" and (
        s.get("block_proposed") != "1"
        or float(s.get("block_linear_true_relative", "inf")) > 0.25
    ):
        raise ValueError("uncertified block Newton direction marked selected")
    if s.get("selected") == "jfnk" and (
        s.get("proposed") != "1"
        or float(s.get("linear_true_relative", "inf")) > 0.25
    ):
        raise ValueError("uncertified Jacobi Newton direction marked selected")
    fields = ("k_returned_nonlinear_raw", "omega_returned_nonlinear_raw",
              "k_local_relative", "omega_local_relative")
    metrics = {k: float(d[k]) for k in fields}
    metrics["merit"] = max(metrics["k_local_relative"],
                           metrics["omega_local_relative"])
    if any(not math.isfinite(v) or v < 0 for v in metrics.values()):
        raise ValueError("invalid physical residual sample")
    return {
        "selected": s["selected"],
        "block_probe": s.get("block_probe", "0") == "1",
        "block_proposed": s.get("block_proposed", "0") == "1",
        "block_linear_true_relative": (
            float(s["block_linear_true_relative"]) if s.get("block_probe") == "1"
            and math.isfinite(float(s.get("block_linear_true_relative", "nan"))) else None
        ),
        "coupled_blocks": int(s.get("coupled_blocks", "0")),
        "jacobi_linear_true_relative": float(s["linear_true_relative"]),
        "tolerance": float(g["tolerance"]),
        "outer_update_ok": g["outer_update_ok"] == "1",
        "returned_equations_ok": g["returned_equations_ok"] == "1",
        **metrics,
    }


def compare(baseline: str, experimental: str) -> dict:
    a, b = first_retry_samples(baseline), first_retry_samples(experimental)
    if a["tolerance"] != b["tolerance"]:
        raise ValueError("non-comparable nonlinear SST tolerances")
    if not b["block_probe"] or b["coupled_blocks"] == 0:
        raise ValueError("V23 blocks not actually exercised")
    result = {
        "schema": "CFDX-BFS-SST-V23",
        "baseline": a,
        "v23": b,
        "v23_nonworse_merit_at_iteration_10": b["merit"] <= a["merit"] + 1e-12,
        "physical_convergence_certified": False,
        "nasa_certified": False,
        "warning": "Ten SST iterations of one outer attempt are diagnostics, not converged BFS",
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        report = compare(args.baseline.read_text(encoding="utf-8", errors="replace"),
                         args.candidate.read_text(encoding="utf-8", errors="replace"))
    except (OSError, KeyError, ValueError) as exc:
        parser.error(str(exc))
    doc = json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(doc, encoding="utf-8")
    else:
        print(doc, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
