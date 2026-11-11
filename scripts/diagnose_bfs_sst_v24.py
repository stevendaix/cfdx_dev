#!/usr/bin/env python3
"""Fail-closed V23/V24 first-outer-attempt SST report; never certify physics."""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

PAIR = re.compile(r"([a-zA-Z_][a-zA-Z_0-9]*)=([^\s]+)")


def first_attempt(text: str) -> dict:
    steps: dict[int, dict] = {}
    diag: dict[int, dict] = {}
    gates: dict[int, dict] = {}
    audits: dict[int, dict] = {}
    for line in text.splitlines():
        if line.startswith("CFDX nonlinear iteration 1 rejected"):
            break
        row = dict(PAIR.findall(line))
        if "iteration" not in row:
            continue
        if line.startswith("SST_JFNK_V22 "):
            steps[int(row["iteration"])] = row
        elif line.startswith("SST_TRANSPORT_DIAGNOSTIC "):
            diag[int(row["iteration"])] = row
        elif line.startswith("SST_ACCEPTANCE_GATES "):
            gates[int(row["iteration"])] = row
        elif line.startswith("SST_JACOBIAN_AUDIT_V24 ") and row.get("iteration") == "10":
            if "max_relative_difference" in row:
                audits[int(row["cell"])] = row
    if 10 not in steps or 10 not in diag or 10 not in gates:
        raise ValueError("missing complete first-attempt physical iteration 10")
    s, d, g = steps[10], diag[10], gates[10]
    for method, prefix in (("jfnk", ""), ("jfnk_block", "block_"),
                           ("jfnk_spatial", "spatial_")):
        if s["selected"] == method:
            key = "proposed" if not prefix else f"{prefix}proposed"
            relative = ("linear_true_relative" if not prefix else
                        f"{prefix}linear_true_relative")
            if s.get(key) != "1" or float(s.get(relative, "inf")) > 0.25:
                raise ValueError(f"uncertified {method} direction selected")
    data = {}
    for field in ("k_returned_nonlinear_raw", "omega_returned_nonlinear_raw",
                  "k_local_relative", "omega_local_relative"):
        data[field] = float(d[field])
    if any(not math.isfinite(v) or v < 0 for v in data.values()):
        raise ValueError("invalid physical residual sample")
    data["merit"] = max(data["k_local_relative"], data["omega_local_relative"])
    errors = {str(cell): float(v["max_relative_difference"])
              for cell, v in sorted(audits.items())}
    if any(not math.isfinite(v) or v < 0 for v in errors.values()):
        raise ValueError("invalid local Jacobian audit")
    return {
        "selected": s["selected"],
        "tolerance": float(g["tolerance"]),
        "outer_equations_ok": g["returned_equations_ok"] == "1",
        "outer_update_ok": g["outer_update_ok"] == "1",
        "spatial_probe": s.get("spatial_probe") == "1",
        "spatial_proposed": s.get("spatial_proposed") == "1",
        "spatial_offcell": int(s.get("spatial_offcell", "0")),
        "spatial_linear_true_relative": (float(s["spatial_linear_true_relative"])
            if s.get("spatial_probe") == "1" else None),
        "jacobi_linear_true_relative": float(s["linear_true_relative"]),
        "block_linear_true_relative": (float(s["block_linear_true_relative"])
            if s.get("block_probe") == "1" else None),
        "jacobian_errors_at_iteration_10": errors,
        **data,
    }


def compare(v23_log: str, v24_log: str) -> dict:
    a, b = first_attempt(v23_log), first_attempt(v24_log)
    if a["tolerance"] != b["tolerance"]:
        raise ValueError("noncomparable physical SST tolerances")
    if not b["spatial_probe"] or b["spatial_offcell"] == 0:
        raise ValueError("V24 spatial preconditioner not actually exercised")
    if not b["spatial_proposed"] and b["selected"] == "jfnk_spatial":
        raise ValueError("V24 selected a failed Newton direction")
    return {
        "schema": "CFDX-BFS-SST-V24-UNQUALIFIED",
        "v23": a, "v24": b,
        "v24_merit_not_worse_at_iteration_10": b["merit"] <= a["merit"] + 1e-12,
        "v24_spatial_linear_better_than_v23_block": (
            b["spatial_linear_true_relative"] is not None and
            b["block_linear_true_relative"] is not None and
            b["spatial_linear_true_relative"] < b["block_linear_true_relative"]),
        "sst_full_convergence_certified": False,
        "bfs_global_convergence_certified": False,
        "nasa_physical_qualification_certified": False,
        "warning": "Only ten SST iterations on one outer attempt. Missing full convergence and validation.",
    }


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline",type=Path,required=True)
    p.add_argument("--candidate",type=Path,required=True)
    p.add_argument("--output",type=Path)
    opt=p.parse_args()
    try:
        report=compare(opt.baseline.read_text(errors="replace"),
                       opt.candidate.read_text(errors="replace"))
    except (OSError,KeyError,ValueError) as exc:
        p.error(str(exc))
    doc=json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n"
    if opt.output:
        opt.output.write_text(doc,encoding="utf-8")
    else:
        print(doc,end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
