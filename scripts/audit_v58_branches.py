"""Fail-closed proof check for V58 face-sided differentiation.

Never treats a local flux derivative as a Newton/BFS physical qualification.
"""
import argparse
import json
import math
import re
from pathlib import Path


def parse_lines(lines):
    groups = {}
    detail = {}
    for line in lines:
        if not line.startswith(("FPHYS_BRANCH_V58 ", "FPHYS_BRANCH_FACE_V58 ")):
            continue
        parts = dict(re.findall(r"([a-z_]+)=([^\s]+)", line))
        required = {"iteration", "component", "native_rc_reclosure_at_offsets",
                    "physical_qualified", "committed"}
        if not required <= parts.keys() or any(parts[key] != "0" for key in
                                                  ("native_rc_reclosure_at_offsets",
                                                   "physical_qualified", "committed")):
            raise ValueError("V58 unsafe or incomplete proof claim")
        iteration = int(parts["iteration"])
        component = int(parts["component"])
        if iteration < 0 or component not in range(3):
            raise ValueError("V58 invalid iteration/component")
        key = iteration, component
        if line.startswith("FPHYS_BRANCH_V58 "):
            fields = ("faces", "ambiguous", "effective_kinks", "unresolved", "inert",
                      "flux_crossings", "roundoff_limited", "sided", "halving", "certifiable")
            if not set(fields) <= parts.keys() or key in groups:
                raise ValueError("V58 duplicate or truncated face group")
            f = {name: float(parts[name]) if name in ("sided", "halving")
                 else int(parts[name]) for name in fields}
            if (f["faces"] <= 0 or f["certifiable"] not in (0, 1) or
                any(not math.isfinite(f[s]) or f[s] < 0 for s in ("sided", "halving")) or
                any(f[s] < 0 or f[s] > f["faces"] for s in
                    ("ambiguous", "effective_kinks", "unresolved", "inert", "flux_crossings", "roundoff_limited")) or
                f["inert"] > f["ambiguous"] or
                (f["certifiable"] and (f["effective_kinks"] or f["unresolved"] or
                                         f["sided"] > 1e-2 or f["halving"] > 1e-2))):
                raise ValueError("V58 inconsistent evidence")
            groups[key] = f
        else:
            fields = ("face", "owner", "neighbour", "phi", "dphi", "right", "left",
                      "sided", "halving", "kink", "unresolved")
            if not set(fields) <= parts.keys():
                raise ValueError("V58 incomplete face sample")
            f = {name: float(parts[name]) if name in
                 ("phi", "dphi", "right", "left", "sided", "halving")
                 else int(parts[name]) for name in fields}
            if not all(math.isfinite(f[s]) for s in
                       ("phi", "dphi", "right", "left", "sided", "halving")):
                raise ValueError("V58 nonfinite face sample")
            detail.setdefault(key, []).append(f)
    if not groups:
        raise ValueError("V58 no native face audits found")
    iterations = sorted({k[0] for k in groups})
    for iteration in iterations:
        if set(k[1] for k in groups if k[0] == iteration) != set(range(3)):
            raise ValueError("V58 missing component in physical chain")
    for key, f in groups.items():
        records = detail.get(key, [])
        if len(records) != min(3, f["faces"]):
            raise ValueError("V58 missing top face evidence")
        for face in records:
            if face["face"] < 0 or face["face"] >= f["faces"] or face["owner"] < 0:
                raise ValueError("V58 invalid face topology")
    return {"status": "BLOCKED_GLOBAL_JV_NOT_AUDITED",
            "iterations": iterations, "audited_components": len(groups),
            "effective_kinks": sum(g["effective_kinks"] for g in groups.values()),
            "unresolved": sum(g["unresolved"] for g in groups.values()),
            "global_newton_qualified": False,
            "native_rc_reclosure_at_offsets": False}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--log", type=Path, required=True)
    p.add_argument("--json", type=Path)
    args = p.parse_args()
    result = parse_lines(args.log.read_text(errors="replace").splitlines())
    message = json.dumps(result, indent=2, sort_keys=True)
    if args.json:
        args.json.write_text(message + "\n")
    print(message)


if __name__ == "__main__":
    main()
