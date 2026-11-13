"""V54 fail-closed forensic gate: frozen 6N Amat versus physical F(X).

This only inspects a recorded diagnostic; it cannot certify a physical BFS or
promote PCD after a fallback. Schur runs must not be called qualified unless
ALL basis checks (and the trial's true residual/physical gates) pass.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

REQUIRED = (
    "v54_frozen_amat", "v54_jacobian_evaluations", "v54_basis_audits",
    "v54_fidelity", "v54_halving", "v54_fidelity_ok", "v54_sided_nonsmooth",
    "v54_centered_columns", "v54_one_sided_columns", "private_accepted",
    "production_committed", "qualified", "linear_true_relative",
    "arnoldi_relative", "reason", "common_state_fingerprint",
)

def audit(text: str, *, expected_cells: int = 152,
          fidelity_limit: float = 0.05) -> dict:
    if not (math.isfinite(fidelity_limit) and 0 < fidelity_limit <= 0.5):
        raise ValueError("invalid new fidelity gate")
    lines = [x for x in text.splitlines() if x.startswith("FPHYS_JFNK_V52 ")]
    if len(lines) != 1:
        raise ValueError("missing or duplicated production V54 trial")
    fields = {}
    header, marker, reason = lines[0].partition(" reason=")
    if not marker or not reason.strip():
        raise ValueError("missing full physical rejection reason")
    fields["reason"] = reason.strip()
    for token in header.split()[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            if key in fields:
                raise ValueError(f"duplicate field {key}")
            fields[key] = value
    for key in REQUIRED:
        if key not in fields:
            raise ValueError(f"missing {key}")
    def integer(name: str) -> int:
        raw=fields[name]
        if not raw.isdigit():
            raise ValueError(f"invalid integer {name}")
        return int(raw)
    def number(name: str, *, finite: bool = True) -> float:
        try:
            value=float(fields[name])
        except ValueError as exc:
            raise ValueError(f"invalid numeric {name}") from exc
        if finite and not math.isfinite(value):
            raise ValueError(f"nonfinite {name}")
        if math.isnan(value):
            raise ValueError(f"NaN {name}")
        return value
    if integer("v54_frozen_amat")!=1 or integer("v54_basis_audits")<1:
        raise ValueError("frozen Jacobian missing or basis not audited")
    if integer("v54_jacobian_evaluations") < 6*expected_cells:
        raise ValueError("not all Jacobian columns evaluated")
    if integer("v54_centered_columns")+integer("v54_one_sided_columns") != 6*expected_cells:
        raise ValueError("incorrect column provenance")
    if integer("production_committed")!=0 or integer("qualified")!=0:
        raise ValueError("unqualified private Newton promoted to production")
    if not integer("common_state_fingerprint"):
        raise ValueError("missing common-state fingerprint")
    fidelity=number("v54_fidelity")
    halving=number("v54_halving")
    sided=number("v54_sided_nonsmooth")
    if min(fidelity,halving,sided)<0:
        raise ValueError("negative error metric")
    valid=fidelity<=fidelity_limit and halving<=fidelity_limit and sided<=fidelity_limit
    if integer("v54_fidelity_ok") != int(valid):
        raise ValueError("fidelity flag contradicts physical derivative checks")
    linear_true=number("linear_true_relative", finite=False)
    arnoldi=number("arnoldi_relative", finite=False)
    if not valid:
        if integer("private_accepted")!=0:
            raise ValueError("Newton step accepted despite unfaithful Jv")
        if "fidelity" not in fields["reason"]:
            raise ValueError("wrong rejection reason for derivative failure")
    if math.isfinite(linear_true) and linear_true<0:
        raise ValueError("negative true linear residual")
    if math.isfinite(arnoldi) and arnoldi<0:
        raise ValueError("negative Arnoldi norm")
    return {
        "status": "BLOCKED_JV_FIDELITY" if not valid else "JACOBIAN_FIDELITY_ONLY",
        "qualified": False,
        "global_newton_committed": False,
        "six_schur_6N_requalified": False,
        "cells": expected_cells,
        "frozen_jacobian_rows": 6*expected_cells,
        "jacobian_evaluations": integer("v54_jacobian_evaluations"),
        "krylov_basis_audits": integer("v54_basis_audits"),
        "fidelity": fidelity,
        "halving": halving,
        "sided_nonsmoothness": sided,
        "new_fidelity_limit": fidelity_limit,
        "linear_true_relative": linear_true if math.isfinite(linear_true) else None,
        "common_state_fingerprint": fields["common_state_fingerprint"],
    }

def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--log",required=True,type=Path)
    p.add_argument("--output",type=Path)
    p.add_argument("--cells",type=int,default=152)
    ns=p.parse_args()
    result=audit(ns.log.read_text(encoding="utf-8",errors="replace"),
                 expected_cells=ns.cells)
    serialized=json.dumps(result,indent=2,allow_nan=False,sort_keys=True)+"\n"
    if ns.output:
        ns.output.parent.mkdir(parents=True,exist_ok=True)
        ns.output.write_text(serialized,encoding="utf-8")
    print(serialized,end="")
if __name__=="__main__":main()
