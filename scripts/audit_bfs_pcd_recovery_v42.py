"""Verify a V42 PCD -> BlockLocal linear recovery from real production logs.

Fail closed on mismatched frozen inputs or an unqualified/fabricated result.
The certificate is ALWAYS linear-only; neither exit status nor SST activity
is interpreted as global physical convergence.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

KEY_VALUE = re.compile(r"([A-Za-z_][A-Za-z_0-9]*)=([^ ]+)")
HASH_KEYS = ("A", "b", "x0", "rows", "nnz")


def _events(path: Path, tag: str) -> list[dict[str, str]]:
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f"missing or empty production log: {path}")
    return [dict(KEY_VALUE.findall(line)) for line in
            path.read_text(encoding="utf-8").splitlines()
            if line.startswith(tag + " ")]


def _number(event: dict[str, str], key: str) -> float:
    value = float(event[key])
    if not math.isfinite(value):
        raise ValueError(f"nonfinite {key} in production log")
    return value


def _fingerprint(path: Path) -> dict[str, str]:
    rows = _events(path, "COUPLED_MATRIX_FINGERPRINT_V40")
    if not rows:
        raise ValueError(f"missing frozen A,b,x0 fingerprint in {path}")
    fingerprint = {key: rows[0][key] for key in HASH_KEYS}
    if len(fingerprint["A"]) != 16 or len(fingerprint["b"]) != 16 or len(fingerprint["x0"]) != 16:
        raise ValueError("malformed V40 FNV64 fingerprint")
    if int(fingerprint["rows"]) <= 0 or int(fingerprint["nnz"]) <= 0:
        raise ValueError("invalid original matrix dimensions")
    return fingerprint


def audit(base_log: Path, recovery_log: Path) -> dict:
    base_log, recovery_log = Path(base_log), Path(recovery_log)
    if base_log.resolve() == recovery_log.resolve():
        raise ValueError("baseline and recovery must be distinct solver logs")
    first = _fingerprint(base_log)
    second = _fingerprint(recovery_log)
    if first != second:
        raise ValueError("different frozen A,b,x0: cannot compare PCD and BlockLocal")
    baseline = _events(base_log, "COUPLED_BLOCK_TRUE_RESIDUAL_V16")
    if not baseline:
        raise ValueError("missing true original-system residual of failed PCD")
    failed = baseline[0]
    if int(failed["krylov_status"]) == 0 or _number(failed, "full_relative") <= 0:
        raise ValueError("baseline must prove a failed nonzero-residual PCD solve")
    retries = _events(recovery_log, "COUPLED_PCD_BLOCKLOCAL_RECOVERY_V42")
    if not retries:
        raise ValueError("V42 recovery never ran in production")
    linear_targets = []
    for r in retries:
        if (r.get("attempted") != "1" or r.get("adopted") != "1" or
            r.get("primary_model") != "PCD" or r.get("recovery_model") != "BlockLocal" or
            r.get("physical_qualified") != "0" or r.get("original_A_unchanged") != "1" or
            int(r["primary_status"]) == 0 or int(r["recovery_status"]) != 0 or
            int(r["recovery_iterations"]) < 1):
            raise ValueError("V42 log does not prove an honest BlockLocal recovery")
        target = _number(r, "unchanged_tolerance")
        measured = _number(r, "true_relative")
        if not (target > 0 and 0 <= measured <= target):
            raise ValueError("V42 retry does not meet the unchanged true-residual gate")
        linear_targets.append(measured)
    def proof(path: Path):
        return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    return {
        "schema": "CFDX_BFS_PCD_BLOCKLOCAL_RECOVERY_V42/1",
        "original_system_fingerprint": first,
        "baseline_pcd_true_relative": _number(failed, "full_relative"),
        "recovery_count": len(retries),
        "first_recovery_true_relative": linear_targets[0],
        "frozen_linear_system_converged": True,
        "pure_pcd_convergence_proven": False,
        "global_physical_convergence_proven": False,
        "qualification": "LINEAR_RECOVERY_ONLY",
        "baseline_log": proof(base_log), "recovery_log": proof(recovery_log),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pcd-log", type=Path, required=True)
    p.add_argument("--recovery-log", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    try:
        result = audit(a.pcd_log, a.recovery_log)
    except (ValueError, KeyError, OverflowError) as exc:
        p.error(str(exc))
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("PCD_BLOCKLOCAL_RECOVERY_V42_CERTIFIED linear_only=1 physical_qualified=0")


if __name__ == "__main__":
    main()
