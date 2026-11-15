"""Reproduce V64 read-only 2x2 momentum/SST attribution on the SAME BFS input.

Runs one native C++ baseline (V51 only) and one V64 C++ diagnostic (V51+V64).
The post-SST nonlinear solver remains unchanged. Non-convergence exit code is
recorded, not turned into a physical PASS. Complete input SHA and original
committed nonlinear evidence must match; partial/timeout runs fail closed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

from analyze_bfs_causal_v64 import parse

AUDITED_PREFIXES = (
    "BFS_PHASE_V10 ",
    "FPHYS_CANDIDATE_V51 ",
    "POST_SST_MOMENTUM ",
    "INCOMPRESSIBLE_ITER ",
    "GLOBAL_NONLINEAR_V8 ",
    "GLOBAL_NONLINEAR_COMMITTED_V13 ",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run(
    solver: Path,
    case: Path,
    folder: Path,
    diagnostic: bool,
    timeout: float,
    iterations: int,
) -> dict:
    folder.mkdir(parents=True, exist_ok=False)
    log = folder / "solver.log"
    env = {
        k: v
        for k, v in os.environ.items()
        if not (
            k.startswith(("CFDX_FPHYS_", "CFDX_COUPLED_SELFP_", "CFDX_POST_SST_", "CFDX_BFS_"))
        )
    }
    env["CFDX_FPHYS_CANDIDATE_V51"] = "1"
    env["CFDX_BFS_STAGE_TRACE"] = "1"
    if diagnostic:
        env["CFDX_FPHYS_CAUSAL_V64"] = "1"
    cmd = [
        str(solver),
        "--mesh",
        str(case),
        "--iterations",
        str(iterations),
        "--experimental-sst",
        "--output-dir",
        str(folder / "run"),
    ]
    with log.open("w") as out:
        completed = subprocess.run(
            cmd,
            env=env,
            stdout=out,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
    return {
        "exit_code": completed.returncode,
        "log": str(log),
        "log_sha256": sha256(log),
        "cmd": cmd,
        "text": log.read_text(errors="replace"),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--solver", type=Path, required=True)
    ap.add_argument("--case", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--iterations", type=int, default=1)
    ap.add_argument("--wall-timeout", type=float, default=90)
    args = ap.parse_args()
    if args.iterations != 1:
        ap.error(
            "bounded common-input parity comparison requires exactly one iteration"
        )
    out = args.out.resolve()
    if out.exists():
        ap.error("refusing existing evidence output directory")
    solver = args.solver.resolve(strict=True)
    case = args.case.resolve(strict=True)
    out.mkdir(parents=True)
    base = run(
        solver, case, out / "baseline", False, args.wall_timeout, args.iterations
    )
    trial = run(solver, case, out / "causal", True, args.wall_timeout, args.iterations)
    blocks = {
        key: [
            line
            for line in data["text"].splitlines()
            if line.startswith(AUDITED_PREFIXES)
        ]
        for key, data in (("baseline", base), ("causal", trial))
    }
    same_input = sha256(case)
    same_trace = blocks["baseline"] == blocks["causal"]
    complete = all(
        any(line.startswith("GLOBAL_NONLINEAR_COMMITTED_V13 ") for line in blocks[key])
        for key in blocks
    )
    has_sst = all(
        any("phase=SST_TRANSPORT_END" in line for line in blocks[key]) for key in blocks
    )
    if (
        not complete
        or not has_sst
        or not same_trace
        or base["exit_code"] != trial["exit_code"]
    ):
        raise RuntimeError("V64 evidence incomplete or baseline physics was changed")
    causal = parse(trial["text"])
    if len(causal["iterations"]) != 1:
        raise RuntimeError("V64 expected one diagnostic iteration")
    manifest = {
        "schema": "CFDX-BFS-CAUSAL-REPLAY/V64",
        "same_physical_case_sha256": same_input,
        "base_binary_sha256": sha256(solver),
        "baseline_physical_trace_identical": same_trace,
        "baseline_solver_exit": base["exit_code"],
        "trial_solver_exit": trial["exit_code"],
        "baseline_log_sha256": base["log_sha256"],
        "trial_log_sha256": trial["log_sha256"],
        "physical_convergence_proven": False,
        "diagnostic": causal,
    }
    (out / "result.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print("V64_REPLAY_PARITY_VERIFIED diagnostic_only=1 physical_converged=0")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
        print("V64_CAUSAL_REPLAY_REJECTED: " + str(exc), file=sys.stderr)
        raise SystemExit(2)
