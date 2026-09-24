#!/usr/bin/env python3
"""Create CI metadata/inventory for the numerical evidence artifact.

The actual solver output is captured by the authoritative verbose CTest step
in cfdx-validation.yml.  This helper never reruns a solver; it records the
test inventory and CI provenance even when the preceding gate failed.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    proc = subprocess.run(
        ["ctest", "--test-dir", str(args.build_dir), "-N"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    (args.output_dir / "ctest-inventory.log").write_text(proc.stdout, encoding="utf-8")

    metadata = {
        "status": "DIAGNOSTIC_ONLY",
        "inventory_return_code": proc.returncode,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "commit": os.environ.get("GITHUB_SHA"),
        "workflow": os.environ.get("GITHUB_WORKFLOW"),
        "run_id": os.environ.get("GITHUB_RUN_ID"),
        "run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT"),
        "ref": os.environ.get("GITHUB_REF"),
        "full_solver_log": "../ctest-full.log",
    }
    (args.output_dir / "diagnostic_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )

    # This step is deliberately non-gating and must succeed after a failed
    # solver test so the evidence upload is not suppressed.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
