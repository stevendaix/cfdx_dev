#!/usr/bin/env python3
"""Run the canonical CFDX benchmark campaign and emit machine-readable evidence.\n\nThe campaign is intended to be replayed on a configured validation build or by CI.

This driver does not invent or infer qualification.  It executes the repository's
existing validation entry points, captures stdout/stderr and records a maturity
state separately from the process exit code.

Default mode runs the short CTest regressions.  --full additionally executes
the available solver-level campaigns and their refinement series.  A successful
quick contract is reported as VERIFIED_CONTRACT, never as QUALIFIED.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path


@dataclass
class CaseResult:
    case_id: str
    title: str
    maturity: str
    command: list[str]
    return_code: int
    log: str
    notes: str


CASES = {
    "couette": (
        "LAM-COUETTE",
        "Couette flow",
        ["ctest", "--test-dir", "{build}", "-R", "^test_couette_quick$", "--output-on-failure"],
        "Short coupled SIMPLE/PISO/COUPLED regression. Full quantitative qualification remains tied to the declared mesh/QoI evidence.",
    ),
    "poiseuille": (
        "LAM-POISEUILLE",
        "Plane Poiseuille flow",
        ["ctest", "--test-dir", "{build}", "-R", "^test_poiseuille_quick$", "--output-on-failure"],
        "Short production-solver regression. The independent refinement/QoI campaign is the qualification evidence.",
    ),
    "ghia": (
        "INC-GHIA",
        "Ghia lid-driven cavity, Re=100",
        ["ctest", "--test-dir", "{build}", "-R", "^test_ghia_cavity_quick$", "--output-on-failure"],
        "Short Re=100 cavity regression. Full mesh/reference study is separate evidence.",
    ),
    "thermal": (
        "THERMAL",
        "Thermal / conduction / CHT verification subset",
        ["ctest", "--test-dir", "{build}", "-R", "^(test_thermal_vv|test_cht_validation)$", "--output-on-failure"],
        "Component/solver V&V entry points; a passing process does not by itself qualify full multi-region CHT.",
    ),
    "radiation": (
        "RADIATION",
        "Radiation and S2S verification subset",
        ["ctest", "--test-dir", "{build}", "-R", "^(test_radiation_vv|test_s2s_radiation_vv)$", "--output-on-failure"],
        "Radiation model verification entry points; coupled thermal-radiation qualification remains distinct.",
    ),
    "bfs": (
        "SEP-BFS",
        "Backward-facing step, Re=200",
        ["ctest", "--test-dir", "{build}", "-R", "^test_bfs_quick$", "--output-on-failure"],
        "Quick mesh/BC/contract smoke only. It is not a reattachment-length qualification.",
    ),
    "naca0012": (
        "EXT-NACA0012",
        "NACA0012 laminar alpha=0 contract",
        ["ctest", "--test-dir", "{build}", "-R", "^test_naca0012_quick$", "--output-on-failure"],
        "Quick geometry/BC contract only. Force polar qualification is still a separate campaign.",
    ),
    "n9_s8_naca0012": (
        "N9-S8-NACA0012",
        "N9-S8 NASA/TMR NACA0012 laminar contract",
        ["ctest", "--test-dir", "{build}", "-R", "^test_n9_s8_naca0012_contract$", "--output-on-failure"],
        "External NASA/TMR mesh provenance and frozen Re=1000 alpha=0 contract. This is not numerical qualification.",
    ),
    "vmfl036": (
        "FORCE-VMFL036",
        "Axisymmetric sphere, Re=100",
        ["ctest", "--test-dir", "{build}", "-R", "^test_vmfl036_axisymmetric$", "--output-on-failure"],
        "Physical axisymmetric campaign. Its result is parsed from the retained log; a failure remains a failure.",
    ),
}


def run(command: list[str], log_path: Path) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as stream:
        proc = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, text=True, check=False)
    return proc.returncode


def full_command(key: str, build: Path) -> list[str] | None:
    if key == "couette":
        return [str(build / "test_phase9_acceptance")]
    if key == "poiseuille":
        return [str(build / "test_poiseuille_diagnostics")]
    if key == "ghia":
        return [str(build / "test_ghia_cavity")]
    if key == "thermal":
        return [str(build / "test_thermal_vv")]
    if key == "radiation":
        return [str(build / "test_radiation_vv")]
    if key == "vmfl036":
        exe = build / "test_vmfl036_axisymmetric"
        return [str(exe)] if exe.exists() else None
    if key == "bfs":
        meshes = [build / "validation_meshes" / f"bfs_re200_n{n}.h5" for n in (16, 32, 64)]
        exe = build / "test_bfs_qualification"
        return [str(exe), *(str(p) for p in meshes)] if exe.exists() and all(p.exists() for p in meshes) else None
    if key == "naca0012":
        meshes = [build / "validation_meshes" / f"naca0012_laminar_n{n}.h5" for n in (64, 128, 256)]
        exe = build / "test_naca0012_qualification"
        return [str(exe), *(str(p) for p in meshes)] if exe.exists() and all(p.exists() for p in meshes) else None
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    parser.add_argument("--output-dir", type=Path, default=Path("build/benchmark-campaign"))
    parser.add_argument("--full", action="store_true", help="run solver-level campaigns instead of short CTest contracts")
    parser.add_argument("--case", choices=sorted(CASES), action="append")
    args = parser.parse_args()

    if not (args.build_dir / "CTestTestfile.cmake").exists():
        print(f"Build directory is not configured: {args.build_dir}", file=sys.stderr)
        return 2

    keys = args.case or list(CASES)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results: list[CaseResult] = []

    for key in keys:
        ident, title, template, notes = CASES[key]
        if args.full:
            command = full_command(key, args.build_dir)
            if command is None:
                results.append(CaseResult(
                    ident, title, "NOT_RUN", [], -1,
                    str(args.output_dir / f"{key}.log"),
                    "Full campaign entry point or required refinement meshes are unavailable.",
                ))
                continue
        else:
            command = [x.format(build=str(args.build_dir)) for x in template]

        log = args.output_dir / f"{key}.log"
        rc = run(command, log)

        if rc != 0:
            maturity = "FAIL"
        elif args.full:
            maturity = "EXECUTED"
        else:
            maturity = "VERIFIED_CONTRACT"

        results.append(CaseResult(
            ident, title, maturity, command, rc, str(log),
            notes,
        ))

    payload = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "mode": "full" if args.full else "quick",
        "repository_rule": "No quick contract is promoted to qualification.",
        "cases": [asdict(r) for r in results],
    }
    (args.output_dir / "results.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )

    print(json.dumps(payload, indent=2))
    return 1 if any(r.return_code != 0 for r in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
