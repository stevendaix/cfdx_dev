#!/usr/bin/env python3
"""Run CFDX numerical-model verification campaigns.

Phase 13 deliberately separates executable evidence from the closure decision:
software tests must pass; hardware-only evidence (CUDA/OOC) is reported as
BLOCKED when the required device is unavailable rather than being counted as
PASS. Use --strict to fail when the documented Phase-13 closure matrix still
contains pending campaigns.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
import xml.etree.ElementTree as ET


P13_LABEL = "phase13"

# Campaigns with executable evidence already present in the repository.
SOFTWARE_CAMPAIGNS = {
    "13.1": ["test_mms_scalar_diffusion"],
    "13.2": ["test_mesh_refinement_order", "test_mms_scalar_diffusion"],
    "13.3": ["test_core_temporal", "test_temporal_physics"],
    "13.4": ["test_m1_m4_validation", "test_m2_m4_solver_validation", "test_cht_validation"],
    "13.5": ["test_matrix_free_fv_operator"],
    "13.6": ["test_analytical_benchmarks", "test_benchmark_matrix"],
    "13.7": ["test_steady_incompressible_solver", "test_phase9_acceptance"],
    "13.8": ["test_cht_validation", "test_level_c_coupled_verification"],
    "13.9": ["test_phase10_turbulence_hardening", "test_thermophysical_models_vv"],
    "13.10": ["test_level_c_coupled_verification", "test_analytical_benchmarks"],
    "13.11": ["test_phase8_mpi_poisson"],
}

# These require capabilities not guaranteed on a normal CPU runner.
HARDWARE_CAMPAIGNS = {
    "13.12": "real CUDA CPU/GPU equivalence for production M1-M4 solvers",
    "13.13": "real CUDA out-of-core CFD equivalence beyond device VRAM",
}


def run(cmd: list[str], cwd: Path | None = None) -> int:
    print("$", " ".join(cmd))
    return subprocess.run(cmd, cwd=cwd, check=False).returncode


def write_junit(path: Path, results: list[tuple[str, int]]) -> None:
    suite = ET.Element(
        "testsuite",
        {
            "name": "cfdx-phase13",
            "tests": str(len(results)),
            "failures": str(sum(rc != 0 for _, rc in results)),
            "errors": "0",
        },
    )
    for name, rc in results:
        case = ET.SubElement(suite, "testcase", {"name": name})
        if rc != 0:
            ET.SubElement(case, "failure", {"message": f"exit code {rc}"})
    path.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(suite).write(path, encoding="utf-8", xml_declaration=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    parser.add_argument("--junit", type=Path)
    parser.add_argument("--json", type=Path)
    parser.add_argument(
        "--strict",
        action="store_true",
        help="also fail while any documented Phase-13 campaign lacks closure evidence",
    )
    args = parser.parse_args()

    if not (args.build_dir / "CTestTestfile.cmake").exists():
        print(f"Build directory is not configured: {args.build_dir}", file=sys.stderr)
        return 2

    results: list[tuple[str, int]] = []
    ctest = ["ctest", "--test-dir", str(args.build_dir), "-L", P13_LABEL, "--output-on-failure"]
    rc = run(ctest)
    results.append(("phase13-labeled-software-campaign", rc))

    print("\nPHASE 13 CAMPAIGN MATRIX")
    print("========================")
    evidence = {}
    for item, tests in SOFTWARE_CAMPAIGNS.items():
        evidence[item] = {"status": "EXECUTABLE", "tests": tests}
        print(f"{item}: EXECUTABLE -> {', '.join(tests)}")

    for item, description in HARDWARE_CAMPAIGNS.items():
        evidence[item] = {
            "status": "BLOCKED_WITHOUT_HARDWARE",
            "requirement": description,
        }
        print(f"{item}: BLOCKED_WITHOUT_HARDWARE -> {description}")

    pending = [
        "13.1 full NS/transient/turbulence/thermal/radiation MMS",
        "13.2 production-scheme systematic order matrix",
        "13.3 transient PDE order for every production integrator",
        "13.4 independent integral conservation across all M1-M4 physics",
        "13.5 production assembled/matrix-free solution equivalence",
        "13.6 production conditioning/scale robustness",
        "13.7 pressure-velocity manufactured solution",
        "13.8 thermal/CHT MMS and interface convergence",
        "13.9 turbulence transport MMS for k-epsilon/SST/SA",
        "13.10 radiation transport and angular convergence campaign",
        "13.11 serial/MPI equivalence for production physics",
        "13.12 real CUDA production CPU/GPU equivalence",
        "13.13 real CUDA beyond-VRAM OOC equivalence",
    ]

    report = {
        "phase": 13,
        "software_exit_code": rc,
        "campaigns": evidence,
        "pending_closure": pending,
        "closure": "VALIDATION_IN_PROGRESS",
    }

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    if args.junit:
        write_junit(args.junit, results)

    print("\nCLOSURE: VALIDATION_IN_PROGRESS")
    print("The driver never promotes skipped/unavailable hardware into PASS.")
    if args.strict:
        print("STRICT: FAIL — Phase 13 still has pending closure campaigns.")
        return 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
