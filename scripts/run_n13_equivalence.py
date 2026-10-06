#!/usr/bin/env python3
"""Run the N13 cross-backend equivalence campaign and retain auditable evidence."""
from __future__ import annotations
import argparse, json, os, re, subprocess
from pathlib import Path

CAMPAIGNS = {
    "serial_mpi": ["test_n13_mpi_equivalence"],
    "assembled_matrix_free": ["test_matrix_free_fv_operator", "test_backend_equivalence"],
    "deterministic_reduction": ["test_n13_mpi_equivalence", "test_mpi_deterministic"],
    "restart": ["test_n13_restart_equivalence", "test_restart_mapping"],
}
N_TO_M_TESTS = ["test_production_mpi_poisson_writer", "test_production_mpi_poisson_reader"]

def run(cmd: list[str], capture: bool = True, timeout: int = 30):
    try:
        return subprocess.run(cmd, check=False, capture_output=capture, text=True, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired):
        return None

def command_text(cmd: list[str]) -> str | None:
    p = run(cmd)
    if p is None:
        return None
    out = (p.stdout + p.stderr).strip()
    return out or None

def cache_value(build: Path, key: str) -> str | None:
    path = build / "CMakeCache.txt"
    if not path.exists():
        return None
    prefix = key + ":"
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(prefix) and "=" in line:
            return line.split("=", 1)[1].strip()
    return None

def provenance(build: Path) -> dict:
    nvcc = command_text(["nvcc", "--version"])
    nvidia = command_text(["nvidia-smi", "--query-gpu=name,driver_version", "--format=csv,noheader"])
    mpi = command_text(["mpirun", "--version"])
    return {
        "commit_sha": command_text(["git", "rev-parse", "HEAD"]),
        "build_type": cache_value(build, "CMAKE_BUILD_TYPE"),
        "compiler": {"path": cache_value(build, "CMAKE_CXX_COMPILER"),
                     "id": cache_value(build, "CMAKE_CXX_COMPILER_ID"),
                     "version": cache_value(build, "CMAKE_CXX_COMPILER_VERSION")},
        "mpi": mpi.splitlines()[0] if mpi else None,
        "cuda": {"device": nvidia, "toolkit": nvcc.splitlines()[-1] if nvcc else None},
        "runner": {k: os.environ.get(k) for k in
                   ("GITHUB_RUN_ID", "GITHUB_WORKFLOW", "RUNNER_NAME", "RUNNER_OS", "RUNNER_ARCH")},
    }

def exists(build: Path, name: str) -> bool:
    p = run(["ctest", "--test-dir", str(build), "-N", "-R", f"^{re.escape(name)}$"])
    return p is not None and p.returncode == 0 and re.search(r"Total Tests:\s*1", p.stdout) is not None

def metrics(text: str) -> dict[str, float]:
    patterns = {"l2_relative": r"L2rel=([0-9.eE+-]+)",
                "linf_relative": r"Linfrel=([0-9.eE+-]+)",
                "max_error": r"max(?:_error| error)[ =]([0-9.eE+-]+)",
                "reduction_delta": r"reduction equivalence failed:\s*([0-9.eE+-]+)"}
    out = {}
    for key, pattern in patterns.items():
        m = re.search(pattern, text, re.I)
        if m:
            out[key] = float(m.group(1))
    return out

def ctest(build: Path, name: str) -> dict:
    if not exists(build, name):
        return {"test": name, "status": "BLOCKED", "exit_code": 2,
                "reason": "test is not registered in this build"}
    p = run(["ctest", "--test-dir", str(build), "-R", f"^{re.escape(name)}$", "--output-on-failure"])
    if p is None:
        return {"test": name, "status": "FAIL", "exit_code": 1, "reason": "ctest could not execute"}
    output = (p.stdout + p.stderr).strip()
    return {"test": name, "status": "PASS" if p.returncode == 0 else "FAIL",
            "exit_code": p.returncode, "metrics": metrics(output), "diagnostics": output[-4000:]}

def cuda_test(build: Path) -> dict:
    name = "test_n13_cuda_equivalence"
    if not exists(build, name) or not (build / name).exists():
        return {"test": name, "status": "BLOCKED", "exit_code": 2,
                "reason": "CUDA executable is not registered/built"}
    p = run([str(build / name)])
    if p is None:
        return {"test": name, "status": "FAIL", "exit_code": 1, "reason": "CUDA executable could not execute"}
    output = (p.stdout + p.stderr).strip()
    return {"test": name, "status": "PASS" if p.returncode == 0 else ("BLOCKED" if p.returncode == 2 else "FAIL"),
            "exit_code": p.returncode, "metrics": metrics(output), "diagnostics": output[-4000:]}

def campaign(build: Path, names: list[str]) -> dict:
    tests = [ctest(build, name) for name in names]
    status = "FAIL" if any(x["status"] == "FAIL" for x in tests) else              "BLOCKED" if any(x["status"] == "BLOCKED" for x in tests) else "PASS"
    return {"status": status, "tests": tests,
            "conservation": {"status": "NOT_APPLICABLE", "reason": "backend/operator equivalence"},
            "true_residual": {"status": "NOT_APPLICABLE", "reason": "backend/operator equivalence"},
            "restart_metadata": {"status": "NOT_APPLICABLE"}}

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", type=Path, default=Path("build"))
    ap.add_argument("--json", type=Path, default=Path("artifacts/n13_equivalence.json"))
    ap.add_argument("--require-cuda", action="store_true")
    ap.add_argument("--require-parallel-hdf5", action="store_true")
    args = ap.parse_args()
    build = args.build_dir
    if not (build / "CTestTestfile.cmake").exists():
        raise SystemExit(f"build directory is not configured: {build}")

    report = {
        "schema_version": "1.1", "package": "N13",
        "name": "Cross-backend numerical equivalence",
        "tolerances": {"operator_l2_relative": 1e-12, "operator_linf_relative": 1e-12,
                       "restart": 0.0, "deterministic_reduction": 1e-14},
        "provenance": provenance(build),
        "problem_definition": {
            "mesh": {"serial_mpi": "canonical two-cell Poisson partition",
                     "assembled_matrix_free": "registered FV operator campaign",
                     "cuda": "two-cell constant-gradient kernel",
                     "restart": "global-ID remapping and temporal history",
                     "n_to_m_mpi_restart": "production Poisson checkpoint, 2 ranks -> 3 ranks"},
            "numerics": "N13 equivalence contracts; no convergence-tolerance substitution"},
        "campaigns": {},
    }

    failures = 0
    for name, tests in CAMPAIGNS.items():
        report["campaigns"][name] = campaign(build, tests)
        if report["campaigns"][name]["status"] == "FAIL":
            failures += 1

    cuda = cuda_test(build)
    report["campaigns"]["cuda"] = {
        "status": cuda["status"], "tests": [cuda],
        "conservation": {"status": "NOT_APPLICABLE", "reason": "kernel equivalence"},
        "true_residual": {"status": "NOT_APPLICABLE", "reason": "kernel equivalence"},
        "restart_metadata": {"status": "NOT_APPLICABLE"},
    }
    if cuda["status"] == "FAIL" or (cuda["status"] == "BLOCKED" and args.require_cuda):
        failures += 1

    n2m = campaign(build, N_TO_M_TESTS)
    n2m["mpi_transition"] = "2->3"
    n2m["restart_metadata"] = {"status": n2m["status"], "mpi_transition": "2->3",
                               "checkpoint": str(build / "production_poisson_checkpoint.h5"),
                               "reference": str(build / "production_poisson_reference.h5")}
    report["campaigns"]["n_to_m_mpi_restart"] = n2m
    if n2m["status"] == "FAIL" or (n2m["status"] == "BLOCKED" and args.require_parallel_hdf5):
        failures += 1

    blocked = [k for k, v in report["campaigns"].items() if v["status"] == "BLOCKED"]
    report["blocked_axes"] = blocked
    report["closure"] = "FAIL" if failures else ("PASS_WITH_BLOCKED_AXES" if blocked else "PASS")
    report["note"] = ("BLOCKED is retained as evidence of unavailable execution capability; "
                      "it is never converted to PASS. Hardware qualification should use "
                      "--require-cuda and --require-parallel-hdf5.")

    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return 1 if failures else 0

if __name__ == "__main__":
    raise SystemExit(main())
