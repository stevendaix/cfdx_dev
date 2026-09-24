#!/usr/bin/env python3
"""Run and summarize the complete CFDX thermal/radiation regression campaign."""
from __future__ import annotations
import argparse, json, re, subprocess
from pathlib import Path

SUITES = [
    ("thermophysical_models", "test_thermophysical_models_vv"),
    ("thermal_solver", "test_thermal_vv"),
    ("radiation_solver", "test_radiation_vv"),
    ("s2s_radiation", "test_s2s_radiation_vv"),
    ("model_matrix", "test_thermal_radiation_model_matrix"),
    ("regression_core", "test_thermal_radiation_regression"),
]
RESULT_RE = re.compile(r"Results:\\s+(\\d+) passed, (\\d+) failed")
PASS_RE = re.compile(r"^\\s*PASS:\\s+(.+)$", re.MULTILINE)
FAIL_RE = re.compile(r"^\\s*FAIL:\\s+(.+)$", re.MULTILINE)
MARKER_RE = re.compile(r"^(?:THERMAL|RADIATION)_RESIDUAL:.*$", re.MULTILINE)
REG_RE = re.compile(r"^THERMAL_RADIATION_REGRESSION:.*$", re.MULTILINE)

REQUIRED = {
    "thermophysical_models": ("THERMAL_RESIDUAL:",),
    "thermal_solver": ("THERMAL_VV:", "THERMAL_RESIDUAL:"),
    "radiation_solver": ("RADIATION_RESIDUAL:",),
    "s2s_radiation": ("RADIATION_RESIDUAL:",),
    "model_matrix": ("THERMAL_RESIDUAL:", "RADIATION_RESIDUAL:"),
    "regression_core": ("THERMAL_RADIATION_REGRESSION:",),
}

def run_suite(build_dir, name, executable, log_dir):
    exe = build_dir / executable
    log = log_dir / (name + ".log")
    if not exe.exists():
        log.write_text("ERROR: missing executable " + str(exe) + "\n", encoding="utf-8")
        return {"suite": name, "returncode": 127, "passed": 0, "failed": 1,
                "markers": [], "log": str(log), "error": "missing executable"}
    p = subprocess.run([str(exe)], cwd=build_dir.parent, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    log.write_text(p.stdout, encoding="utf-8")
    m = RESULT_RE.search(p.stdout)
    passed = int(m.group(1)) if m else 0
    failed = int(m.group(2)) if m else (1 if p.returncode else 0)
    return {"suite": name, "returncode": p.returncode, "passed": passed,
            "failed": failed, "cases": PASS_RE.findall(p.stdout) +
            ["FAIL: " + x for x in FAIL_RE.findall(p.stdout)],
            "markers": MARKER_RE.findall(p.stdout) + REG_RE.findall(p.stdout),
            "log": str(log)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)
    results = [run_suite(a.build_dir, n, e, a.output_dir) for n, e in SUITES]
    for r in results:
        required = REQUIRED[r["suite"]]
        r["missing_markers"] = [x for x in required if not any(m.startswith(x) for m in r["markers"])]
        r["status"] = "PASS" if r["returncode"] == 0 and r["failed"] == 0 and not r["missing_markers"] else "FAIL"
    status = "PASS" if all(r["status"] == "PASS" for r in results) else "FAIL"
    summary = {"status": status, "suite_count": len(results),
               "total_passed": sum(r["passed"] for r in results),
               "total_failed": sum(r["failed"] for r in results), "suites": results,
               "policy": {"analytical_oracles_are_executable_gates": True,
                          "solver_returncode_must_be_zero": True,
                          "required_diagnostics_must_be_present": True,
                          "historical_baseline_not_assumed": True}}
    (a.output_dir / "thermal-radiation-regression.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    lines = ["# CFDX Thermal/Radiation Regression Report", "",
             "**Status:** " + status,
             "**Suites:** " + str(len(results)),
             "**Cases passed:** " + str(summary["total_passed"]),
             "**Cases failed:** " + str(summary["total_failed"]), "",
             "| Suite | Status | Passed | Failed | Missing diagnostics |",
             "|---|---|---:|---:|---|"]
    for r in results:
        lines.append("| " + r["suite"] + " | " + r["status"] + " | " +
                     str(r["passed"]) + " | " + str(r["failed"]) + " | " +
                     (", ".join(r["missing_markers"]) or "none") + " |")
    lines += ["", "## Residual / balance diagnostics", ""]
    for r in results:
        lines.append("### " + r["suite"])
        if r["markers"]:
            lines.extend("- " + m for m in r["markers"])
        else:
            lines.append("- none")
    lines += ["", "## Gate", "",
              "FAIL = non-zero executable, failed test case, or missing required diagnostics.",
              "Analytical expectations remain inside the executable tests; no historical numerical baseline is assumed."]
    (a.output_dir / "thermal-radiation-regression.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("THERMAL_RADIATION_REGRESSION: status=" + status +
          " suites=" + str(len(results)) +
          " passed=" + str(summary["total_passed"]) +
          " failed=" + str(summary["total_failed"]))
    for r in results:
        print("THERMAL_RADIATION_REGRESSION: suite=" + r["suite"] +
              " status=" + r["status"] + " passed=" + str(r["passed"]) +
              " failed=" + str(r["failed"]))
    return 0 if status == "PASS" else 1

if __name__ == "__main__":
    raise SystemExit(main())
