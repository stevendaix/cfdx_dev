"""Report the longest executed tests from CTest JUnit without hiding failures."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def summarize(paths: list[Path]) -> dict:
    records = []
    for path in paths:
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"missing or empty CTest JUnit evidence: {path}")
        root = ET.parse(path).getroot()
        cases = root.findall(".//testcase")
        if not cases:
            raise ValueError(f"no testcase in CTest JUnit: {path}")
        for case in cases:
            status = "PASS"
            for result in ("failure", "error", "skipped"):
                if case.find(result) is not None:
                    status = result.upper()
                    break
            try:
                seconds = float(case.get("time", "0"))
            except ValueError as exc:
                raise ValueError(f"invalid elapsed time in {path}") from exc
            if not math.isfinite(seconds) or seconds < 0:
                raise ValueError(f"invalid elapsed time in {path}")
            if not case.get("name"):
                raise ValueError(f"unnamed test case in {path}")
            records.append({
                "name": case.get("name", ""),
                "seconds": seconds,
                "status": status,
                "file": str(path),
            })
    return {
        "executed": len(records),
        "failed_or_skipped": [r for r in records if r["status"] != "PASS"],
        "sum_test_seconds": round(sum(r["seconds"] for r in records), 3),
        "slowest": sorted(records, key=lambda r: r["seconds"], reverse=True)[:20],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("junit", nargs="+", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-names-file", type=Path)
    args = parser.parse_args()
    try:
        result = summarize(args.junit)
    except (ValueError, OSError, ET.ParseError) as exc:
        print(f"Cannot report CTest timings: {exc}", file=sys.stderr)
        return 2
    if args.require_names_file:
        expected = set(args.require_names_file.read_text(encoding="utf-8").splitlines())
        all_cases = [
            case
            for path in args.junit
            for case in ET.parse(path).getroot().findall(".//testcase")
        ]
        counts = {name: 0 for name in expected}
        incomplete = []
        for case in all_cases:
            name = case.get("name", "")
            if name in counts:
                counts[name] += 1
                if any(case.find(kind) is not None for kind in ("failure", "error", "skipped")):
                    incomplete.append(name)
        missing_or_duplicate = [name for name, n in counts.items() if n != 1]
        if not expected or missing_or_duplicate or incomplete:
            print(
                "Fail-closed CTest JUnit coverage violation: "
                f"missing/duplicate={missing_or_duplicate}, failing/skipped={incomplete}",
                file=sys.stderr,
            )
            return 2
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        f"CTest recorded {result['executed']} tests, "
        f"sum of individual test wall times={result['sum_test_seconds']:.1f}s"
    )
    for test in result["slowest"][:12]:
        print(f"  {test['seconds']:8.1f}s {test['status']:8s} {test['name']}")
    # GitHub automatically displays this in the Job summary. It is only a
    # diagnostic; the CTest exit status remains the source of pass/fail truth.
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with Path(summary).open("a", encoding="utf-8") as stream:
            stream.write("\n### CFDX CTest duration breakdown\n\n")
            stream.write("| Test | Wall time (s) | Status |\n|---|---:|---|\n")
            for test in result["slowest"][:12]:
                name = test["name"].replace("|", "\\|")
                stream.write(f"| `{name}` | {test['seconds']:.3f} | {test['status']} |\n")
            stream.write(
                "\n*Sum of test durations is not elapsed CI job time.*\n"
            )
    return 0  # the actual test run's exit code is the authoritative gate


if __name__ == "__main__":
    raise SystemExit(main())