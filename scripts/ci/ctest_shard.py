"""Deterministic, fail-closed CTest lane planner (stdlib only).

The CTest JSON manifest is the source of truth; no test names are hand-copied
into lane configuration. Each selected test is assigned to exactly one lane.
CTest DEPENDS links between selected tests remain in the same lane. Required
fixture setup/cleanup tests are added automatically by CTest when selected
with -R, exactly as in an unsharded invocation.

Weights are *planning hints*, not elapsed-time measurements: prefer the CTest
cost file if supplied, otherwise fall back to configured TIMEOUT upper bounds.
This script never reduces solver budgets, physical acceptance gates or test
TIMEOUT values.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def properties(test: dict) -> dict:
    return {item["name"]: item["value"] for item in test.get("properties", [])}


def array(value: object) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(x) for x in value]
    return [str(x) for x in str(value).split(";") if x]


def load_manifest(build_dir: Path) -> list[dict]:
    if not (build_dir / "CTestTestfile.cmake").is_file():
        raise ValueError(f"not a CTest build directory: {build_dir}")
    result = subprocess.run(
        ["ctest", "--test-dir", str(build_dir), "--show-only=json-v1"],
        text=True,
        capture_output=True,
        check=True,
    )
    tests = json.loads(result.stdout)["tests"]
    names = [t["name"] for t in tests]
    if len(names) != len(set(names)):
        raise ValueError("CTest has duplicate test names")
    return tests


def costs_from_ctest(path: Path | None) -> dict[str, float]:
    """Read CTestCostData.txt (name, previous run count, mean elapsed seconds)."""
    if path is None or not path.exists():
        return {}
    costs: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("---"):
            break
        parts = line.split()
        if len(parts) != 3:
            continue
        name, count, elapsed = parts
        try:
            if int(count) > 0 and math.isfinite(float(elapsed)) and float(elapsed) > 0:
                costs[name] = float(elapsed)
        except ValueError:
            continue
    return costs


def build_plan(
    tests: list[dict],
    label: str,
    lanes: int,
    required_tests: tuple[str, ...] = (),
    measured_costs: dict[str, float] | None = None,
) -> dict:
    if lanes < 1:
        raise ValueError("lanes must be positive")
    if not re.fullmatch(r"[A-Za-z0-9_-]+", label):
        raise ValueError("invalid CTest label")
    measured_costs = measured_costs or {}
    by_name = {t["name"]: t for t in tests}
    selected = {
        name: test
        for name, test in by_name.items()
        if label in array(properties(test).get("LABELS"))
    }
    if len(selected) < lanes:
        raise ValueError(
            f"label {label!r} selects only {len(selected)} tests for {lanes} lanes"
        )
    for name in required_tests:
        if name not in selected:
            raise ValueError(f"required {label} test missing: {name}")
    for name, test in selected.items():
        props = properties(test)
        if str(props.get("DISABLED", "False")).lower() in ("true", "1", "on"):
            raise ValueError(f"disabled {label} test would be silently skipped: {name}")
        # CTest JSON omits command when a target has not been built yet.
        # Execution checks belong to CTest itself, not the planner.

    parent = {name: name for name in selected}

    def root(name: str) -> str:
        while parent[name] != name:
            parent[name] = parent[parent[name]]
            name = parent[name]
        return name

    def union(a: str, b: str) -> None:
        pa, pb = root(a), root(b)
        if pa != pb:
            parent[max(pa, pb)] = min(pa, pb)

    for name, test in selected.items():
        for dependency in array(properties(test).get("DEPENDS")):
            if dependency not in by_name:
                raise ValueError(f"{name} depends on missing test {dependency}")
            if dependency in selected:
                union(name, dependency)
            else:
                # CTest auto-includes a fixture producer when a requiring test
                # is selected. An arbitrary cross-suite DEPENDS is not run by
                # CTest -R; fail rather than silently lose that prerequisite.
                p = properties(test)
                producers = properties(by_name[dependency])
                required = set(array(p.get("FIXTURES_REQUIRED")))
                setup = set(array(producers.get("FIXTURES_SETUP")))
                if not required.intersection(setup):
                    raise ValueError(
                        f"{name} has dependency outside {label} suite: {dependency}"
                    )

    groups: dict[str, list[str]] = defaultdict(list)
    for name in selected:
        groups[root(name)].append(name)

    def weight(name: str) -> float:
        measured = measured_costs.get(name)
        if measured is not None and math.isfinite(measured) and measured > 0:
            return measured
        timeout = properties(selected[name]).get("TIMEOUT")
        if timeout is not None:
            try:
                v = float(timeout)
                if math.isfinite(v) and v > 0:
                    return v
            except (ValueError, TypeError):
                pass
        return 180.0  # heuristic only, not a measured duration

    grouped = [
        (sum(weight(n) for n in names), sorted(names)) for names in groups.values()
    ]
    grouped.sort(key=lambda pair: (-pair[0], pair[1]))
    assigned: list[list[str]] = [[] for _ in range(lanes)]
    totals = [0.0] * lanes
    for cost, names in grouped:
        lane = min(range(lanes), key=lambda i: (totals[i], len(assigned[i]), i))
        assigned[lane].extend(names)
        totals[lane] += cost

    flat = [name for lane in assigned for name in lane]
    if set(flat) != set(selected) or len(flat) != len(selected):
        raise ValueError("FAIL-CLOSED: lane coverage missing or duplicated tests")
    if any(not names for names in assigned):
        raise ValueError("FAIL-CLOSED: empty lane would produce a false-green job")
    return {
        "suite_label": label,
        "total_selected": len(selected),
        "number_of_lanes": lanes,
        "cost_policy": "observed CTestCostData when present, otherwise TIMEOUT upper bound / 180s heuristic; not a runtime estimate",
        "lanes": [
            {
                "index": i,
                "test_count": len(names),
                "planning_weight": round(totals[i], 3),
                "tests": sorted(names),
            }
            for i, names in enumerate(assigned)
        ],
    }


def write_plan(plan: dict, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    for lane in plan["lanes"]:
        prefix = output / f"{plan['suite_label']}-shard-{lane['index']}"
        names = lane["tests"]
        prefix.with_suffix(".txt").write_text("\n".join(names) + "\n", encoding="utf-8")
        # Compatible with CTest 3.28 on ubuntu-24.04; --tests-from-file needs
        # CTest 3.29+, so we use an anchored, escaped -R regex instead.
        prefix.with_suffix(".regex").write_text(
            "^(" + "|".join(re.escape(name) for name in names) + ")$\n",
            encoding="utf-8",
        )
    (output / f"{plan['suite_label']}-manifest.json").write_text(
        json.dumps(plan, indent=2, ensure_ascii=True) + "\n", encoding="utf-8"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--label", default="long")
    parser.add_argument("--lanes", type=int, default=4)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--required-test", action="append", default=[])
    parser.add_argument("--cost-data", type=Path, default=None)
    args = parser.parse_args()
    try:
        tests = load_manifest(args.build_dir.resolve())
        costs = costs_from_ctest(args.cost_data)
        plan = build_plan(
            tests,
            args.label,
            args.lanes,
            required_tests=tuple(args.required_test),
            measured_costs=costs,
        )
        write_plan(plan, args.output_dir)
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        print(f"CTest shard planning FAILED: {error}", file=sys.stderr)
        return 2
    print(
        f"Fail-closed {args.label} coverage: "
        f"{plan['total_selected']} tests in {args.lanes} lanes"
    )
    for lane in plan["lanes"]:
        print(
            f"  shard {lane['index']}: {lane['test_count']} tests; "
            f"planning weight={lane['planning_weight']:.1f} (NOT measured seconds)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())