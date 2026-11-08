"""Fail-closed aggregate evidence check for the split total-validation CTest suite.

This examines downloaded GitHub Actions artifacts, *not* a fresh local CTest
manifest: every shard must report the same deterministic plan and run every
assigned acceptance test exactly once with a passing (non-skipped) result.
CTest fixture setup tests outside the planned label may run on several workers.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path


def exactly_one(folder: Path, name: str) -> Path:
    matches = sorted(folder.rglob(name))
    if len(matches) != 1 or not matches[0].is_file():
        raise ValueError(
            f"expected exactly one {name} in {folder}; found {len(matches)}"
        )
    return matches[0]


def verify(root: Path, label: str, lanes: int) -> dict:
    if lanes < 1:
        raise ValueError("shard count must be positive")
    if not root.is_dir():
        raise ValueError(f"missing artifact directory: {root}")

    canonical = None
    executed = Counter()
    expected = Counter()
    extra = set()
    durations: dict[str, float] = {}
    for index in range(lanes):
        artifact = root / f"cfdx-total-validation-long-evidence-{index}"
        if not artifact.is_dir():
            raise ValueError(f"missing long-test artifact for shard {index}")
        plan_file = exactly_one(artifact, f"{label}-manifest.json")
        plan = json.loads(plan_file.read_text(encoding="utf-8"))
        if plan.get("suite_label") != label or plan.get("number_of_lanes") != lanes:
            raise ValueError(f"invalid plan metadata in shard {index}")
        if canonical is None:
            canonical = plan
        elif canonical != plan:
            raise ValueError("shards did not use the same test manifest")

        assigned = plan["lanes"][index]
        if assigned["index"] != index or assigned["test_count"] != len(
            assigned["tests"]
        ):
            raise ValueError(f"invalid shard {index} test assignment")
        test_names = exactly_one(artifact, f"{label}-shard-{index}.txt")
        if test_names.read_text(encoding="utf-8").splitlines() != assigned["tests"]:
            raise ValueError(f"shard {index} plan and selected-test file differ")

        xml = exactly_one(artifact, f"ctest-long-{index}.xml")
        cases = ET.parse(xml).getroot().findall(".//testcase")
        if not cases:
            raise ValueError(f"no CTest JUnit test cases for shard {index}")
        own = Counter()
        for case in cases:
            name = case.get("name")
            if not name:
                raise ValueError(f"empty JUnit test name in shard {index}")
            own[name] += 1
            for result in ("failure", "error", "skipped"):
                if case.find(result) is not None:
                    raise ValueError(
                        f"non-passing test: {name} on shard {index}: {result}"
                    )
            duration = float(case.get("time", "0"))
            if not math.isfinite(duration) or duration < 0:
                raise ValueError(f"invalid elapsed time for {name}: {duration}")
            durations[name] = duration
        for name in assigned["tests"]:
            if own[name] != 1:
                raise ValueError(
                    f"shard {index} expected one execution of {name}, got {own[name]}"
                )
            executed[name] += 1
        unexpected_long = set(own) & (
            set().union(*(set(l["tests"]) for l in plan["lanes"]))
            - set(assigned["tests"])
        )
        if unexpected_long:
            raise ValueError(
                f"long tests from wrong shard {index}: {sorted(unexpected_long)}"
            )
        extra.update(set(own) - set(assigned["tests"]))

    assert canonical is not None
    for shard in canonical["lanes"]:
        for name in shard["tests"]:
            expected[name] += 1
    if (
        not expected
        or len(expected) != canonical["total_selected"]
        or set(expected.values()) != {1}
    ):
        raise ValueError("duplicate or incomplete assignment in plan")
    if executed != expected:
        raise ValueError("aggregate long-test coverage incomplete or duplicated")
    return {
        "status": "PASS",
        "label": label,
        "shards": lanes,
        "long_tests_completed": len(expected),
        "fixture_or_auxiliary_tests": sorted(extra),
        "top_long_durations_seconds": [
            {"name": n, "seconds": durations[n]}
            for n in sorted(expected, key=lambda n: -durations[n])[:15]
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--lanes", required=True, type=int)
    parser.add_argument("--label", default="long")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        report = verify(args.root, args.label, args.lanes)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(
                json.dumps(report, indent=2) + "\n", encoding="utf-8"
            )
    except (
        ValueError,
        KeyError,
        IndexError,
        TypeError,
        OSError,
        ET.ParseError,
        json.JSONDecodeError,
    ) as exc:
        print(f"FAIL-CLOSED CTest shard verification: {exc}", file=sys.stderr)
        return 2
    print(
        f"Verified all {report['long_tests_completed']} long CTest checks "
        f"in {report['shards']} shards"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
