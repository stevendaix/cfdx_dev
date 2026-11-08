"""Basic workflow topology and YAML safety checks (not a GitHub Actions emulator).

Run before a PR to catch accidental duplicate keys, duplicate N9 triggers,
missing total-validation partition gates, and missing long-test scheduling.
Requires PyYAML; authoritative GitHub execution still happens on runners.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml


class UniqueKeysLoader(yaml.BaseLoader):
    """YAML 1.2-like scalar behavior for the GitHub `on` trigger key."""


def unique_mapping(loader: yaml.BaseLoader, node: yaml.MappingNode) -> dict:
    result = {}
    for k_node, v_node in node.value:
        key = loader.construct_object(k_node, deep=True)
        if key in result:
            raise ValueError(f"duplicate YAML key: {key!r} at line {k_node.start_mark.line+1}")
        result[key] = loader.construct_object(v_node, deep=True)
    return result


UniqueKeysLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, unique_mapping)


def load(path: Path) -> dict:
    try:
        result = yaml.load(path.read_text(encoding="utf-8"), Loader=UniqueKeysLoader)
    except yaml.YAMLError as exc:
        raise ValueError(f"invalid YAML {path}: {exc}") from exc
    if not isinstance(result, dict) or not isinstance(result.get("jobs"), dict):
        raise ValueError(f"invalid GitHub workflow structure: {path}")
    if not isinstance(result.get("on"), dict):
        raise ValueError(f"missing workflow event map: {path}")
    return result


def audit(workflows: Path) -> dict:
    definitions = {}
    for path in sorted(workflows.glob("*.yml")):
        definitions[path.name] = load(path)
    required = {
        "cfdx-ci.yml", "cfdx-numerical.yml", "cfdx-validation.yml",
        "n9-campaign.yml", "n9-physical-validation.yml", "precommit.yml",
    }
    if not required.issubset(definitions):
        raise ValueError(f"missing mandatory workflow definitions: {sorted(required - set(definitions))}")
    full = definitions["cfdx-validation.yml"]
    jobs = full["jobs"]
    partitions = {
        "prepare", "fast-tests", "long-tests", "wall-distance", "wall-distance-fine",
        "thermal-radiation", "application-e2e", "validation-report",
    }
    if not partitions.issubset(jobs):
        raise ValueError(f"missing total-validation jobs: {sorted(partitions - set(jobs))}")
    final = jobs.get("final-gate")
    if not final or set(final.get("needs", [])) != partitions:
        raise ValueError("final-gate must depend on every total-validation partition")
    if "schedule" not in full["on"] or "workflow_dispatch" not in full["on"]:
        raise ValueError("missing weekly/manual total-validation triggers")
    long = jobs["long-tests"]
    matrix = long.get("strategy", {}).get("matrix", {})
    if matrix.get("shard") != ["0", "1", "2", "3"] or long.get("strategy", {}).get("fail-fast") != "false":
        raise ValueError("long validation must have 4 independent, non-fail-fast shards")
    if "long-tests" not in jobs["validation-report"].get("needs", []):
        raise ValueError("report must wait for all long-test evidence")
    for name in ("n9-campaign.yml", "n9-physical-validation.yml"):
        triggers = set(definitions[name]["on"])
        if triggers != {"workflow_dispatch"}:
            raise ValueError(f"{name} repeats a N9 PR/scheduled campaign: {sorted(triggers)}")
    if "pull_request" not in definitions["cfdx-ci.yml"]["on"]:
        raise ValueError("mandatory fast PR workflow is not triggered on PRs")
    return {"workflows_checked": len(definitions), "mandatory_partitions": len(partitions)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workflows", type=Path, default=Path(".github/workflows"))
    args = parser.parse_args()
    try:
        summary = audit(args.workflows)
    except (ValueError, OSError) as exc:
        print(f"CFDX CI topology FAILED: {exc}", file=sys.stderr)
        return 2
    print(f"CI topology PASS: {summary['workflows_checked']} workflows; {summary['mandatory_partitions']} gated full-validation partitions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())