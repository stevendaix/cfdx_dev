"""Single fail-closed entrypoint for CFDX CTest qualification profiles.

Only orchestrates existing tests; never edits physics thresholds or promotes
physical qualification from the passing status of a contract-only test.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
CATALOG = HERE / "qualification_profiles.json"


class QualificationError(ValueError):
    pass


def read_catalog(path: Path = CATALOG) -> dict[str, dict[str, Any]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema") != 1 or not isinstance(data.get("profiles"), list):
        raise QualificationError("invalid qualification catalog schema")
    result: dict[str, dict[str, Any]] = {}
    labels: set[str] = set()
    for profile in data["profiles"]:
        name = profile.get("id")
        label = profile.get("label")
        if not isinstance(name, str) or not re.fullmatch(r"[a-z][a-z0-9-]*", name):
            raise QualificationError(f"invalid profile id: {name!r}")
        if name in result or label in labels or label != f"qualification:{name}":
            raise QualificationError(f"duplicate/invalid label for profile: {name}")
        for key in ("quick", "quick_targets", "required_extended"):
            items = profile.get(key)
            if not isinstance(items, list) or not all(
                isinstance(i, str) and i for i in items
            ):
                raise QualificationError(f"invalid {key} in profile {name}")
            if len(items) != len(set(items)):
                raise QualificationError(f"duplicate {key} in profile {name}")
        if not profile["quick"] or not profile["quick_targets"]:
            raise QualificationError(f"empty quick profile: {name}")
        if not isinstance(profile.get("description"), str):
            raise QualificationError(f"missing description in profile: {name}")
        pattern = profile.get("extended_regex")
        if not isinstance(pattern, str) or not pattern.startswith("^"):
            raise QualificationError(f"invalid extended_regex in profile {name}")
        re.compile(pattern)
        result[name] = profile
        labels.add(label)
    if not result:
        raise QualificationError("no profiles configured")
    return result


def props(test: dict[str, Any]) -> dict[str, Any]:
    return {p["name"]: p.get("value") for p in test.get("properties", [])}


def values(obj: Any) -> list[str]:
    if obj is None:
        return []
    if isinstance(obj, str):
        return [obj]
    if isinstance(obj, list) and all(isinstance(v, str) for v in obj):
        return obj
    raise QualificationError(f"invalid CTest property shape: {obj!r}")


def read_manifest(
    build_dir: Path, manifest_file: Path | None = None
) -> list[dict[str, Any]]:
    if manifest_file is not None:
        data = json.loads(manifest_file.read_text(encoding="utf-8"))
    else:
        output = subprocess.run(
            ["ctest", "--test-dir", str(build_dir), "--show-only=json-v1"],
            capture_output=True,
            text=True,
            check=True,
        )
        data = json.loads(output.stdout)
    tests = data.get("tests")
    if not isinstance(tests, list) or not tests:
        raise QualificationError(
            "CTest manifest is absent/empty; configure with CFDX_BUILD_TESTS=ON"
        )
    names = [t.get("name") for t in tests]
    if any(not isinstance(n, str) or not n for n in names) or len(set(names)) != len(
        names
    ):
        raise QualificationError("CTest manifest contains missing/duplicate test names")
    return tests


def plan_profile(
    tests: list[dict[str, Any]],
    name: str,
    scope: str,
    catalog: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    if name not in catalog:
        raise QualificationError(f"unknown qualification profile {name!r}")
    if scope not in ("quick", "extended"):
        raise QualificationError(f"unknown qualification scope {scope!r}")
    profile = catalog[name]
    by_name = {t["name"]: t for t in tests}
    anchor_names = list(profile["quick"])
    if scope == "extended":
        expression = re.compile(profile["extended_regex"])
        primary = set(anchor_names) | {n for n in by_name if expression.search(n)}
        for needed in profile["required_extended"]:
            if needed not in primary:
                raise QualificationError(
                    f"extended profile {name} does not select mandatory {needed}"
                )
        anchor_names.extend(profile["required_extended"])
    else:
        primary = set(anchor_names)

    missing = [n for n in anchor_names if n not in by_name]
    if missing:
        raise QualificationError(
            f"required qualification tests missing from CTest: {missing}"
        )
    if not primary:
        raise QualificationError(f"profile {name} selects no tests")

    setup: dict[str, set[str]] = {}
    for t in tests:
        for fixture in values(props(t).get("FIXTURES_SETUP")):
            setup.setdefault(fixture, set()).add(t["name"])

    selected = set(primary)
    pending = list(primary)
    while pending:
        current = pending.pop()
        if current not in by_name:
            raise QualificationError(
                f"qualification prerequisite not registered: {current}"
            )
        p = props(by_name[current])
        if p.get("DISABLED") in (True, "TRUE", "ON", "1", 1):
            raise QualificationError(
                f"disabled qualification must not be counted as passing: {current}"
            )
        if scope == "quick" and "long" in values(p.get("LABELS")):
            raise QualificationError(
                f"quick profile {name} includes a long test: {current}"
            )
        required = set(values(p.get("DEPENDS")))
        for fixture in values(p.get("FIXTURES_REQUIRED")):
            if fixture not in setup:
                raise QualificationError(
                    f"{current} requires missing fixture producer: {fixture}"
                )
            required.update(setup[fixture])
        for requirement in required:
            if requirement not in by_name:
                raise QualificationError(
                    f"{current} depends on missing test: {requirement}"
                )
            if requirement not in selected:
                selected.add(requirement)
                pending.append(requirement)

    names = sorted(selected)
    return {
        "schema": 1,
        "profile": name,
        "scope": scope,
        "description": profile["description"],
        "primary_tests": sorted(primary),
        "prerequisites": sorted(selected - primary),
        "selected_tests": names,
        "test_count": len(names),
        "includes_long": any(
            "long" in values(props(by_name[n]).get("LABELS")) for n in names
        ),
        "scope_warning": (
            "CTest profile only; passing does not certify reference physical accuracy. "
            "The total workflow additionally runs wall-distance N128, thermal/E2E, "
            "cross-shard coverage and physical-report gates."
        ),
    }


def save_plan(plan: dict[str, Any], directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    names = plan["selected_tests"]
    (directory / "plan.json").write_text(
        json.dumps(plan, indent=2) + "\n", encoding="utf-8"
    )
    (directory / "tests.txt").write_text("\n".join(names) + "\n", encoding="utf-8")
    (directory / "tests.regex").write_text(
        "^(" + "|".join(re.escape(n) for n in names) + ")$\n",
        encoding="utf-8",
    )


def result_summary(plan: dict[str, Any], mode: str) -> None:
    print(
        f"CFDX {mode}: {plan['profile']} ({plan['scope']}) - {plan['test_count']} CTest tests"
    )
    print(
        f"  primary: {len(plan['primary_tests'])}; prerequisites: {len(plan['prerequisites'])}"
    )
    print(f"  includes long tests: {plan['includes_long']}")
    print(f"  {plan['scope_warning']}")
    for n in plan["selected_tests"]:
        print(f"  - {n}")
    path = os.environ.get("GITHUB_STEP_SUMMARY")
    if path:
        with Path(path).open("a", encoding="utf-8") as stream:
            stream.write(
                f"\n### CFDX qualification: `{plan['profile']}` / `{plan['scope']}`\n\n"
                f"Selected CTest: **{plan['test_count']}** (prerequisites: "
                f"{len(plan['prerequisites'])}, long: "
                f"{'yes' if plan['includes_long'] else 'no'}).\n\n"
                f"> {plan['scope_warning']}\n\n"
            )


def resolve(
    event_name: str,
    event: dict[str, Any],
    manual_profile: str,
    manual_scope: str,
    catalog: dict[str, dict[str, Any]],
) -> tuple[str, str, str]:
    if event_name == "workflow_dispatch":
        if manual_profile not in catalog and manual_profile != "total":
            raise QualificationError(f"unknown manual profile: {manual_profile}")
        if manual_scope not in ("quick", "extended"):
            raise QualificationError(f"invalid manual scope: {manual_scope}")
        return manual_profile, manual_scope, "manual selection"
    if event_name == "pull_request":
        pr = event.get("pull_request") or {}
        label_items = pr.get("labels") or []
        labels = {item.get("name") for item in label_items if isinstance(item, dict)}
        active = sorted(labels & {p["label"] for p in catalog.values()})
        if "validation-total" in labels:
            return "none", "none", "full validation owns this PR (validation-total)"
        if len(active) > 1:
            raise QualificationError(
                f"multiple qualification labels; keep one: {active}"
            )
        if not active:
            return "none", "none", "no qualification label on PR"
        return active[0].removeprefix("qualification:"), "extended", "opt-in PR label"
    raise QualificationError(f"unsupported event: {event_name}")


def write_output(name: str, scope: str, reason: str, output: Path | None) -> None:
    print(f"CFDX profile resolution: {name} / {scope} ({reason})")
    if output:
        with output.open("a", encoding="utf-8") as f:
            f.write(f"profile={name}\nscope={scope}\n")
    path = os.environ.get("GITHUB_STEP_SUMMARY")
    if path:
        with Path(path).open("a", encoding="utf-8") as f:
            f.write(
                f"\nCFDX qualification activation: **{name}** ({scope}) - {reason}.\n"
            )


def cli() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="list available qualification profiles and PR labels")
    labels = sub.add_parser("labels", help="emit labels for one-time GitHub setup")
    labels.add_argument(
        "--tsv", action="store_true", help="name, description TSV for gh label create"
    )
    res = sub.add_parser(
        "resolve", help="resolve manual dispatch or PR qualification label"
    )
    res.add_argument("--event-name", required=True)
    res.add_argument("--event-file", type=Path)
    res.add_argument("--manual-profile", default="naca")
    res.add_argument("--manual-scope", default="quick")
    res.add_argument("--github-output", type=Path)
    for kind in ("plan", "build", "run"):
        p = sub.add_parser(kind)
        p.add_argument("--profile", required=True)
        p.add_argument("--scope", choices=("quick", "extended"), default="quick")
        p.add_argument("--build-dir", type=Path, required=True)
        p.add_argument("--output-dir", type=Path)
        p.add_argument("--parallel", type=int, default=3)
        p.add_argument("--manifest", type=Path, help="offline planning/tests only")
    args = parser.parse_args()
    try:
        catalog = read_catalog(args.catalog)
        if args.command == "list":
            print("CFDX qualification profiles (quick / extended CTest):")
            for name, profile in catalog.items():
                print(f"  {name:14s} {profile['label']:30s} {profile['description']}")
            print(
                "  total          validation-total               entire canonical CI workflow"
            )
            print("Note: extended CTest is not the total physical qualification.")
            return 0
        if args.command == "labels":
            for profile in catalog.values():
                if args.tsv:
                    print(profile["label"] + "\t" + profile["description"])
                else:
                    print(profile["label"])
            return 0
        if args.command == "resolve":
            event = json.loads(args.event_file.read_text()) if args.event_file else {}
            name, scope, reason = resolve(
                args.event_name,
                event,
                args.manual_profile,
                args.manual_scope,
                catalog,
            )
            write_output(name, scope, reason, args.github_output)
            return 0
        if args.profile == "total":
            raise QualificationError(
                "total must run via cfdx-validation.yml, not a partial CTest selection"
            )
        if args.parallel < 1:
            raise QualificationError("parallel must be positive")
        build_dir = args.build_dir.resolve()
        out = (
            args.output_dir
            or build_dir / "qualification" / f"{args.profile}-{args.scope}"
        ).resolve()
        if args.command == "build":
            # Quick mode intentionally builds only listed executables; extended
            # mode builds everything because some Python tests launch other tools.
            command = [
                "cmake",
                "--build",
                str(build_dir),
                "--parallel",
                str(args.parallel),
            ]
            if args.scope == "quick":
                command += ["--target", *catalog[args.profile]["quick_targets"]]
            print("Build: " + " ".join(command), flush=True)
            return subprocess.call(command)
        tests = read_manifest(build_dir, args.manifest)
        plan = plan_profile(tests, args.profile, args.scope, catalog)
        save_plan(plan, out)
        result_summary(plan, args.command)
        if args.command == "plan":
            return 0
        if args.manifest:
            raise QualificationError(
                "--manifest is for offline planning only, not execution"
            )
        junit = out / "ctest.xml"
        command = [
            "ctest",
            "--test-dir",
            str(build_dir),
            "-R",
            (out / "tests.regex").read_text(encoding="utf-8").strip(),
            "--no-tests=error",
            "--parallel",
            str(args.parallel),
            "--output-on-failure",
            "--output-junit",
            str(junit),
        ]
        print("Execute selected CTest tests with JUnit evidence", flush=True)
        status = subprocess.call(command)
        # Always validate evidence even when CTest failed or crashed; the
        # original exit status remains authoritative and is never hidden.
        coverage = subprocess.call(
            [
                sys.executable,
                str(HERE / "ctest_timings.py"),
                str(junit),
                "--require-names-file",
                str(out / "tests.txt"),
                "--output",
                str(out / "timings.json"),
            ]
        )
        if status != 0 or coverage != 0:
            print(
                f"QUALIFICATION FAIL: ctest={status} coverage={coverage}",
                file=sys.stderr,
            )
            return 1
        print(f"QUALIFICATION PASS (CTest scope only): {args.profile} / {args.scope}")
        return 0
    except (
        QualificationError,
        OSError,
        subprocess.CalledProcessError,
        json.JSONDecodeError,
        re.error,
    ) as exc:
        print(f"CFDX qualification FAIL-CLOSED: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(cli())