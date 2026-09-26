#!/usr/bin/env python3
"""Validate the CFDX engineering qualification registry.

This is a metadata/contract check, not a numerical validation. Numerical PASS
requires the case executable and its retained calculation evidence.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REQUIRED = {
    "reference",
    "parameters",
    "mesh",
    "boundary_conditions",
    "numerical_schemes",
    "solver_settings",
    "convergence_criteria",
    "conservation_criteria",
    "qoi",
    "reference_value",
    "error",
    "observed_order",
    "runtime",
    "regression_status",
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "registry",
        nargs="?",
        default="docs/validation/CFDX_QUALIFICATION_REGISTRY.json",
    )
    args = parser.parse_args()

    path = Path(args.registry)
    data = json.loads(path.read_text(encoding="utf-8"))

    if data.get("schema_version") != "1.0":
        raise SystemExit("invalid qualification registry schema_version")

    cases = data.get("cases")
    if not isinstance(cases, list) or not cases:
        raise SystemExit("qualification registry contains no cases")

    ids: set[str] = set()
    failures: list[str] = []

    for case in cases:
        ident = case.get("id", "<missing-id>")
        if ident in ids:
            failures.append(f"{ident}: duplicate case id")
        ids.add(ident)

        missing = sorted(REQUIRED - set(case))
        if missing:
            failures.append(f"{ident}: missing fields: {', '.join(missing)}")

        for field in REQUIRED:
            value = case.get(field)
            if value is None or (isinstance(value, str) and not value.strip()):
                failures.append(f"{ident}: empty required field '{field}'")

        status = case.get("status")
        if status not in {"PLANNED", "READY", "RUNNING", "DIAGNOSTIC", "PASS", "BLOCKED"}:
            failures.append(f"{ident}: invalid status '{status}'")

        if status == "PASS":
            text = " ".join(str(case.get(k, "")) for k in REQUIRED)
            if "PASS" not in text and "regression" not in text.lower():
                failures.append(f"{ident}: PASS case lacks regression evidence declaration")

    if failures:
        print("CFDX_QUALIFICATION_REGISTRY: FAIL")
        for failure in failures:
            print(f"  - {failure}")
        return 1

    print(
        "CFDX_QUALIFICATION_REGISTRY: PASS "
        f"cases={len(cases)} "
        f"required_fields={len(REQUIRED)}"
    )
    for case in cases:
        print(
            f"QUALIFICATION_CASE id={case['id']} "
            f"domain={case['domain']} status={case['status']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
