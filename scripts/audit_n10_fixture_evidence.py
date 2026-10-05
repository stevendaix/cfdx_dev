#!/usr/bin/env python3
"""Audit the verified N10 public-fixture evidence chain.

The gate combines independent byte verification with the production importer
qualification. It does not claim numerical solver support: that remains a
separate N10 milestone requiring solver iterations/residual evidence.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: Path) -> dict[str, object]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"error: cannot read valid JSON report {path}: {exc}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verification", type=Path, required=True)
    parser.add_argument("--qualification", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    verification = load(args.verification)
    qualification = load(args.qualification)

    verification_rows = {
        str(row["id"]): row
        for row in verification.get("fixtures", [])
        if isinstance(row, dict) and "id" in row
    }
    qualification_rows = {
        str(row["id"]): row
        for row in qualification.get("fixtures", [])
        if isinstance(row, dict) and "id" in row
    }

    ids = sorted(set(verification_rows) | set(qualification_rows))
    mismatches: list[str] = []
    rows: list[dict[str, object]] = []

    for fixture_id in ids:
        verified = verification_rows.get(fixture_id)
        qualified = qualification_rows.get(fixture_id)
        if verified is None:
            mismatches.append(f"{fixture_id}: missing from independent verification")
            continue
        if qualified is None:
            mismatches.append(f"{fixture_id}: missing from production qualification")
            continue

        byte_status = verified.get("status")
        imported = qualified.get("imported")
        topology = qualified.get("topology_valid")
        geometry = qualified.get("geometry_quality_valid")
        counts = qualified.get("counts", {})

        if byte_status != "VERIFIED":
            mismatches.append(f"{fixture_id}: byte verification status={byte_status}")
        if imported is not True:
            mismatches.append(f"{fixture_id}: imported={imported}")
        if topology is not True:
            mismatches.append(f"{fixture_id}: topology_valid={topology}")
        if geometry is not True:
            mismatches.append(f"{fixture_id}: geometry_quality_valid={geometry}")
        if not isinstance(counts, dict) or any(
            not isinstance(counts.get(key), int) or counts.get(key, 0) <= 0
            for key in ("points", "faces", "cells")
        ):
            mismatches.append(f"{fixture_id}: invalid positive mesh counts")

        rows.append(
            {
                "id": fixture_id,
                "byte_verified": byte_status == "VERIFIED",
                "imported": imported is True,
                "topology_valid": topology is True,
                "geometry_quality_valid": geometry is True,
                "counts": counts,
            }
        )

    report = {
        "campaign": "N10 verified public-fixture evidence gate",
        "status": "PASS" if ids and not mismatches else "FAIL",
        "fixture_count": len(ids),
        "fixtures": rows,
        "mismatches": mismatches,
        "scope": {
            "byte_integrity": "independent SHA-256 verification",
            "production_import": "import_mesh dispatcher",
            "topology": "Mesh::topo_validate",
            "geometry": "existing 3D validator or explicit 2D edge-mesh validation",
            "numerical_solver_support": "NOT_CLAIMED",
        },
        "policy": {
            "changes_numerical_tolerances": False,
            "disables_validation": False,
            "silent_fallbacks": False,
            "solver_qualification_deferred": True,
        },
    }

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
