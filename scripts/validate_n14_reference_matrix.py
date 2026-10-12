#!/usr/bin/env python3
"""Validate the machine-readable N14 reference-code comparison matrix."""
from __future__ import annotations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
MATRIX=ROOT/"docs/validation/N14_REFERENCE_CODE_MATRIX.json"
CATEGORIES={
    "same_mathematical_method",
    "similar_method",
    "different_method_same_purpose",
    "unavailable_in_cfdx",
}
STATUSES={"planned","missing","partial","implemented","verified","validated","qualified"}
REQUIRED={
    "cfdx_method","mesh_assumptions","implementation_status","verification_status",
    "validation_status","qualification_status","known_differences",
    "benchmark_evidence","parity_claim",
}

def evidence_paths_exist(paths, label):
    for path in paths:
        if not isinstance(path,str):
            raise SystemExit(f"{label}: evidence path is not a string: {path!r}")
        if not (ROOT/path).exists():
            raise SystemExit(f"{label}: missing evidence path: {path}")

def main() -> int:
    data=json.loads(MATRIX.read_text(encoding="utf-8"))
    if data.get("issue") != 461 or data.get("package") != "N14":
        raise SystemExit("N14 matrix must identify issue 461 and package N14")
    refs={item["id"] for item in data["reference_codes"]}
    expected={"openfoam","su2","fluent","cfx","starccm"}
    if refs != expected:
        raise SystemExit(f"reference-code set mismatch: {sorted(refs)}")
    if set(data["comparison_categories"]) != CATEGORIES:
        raise SystemExit("comparison categories do not match the N14 contract")
    for record in data["records"]:
        missing=REQUIRED-set(record)
        if missing:
            raise SystemExit(f"{record.get('id','<unknown>')}: missing {sorted(missing)}")
        for key in ("implementation_status","verification_status","validation_status","qualification_status"):
            if record[key] not in STATUSES:
                raise SystemExit(f"{record['id']}: invalid {key}={record[key]!r}")
        if record["parity_claim"] is not False:
            raise SystemExit(f"{record['id']}: N14 parity_claim must remain false")
        evidence=record.get("benchmark_evidence",[])
        evidence_paths_exist(evidence, record["id"])
        for ref_id,ref in record["references"].items():
            if ref_id not in expected:
                raise SystemExit(f"{record['id']}: unknown reference {ref_id}")
            if ref.get("category") not in CATEGORIES:
                raise SystemExit(f"{record['id']}/{ref_id}: invalid comparison category")
            if not ref.get("source_evidence"):
                raise SystemExit(f"{record['id']}/{ref_id}: missing source evidence")
    for ref_id,path in data["source_dossiers"].items():
        if ref_id not in expected:
            raise SystemExit(f"unknown source dossier {ref_id}")
        if not (ROOT/path).exists():
            raise SystemExit(f"missing source dossier: {path}")
    for item in data["unavailable_reference_capabilities"]:
        if item.get("category") != "unavailable_in_cfdx":
            raise SystemExit("unavailable reference capability has wrong category")
        if item.get("reference") not in expected:
            raise SystemExit("unavailable reference capability has unknown reference")
    print(f"N14 matrix OK: {len(data['records'])} method records, {len(data['unavailable_reference_capabilities'])} unavailable capabilities")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
