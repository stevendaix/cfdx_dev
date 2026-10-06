#!/usr/bin/env python3
"""Generic CFDX validation-result comparator.

Reads CSV files containing a case identifier and one or more quantities of
interest (QoIs), compares them with a reference CSV, and emits a compact
JSON report. No CFD solver assumptions are encoded here: the tool is intended
to be reusable for N9, later physical validation, and external benchmarks.

Required columns in both files: case
Reference and CFDX QoI columns must have identical names.

Usage:
  python scripts/compare_validation_results.py --computed computed.csv --reference reference.csv --qoi Cd Cl
"""
from __future__ import annotations
import argparse, csv, json, math, pathlib

def read_csv(path: pathlib.Path) -> dict[str, dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as f:
        rows = csv.DictReader(f)
        if not rows.fieldnames or "case" not in rows.fieldnames:
            raise ValueError(f"{path}: CSV must contain a 'case' column")
        out = {}
        for row in rows:
            case = row["case"]
            out[case] = {k: float(v) for k, v in row.items() if k != "case" and v not in ("", None)}
        return out

def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--computed", type=pathlib.Path, required=True)
    p.add_argument("--reference", type=pathlib.Path, required=True)
    p.add_argument("--qoi", nargs="+", required=True)
    p.add_argument("--output", type=pathlib.Path)
    a = p.parse_args()
    c, r = read_csv(a.computed), read_csv(a.reference)
    rows = []
    for case in sorted(set(c) & set(r)):
        for q in a.qoi:
            if q not in c[case] or q not in r[case]:
                continue
            cv, rv = c[case][q], r[case][q]
            abs_err = abs(cv-rv)
            rel_err = abs_err/abs(rv) if rv != 0 else None
            rows.append({"case":case,"qoi":q,"computed":cv,"reference":rv,
                         "absolute_error":abs_err,"relative_error":rel_err})
    report={"computed":str(a.computed),"reference":str(a.reference),"results":rows}
    data=json.dumps(report, indent=2)+"\n"
    if a.output: a.output.write_text(data, encoding="utf-8")
    else: print(data, end="")
    return 0
if __name__ == "__main__":
    raise SystemExit(main())
