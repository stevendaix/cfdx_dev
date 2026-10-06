#!/usr/bin/env python3
"""Aggregate generic CFDX validation run CSVs into one deterministic table."""
from __future__ import annotations
import argparse, csv, pathlib

def main() -> int:
    p=argparse.ArgumentParser()
    p.add_argument("--input-dir",type=pathlib.Path,required=True)
    p.add_argument("--pattern",default="*.csv")
    p.add_argument("--output",type=pathlib.Path,required=True)
    a=p.parse_args()
    files=sorted(a.input_dir.glob(a.pattern))
    if not files: raise SystemExit(f"No files matching {a.pattern} in {a.input_dir}")
    rows=[]; fields=[]
    for f in files:
        with f.open(newline="",encoding="utf-8") as h:
            r=csv.DictReader(h)
            if r.fieldnames:
                fields += [x for x in r.fieldnames if x not in fields]
                rows.extend(dict(x) for x in r)
    with a.output.open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=fields,extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    return 0
if __name__=="__main__": raise SystemExit(main())
