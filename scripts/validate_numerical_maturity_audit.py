#!/usr/bin/env python3
"""Validate the machine-readable numerical-method audits for #461.

Handles both audit schemas: the package-level maturity audit (``packages``, keyed
by N1-N17) and the finer-grained requirement matrix (``requirements``, keyed by
``package_id``). Every evidence entry must name a file that actually exists in
the repository: an audit whose evidence cannot be opened does not document
maturity, it only asserts it, which is exactly what rule R001 forbids.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
ALLOWED={"planned","missing","partial","implemented","verified","validated","qualified"}
REQUIRED={"id","name","status","evidence","gaps"}
REQUIRED_MATRIX={"package_id","requirement","status","evidence"}
PROMOTED={"verified","validated","qualified"}
REPO_ROOT=Path(__file__).resolve().parents[1]

def check_evidence(label,evidence,dangling):
    if not isinstance(evidence,list): raise SystemExit(f"{label}: evidence must be a list")
    for path in evidence:
        if not isinstance(path,str): raise SystemExit(f"{label}: evidence entries must be strings")
        if not (REPO_ROOT/path).exists(): dangling.append((label,path))

def validate_packages(packages,dangling):
    ids=set()
    for item in packages:
        missing=REQUIRED-set(item)
        if missing: raise SystemExit(f"{item.get('id','<unknown>')}: missing keys {sorted(missing)}")
        ident=item["id"]
        if ident in ids: raise SystemExit(f"duplicate package id: {ident}")
        ids.add(ident)
        if item["status"] not in ALLOWED: raise SystemExit(f"{ident}: invalid status {item['status']!r}")
        if not isinstance(item["gaps"],list): raise SystemExit(f"{ident}: gaps must be a list")
        check_evidence(ident,item["evidence"],dangling)
        if item["status"] in PROMOTED and not item["evidence"]:
            raise SystemExit(f"{ident}: promoted status requires evidence")
    expected={f"N{i}" for i in range(1,18)}
    if ids!=expected: raise SystemExit(f"expected N1-N17 exactly, got {sorted(ids)}")
    counts={s:sum(x["status"]==s for x in packages) for s in sorted(ALLOWED)}
    print(f"  {len(packages)} packages, status counts: {counts}")

def validate_requirements(requirements,dangling):
    seen=set()
    for item in requirements:
        missing=REQUIRED_MATRIX-set(item)
        if missing:
            raise SystemExit(f"{item.get('requirement','<unknown>')}: missing keys {sorted(missing)}")
        label=f"{item['package_id']}/{item['requirement']}"
        if label in seen: raise SystemExit(f"duplicate requirement: {label}")
        seen.add(label)
        if item["package_id"] not in {f"N{i}" for i in range(1,18)}:
            raise SystemExit(f"{label}: package_id outside N1-N17")
        if item["status"] not in ALLOWED: raise SystemExit(f"{label}: invalid status {item['status']!r}")
        check_evidence(label,item["evidence"],dangling)
        if item["status"] in PROMOTED and not item["evidence"]:
            raise SystemExit(f"{label}: promoted status requires evidence")
    counts={s:sum(x["status"]==s for x in requirements) for s in sorted(ALLOWED)}
    print(f"  {len(requirements)} requirements, status counts: {counts}")

def main():
    if len(sys.argv)<2:
        print("usage: validate_numerical_maturity_audit.py AUDIT.json [AUDIT.json ...]",
              file=sys.stderr); return 2
    dangling=[]
    for arg in sys.argv[1:]:
        data=json.loads(Path(arg).read_text(encoding="utf-8"))
        if data.get("issue")!=461: raise SystemExit(f"{arg}: audit issue must be 461")
        print(f"validating #461 audit: {arg}")
        if "packages" in data:
            packages=data["packages"]
            if not isinstance(packages,list) or not packages:
                raise SystemExit(f"{arg}: packages must be non-empty")
            validate_packages(packages,dangling)
        elif "requirements" in data:
            requirements=data["requirements"]
            if not isinstance(requirements,list) or not requirements:
                raise SystemExit(f"{arg}: requirements must be non-empty")
            validate_requirements(requirements,dangling)
        else:
            raise SystemExit(f"{arg}: audit has neither 'packages' nor 'requirements'")
    if dangling:
        for label,path in dangling:
            print(f"dangling evidence: {label} cites {path}, which does not exist",file=sys.stderr)
        raise SystemExit(f"{len(dangling)} evidence path(s) do not exist")
    return 0
if __name__=="__main__": raise SystemExit(main())
