#!/usr/bin/env python3
"""Validate the machine-readable numerical-method maturity audit for #461."""
from __future__ import annotations
import json, sys
from pathlib import Path
ALLOWED={"planned","missing","partial","implemented","verified","validated","qualified"}
REQUIRED={"id","name","status","evidence","gaps"}
def main():
    if len(sys.argv)!=2:
        print("usage: validate_numerical_maturity_audit.py AUDIT.json",file=sys.stderr); return 2
    data=json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    if data.get("issue")!=461: raise SystemExit("audit issue must be 461")
    packages=data.get("packages")
    if not isinstance(packages,list) or not packages: raise SystemExit("packages must be non-empty")
    ids=set()
    for item in packages:
        missing=REQUIRED-set(item)
        if missing: raise SystemExit(f"{item.get('id','<unknown>')}: missing keys {sorted(missing)}")
        ident=item["id"]
        if ident in ids: raise SystemExit(f"duplicate package id: {ident}")
        ids.add(ident)
        if item["status"] not in ALLOWED: raise SystemExit(f"{ident}: invalid status {item['status']!r}")
        if not isinstance(item["evidence"],list) or not isinstance(item["gaps"],list):
            raise SystemExit(f"{ident}: evidence and gaps must be lists")
        if item["status"] in {"verified","validated","qualified"} and not item["evidence"]:
            raise SystemExit(f"{ident}: promoted status requires evidence")
    expected={f"N{i}" for i in range(1,18)}
    if ids!=expected: raise SystemExit(f"expected N1-N17 exactly, got {sorted(ids)}")
    counts={s:sum(x["status"]==s for x in packages) for s in sorted(ALLOWED)}
    print(f"validated #461 maturity audit: {len(packages)} packages")
    print("status counts:",counts)
    return 0
if __name__=="__main__": raise SystemExit(main())
