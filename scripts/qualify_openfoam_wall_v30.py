"""Strict offline CFDX/OpenFOAM Foundation omega-wall parity checker.

Accepts independent exports recorded at the same mesh and state. It NEVER
constructs an OpenFOAM reference from CFDX fields and cannot certify without
an attested OpenFOAM-native export. No OpenFOAM binaries are invoked here.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path

_HEX64 = re.compile(r"^[0-9a-f]{64}$")
_GROUPS = {"cells": ("cell", ("omega", "G", "nut", "corner_faces")),
           "faces": ("face", ("owner", "coverage", "weight", "omega", "G", "nut")),
           "matrix_rows": ("row", ("diagonal", "rhs"))}


def _load(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("schema") != "CFDX_OPENFOAM_OMEGA_V30":
        raise ValueError("missing/incorrect V30 oracle schema")
    for name in ("mesh_sha256", "state_sha256"):
        if not isinstance(data.get(name), str) or not _HEX64.fullmatch(data[name]):
            raise ValueError(f"unverified {name}")
    if not isinstance(data.get("provenance"), dict):
        raise ValueError("missing independent provenance")
    # A pair of equally truncated exports is NOT a complete wall certificate.
    # Every expected owner, face, and affected matrix row must be declared.
    scope=data.get("scope")
    if not isinstance(scope, dict) or set(scope) != set(_GROUPS):
        raise ValueError("missing explicit wall-oracle scope")
    for group in _GROUPS:
        ids=scope[group]
        if (not isinstance(ids,list) or not ids or
            any(type(x) is not int or x<0 for x in ids) or
            ids!=sorted(set(ids))):
            raise ValueError(f"invalid/incomplete scope for {group}")
    for group, (id_name, keys) in _GROUPS.items():
        entries = data.get(group)
        if not isinstance(entries, list) or not entries:
            raise ValueError(f"missing {group} data")
        seen = set()
        for entry in entries:
            if not isinstance(entry, dict) or not isinstance(entry.get(id_name), int):
                raise ValueError(f"malformed {group} identifier")
            key = entry[id_name]
            if key in seen or key < 0:
                raise ValueError(f"duplicated/negative {group} identifier")
            seen.add(key)
            for field in keys:
                value = entry.get(field)
                if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
                    raise ValueError(f"nonfinite/missing {group}.{field}")
            if group == "faces":
                if entry["owner"] not in scope["cells"] or entry["owner"] not in scope["matrix_rows"]:
                    raise ValueError("face owner omitted from wall cell/matrix scope")
                if not 0 <= entry["coverage"] <= 1:
                    raise ValueError("invalid patch coverage")
                expected = 0. if entry["coverage"] <= 0.1 else (entry["coverage"] - 0.1)/0.9
                if abs(entry["weight"] - expected) > 2e-12:
                    raise ValueError("patch weight inconsistent with OpenFOAM Foundation v11")
            if group == "matrix_rows":
                cols = entry.get("entries")
                if not isinstance(cols, list) or not cols:
                    raise ValueError("missing complete matrix row")
                col_ids = set()
                for val in cols:
                    if not (isinstance(val, list) and len(val) == 2 and
                            isinstance(val[0], int) and val[0] >= 0 and
                            isinstance(val[1], (float, int)) and not isinstance(val[1], bool) and
                            math.isfinite(val[1]) and val[0] not in col_ids):
                        raise ValueError("invalid/missing/duplicate matrix coefficient")
                    col_ids.add(val[0])
                diag = next((v for j, v in cols if j == key), None)
                if diag is None or not math.isclose(diag, entry["diagonal"], rel_tol=1e-13, abs_tol=1e-13):
                    raise ValueError("inconsistent matrix diagonal")
        if seen != set(scope[group]):
            raise ValueError(f"missing/extra wall-oracle entries in {group}")
    return data


def compare(a: dict, b: dict, *, rtol: float=1e-8, atol: float=1e-10,
            verified_native_run: bool=False) -> dict:
    if not (math.isfinite(rtol) and math.isfinite(atol) and rtol >= 0 and atol >= 0):
        raise ValueError("invalid comparison tolerances")
    for field in ("mesh_sha256", "state_sha256"):
        if a[field] != b[field]:
            raise ValueError(f"input mismatch: {field}; cannot compare different states")
    if (a["provenance"].get("engine") != "OpenFOAM Foundation" or
            b["provenance"].get("engine") != "CFDX"):
        raise ValueError("reference must be independent OpenFOAM, candidate CFDX")
    # The known OpenFOAM Foundation v11 codepath is distinct from OpenCFD forks.
    if a["provenance"].get("version") != "11":
        raise ValueError("unsupported OpenFOAM Foundation version")
    if a["scope"] != b["scope"]:
        raise ValueError("wall-oracle declared scope differs between codes")
    native_claim = a["provenance"].get("native_export") is True and all(
        isinstance(a["provenance"].get(key), str) and a["provenance"][key]
        for key in ("executable", "case", "source_commit"))
    errors=[]
    max_relative=0.
    compared=0
    def check(path: str,x: float,y: float):
        nonlocal max_relative,compared
        compared+=1
        scale=atol+rtol*max(abs(x),abs(y))
        diff=abs(x-y)
        if diff>scale and len(errors)<30:
            errors.append({"at":path,"reference":x,"candidate":y,
                           "error":diff,"allowed":scale})
        max_relative=max(max_relative,diff/max(scale,1e-300))
    for group,(id_field,fields) in _GROUPS.items():
        ref={e[id_field]:e for e in a[group]}
        got={e[id_field]:e for e in b[group]}
        if ref.keys()!=got.keys():
            raise ValueError(f"different owner/face/matrix IDs in {group}")
        for index,e in ref.items():
            f=got[index]
            for field in fields:
                check(f"{group}[{index}].{field}",e[field],f[field])
            if group == "matrix_rows":
                r_entries={k:v for k,v in e["entries"]}
                c_entries={k:v for k,v in f["entries"]}
                if r_entries.keys()!=c_entries.keys():
                    raise ValueError(f"different matrix sparsity in row {index}")
                for j,x in r_entries.items():
                    check(f"matrix_rows[{index}].entries[{j}]",x,c_entries[j])
    # A self-declared native_export=true in JSON is NOT enough. The CLI must
    # verify the actual reference run log against the reference's SHA-256.
    native = native_claim and verified_native_run
    status = "FAIL_NUMERICAL_MISMATCH" if errors else (
        "PASS_NATIVE" if native else "INCOMPLETE_NO_NATIVE_OPENFOAM_EXPORT")
    return {"status":status,"native_reference":bool(native),"compared_scalars":compared,
            "maximum_normalized_error":max_relative,"first_errors":errors,
            "qualified":status=="PASS_NATIVE"}


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("openfoam_json",type=Path)
    parser.add_argument("cfdx_json",type=Path)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--rtol",type=float,default=1e-8)
    parser.add_argument("--atol",type=float,default=1e-10)
    parser.add_argument("--native-run-log",type=Path,
        help="Mandatory independent OpenFOAM runtime log for a positive certificate")
    args=parser.parse_args()
    try:
        ref=_load(args.openfoam_json)
        got=_load(args.cfdx_json)
        verified=False
        if args.native_run_log:
            expected=ref["provenance"].get("run_log_sha256")
            if not isinstance(expected,str) or not _HEX64.fullmatch(expected):
                raise ValueError("missing native run-log SHA256 in provenance")
            observed=hashlib.sha256(args.native_run_log.read_bytes()).hexdigest()
            if observed != expected:
                raise ValueError("native OpenFOAM runtime log SHA256 mismatch")
            verified=True
        out=compare(ref,got,rtol=args.rtol,atol=args.atol,
                    verified_native_run=verified)
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        out={"status":"FAIL_INVALID_EVIDENCE","qualified":False,"reason":str(exc)}
    if args.output:
        args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(out,sort_keys=True))
    return 0 if out.get("qualified") else 3

if __name__ == "__main__":
    sys.exit(main())
