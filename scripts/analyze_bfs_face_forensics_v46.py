"""Fail-closed, read-only V46 BFS post-SST face audit, NOT physical validation.

Consumes POST_SST_FACE_CELL_V46 and POST_SST_FACE_V46 emitted by the production
solver. Never turns linear/SST certificates into a CFD qualification.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
from typing import Any

CELL_PREFIX = "POST_SST_FACE_CELL_V46 "
FACE_PREFIX = "POST_SST_FACE_V46 "
ERROR_PREFIX = "POST_SST_FACE_ERROR_V46 "
NUMERICS = ("u", "conv", "visc", "transpose", "pressure", "body", "mass_div",
            "bounded", "collocated_delta", "physical", "face_operator", "csr",
            "closure", "closure_scale", "nu_t", "k", "omega")
FACE_VALUES = ("mass", "conv", "visc", "transpose", "collocated_delta")


def _line_fields(line: str, prefix: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in line[len(prefix):].split():
        if "=" not in item:
            raise ValueError("malformed V46 key-value token")
        key, value = item.split("=", 1)
        if not key or not value or key in result:
            raise ValueError("missing/duplicate V46 field")
        result[key] = value
    return result


def _int(d: dict[str, str], key: str) -> int:
    value = d[key]
    if not value or value.startswith("+") or value.strip("-0123456789"):
        raise ValueError("bad V46 integer " + key)
    return int(value)


def _num(d: dict[str, str], key: str) -> float:
    x = float(d[key])
    if not math.isfinite(x):
        raise ValueError("nonfinite V46 field " + key)
    return x


def _close(a: float, b: float, scale: float) -> bool:
    return abs(a-b) <= 2e-10 * max(1e-15, scale)


def certify(log: str, expected_cells: set[int] | None = None) -> dict[str, Any]:
    summaries: dict[tuple[int,int,int], dict[str,str]] = {}
    faces: dict[tuple[int,int,int], list[dict[str,str]]] = {}
    for line in log.splitlines():
        if line.startswith(ERROR_PREFIX):
            raise ValueError("production face diagnostic failed: " + line)
        if line.startswith(CELL_PREFIX):
            d = _line_fields(line,CELL_PREFIX)
            key = tuple(_int(d,k) for k in ("iteration","component","cell"))
            if key in summaries:
                raise ValueError("duplicate V46 cell report")
            summaries[key] = d
        if line.startswith(FACE_PREFIX):
            d = _line_fields(line,FACE_PREFIX)
            key = tuple(_int(d,k) for k in ("iteration","component","cell"))
            faces.setdefault(key,[]).append(d)
    if not summaries or set(faces) != set(summaries):
        raise ValueError("missing V46 face/cell evidence")
    by_phase: dict[str, dict[str, float]] = {}
    count=0
    max_gap=0.0
    face_counts=[]
    for key in sorted(summaries):
        iteration, component, cell = key
        s=summaries[key]
        if iteration < 1 or component not in (0,1,2) or cell < 0:
            raise ValueError("invalid V46 cell identifiers")
        if expected_cells is not None and cell not in expected_cells:
            raise ValueError("unexpected V46 cell")
        if any(_int(s,field)!=0 for field in ("committed","qualified")):
            raise ValueError("false V46 physical qualification claim")
        if any(_int(s,field)!=1 for field in ("face_ok","csr_ok","all_csr_ok")):
            raise ValueError("face/CSR closure not certified")
        for name in NUMERICS:
            _num(s,name)
        bounded = _int(s,"bounded_mode")
        if bounded not in (0,1):
            raise ValueError("invalid boundedness switch")
        face_rows=faces[key]
        expected=_int(s,"faces")
        if expected != len(face_rows) or expected < 1:
            raise ValueError("missing/extra V46 faces")
        face_counts.append(expected)
        seen=set()
        totals={name:0.0 for name in FACE_VALUES}
        for f in face_rows:
            fid=_int(f,"face")
            if fid in seen or fid < 0:
                raise ValueError("duplicate/bad V46 face id")
            seen.add(fid)
            owner, neighbour, sign = (_int(f,name) for name in ("owner","neighbour","sign"))
            if owner < 0 or neighbour < -1 or sign not in (-1,1):
                raise ValueError("invalid V46 orientation")
            if (sign==1 and cell!=owner) or (sign==-1 and cell!=neighbour):
                raise ValueError("inconsistent V46 owner/neighbour sign")
            if neighbour==-1 and f.get("patch") in (None,"internal"):
                raise ValueError("boundary patch identity missing")
            if neighbour>=0 and f.get("patch")!="internal":
                raise ValueError("internal face misidentified")
            if any(_int(f,name)!=0 for name in ("committed","qualified")):
                raise ValueError("false V46 face qualification")
            for name in FACE_VALUES:
                totals[name] += _num(f,name)
        scale=max(abs(_num(s,name)) for name in ("conv","visc","transpose","mass_div","collocated_delta"))
        mapping=("conv","visc","transpose","mass_div","collocated_delta")
        for name in mapping:
            source="mass" if name=="mass_div" else name
            if not _close(totals[source],_num(s,name),scale):
                raise ValueError("face sum does not reproduce cell " + name)
        physical=(_num(s,"conv")+_num(s,"visc")+_num(s,"transpose")+
                  _num(s,"pressure")+_num(s,"body"))
        if not _close(physical,_num(s,"physical"),max(scale,abs(physical))):
            raise ValueError("physical cell decomposition does not close")
        if not _close(_num(s,"bounded"),-bounded*_num(s,"u")*_num(s,"mass_div"),
                      max(scale,abs(_num(s,"bounded")))):
            raise ValueError("bounded divergence incorrect")
        if not _close(physical+_num(s,"bounded"),_num(s,"face_operator"),
                      max(scale,abs(physical))):
            raise ValueError("physical vs bounded operator mismatch")
        gap=_num(s,"csr")-_num(s,"face_operator")
        closure=_num(s,"closure")
        arithmetic_scale=_num(s,"closure_scale")
        if arithmetic_scale<0 or not _close(gap,closure,arithmetic_scale) or \
                not _close(closure,0.0,arithmetic_scale):
            raise ValueError("CSR/face evidence inconsistent")
        max_gap=max(max_gap,abs(gap))
        phase=f"{iteration}:{component}:{cell}"
        by_phase[phase]={k:_num(s,k) for k in ("physical","face_operator","csr",
                                                 "pressure","visc","conv","transpose","nu_t")}
        count += 1
    # Compare exactly consecutive physical iterates. A reset/restart is a
    # separate trajectory: never invent a comparison across a missing state.
    # The bounded term and pressure/source deltas must be included; face fluxes
    # alone are not the full nonlinear FV momentum residual.
    term_names=("conv","visc","transpose","pressure","body","bounded")
    changes=[]
    for iteration,component,cell in sorted(summaries):
        key=(iteration,component,cell)
        previous=(iteration-1,component,cell)
        if previous not in summaries:
            continue
        before,after=summaries[previous],summaries[key]
        deltas={name:_num(after,name)-_num(before,name) for name in term_names}
        residual_change=_num(after,"csr")-_num(before,"csr")
        residual_scale=max(abs(_num(after,"csr")),abs(_num(before,"csr")),
                           *(abs(x) for x in deltas.values()))
        if not _close(sum(deltas.values()),residual_change,residual_scale):
            raise ValueError("nonlinear phase delta fails FV/CSR recomposition")
        prior_face={_int(face,"face"):face for face in faces[previous]}
        latest_face={_int(face,"face"):face for face in faces[key]}
        if set(prior_face)!=set(latest_face):
            raise ValueError("face identity changed across consecutive phases")
        changed_faces=[]
        for fid in sorted(prior_face):
            a,b=prior_face[fid],latest_face[fid]
            if any(a[name]!=b[name] for name in ("owner","neighbour","sign","patch")):
                raise ValueError("face topology changed across consecutive phases")
            increments={name:_num(b,name)-_num(a,name) for name in FACE_VALUES}
            previous_mass=_num(a,"mass")
            current_mass=_num(b,"mass")
            flipped=(previous_mass*current_mass<0.0)
            changed_faces.append({"face":fid,"patch":b["patch"],
                                  "mass_direction_reversed":flipped,
                                  "delta_mass":increments["mass"],
                                  "delta_convection":increments["conv"],
                                  "delta_viscosity":increments["visc"],
                                  "delta_transpose":increments["transpose"],
                                  "delta_collocated":increments["collocated_delta"],
                                  "delta_fv_face":increments["conv"]+increments["visc"]+
                                                  increments["transpose"]})
        changed_faces.sort(key=lambda row:(-abs(row["delta_fv_face"]),row["face"]))
        changes.append({"from":iteration-1,"to":iteration,"component":component,
                        "cell":cell,"delta_csr":residual_change,
                        "delta_terms":deltas,
                        "delta_nu_t":_num(after,"nu_t")-_num(before,"nu_t"),
                        "dominant_terms":sorted(deltas,key=lambda n:-abs(deltas[n])),
                        "mass_flux_direction_reversals":sum(
                            1 for row in changed_faces if row["mass_direction_reversed"]),
                        "most_changed_faces":changed_faces[:6]})
    return {
        "status": "DIAGNOSTIC_ONLY_UNQUALIFIED",
        "physically_converged": False,
        "cells": count,
        "faces": sum(face_counts),
        "iterations": sorted({k[0] for k in summaries}),
        "components": sorted({k[1] for k in summaries}),
        "selected_cells": sorted({k[2] for k in summaries}),
        "max_matrix_face_gap": max_gap,
        "per_phase": by_phase,
        "consecutive_phase_deltas": changes,
    }


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--cells",default="1,2,12,30")
    args=parser.parse_args()
    expected={int(x) for x in args.cells.split(",")}
    result=certify(args.log.read_text(),expected_cells=expected)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(f"V46 diagnostic-only: {result['cells']} rows, {result['faces']} faces, physical qualification=NO")

if __name__ == "__main__":
    main()
