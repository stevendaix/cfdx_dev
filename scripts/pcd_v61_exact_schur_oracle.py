#!/usr/bin/env python3
"""V61 independent frozen 4N BFS/PCD oracle (optional SciPy; never CFD qualification).

Reads the six ORIGINAL MatrixMarket A,b,x0 captures, verifies their SHA-256
identity before doing any algebra, forms S=C-D Auu^-1 G, and directly solves
A x=b. This diagnoses a PCD model failure without renaming an LU solution
as a production PCD or as a six-field Newton convergence certificate.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.io import mmread
from scipy.sparse import csc_matrix, issparse
from scipy.sparse.linalg import splu

METHODS=("pcd", "lsc", "bfbt", "block_local", "simple", "simplec")
SCHEMA="cfdx-bfs-pcd-exact-saddle-oracle/v61"


def sha256(path: Path) -> str:
    if not path.is_file():
        raise ValueError(f"missing original Schur capture: {path}")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read(path: Path):
    if not path.is_file():
        raise ValueError(f"missing captured original MatrixMarket input: {path}")
    value=mmread(path)
    if issparse(value):
        if not np.isfinite(value.data).all():
            raise ValueError(f"nonfinite captured sparse matrix: {path}")
        return value.tocsc()
    value=np.asarray(value,dtype=np.float64)
    if value.ndim!=2 or not np.isfinite(value).all():
        raise ValueError(f"nonfinite or malformed capture: {path}")
    return value


def analyze(directory: Path, *, maximum_cells: int=256) -> dict:
    if maximum_cells<1 or maximum_cells>1024:
        raise ValueError("max cells must be 1..1024, no unbounded dense Schur")
    fingerprints={}
    for m in METHODS:
        fingerprints[m]={kind:sha256(directory/f"{m}_{kind}.mtx")
                         for kind in ("A","b","x0")}
    expected=fingerprints[METHODS[0]]
    for model,digest in fingerprints.items():
        if digest!=expected:
            raise ValueError(f"six-method input mismatch: {model} has changed A,b or x0")
    A=_read(directory/"pcd_A.mtx")
    b=_read(directory/"pcd_b.mtx")
    x0=_read(directory/"pcd_x0.mtx")
    if not issparse(A):
        A=csc_matrix(A)
    if A.shape[0]!=A.shape[1] or A.shape[0]%4:
        raise ValueError("expected square frozen coupled 4N matrix")
    n=A.shape[0]//4
    if n>maximum_cells:
        raise ValueError("dense reference Schur disabled on unbounded case")
    if b.shape!=(4*n,1) or x0.shape!=(4*n,1):
        raise ValueError("wrong original RHS or warm-start dimensions")
    b=b[:,0];x0=x0[:,0]
    rhs_norm=float(np.linalg.norm(b))
    if not math.isfinite(rhs_norm) or rhs_norm==0.:
        raise ValueError("zero or nonfinite original coupled RHS")
    Mu=A[:3*n,:3*n]
    G=A[:3*n,3*n:]
    D=A[3*n:,:3*n]
    C=A[3*n:,3*n:]
    try:
        lu_uu=splu(Mu)
        S=C.toarray()-D.dot(lu_uu.solve(G.toarray()))
        reduced_b=b[3*n:]-D.dot(lu_uu.solve(b[:3*n]))
        sv=np.linalg.svd(S,compute_uv=False)
        lu=splu(A)
        solution=lu.solve(b)
    except (RuntimeError,ValueError,np.linalg.LinAlgError) as exc:
        raise ValueError(f"singular or undecidable original 4N saddle system: {exc}") from exc
    if not np.all(np.isfinite(sv)) or sv[-1]<=1e-14*sv[0]:
        raise ValueError("pressure Schur singular at the stated numerical-rank cutoff")
    full_delta=b-A.dot(solution)
    initial_defect=b-A.dot(x0)
    pressure_initial=np.linalg.norm(initial_defect[3*n:])
    momentum_initial=np.linalg.norm(initial_defect[:3*n])
    norm_initial=float(np.linalg.norm(initial_defect))
    rel_direct=float(np.linalg.norm(full_delta)/rhs_norm)
    if not math.isfinite(rel_direct):
        raise ValueError("nonfinite independent direct solution certificate")
    # This is an independent ALGEBRAIC evidence gate, never a physical one.
    if rel_direct>1e-9:
        raise ValueError("sparse LU failed to close the original coupled true residual")
    worst_pressure=int(np.argmax(np.abs(initial_defect[3*n:])))
    schur_residual=S@solution[3*n:]-reduced_b
    result={
        "schema":SCHEMA,
        "source_kind":"frozen-original-4N-production-capture",
        "methods":list(METHODS),
        "same_A_b_x0_sha256":True,
        "sha256":expected,
        "rows":int(A.shape[0]),"cells":n,"nnz":int(A.nnz),
        "original_initial_relative_residual":norm_initial/rhs_norm,
        "initial_momentum_residual_l2":float(momentum_initial),
        "initial_continuity_residual_l2":float(pressure_initial),
        "initial_continuity_squared_fraction":float((pressure_initial/norm_initial)**2),
        "initial_worst_continuity_cell":worst_pressure,
        "pressure_schur_sigma_max":float(sv[0]),
        "pressure_schur_sigma_min":float(sv[-1]),
        "pressure_schur_condition_number":float(sv[0]/sv[-1]),
        "pressure_schur_rhs_l2":float(np.linalg.norm(reduced_b)),
        "direct_LU_original_true_relative":rel_direct,
        "direct_LU_momentum_residual_l2":float(np.linalg.norm(full_delta[:3*n])),
        "direct_LU_continuity_residual_l2":float(np.linalg.norm(full_delta[3*n:])),
        "direct_LU_schur_residual_l2":float(np.linalg.norm(schur_residual)),
        "pcd_converged":False,
        "six_n_newton_accepted":False,
        "physical_BFS_converged":False,
        "description":"Dense Schur and sparse LU are independent reference algebra, not PCD actions",
    }
    json.dumps(result,allow_nan=False)
    return result


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--capture-dir",type=Path,required=True)
    ap.add_argument("--max-cells",type=int,default=256)
    ap.add_argument("--out",type=Path)
    args=ap.parse_args()
    r=analyze(args.capture_dir,maximum_cells=args.max_cells)
    serialized=json.dumps(r,indent=2,sort_keys=True,allow_nan=False)+"\n"
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(serialized,encoding="utf8")
    print(f"V61_PCD_4N_ORACLE rows={r['rows']} nnz={r['nnz']} "
          f"LU_true={r['direct_LU_original_true_relative']:.12g} "
          f"S_condition={r['pressure_schur_condition_number']:.12g} "
          "pcd_qualified=0 physical_qualified=0")
    return 0


if __name__=="__main__":
    raise SystemExit(main())