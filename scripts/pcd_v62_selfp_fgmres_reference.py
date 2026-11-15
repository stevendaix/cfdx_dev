@@ -0,0 +1,174 @@
#!/usr/bin/env python3
"""V62: bounded, read-only PETSc-selfp Schur diagnostic on ORIGINAL six captures.

This is NOT a PCD implementation or a CFDX production solve. Both the full
coupled residual b-Ax and the independent Schur residual are re-evaluated.
SciPy GMRES applies a LEFT preconditioner with block-upper/lower ILU velocity and a sparse
approximate pressure Schur Sp=C-D diag(Auu)^-1 G. All tests use IDENTICAL A,b,x0.
No tolerance from CFDX's physics or PCD production is changed.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.io import mmread
from scipy.sparse import diags, issparse
from scipy.sparse.linalg import LinearOperator, gmres, splu, spilu

METHODS=("pcd","lsc","bfbt","block_local","simple","simplec")
SCHEMA="cfdx-bfs-pcd-selfp-fgmres-reference/v62"


def _sha(path:Path)->str:
    if not path.is_file():
        raise ValueError(f"missing captured BFS matrix: {path}")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read(path:Path):
    matrix=mmread(path)
    if issparse(matrix):
        matrix=matrix.tocsc()
        if not np.isfinite(matrix.data).all():
            raise ValueError(f"nonfinite sparse matrix in {path}")
    else:
        matrix=np.asarray(matrix,dtype=np.float64)
        if matrix.ndim!=2 or not np.isfinite(matrix).all():
            raise ValueError(f"nonfinite or incorrect dense input in {path}")
    return matrix


def analyze(folder:Path, *, max_cells:int=256,
            ilu_drop_tol:float=1e-3, ilu_fill_factor:float=5.0,
            relative_target:float=1e-9)->dict:
    if not (1<=max_cells<=256):
        raise ValueError("V62 reference dense/ILU audit is bounded to 256 pressure cells")
    if not (math.isfinite(ilu_drop_tol) and 0.<=ilu_drop_tol<=0.1 and
            math.isfinite(ilu_fill_factor) and 1.<=ilu_fill_factor<=20. and
            math.isfinite(relative_target) and 1e-12<=relative_target<=1e-5):
        raise ValueError("unbounded ILU or diagnostic tolerance setting")
    hashes={name:{key:_sha(folder/f"{name}_{key}.mtx")
                  for key in ("A","b","x0")} for name in METHODS}
    if any(value!=hashes["pcd"] for value in hashes.values()):
        raise ValueError("six Schur captures are not identical A,b,x0; comparison refused")
    A=_read(folder/"pcd_A.mtx")
    b=_read(folder/"pcd_b.mtx")
    x0=_read(folder/"pcd_x0.mtx")
    if not issparse(A) or A.shape[0]!=A.shape[1] or A.shape[0]%4:
        raise ValueError("expected captured square sparse 4N coupled system")
    n=A.shape[0]//4
    if n>max_cells or A.nnz<4*n:
        raise ValueError("invalid or excessive original captured system")
    if b.shape!=(4*n,1) or x0.shape!=(4*n,1):
        raise ValueError("captured RHS or initial iterate dimensions changed")
    b=np.asarray(b).ravel();x0=np.asarray(x0).ravel()
    bnorm=float(np.linalg.norm(b))
    if not math.isfinite(bnorm) or bnorm<=0.:
        raise ValueError("bad original coupled right hand side")
    U=A[:3*n,:3*n];G=A[:3*n,3*n:]
    D=A[3*n:,:3*n];C=A[3*n:,3*n:]
    diagonal=U.diagonal()
    if (len(diagonal)!=3*n or not np.all(np.isfinite(diagonal)) or
            np.any(diagonal<=0.)):
        raise ValueError("velocity diagonal must be strictly positive and finite")
    # Original C is retained: in particular, the pressure gauge is not
    # replaced by an arbitrary new reference equation.
    Sp=(C-D@diags(1./diagonal)@G).tocsc()
    sv=np.linalg.svd(Sp.toarray(),compute_uv=False)
    if not np.isfinite(sv).all() or sv[0]<=0. or sv[-1]<=sv[0]*1e-12:
        raise ValueError("captured selfp pressure Schur is singular/undecidable")
    try:
        lu_U=splu(U);lu_S=splu(Sp)
        ilu_U=spilu(U,drop_tol=ilu_drop_tol,fill_factor=ilu_fill_factor)
        ilu_S=spilu(Sp,drop_tol=ilu_drop_tol,fill_factor=ilu_fill_factor)
    except RuntimeError as e:
        raise ValueError(f"selfp LU/ILU setup failed: {e}") from e
    b_initial=b-A@x0
    initial_rel=float(np.linalg.norm(b_initial)/bnorm)
    if not np.isfinite(initial_rel):
        raise ValueError("initial coupled residual invalid")
    def solve_family(u_solve,s_solve):
        def action(residual):
            u=u_solve(residual[:3*n]);pr=residual[3*n:]-D@u
            pressure=s_solve(pr)
            return np.concatenate((u-u_solve(G@pressure),pressure))
        count=0
        history=[]
        def cb(value):
            nonlocal count
            count+=1
            if len(history)<5000:history.append(float(value))
        matrix=LinearOperator(A.shape,matvec=action,dtype=np.float64)
        start=perf_counter()
        x,info=gmres(A,b,x0=x0,M=matrix,rtol=relative_target,atol=0.,
                     restart=40,maxiter=40,callback=cb,callback_type="pr_norm")
        elapsed=perf_counter()-start
        # This is the actual original A, never an Arnoldi/Hessenberg estimate.
        exact=b-A@x
        true_rel=float(np.linalg.norm(exact)/bnorm)
        momentum=float(np.linalg.norm(exact[:3*n]))
        continuity=float(np.linalg.norm(exact[3*n:]))
        if not all(map(math.isfinite,(true_rel,momentum,continuity,elapsed))):
            raise ValueError("nonfinite GMRES result on original coupled system")
        return {
            "solver":"scipy_gmres_left_preconditioned_block_schur", "iterations":count,
            "status_code":int(info),"true_relative_residual":true_rel,
            "true_momentum_l2":momentum,"true_continuity_l2":continuity,
            "meets_reference_target":bool(info==0 and true_rel<=relative_target),
            "elapsed_seconds_diagnostic_only":elapsed,
            "preconditioned_residual_history":history,
        }
    reference=solve_family(lu_U.solve,lu_S.solve)
    approximate=solve_family(ilu_U.solve,ilu_S.solve)
    result={
        "schema":SCHEMA,"purpose":"independent selfp baseline, NOT production PCD",
        "original_A_b_x0_sha256":hashes["pcd"],
        "six_original_captures_identical":True,
        "methods_present":list(METHODS),
        "rows":int(A.shape[0]),"cells":n,"nnz_A":int(A.nnz),
        "nnz_Sp":int(Sp.nnz),
        "selfp_condition_number":float(sv[0]/sv[-1]),
        "selfp_rank":int((sv>sv[0]*1e-12).sum()),
        "original_initial_true_relative":initial_rel,
        "relative_target_for_algebraic_reference_only":relative_target,
        "ilu_drop_tol":ilu_drop_tol,"ilu_fill_factor":ilu_fill_factor,
        "reference_direct_blocks":reference,"approximate_ilu_blocks":approximate,
        "PCD_production_converged":False,"six_field_newton_accepted":False,
        "BFS_physical_qualified":False,
        "limitations":"Only the original frozen 4N saddle system, not a candidate 6N Jacobian or full CFD production run",
    }
    json.dumps(result,allow_nan=False)
    return result


def main()->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--captures",type=Path,required=True)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--maximum-cells",type=int,default=256)
    parser.add_argument("--ilu-drop-tol",type=float,default=1e-3)
    parser.add_argument("--ilu-fill-factor",type=float,default=5.)
    parser.add_argument("--relative-target",type=float,default=1e-9)
    args=parser.parse_args()
    result=analyze(args.captures,max_cells=args.maximum_cells,
                   ilu_drop_tol=args.ilu_drop_tol,
                   ilu_fill_factor=args.ilu_fill_factor,
                   relative_target=args.relative_target)
    serialized=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(serialized,encoding="utf-8")
    approx=result["approximate_ilu_blocks"]
    print(f"V62_SELFP_REFERENCE cells={result['cells']} iter={approx['iterations']} "
          f"true={approx['true_relative_residual']:.12g} "
          f"reference_ok={int(approx['meets_reference_target'])} "
          "pcd_qualified=0 bfs_qualified=0")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
