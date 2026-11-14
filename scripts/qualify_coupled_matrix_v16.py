"""Independent linear ALGORITHMIC oracle for an exact CFDX coupled KSP snapshot.

Requires optional SciPy; reads three MatrixMarket files emitted at the exact
Krylov call point.  Never certifies a nonlinear or experimental CFD result.
"""
import argparse
import json
import math
from pathlib import Path


def qualify(prefix: Path) -> dict:
    import numpy as np
    from scipy.io import mmread
    from scipy.sparse.linalg import splu, spilu, LinearOperator, gmres
    try:
        A=mmread(str(prefix)+'_A.mtx').tocsr()
        b=np.asarray(mmread(str(prefix)+'_b.mtx')).reshape(-1)
        x0=np.asarray(mmread(str(prefix)+'_x0.mtx')).reshape(-1)
    except (OSError,ValueError) as exc:
        raise ValueError(f'invalid/incomplete original matrix capture: {exc}') from exc
    n=A.shape[0]
    if A.shape[0]!=A.shape[1] or n%4 or not (4<=n<=6000) or b.shape!=(n,) or x0.shape!=(n,) or not (
            np.isfinite(A.data).all() and np.isfinite(b).all() and np.isfinite(x0).all()):
        raise ValueError('invalid coupled CSR or U/p field dimensions/nonfinite data')
    nv=3*(n//4)
    bnorm=float(np.linalg.norm(b))
    if not bnorm>0 or not math.isfinite(bnorm):
        raise ValueError('unusable original coupled right-hand side')
    def defects(x):
        r=b-A@x
        return {'full_relative':float(np.linalg.norm(r)/bnorm),
                'momentum_relative':float(np.linalg.norm(r[:nv])/max(np.linalg.norm(b[:nv]),1e-300)),
                'continuity_relative':float(np.linalg.norm(r[nv:])/max(np.linalg.norm(b[nv:]),1e-300)),
                'solution_norm':float(np.linalg.norm(x))}
    report={'format':'CFDX-COUPLED-ORIGINAL-CSR-FORENSICS/v16','unknowns':n,
            'nnz':int(A.nnz),'initial':defects(x0),
            'physical_validation':'NOT_QUALIFIED','no_schur_substitution_in_original_operator':True}
    try:
        lu=splu(A.tocsc())
        report['direct_lu']=defects(lu.solve(b))
    except (RuntimeError,ValueError) as exc:
        report['direct_lu_error']=str(exc)
    try:
        ilu=spilu(A.tocsc(),drop_tol=1e-5,fill_factor=30)
        preconditioner=LinearOperator(A.shape,matvec=ilu.solve)
        history=[]
        sol,flag=gmres(A,b,x0=x0,M=preconditioner,atol=0,rtol=1e-12,
                      restart=160,maxiter=10,callback=history.append,
                      callback_type='pr_norm')
        report['ilut_gmres']={**defects(sol),'status':int(flag),'inner_iterations':len(history),
                              'true_target':1e-12,
                              'original_system_converged':flag==0 and defects(sol)['full_relative']<=1e-12}
    except (RuntimeError,ValueError) as exc:
        report['ilut_gmres_error']=str(exc)
    return report


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--prefix',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    try:
        result=qualify(args.prefix)
        rc=0 if result.get('direct_lu',{}).get('full_relative',math.inf)<1e-10 else 3
    except (ImportError,ValueError,RuntimeError) as exc:
        rc=3
        result={'format':'CFDX-COUPLED-ORIGINAL-CSR-FORENSICS/v16',
                'physical_validation':'NOT_QUALIFIED','error':str(exc)}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return rc

if __name__=='__main__':
    raise SystemExit(main())
