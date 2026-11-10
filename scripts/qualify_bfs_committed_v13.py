"""Fail-closed V13 committed-state evidence. Trials and inner SST solves don't count.

The production solver emits GLOBAL_NONLINEAR_V8 on a *tentative* state and
GLOBAL_NONLINEAR_COMMITTED_V13 only after U,p,phi and SST transaction.commit().
A missing or truncated marker means the outer step is not proven accepted.
"""
from __future__ import annotations
import argparse
import json
import math
import sys
from pathlib import Path

FIELDS=("U","p_continuity","k","omega","nut_closure","merit")

def parse(log: str) -> dict:
    pending: dict[int,list[dict]]={}
    committed: list[dict]=[]
    rejection_count=0
    for line in log.splitlines():
        if line.startswith("CFDX nonlinear iteration") and "rejected" in line:
            rejection_count+=1
        if not line.startswith(("GLOBAL_NONLINEAR_V8 ","GLOBAL_NONLINEAR_COMMITTED_V13 ")):
            continue
        pieces=line.split()
        fields={}
        for piece in pieces[1:]:
            if "=" in piece:
                k,v=piece.split("=",1)
                if k in fields:
                    raise ValueError("duplicate global certificate field")
                fields[k]=v
        try:
            iteration=int(fields["iteration"])
            values={name:float(fields[name]) for name in FIELDS}
        except (KeyError,ValueError) as exc:
            raise ValueError("incomplete global nonlinear certificate") from exc
        if iteration<=0 or any(not math.isfinite(v) or v<0 for v in values.values()):
            raise ValueError("nonfinite/invalid global physical residual")
        if values["merit"]+1e-12 < (1-1e-8)*max(values[k] for k in FIELDS if k!="merit"):
            raise ValueError("global merit hides a physical component")
        if pieces[0]=="GLOBAL_NONLINEAR_V8":
            if fields.get("qualified")!="0":
                raise ValueError("tentative iteration cannot qualify BFS")
            pending.setdefault(iteration,[]).append(values)
            continue
        if fields.get("converged") not in ("0","1") or fields.get("qualified")!="0":
            raise ValueError("invalid committed/qualified flags")
        if not pending.get(iteration):
            raise ValueError("committed record without a preceding same-iteration trial")
        matching=[i for i,v in enumerate(pending[iteration]) if all(
            math.isclose(v[k], values[k], rel_tol=1e-8,abs_tol=1e-12) for k in FIELDS)]
        if not matching:
            raise ValueError("committed state differs from every tentative physical certificate")
        del pending[iteration][matching[-1]]
        if committed and iteration<=committed[-1]["iteration"]:
            raise ValueError("duplicate/reordered committed nonlinear iteration")
        committed.append({"iteration":iteration,"fields":values,
            "converged":fields["converged"]=="1","physical_qualified":False})
    return {
        "schema":"cfdx-bfs-committed-v13",
        "trial_count":sum(len(v) for v in pending.values())+len(committed),
        "committed_count":len(committed),
        "rejected_retry_markers":rejection_count,
        "uncommitted_trial_count":sum(len(v) for v in pending.values()),
        "committed_states":committed,
        "all_committed_converged":bool(committed) and all(r["converged"] for r in committed),
        "physical_qualification":False,
    }

def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--log',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--require-committed',action='store_true')
    args=ap.parse_args(argv)
    try:
        result=parse(args.log.read_text(errors='replace'))
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print('BFS_COMMITTED_V13 committed=',result['committed_count'],
              'uncommitted=',result['uncommitted_trial_count'],
              'physical_qualified=0')
        return 0 if (not args.require_committed or result['committed_count']) else 3
    except (ValueError,OSError) as exc:
        print(f'BFS_COMMITTED_V13_ERROR: {exc}',file=sys.stderr)
        return 2
if __name__=='__main__':
    raise SystemExit(main())
