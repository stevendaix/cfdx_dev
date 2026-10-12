#!/usr/bin/env python3
"""Reproducible CFDX N15 performance benchmark runner.

This tool records real measurements only. Missing application-internal timings
remain missing; they are never inferred from iteration counts.
"""
from __future__ import annotations
import argparse, datetime as dt, json, os, platform, resource, shutil, socket, subprocess, sys, time
from pathlib import Path
from typing import Any
SCHEMA_VERSION = "1.0"

def git_sha(repo: Path) -> str | None:
    try:
        return subprocess.check_output(["git","-C",str(repo),"rev-parse","HEAD"], text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError): return None

def compiler_info() -> dict[str, str | None]:
    for name in ("c++","g++","clang++","icpx"):
        exe=shutil.which(name)
        if exe:
            try: version=subprocess.check_output([exe,"--version"],text=True,stderr=subprocess.STDOUT).splitlines()[0]
            except (OSError,subprocess.CalledProcessError): version=None
            return {"name":name,"path":exe,"version":version}
    return {"name":None,"path":None,"version":None}

def gpu_info() -> dict[str,Any]:
    exe=shutil.which("nvidia-smi")
    if not exe: return {"available":False,"reason":"nvidia-smi-not-found"}
    try:
        out=subprocess.check_output([exe,"--query-gpu=name,driver_version,memory.total","--format=csv,noheader"],text=True,stderr=subprocess.STDOUT)
        devices=[]
        for line in out.splitlines():
            f=[x.strip() for x in line.split(",")]
            devices.append({"name":f[0],"driver":f[1] if len(f)>1 else None,"memory":f[2] if len(f)>2 else None})
        return {"available":bool(devices),"devices":devices}
    except (OSError,subprocess.CalledProcessError) as exc: return {"available":False,"reason":str(exc)}

def hardware_metadata() -> dict[str,Any]:
    memory_bytes=None
    p=Path("/proc/meminfo")
    if p.exists():
        for line in p.read_text().splitlines():
            if line.startswith("MemTotal:"): memory_bytes=int(line.split()[1])*1024; break
    return {"hostname":socket.gethostname(),"platform":platform.platform(),"system":platform.system(),"machine":platform.machine(),"processor":platform.processor(),"logical_cpus":os.cpu_count(),"memory_bytes":memory_bytes,"gpu":gpu_info()}

def peak_rss_bytes() -> int | None:
    value=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if value<=0: return None
    return value*1024 if sys.platform!="darwin" else value

def load_json(path:Path)->dict[str,Any]:
    data=json.loads(path.read_text())
    if not isinstance(data,dict): raise ValueError(f"{path} must contain a JSON object")
    return data

def run_once(args:argparse.Namespace,case:dict[str,Any],rank:int|None=None)->dict[str,Any]:
    command=list(args.command); env=os.environ.copy()
    for key,value in args.env: env[key]=value
    if rank is not None: env["CFDX_N15_MPI_RANKS"]=str(rank)
    started=dt.datetime.now(dt.timezone.utc); t0=time.perf_counter()
    completed=subprocess.run(command,cwd=args.cwd,env=env,text=True,capture_output=not args.inherit_output)
    elapsed=time.perf_counter()-t0; finished=dt.datetime.now(dt.timezone.utc)
    metrics=load_json(Path(args.metrics_file)) if args.metrics_file else {}
    metrics.setdefault("total",elapsed); metrics.setdefault("peak_rss_bytes",peak_rss_bytes())
    return {"schema_version":SCHEMA_VERSION,"benchmark":case,
            "provenance":{"git_sha":git_sha(Path(args.repo)),"timestamp_start_utc":started.isoformat(),"timestamp_end_utc":finished.isoformat(),"command":command,"cwd":str(Path(args.cwd).resolve()),"returncode":completed.returncode},
            "software":{"python":platform.python_version(),"compiler":compiler_info(),"build_type":args.build_type,"mpi":args.mpi,"cuda":args.cuda},
            "hardware":hardware_metadata(),"mesh":case.get("mesh",{}),"numerics":case.get("numerics",{}),"performance":metrics,
            "status":"PASS" if completed.returncode==0 else "FAIL",
            "stdout":completed.stdout if not args.inherit_output else None,"stderr":completed.stderr if not args.inherit_output else None}

def write_json(data:dict[str,Any],path:Path)->None:
    path.parent.mkdir(parents=True,exist_ok=True); path.write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")

def scaling(args:argparse.Namespace)->int:
    case=load_json(Path(args.case)); values=[int(x) for x in args.values.split(",") if x]; runs=[]
    for value in values:
        run_case=json.loads(json.dumps(case)); run_case.setdefault("parallel",{})["ranks"]=value
        if args.kind=="weak":
            cpr=int(run_case.get("mesh",{}).get("cells_per_rank",0))
            if cpr<=0: raise ValueError("weak scaling requires mesh.cells_per_rank")
            run_case.setdefault("mesh",{})["cells"]=cpr*value
        runs.append(run_once(args,run_case,value))
        if runs[-1]["status"]!="PASS": break
    write_json({"schema_version":SCHEMA_VERSION,"campaign":{"kind":args.kind,"values":values},"runs":runs},Path(args.output))
    return 0 if all(x["status"]=="PASS" for x in runs) else 1

def analyze_scaling(data:dict[str,Any])->dict[str,Any]:
    runs=data["runs"]
    if not runs: raise ValueError("campaign contains no runs")
    t0=float(runs[0]["performance"]["total"]); p0=int(runs[0]["benchmark"].get("parallel",{}).get("ranks",1)); rows=[]
    for run in runs:
        p=int(run["benchmark"].get("parallel",{}).get("ranks",1)); t=float(run["performance"]["total"]); s=t0/t if t>0 else None
        rows.append({"ranks":p,"cells":run.get("mesh",{}).get("cells"),"time_s":t,"speedup":s,"efficiency":s/(p/p0) if t>0 and p else None,"peak_rss_bytes":run["performance"].get("peak_rss_bytes")})
    return {"kind":data["campaign"]["kind"],"rows":rows}

def regression(args:argparse.Namespace)->int:
    b=float(load_json(Path(args.baseline))["performance"]["total"]); candidate=load_json(Path(args.candidate)); c=float(candidate["performance"]["total"]); base=load_json(Path(args.baseline))
    change=(c-b)/b if b else float("inf")
    result={"baseline_s":b,"candidate_s":c,"relative_change":change,"relative_change_percent":100*change,"threshold_percent":args.threshold_percent,"regression":change>args.threshold_percent/100,"baseline_git_sha":base.get("provenance",{}).get("git_sha"),"candidate_git_sha":candidate.get("provenance",{}).get("git_sha")}
    write_json(result,Path(args.output)); print(json.dumps(result,indent=2)); return 1 if result["regression"] else 0

def make_parser()->argparse.ArgumentParser:
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest="mode",required=True); common=argparse.ArgumentParser(add_help=False)
    common.add_argument("--case",required=True); common.add_argument("--repo",default="."); common.add_argument("--cwd",default="."); common.add_argument("--build-type",default="Release"); common.add_argument("--mpi",default="unknown"); common.add_argument("--cuda",default="unknown"); common.add_argument("--metrics-file"); common.add_argument("--inherit-output",action="store_true"); common.add_argument("--env",action="append",nargs=2,default=[]); common.add_argument("command",nargs=argparse.REMAINDER)
    run=sub.add_parser("run",parents=[common]); run.add_argument("--output",required=True)
    scale=sub.add_parser("scale",parents=[common]); scale.add_argument("--kind",choices=("strong","weak"),required=True); scale.add_argument("--values",required=True); scale.add_argument("--output",required=True)
    analyze=sub.add_parser("analyze-scaling"); analyze.add_argument("--input",required=True); analyze.add_argument("--output",required=True)
    reg=sub.add_parser("regression"); reg.add_argument("--baseline",required=True); reg.add_argument("--candidate",required=True); reg.add_argument("--threshold-percent",type=float,default=10.0); reg.add_argument("--output",required=True)
    return p

def main()->int:
    args=make_parser().parse_args()
    if getattr(args, "command", None) and args.command[0] == "--":
        args.command = args.command[1:]
    if args.mode=="run":
        if not args.command: raise SystemExit("run requires a command after '--'")
        result=run_once(args,load_json(Path(args.case))); write_json(result,Path(args.output)); return 0 if result["status"]=="PASS" else 1
    if args.mode=="scale":
        if not args.command: raise SystemExit("scale requires a command after '--'")
        return scaling(args)
    if args.mode=="analyze-scaling":
        write_json(analyze_scaling(load_json(Path(args.input))),Path(args.output)); return 0
    return regression(args)
if __name__=="__main__": raise SystemExit(main())
