"""Make an opt-in BFS20540 V83 experiment without changing mesh or physics.

This does NOT enable final SST convergence or physical certification. It only
allows *bounded provisional* SST states to enter the full nonlinear residual
progress ledger V8/V13/V14; final thresholds and mesh/operator are unchanged.
"""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

try:
    import h5py
except ImportError as exc:
    raise SystemExit("h5py is required for CFDX HDF5 cases") from exc

REQUIRED={
    "v7_post_sst_momentum_gate": "true",
    "v8_global_nonlinear_residual": "true",
    "v13_joint_globalization": "true",
    "v14_monotone_global_merit": "true",
    "pressure_high_low_production": "true",
}

def make_progress_case(source: Path, destination: Path) -> dict:
    source, destination = Path(source), Path(destination)
    if source.resolve()==destination.resolve():
        raise ValueError("V83 never overwrites the reference case")
    if not source.is_file():
        raise FileNotFoundError(source)
    with h5py.File(source,"r") as h5:
        n=int(h5.attrs["n_cells"])
        if n!=20540 or str(h5.attrs.get("source_case_name",""))!="backwardFacingStep2D":
            raise ValueError("V83 experiment requires the actual six-block BFS20540")
        setup=json.loads(h5.attrs["case_setup_json"])
        numerics=setup["numerics"]
        if numerics["coupled_solver"]!="SIMPLE" or setup.get("transient",False):
            raise ValueError("V83 supports steady SIMPLE only")
        raw=numerics["raw_settings"]
        for key,value in REQUIRED.items():
            if raw.get(key)!=value:
                raise ValueError(f"V83 requires {key}={value}")
        if float(raw.get("pressure_linear_tolerance","1"))>1e-11:
            raise ValueError("V83 must not relax 1e-11 pressure gate")
        if int(raw.get("sst_transport_max_iterations","0"))<10:
            raise ValueError("V83 must not shorten existing SST diagnostic budget")
    digest=hashlib.sha256(source.read_bytes()).hexdigest()
    destination.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(source,destination)
    with h5py.File(destination,"r+") as h5:
        state=json.loads(h5.attrs["case_setup_json"])
        cfg=state["numerics"]["raw_settings"]
        cfg["sst_outer_provisional"]="enabled"
        cfg["v83_provisional_outer_merit_progress"]="enabled"
        h5.attrs.modify("case_setup_json",json.dumps(state,separators=(",", ":")))
        # The only changed metadata is the serialized opt-in settings. All
        # mesh datasets and other attributes remain exactly from the source.
    with h5py.File(source,"r") as a,h5py.File(destination,"r") as b:
        assert a.attrs["mesh_hash"]==b.attrs["mesh_hash"]
        assert a.attrs["topology_hash"]==b.attrs["topology_hash"]
        assert tuple(a.keys())==tuple(b.keys())
        for name in a:
            if hasattr(a[name],"shape") and not (a[name][...]==b[name][...]).all():
                raise RuntimeError("V83 changed a mesh dataset: "+name)
    return {"source":str(source),"output":str(destination),
        "source_sha256":digest,"cells":n,"physics_qualified":False,
        "final_sst_tolerance":1e-8,"final_pressure_tolerance":1e-11,
        "changed_only": ["numerics.raw_settings.sst_outer_provisional",
                          "numerics.raw_settings.v83_provisional_outer_merit_progress"]}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("source",type=Path);p.add_argument("destination",type=Path)
    args=p.parse_args();print(json.dumps(make_progress_case(args.source,args.destination),indent=2))
if __name__=="__main__":main()