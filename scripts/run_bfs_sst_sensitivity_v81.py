@@ -0,0 +1,104 @@
#!/usr/bin/env python3
"""Bounded, first-private-SST-cycle sensitivity runner for an existing CFDX case.

No change to convergence tolerances or physical models. A truncated run is NOT
certification. Designed to read only one rejected private SST candidate.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import selectors
import subprocess
import time


def hexdigest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def change_case(source, destination, iterations, molecular, newton, omega_first):
    import h5py
    if not 1 <= iterations <= 800:
        raise ValueError('diagnostic SST budget must be in [1,800]')
    with h5py.File(source,'r') as old,h5py.File(destination,'w') as new:
        for key in old.keys():old.copy(key,new)
        for key,value in old.attrs.items():new.attrs[key]=value
        settings=json.loads(new.attrs['case_setup_json'])
        raw=settings['numerics']['raw_settings']
        raw.update({'sst_transport_max_iterations':str(iterations),
                    'sst_molecular_wall_diffusion':str(molecular).lower(),
                    'sst_omega_newton_sink':str(newton).lower(),
                    'sst_sequential_omega_first':str(omega_first).lower()})
        # Do not alter pressure, momentum, SST acceptance thresholds, wall
        # omega factor, Rhie-Chow, or source physics coefficients.
        new.attrs['case_setup_json']=json.dumps(settings)


def launch(args):
    root=Path(args.output_dir);root.mkdir(parents=True,exist_ok=True)
    source=Path(args.case);target=root/'sensitivity.cfdx.h5'
    change_case(source,target,args.iterations,args.molecular,args.omega_newton,args.omega_first)
    env=os.environ.copy()
    env.update({'CFDX_V80_NATIVE_FACE_DEFECT':'1',
                'CFDX_SST_DIAGNOSTICS':'1',
                'CFDX_BFS_STAGE_TRACE':'1',
                'OMP_NUM_THREADS':'2'})
    cmd=[str(Path(args.solver).resolve()),'--mesh',str(target.resolve()),
         '--experimental-sst','--iterations','1','--output-dir',str(root/'outputs')]
    start=time.monotonic();observed_end=False;timeout=False
    with (root/'native.stdout').open('wb') as stdout, (root/'native.stderr').open('wb') as stderr:
        child=subprocess.Popen(cmd,stdout=stdout,stderr=subprocess.PIPE,env=env)
        selector=selectors.DefaultSelector()
        selector.register(child.stderr,selectors.EVENT_READ)
        buffer=b''
        try:
            while time.monotonic()-start<args.seconds:
                events=selector.select(timeout=.25)
                if not events and child.poll() is not None:break
                for key,_ in events:
                    chunk=os.read(key.fd,32768)
                    if not chunk:
                        selector.unregister(key.fileobj)
                        break
                    stderr.write(chunk);stderr.flush();buffer+=chunk
                    while b'\n' in buffer:
                        line,buffer=buffer.split(b'\n',1)
                        if line.startswith(b'BFS_PHASE_V10') and b'phase=SST_TRANSPORT_END' in line:
                            observed_end=True;break
                    if observed_end:break
                if observed_end:break
            timeout=not observed_end
        finally:
            if child.poll() is None:
                child.terminate()
                try:child.wait(timeout=2.)
                except subprocess.TimeoutExpired:child.kill();child.wait(timeout=5.)
            selector.close();child.stderr.close()
    verdict={'provenance':'native CFDX private SST cycle; no accepted flow state',
             'solver_sha256':hexdigest(args.solver),'source_h5_sha256':hexdigest(source),
             'variant_h5_sha256':hexdigest(target),'variant':{'molecular':args.molecular,
             'omega_newton_sink':args.omega_newton,'omega_first':args.omega_first,
             'sst_budget':args.iterations},'full_SST_TRANSPORT_END_seen':observed_end,
             'timeout_or_incomplete':timeout,'outer_accepted':False,
             'elapsed_seconds':time.monotonic()-start,
             'native_returncode_after_stopping':child.returncode}
    (root/'verdict.json').write_text(json.dumps(verdict,indent=2))
    print(json.dumps(verdict,indent=2))
    if not observed_end:raise SystemExit(2)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--case',required=True);p.add_argument('--solver',required=True)
    p.add_argument('--output-dir',required=True);p.add_argument('--iterations',type=int,default=40)
    p.add_argument('--seconds',type=int,default=90)
    p.add_argument('--molecular',action='store_true')
    p.add_argument('--omega-newton',action='store_true')
    p.add_argument('--omega-first',action='store_true')
    launch(p.parse_args())
if __name__=='__main__':main()
