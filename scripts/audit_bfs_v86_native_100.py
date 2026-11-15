#!/usr/bin/env python3
"""Audit independent of a real CFDX V85 100-outer campaign.

Evidence requirements: committed native U/p/SST transaction AND per-state field
marker, case-specific HDF5 provenance, and VTU consistency when a snapshot exists.
No convergence or NASA validation is inferred from a provisional state.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re
from pathlib import Path
import h5py
import numpy as np

COMMIT='GLOBAL_PROGRESS_COMMITTED_V85 '
FIELD='BFS_FIELD_COMMITTED_V85 '
TRIAL='GLOBAL_PROGRESS_TRIAL_V85 '
KEYS=('U','U_global','p_continuity','k','omega','merit')
FIELD_KEYS=('reverse_ux_cells','ux_min','ux_max','pressure_min_pa','pressure_max_pa','pressure_span_pa')

def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def kv(line: str) -> dict[str,str]:
    return dict(re.findall(r'([A-Za-z_][A-Za-z_0-9]*)=([^\s]+)',line))

def read_evidence(text: str, max_outer:int=100):
    commit={}; fields={}; trial=[]
    for line in text.splitlines():
        if line.startswith(COMMIT):
            r=kv(line);i=int(r['iteration'])
            if i in commit or not 1<=i<=max_outer:raise ValueError('Duplicate/out-of-range committed iteration')
            if r.get('qualified')!='0' or r.get('converged')!='0':
                raise ValueError('Provisional native state falsely claims final convergence')
            if int(r.get('committed','-1'))!=len(commit)+1:
                raise ValueError('Commit ledger count inconsistent with native chronology')
            for k in KEYS:
                v=float(r[k]);
                if not math.isfinite(v) or v<0:raise ValueError('Bad native physical residual '+k)
            commit[i]=r
        elif line.startswith(FIELD):
            r=kv(line);i=int(r['iteration'])
            if i in fields or not 1<=i<=max_outer:raise ValueError('Duplicate/out-of-range native field state')
            if int(r['cells'])!=20540:raise ValueError('Field not from BFS20540')
            for key in FIELD_KEYS:
                val=float(r[key]);
                if not math.isfinite(val):raise ValueError('Bad field diagnostic '+key)
            if int(r['reverse_ux_cells'])<0 or int(r['reverse_ux_cells'])>20540:
                raise ValueError('Invalid negative-velocity-cell count')
            if 'ux_min_cell' in r and not 0<=int(r['ux_min_cell'])<20540:
                raise ValueError('Invalid worst-velocity cell index')
            if 'speed_max' in r and (not math.isfinite(float(r['speed_max'])) or float(r['speed_max'])<0):
                raise ValueError('Invalid maximum speed')
            if float(r['pressure_span_pa'])<0:raise ValueError('Negative pressure span')
            if not math.isclose(float(r['pressure_max_pa'])-float(r['pressure_min_pa']),
                                float(r['pressure_span_pa']),rel_tol=1e-8,abs_tol=1e-5):
                raise ValueError('Pressure field statistics inconsistent')
            fields[i]=r
        elif line.startswith(TRIAL):trial.append(kv(line))
    if list(commit)!=sorted(commit):raise ValueError('Nonchronological commit stream')
    if set(commit)!=set(fields):raise ValueError('A committed state lacks its matching field evidence')
    for i in commit:
        if not any(int(t.get('iteration','-1'))==i and t.get('accepted')=='1' for t in trial):
            raise ValueError('Native committed state lacks an accepted evaluated trial')
    return commit,fields,trial

def read_case(case:Path):
    with h5py.File(case,'r') as h:
        if int(h.attrs.get('n_cells',0))!=20540:raise ValueError('Not BFS20540')
        setup=json.loads(h.attrs['case_setup_json'])
        if int(setup['numerics']['max_iterations'])!=100:raise ValueError('Case not configured for 100 outer iterations')
        raw=setup['numerics']['raw_settings']
        if raw.get('v85_bounded_nonmonotone_outer_progress')!='enabled' or raw.get('sst_outer_provisional')!='enabled':
            raise ValueError('Not a V85 provisional outer case')
        provenance=str(h.attrs.get('source_case_path',''))
        # An HDF5 source_format='openfoam' alone is NOT native geometry proof.
        native_claim=provenance.startswith('OpenFOAM polyMesh, ')
        if native_claim and not h.attrs.get('boundary_patch_hash'):
            raise ValueError('Native geometry claimed but boundary ownership hash missing')
        rho=float(setup['materials'][0]['density'])
        inlet=[x for x in setup['boundary_conditions'] if x.get('type')=='inlet']
        if len(inlet)!=1:raise ValueError('Unique inlet boundary condition required for normalisation')
        vec=inlet[0].get('velocity_vector',[])
        if len(vec)!=3 or not all(math.isfinite(float(x)) for x in vec):
            raise ValueError('Invalid inlet velocity vector')
        uref=math.sqrt(sum(float(x)*float(x) for x in vec))
        if not (math.isfinite(rho) and math.isfinite(uref) and rho>0 and uref>0):
            raise ValueError('Invalid dynamic-pressure reference')
        return {'mesh_provenance':'NATIVE_OPENFOAM_REQUIRES_SEPARATE_REIMPORT_GATE' if native_claim
                else 'RECONSTRUCTED_SIX_BLOCK_GEOMETRY_ONLY',
                'source_case_path':provenance,'mesh_hash':str(h.attrs.get('mesh_hash','')),
                'rho_kg_m3':rho,'u_inlet_bc_m_s':uref,
                'case_setup_ref_velocity_m_s':float(setup['ref_velocity']),
                'q_reference_pa':0.5*rho*uref*uref,
                'physical_reference_qualified':False}

def verify_vtu(path:Path, field:dict):
    try:
        import vtk
        from vtk.util.numpy_support import vtk_to_numpy
    except ImportError as exc:raise RuntimeError('VTK required for independent VTU proof') from exc
    reader=vtk.vtkXMLUnstructuredGridReader();reader.SetFileName(str(path));reader.Update();grid=reader.GetOutput()
    if grid.GetNumberOfCells()!=20540:raise ValueError('VTU has wrong cell count')
    for name,lo,hi in [('u_x','ux_min','ux_max'),('p','pressure_min_pa','pressure_max_pa')]:
        a=grid.GetCellData().GetArray(name)
        if a is None:raise ValueError('Missing physical cell array '+name)
        vals=vtk_to_numpy(a)
        if not np.isfinite(vals).all():raise ValueError('Nonfinite VTU cell field '+name)
        for key,v in ((lo,float(vals.min())),(hi,float(vals.max()))):
            if not math.isclose(v,float(field[key]),rel_tol=2e-7,abs_tol=1e-5):
                raise ValueError(f'VTU/committed field mismatch {name}/{key}: {v} != {field[key]}')
        if name=='u_x' and int(np.count_nonzero(vals<0))!=int(field['reverse_ux_cells']):
            raise ValueError('VTU reverse cell count mismatches native field marker')
    return {'file':path.name,'sha256':digest(path),'verified_cells':20540}

def audit(log:Path,case:Path,out:Path,states:Path | None=None):
    case_info=read_case(case)
    commit,field,trials=read_evidence(log.read_text(errors='replace'))
    out.mkdir(parents=True,exist_ok=True)
    rows=[]; q=case_info['q_reference_pa']
    for i,r in commit.items():
        f=field[i]
        row={'iteration':i,**{k:float(r[k]) for k in KEYS},
             **{k:float(f[k]) for k in FIELD_KEYS},
             'k_debt':int(r['k_debt']),'omega_debt':int(r['omega_debt']),
             'p_span_over_q':float(f['pressure_span_pa'])/q,
             'u_x_min_over_u_inlet':float(f['ux_min'])/case_info['u_inlet_bc_m_s'],
             'has_reverse_flow':int(f['reverse_ux_cells'])>0,'provisional':True}
        if 'ux_min_cell' in f:
            row['ux_min_cell']=int(f['ux_min_cell'])
            row['speed_max_over_u_inlet']=float(f['speed_max'])/case_info['u_inlet_bc_m_s']
        rows.append(row)
    if rows and all('ux_min_cell' in row for row in rows):
        with h5py.File(case,'r') as h:
            points=h['points'][...];co=h['cell_offsets'][...];cf=h['cell_faces'][...]
            fo=h['face_offsets'][...];fv=h['face_vertices'][...]
            for row in rows:
                ci=row['ux_min_cell'];fs=cf[co[ci]:co[ci+1]]
                ids=np.unique(np.concatenate([fv[fo[j]:fo[j+1]] for j in fs]))
                centre=np.mean(points[ids],axis=0)
                row['ux_min_x_over_H']=float(centre[0]/.0127)
                row['ux_min_y_over_H']=float(centre[1]/.0127)
    snapshots=[]
    if states is not None:
        for path in sorted(states.glob('result_*.vtu')):
            try:i=int(path.stem.rsplit('_',1)[1])
            except (ValueError,IndexError):raise ValueError('Unexpected VTU filename')
            if i not in commit or not(i<=10 or i%10==0 or i==100):
                raise ValueError('VTU cannot be linked to a permitted committed state')
            snapshots.append(verify_vtu(path,field[i]))
    with (out/'native_progress_v86.csv').open('w',newline='') as fd:
        cols=list(rows[0]) if rows else ['iteration']
        writer=csv.DictWriter(fd,fieldnames=cols);writer.writeheader();writer.writerows(rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if rows:
        xx=[r['iteration'] for r in rows]
        for plotname,keys,title,yscale in [
          ('01_global_residues.png',('U','U_global','k','omega'),'CFDX natif / provisoire : residus U, k, omega','log'),
          ('02_qualite_physique.png',('p_span_over_q','reverse_ux_cells'),'Pression normalisee et cellules en contre-courant','linear')]:
            fig,ax=plt.subplots(figsize=(9,5))
            for k in keys:ax.plot(xx,[r[k] for r in rows],marker='o',label=k)
            if yscale=='log':ax.set_yscale('log')
            ax.set_xlabel('Iteration exterieure reellement acceptee')
            ax.set_title(title+' — NON CONVERGE');ax.grid(alpha=.3);ax.legend();fig.tight_layout()
            fig.savefig(out/plotname,dpi=155);plt.close(fig)
    result={'schema':'CFDX-V86-NATIVE-EVIDENCE/1','case':str(case),'case_sha256':digest(case),
       'stderr':str(log),'stderr_sha256':digest(log),'max_outer_requested':100,
       'committed_states':len(rows),'last_committed':rows[-1]['iteration'] if rows else 0,
       'target_100_reached':len(rows)>=100,'physically_qualified':False,
       'trials':len(trials),'rejected_trials':sum(t.get('accepted')=='0' for t in trials),
       'snapshots_independently_verified':snapshots,'case_metadata':case_info,
       'last_physical_metrics':rows[-1] if rows else None,
       'physics_warnings':(
           (["Pressure span > 10 inlet dynamic pressures; check gauge, walls and momentum/pressure coupling"]
            if rows and rows[-1]['p_span_over_q']>10 else [])+
           (["Reverse speed larger than inlet speed: NOT certified physical recirculation"]
            if rows and rows[-1]['u_x_min_over_u_inlet']<-1 else [])),
       'uncertainty':'Pressure span/kinetic pressure is a diagnostic only. No final CFDX Cf/Cp/U/uv validation. Geometry has no native blockMesh identity unless separately verified.'}
    (out/'evidence_v86.json').write_text(json.dumps(result,indent=2,ensure_ascii=False))
    return result

if __name__=='__main__':
    pa=argparse.ArgumentParser(description=__doc__)
    pa.add_argument('--native-log',type=Path,required=True);pa.add_argument('--case',type=Path,required=True)
    pa.add_argument('--states',type=Path);pa.add_argument('--output',type=Path,required=True)
    a=pa.parse_args();print(json.dumps(audit(a.native_log,a.case,a.output,a.states),indent=2,ensure_ascii=False))