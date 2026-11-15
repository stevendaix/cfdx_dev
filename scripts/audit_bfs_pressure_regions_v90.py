"""Independent, fail-closed V90 pressure/spatial audit on accepted native VTU.

This read-only tool never certifies convergence.  It requires an explicit native
GLOBAL_PROGRESS_COMMITTED_V85 marker for every plotted provisional state; the
HDF5 inlet velocity, NOT 'ref_velocity', defines the pressure-dynamic scale.
"""
from __future__ import annotations
import argparse, csv, json, math, re
from pathlib import Path
import h5py
import numpy as np

H = 0.0127

def inlet_physical_data(case):
    with h5py.File(case, 'r') as f:
        js = json.loads(f.attrs['case_setup_json'])
    inlets = [b for b in js['boundary_conditions'] if b['patch_name']=='inlet']
    if len(inlets)!=1 or inlets[0]['type']!='inlet':
        raise ValueError('A unique fixed-velocity inlet must be present')
    v = np.asarray(inlets[0]['velocity_vector'],dtype=float)
    mat = js['materials']
    if v.shape!=(3,) or len(mat)!=1:raise ValueError('invalid inlet velocity/materials')
    rho = float(mat[0]['density']);speed=float(np.linalg.norm(v))
    if not np.isfinite(v).all() or not math.isfinite(rho) or speed<=0 or rho<=0:
        raise ValueError('nonfinite physical inlet velocity or density')
    return rho, speed, 0.5*rho*speed**2

def accepted_iterations(log):
    rows=[]
    for line in Path(log).read_text(errors='replace').splitlines():
        if not line.startswith('GLOBAL_PROGRESS_COMMITTED_V85 '):continue
        fields=dict(re.findall(r'([A-Za-z_]+)=([^\s]+)',line))
        if fields.get('converged')!='0' or fields.get('qualified')!='0':
            raise ValueError('V90 reads provisional unqualified native states only')
        n=int(fields['committed'])
        if n<=0 or (rows and rows[-1]>=n):raise ValueError('nonincreasing accepted state indices')
        rows.append(n)
    if not rows:raise ValueError('no native accepted state found')
    return rows

def centres_from_native_hdf5(path):
    with h5py.File(path,'r') as f:
        pts=f['points'][...];co=f['cell_offsets'][...];cf=f['cell_faces'][...]
        fo=f['face_offsets'][...];fv=f['face_vertices'][...]
    n=len(co)-1
    if n!=20540:raise ValueError('V90 pressure region audit expects 20540-cell BFS')
    c=np.empty((n,3),dtype=np.float64)
    for i in range(n):
        faces=cf[co[i]:co[i+1]]
        vertices=np.unique(np.concatenate([fv[fo[j]:fo[j+1]] for j in faces]))
        # Reconstruction by mean vertex position is valid for hexahedra used
        # on this case, NOT for general skew polyhedra.
        if len(vertices)!=8:raise ValueError('non-hexahedral geometry: no vertex-centre shortcut')
        c[i]=pts[vertices].mean(axis=0)
    if not np.isfinite(c).all():raise ValueError('nonfinite cell centre')
    return c

def load_vtu(path,n):
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    r=vtk.vtkXMLUnstructuredGridReader();r.SetFileName(str(path));r.Update()
    grid=r.GetOutput()
    if grid.GetNumberOfCells()!=n:raise ValueError('VTU cell count/mesh mismatch')
    arr={}
    for key in ('u_x','u_y','p','k','omega'):
        raw=grid.GetCellData().GetArray(key)
        if raw is None:raise ValueError(f'Missing CFDX cell field {key}')
        v=vtk_to_numpy(raw)
        if v.ndim!=1 or len(v)!=n or not np.isfinite(v).all():raise ValueError('invalid CFDX scalar field '+key)
        arr[key]=v
    return arr

def make_audit(case,log,states,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    rho,u,q=inlet_physical_data(case)
    centre=centres_from_native_hdf5(case);x=centre[:,0]/H;y=centre[:,1]/H
    regions={'all':np.ones(len(x),dtype=bool),'upstream':x<-1,
      'step':(x>0)&(x<20)&(y<3),'outlet':x>40,
      'upper_outlet':(x>40)&(y>8.9)}
    if not all(mask.any() for mask in regions.values()):raise ValueError('incomplete BFS regions')
    rows=[]
    accepted=accepted_iterations(log)
    for i in accepted:
        file=Path(states)/f'result_{i}.vtu'
        if not file.exists():raise FileNotFoundError(f'Native committed state {i} has no VTU')
        a=load_vtu(file,len(x));p=a['p'];j=int(p.argmin());k=int(p.argmax())
        entry={'iteration':i,'p_span_over_q':float(np.ptp(p)/q),
         'p_min_pa':float(p[j]),'p_min_xH':float(x[j]),'p_min_yH':float(y[j]),
         'p_max_pa':float(p[k]),'p_max_xH':float(x[k]),'p_max_yH':float(y[k]),
         'ux_min_ms':float(a['u_x'].min()),'k_max_m2s2':float(a['k'].max()),
         'reverse_ux_cells':int(np.count_nonzero(a['u_x']<0)),
         'p_cell20539_pa':float(p[20539]),'uy_cell20539_ms':float(a['u_y'][20539]),
         'p_cell7406_pa':float(p[7406]), 'k_cell7406_m2s2':float(a['k'][7406]),
         'qualification':'PROVISIONAL_NOT_CONVERGED'}
        for name,mask in regions.items():entry[f'{name}_span_pa']=float(np.ptp(p[mask]))
        rows.append(entry)
    with (out/'regions_v90.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'audit_v90.json').write_text(json.dumps({
      'case':str(case),'log':str(log),'rho':rho,'physical_uin_ms':u,'q_in_pa':q,
      'source':'native HDF5 + native COMMITTED stderr + native per-cell VTU',
      'certified':False,'states':rows},indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(9,4.6))
    t=np.array([r['iteration'] for r in rows])
    for name in ('all','upstream','step','outlet'):
        ax.plot(t,[r[f'{name}_span_pa']/q for r in rows],marker='o',label=name)
    ax.set(xlabel='Accepted external state (provisional)',ylabel='Pressure span / q_in',
           title='CFDX BFS20540: pressure span by physical region (NOT converged)')
    ax.grid(alpha=.3);ax.legend();fig.tight_layout();fig.savefig(out/'pressure_region_drift_v90.png',dpi=150)
    plt.close(fig)
    return rows

def main():
    p=argparse.ArgumentParser();p.add_argument('--case',required=True,type=Path)
    p.add_argument('--log',required=True,type=Path);p.add_argument('--states',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    a=p.parse_args();rows=make_audit(a.case,a.log,a.states,a.output)
    print('ACCEPTED_PROVISIONAL',len(rows),'max_iteration',rows[-1]['iteration'],'CERTIFIED',False)
if __name__=='__main__':main()