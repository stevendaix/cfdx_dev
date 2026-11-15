"""Render gauge-invariant and per-cell step changes from REAL committed CFDX VTUs.

No vtkCellCenters: polyhedral coordinates are computed from the HDF5 mesh
vertices, matching the geometry's cell indexing. Never plots uncommitted
nonlinear trials as a physical solution. V87 can be supplied as a baseline.
"""
import argparse, json, hashlib
from pathlib import Path
import h5py
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def centers(path):
    with h5py.File(path) as f:
        pts=f['points'][...]; offsets=f['cell_offsets'][...]
        cf=f['cell_faces'][...]; fo=f['face_offsets'][...]; fv=f['face_vertices'][...]
        if len(offsets)!=20541: raise ValueError('Requires real BFS20540 mesh')
        out=np.empty((20540,3))
        for c in range(20540):
            faces=cf[offsets[c]:offsets[c+1]]
            verts=np.unique(np.concatenate([fv[fo[face]:fo[face+1]] for face in faces]))
            out[c]=pts[verts].mean(axis=0)
    return out

def read(path):
    rdr=vtk.vtkXMLUnstructuredGridReader();rdr.SetFileName(str(path));rdr.Update()
    grid=rdr.GetOutput()
    if grid.GetNumberOfCells()!=20540:raise ValueError('Wrong VTU cell count')
    required=('u_x','u_y','u_z','p','k','omega')
    r={}
    for key in required:
        array=grid.GetCellData().GetArray(key)
        if array is None:raise ValueError('Missing native field '+key)
        r[key]=vtk_to_numpy(array)
        if not np.all(np.isfinite(r[key])):raise ValueError('Nonfinite CFDX field')
    return r

def step(before,after,inlet=44.2,rho=1.2):
    delta=np.linalg.norm(np.column_stack([after['u_'+k]-before['u_'+k] for k in 'xyz']),axis=1)
    dp=after['p']-before['p']
    q=0.5*rho*inlet**2
    return {'delta_u':delta, 'dp':dp, 'q':q,
            'max_delta_u_over_uin':float(delta.max()/inlet),
            'worst_delta_u_cell':int(delta.argmax()),
            'pressure_increment_span_over_qin':float(np.ptp(dp)/q),
            'max_k':float(after['k'].max()),
            'negative_ux_count':int((after['u_x']<0).sum()),
            'pressure_span_over_qin':float(np.ptp(after['p'])/q)}

def summarize(h5,series,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    ctr=centers(h5); x=ctr[:,0]/0.0127;y=ctr[:,1]/0.0127
    entries=[]
    for label,directory in series:
        directory=Path(directory)
        states=sorted(directory.glob('result_*.vtu'),key=lambda p:int(p.stem.split('_')[-1]))
        if not states:continue
        values={int(file.stem.split('_')[-1]):read(file) for file in states}
        for prev,it in zip(sorted(values)[:-1],sorted(values)[1:]):
            if it!=prev+1:raise ValueError('VTU states must be contiguous; no fictitious interpolation')
            d=step(values[prev],values[it]);entries.append((label,prev,it,d))
    if not entries:raise RuntimeError('Need at least two actually exported accepted states')
    # One scale shared across all candidate methods and across all states.
    max_du=max(float(e[3]['delta_u'].max()/44.2) for e in entries)
    max_dp=max(float(np.max(np.abs(e[3]['dp']))/e[3]['q']) for e in entries)
    rows=[]
    for label,prev,it,d in entries:
        row={key:value for key,value in d.items() if key not in ('delta_u','dp')}
        row.update({'version':label,'previous_accepted':prev,'accepted_iteration':it})
        rows.append(row)
        mask=(x>=-1)&(x<=8)&(y>=0)&(y<=2)
        fig,axs=plt.subplots(2,1,figsize=(13,7),layout='constrained')
        a=d['delta_u']/44.2
        c1=axs[0].scatter(x[mask],y[mask],c=a[mask],s=10,cmap='inferno',vmin=0,vmax=max_du,edgecolor='none')
        fig.colorbar(c1,ax=axs[0],label=r'$|\Delta U|/U_{in}$')
        b=d['dp']/d['q']
        c2=axs[1].scatter(x[mask],y[mask],c=b[mask],s=10,cmap='coolwarm',vmin=-max_dp,vmax=max_dp,edgecolor='none')
        fig.colorbar(c2,ax=axs[1],label=r'$\Delta p/q_{in}$ (gauge increment only)')
        axs[0].set_title(f'{label} accepted {prev}->{it}: max |dU|/Uin={d["max_delta_u_over_uin"]:.3f}, worst cell {d["worst_delta_u_cell"]}')
        axs[1].set_title(f'{label}: pressure increment span / q = {d["pressure_increment_span_over_qin"]:.2f} (NOT Cp)')
        for ax in axs:ax.set_xlim(-1,8);ax.set_ylim(0,2);ax.set_xlabel('x/H');ax.set_ylabel('y/H')
        fig.suptitle('CFDX BFS20540 - real committed fields, NOT CONVERGED',fontsize=15)
        fig.savefig(out/f'local_step_{label}_{prev}_{it}.png',dpi=160)
        plt.close(fig)
    (out/'local_step_metrics.json').write_text(json.dumps(rows,indent=2))
    return rows
if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--h5',type=Path,required=True)
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--trial',type=Path)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    series=[('V87',a.baseline)]
    if a.trial:series.append(('V88',a.trial))
    print(json.dumps(summarize(a.h5,series,a.output),indent=2))