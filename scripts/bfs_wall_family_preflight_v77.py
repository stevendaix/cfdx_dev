@@ -0,0 +1,98 @@
#!/usr/bin/env python3
"""CFDX BFS wall-family preflight: descriptive, fail-closed; NOT a wall-law validation."""
import argparse
import json
import re
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_yplus_counts(logfile):
    text=Path(logfile).read_text(errors='replace')
    matches=re.findall(r'^BFS_WALL_REGIME_V74\s+([^\n]*)$',text,re.M)
    if not matches:raise ValueError('missing native BFS_WALL_REGIME_V74 evidence')
    d={}
    for token in matches[-1].split():
        if '=' in token:
            k,v=token.split('=',1)
            if k in d:raise ValueError('duplicate wall-key '+k)
            d[k]=v
    for k in ['faces','wall_resolved','buffer','log_layer','out_of_range',
              'min_yplus_spalding','max_yplus_spalding']:
        if k not in d:raise ValueError('incomplete yplus evidence: '+k)
    vals={k:int(d[k]) for k in ['faces','wall_resolved','buffer','log_layer','out_of_range']}
    vals['min_yplus']=float(d['min_yplus_spalding'])
    vals['max_yplus']=float(d['max_yplus_spalding'])
    if not np.isfinite([vals['min_yplus'],vals['max_yplus']]).all():
        raise ValueError('nonfinite wall yplus')
    if vals['faces']!=sum(vals[x] for x in ['wall_resolved','buffer','log_layer','out_of_range']):
        raise ValueError('wall yplus classification does not close')
    vals['face_fraction_out_of_range']=vals['out_of_range']/vals['faces']
    return vals


def wall_distances(h5path):
    with h5py.File(h5path) as f:
        p=f['points'][:];fo=f['face_offsets'][:]; fv=f['face_vertices'][:]
        offs=f['patch_face_offsets'][:];fids=f['patch_face_ids'][:]
        owners=f['owner'][:];co=f['cell_offsets'][:];cf=f['cell_faces'][:]
        patch=[s.split(':')[0] for s in f.attrs['boundary_patches'].split(';')]
        def center(fid):return p[fv[fo[fid]:fo[fid+1]]].mean(axis=0)
        def cell_center(cell):
            vids=set()
            for k in cf[co[cell]:co[cell+1]]:
                vids.update(map(int,fv[fo[k]:fo[k+1]]))
            return p[sorted(vids)].mean(axis=0)
        lookup={name:fids[offs[i]:offs[i+1]] for i,name in enumerate(patch)}
        h=abs(center(int(lookup['upstream_wall'][0]))[1]-center(int(lookup['downstream_wall'][0]))[1])
        if h<=0:raise ValueError('bad step height')
        samples=[]
        for name in ['upstream_wall','step_wall','downstream_wall','top']:
            for fid in lookup[name]:
                fid=int(fid);c=cell_center(int(owners[fid]));fcenter=center(fid)
                vert=p[fv[fo[fid]:fo[fid+1]]]
                normal=np.cross(vert[1]-vert[0],vert[2]-vert[0]);length=np.linalg.norm(normal)
                if not length>0:raise ValueError('degenerate wall face')
                distance=abs(float(np.dot(c-fcenter,normal/length)))
                if distance<=0:raise ValueError('degenerate owner-to-wall distance')
                samples.append((name, fid, float(distance/h)))
        return h,samples


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--case-h5',required=True)
    ap.add_argument('--log',required=True)
    ap.add_argument('--output-dir',required=True)
    args=ap.parse_args()
    y=read_yplus_counts(args.log)
    h,samples=wall_distances(args.case_h5)
    out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
    groups={name:[v for p,f,v in samples if p==name] for name in sorted(set(p for p,_,_ in samples))}
    fig,axes=plt.subplots(1,2,figsize=(12,4.6))
    names=['y+<=1','1<y+<30','30<=y+<=300','y+>300']
    axes[0].bar(names,[y[x] for x in ['wall_resolved','buffer','log_layer','out_of_range']])
    axes[0].set(ylabel='Nombre de faces',title='BFS152 : V74 y+ selon Spalding, non converge')
    for i,n in enumerate([y[x] for x in ['wall_resolved','buffer','log_layer','out_of_range']]):
        axes[0].text(i,n+0.3,str(n),ha='center')
    axes[0].set_ylim(0,max(27,y['out_of_range']+3));axes[0].tick_params(axis='x',rotation=18)
    axes[1].boxplot([groups[n] for n in groups],tick_labels=list(groups))
    axes[1].set(ylabel='Distance géométrique centre-face / H',
        title='BFS152 : distance 1re cellule, indépendante de y+')
    axes[1].tick_params(axis='x',rotation=15)
    fig.tight_layout();fig.savefig(out/'08_cfdx_wall_family_preflight.png',dpi=180);plt.close(fig)
    report={'schema':'CFDX-WALL-FAMILY-PREFLIGHT/1','mesh':'BFS152 pilot',
        'nfaces_wall':len(samples),'step_height_m':h,'yplus_native_report':y,
        'geometric_first_cell_d_over_h_by_patch':{n:{'min':min(v),'median':float(np.median(v)),
           'max':max(v),'n':len(v)} for n,v in groups.items()},
        'wall_resolved_qualified':False,'wall_function_qualified':False,
        'reason':'No y+<=1 faces, 22/34 faces above y+=300; final converged field and matched wall laws unavailable',
        'remark':'low shear y+ near reattachment does not certify geometric nearwall resolution'}
    (out/'wall_family_preflight_v77.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
