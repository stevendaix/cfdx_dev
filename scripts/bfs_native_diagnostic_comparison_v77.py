@@ -0,0 +1,378 @@
#!/usr/bin/env python3
"""CFDX BFS152 V78 diagnostic-only comparison using ACTUAL accepted VTU fields.

Strict provenance: the native CFDX solver emits the VTU; these plots are NOT
converged Cf/Cp/uv, nor wall traction from CFDX's own wall operator.
"""
import argparse
import csv
import json
import math
import re
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from bfs_comparison_plots_v76 import CF_EXP, CF_NASA, CP_EXP_BOTTOM, CP_NASA
from bfs_geometry_comparability_v78 import read_case, validate_one


def require(cond, message):
    if not cond:
        raise ValueError(message)


def get_mesh_and_case(path):
    with h5py.File(path) as f:
        cfg = json.loads(f.attrs['case_setup_json'])
        coords = f['points'][:]
        faceoffset = f['face_offsets'][:]
        verts = f['face_vertices'][:]
        owners = f['owner'][:]
        poff = f['patch_face_offsets'][:]
        pfaces = f['patch_face_ids'][:]
        names = [x.split(':')[0] for x in f.attrs['boundary_patches'].split(';')]
        ncell = int(f.attrs['n_cells'])
        cell_faces=f['cell_faces'][:]
        cell_offsets=f['cell_offsets'][:]
        patch = {}
        for i, name in enumerate(names):
            patch[name] = pfaces[poff[i]:poff[i+1]]
        require({'downstream_wall', 'upstream_wall', 'step_wall', 'inlet', 'outlet'}.issubset(patch),
                'missing required wall and inlet/outlet patches')
        def center(fid):
            return coords[verts[faceoffset[fid]:faceoffset[fid+1]]].mean(axis=0)
        wall_lower = np.array([center(int(fid))[1] for fid in patch['downstream_wall']])
        wall_upper = np.array([center(int(fid))[1] for fid in patch['upstream_wall']])
        h = float(np.median(wall_upper)-np.median(wall_lower))
        require(h>0 and np.max(np.abs(wall_lower-np.median(wall_lower)))<1e-8,
                'invalid BFS step height/wall flatness')
        require(np.max(np.abs(wall_upper-np.median(wall_upper)))<1e-8,
                'inconsistent upstream wall height')
        mat=cfg['materials'][0]
        rho=float(mat['density']);mu=float(mat['dynamic_viscosity'])
        inlet=next(x for x in cfg['boundary_conditions'] if x['patch_name']=='inlet')
        u_ref=float(inlet['velocity_magnitude'])
        require(rho>0 and mu>0 and u_ref>0, 'invalid fluid normalization')
        return dict(config=cfg, points=coords, faceoffset=faceoffset, vertices=verts,
                    owners=owners, patch=patch, ncell=ncell, height=h,
                    cell_faces=cell_faces,cell_offsets=cell_offsets,
                    rho=rho, mu=mu, u_ref=u_ref)


def audited_cell_centers(grid,mesh):
    """Use original HDF5 cell connectivity, not vtkCellCenters for VTK_POLYHEDRON.

    This vertex-average is an EXACT centroid only for the rectangular, affine
    BFS152 hexahedra; reject any skewed or warped cell rather than silently
    changing the sampling scheme. The CFDX 20,540-cell case must use a
    verified FVM volume centroid before any field comparison is certified.
    """
    pts=mesh['points']
    fo=mesh['faceoffset'];fv=mesh['vertices']
    cf=mesh['cell_faces'];co=mesh['cell_offsets']
    require(grid.GetNumberOfCells()==mesh['ncell'], 'VTU/HDF5 cell count differs')
    centers=np.empty((mesh['ncell'],3),dtype=float)
    for i in range(mesh['ncell']):
        faces=cf[co[i]:co[i+1]]
        ids=np.unique(np.concatenate([fv[fo[f]:fo[f+1]] for f in faces]))
        vtkids=np.array([grid.GetCell(i).GetPointId(j)
                           for j in range(grid.GetCell(i).GetNumberOfPoints())],dtype=int)
        require(len(ids)==8 and np.array_equal(np.sort(vtkids),ids),
                'VTU/HDF5 cell vertex IDs inconsistent or not a hexahedron')
        p=pts[ids]
        for axis in range(3):
            require(len(np.unique(np.round(p[:,axis],decimals=12)))==2,
                    'non-rectilinear cell: do not use a vertex-average as volume centroid')
        centers[i]=p.mean(axis=0)
    require(np.isfinite(centers).all(),'nonfinite verified cell center')
    return centers


def read_native_vtu(path, mesh):
    reader=vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path));reader.Update()
    grid=reader.GetOutput()
    require(grid.GetNumberOfCells()==mesh['ncell'],
            'VTU cell count differs from original HDF5 mesh')
    require(np.allclose(mesh['points'],vtk_to_numpy(grid.GetPoints().GetData()),
                        rtol=0,atol=1e-12),'VTU/HDF5 point geometry differs')
    xyz=audited_cell_centers(grid,mesh)
    fields={}
    for name in ['u_x','u_y','p','k','omega','nu_t','wall_distance']:
        arr=grid.GetCellData().GetArray(name)
        require(arr is not None, 'missing native VTU field '+name)
        fields[name]=vtk_to_numpy(arr).copy()
        require(len(fields[name])==len(xyz) and np.isfinite(fields[name]).all(),
                'nonfinite or invalid native '+name)
    require((fields['k']>=0).all() and (fields['omega']>0).all() and
            (fields['nu_t']>=0).all(), 'inadmissible native SST state')
    return grid,xyz,fields


def wall_faces(mesh, xyz, fields, iteration):
    """One-point viscous Cf and owner pressure; DO NOT call native wall outputs."""
    row=[]
    q=.5*mesh['rho']*mesh['u_ref']**2
    for patch in ['upstream_wall', 'downstream_wall']:
        for fid in mesh['patch'][patch]:
            fid=int(fid);owner=int(mesh['owners'][fid])
            vind=mesh['vertices'][mesh['faceoffset'][fid]:mesh['faceoffset'][fid+1]]
            vert=mesh['points'][vind];fc=vert.mean(axis=0)
            # Both target wall patches are horizontal; normal distance is y.
            require(np.ptp(vert[:,1])<1e-9, 'one-point wall proxy requires flat y-wall')
            d=abs(float(xyz[owner,1]-fc[1]))
            require(d>0 and math.isfinite(d), 'invalid wall-normal distance')
            cx=float(fc[0]/mesh['height'])
            # rho cancels in Cf = mu du/dn / (0.5 rho Uref^2).
            cf_proxy=(mesh['mu']*float(fields['u_x'][owner])/d)/q
            row.append(dict(patch=patch,face=fid,owner=owner,x_over_h=cx,
                 y=float(fc[1]),u_x_owner=float(fields['u_x'][owner]),
                 p_owner=float(fields['p'][owner]),d_normal=d,
                 cf_one_point_viscous_proxy=cf_proxy,iteration=iteration))
    downstream=sorted((a for a in row if a['patch']=='downstream_wall'),
                      key=lambda a:a['x_over_h'])
    require(len(downstream)>=4, 'too few downstream wall faces')
    # Native p has arbitrary gauge; offset Cp from last downstream wall OWNER.
    # This is NOT an exact wall-face pressure nor an experimental U_ref match.
    p_far=downstream[-1]['p_owner']
    for x in row:
        x['cp_owner_proxy_shifted']=(x['p_owner']-p_far)/q
        x['pressure_reference_x_over_h']=downstream[-1]['x_over_h']
    return sorted(row,key=lambda a:a['x_over_h'])


def wls_uv_proxy(xyz, fields, u_ref, neighbours=16):
    """External 2-D WLS gradient of velocity; modeled uv, NOT a native flux."""
    xy=xyz[:,:2]
    tree=cKDTree(xy)
    neighbors=min(neighbours+1,len(xy))
    uv=np.full(len(xy),np.nan)
    condition=np.full(len(xy),np.nan)
    for i, x in enumerate(xy):
        dist,ids=tree.query(x,k=neighbors)
        dist=np.atleast_1d(dist);ids=np.atleast_1d(ids)
        selected=(dist>1e-12)
        r=dist[selected];ix=ids[selected]
        if len(ix)<4:continue
        delta=xy[ix]-x
        weight=1.0/np.maximum(r,1e-7)
        A=delta*weight[:,None]
        cond=float(np.linalg.cond(A))
        condition[i]=cond
        if not math.isfinite(cond) or cond>1e5:
            continue
        diff_ux=(fields['u_x'][ix]-fields['u_x'][i])*weight
        diff_uy=(fields['u_y'][ix]-fields['u_y'][i])*weight
        dx=np.linalg.lstsq(A,diff_ux,rcond=None)[0]
        dy=np.linalg.lstsq(A,diff_uy,rcond=None)[0]
        uv[i]=-fields['nu_t'][i]*(dx[1]+dy[0])/(u_ref*u_ref)*1000.
    return uv,condition


def interpolated_profiles(xyz, fields, h, u_ref, samples, uv_proxy):
    # Use ONLY downstream cell centers: do not interpolate through the solid step.
    mask=(xyz[:,0]>0) & (xyz[:,1]>0)
    coord=xyz[mask,:2]/h
    # Keep genuine convex-hull extrapolation as NaN: no fabricated near-wall layer.
    u_interp=LinearNDInterpolator(coord,fields['u_x'][mask]/u_ref,fill_value=np.nan)
    uv_interp=LinearNDInterpolator(coord,uv_proxy[mask],fill_value=np.nan)
    rows=[]
    for entry in samples:
        x=float(entry['xH']);y=float(entry['yH']);u=float(entry['uUref'])
        uv=float(entry['uv_times1000'])
        est_u=float(u_interp(np.array([[x,y]]))[0])
        est_uv=float(uv_interp(np.array([[x,y]]))[0])
        rows.append(dict(x_over_h=x,y_over_h=y,exp_u_u_ref=u,
             cfdx_u_interp=est_u if math.isfinite(est_u) else None,
             exp_uv_times1000=uv,
             cfdx_uv_wls_proxy=est_uv if math.isfinite(est_uv) else None))
    return rows


def read_experimental_samples(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path,rows,cols=None):
    if cols is None: cols=list(rows[0])
    with Path(path).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=cols)
        writer.writeheader();writer.writerows(rows)


def compare_metrics(wall,profile):
    d=[x for x in wall if x['patch']=='downstream_wall']
    cp_x=np.array([x['x_over_h'] for x in d])
    cp=np.array([x['cp_owner_proxy_shifted'] for x in d])
    mask=(cp_x>=CP_EXP_BOTTOM[0,0])&(cp_x<=CP_EXP_BOTTOM[-1,0])
    cp_rmse=float(np.sqrt(np.mean((cp[mask]-np.interp(cp_x[mask],CP_EXP_BOTTOM[:,0],CP_EXP_BOTTOM[:,1]))**2))) if mask.sum() else None
    result={'cp_owner_proxy_rmse_vs_exp_on_CFDX_wall_face_x':cp_rmse,
            'cp_faces_evaluated':int(mask.sum())}
    for station in (1,4,6,10):
        points=[r for r in profile if r['x_over_h']==station and r['cfdx_u_interp'] is not None]
        result[f'u_station_{station}']={
            'experimental_stations_interpolable':len(points),
            'rmse_u_u_ref':float(np.sqrt(np.mean([(r['cfdx_u_interp']-r['exp_u_u_ref'])**2 for r in points]))) if points else None}
    return result


def plots(out,mesh,grid,xyz,fields,wall,profile,iteration,mid_native=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    # Native cell map preserving polygons, not a smoothed/digitized contour.
    from matplotlib.collections import PolyCollection
    polys=[];val=[]
    for i in range(grid.GetNumberOfCells()):
        cell=grid.GetCell(i)
        vertices=np.array([grid.GetPoint(cell.GetPointId(j)) for j in range(cell.GetNumberOfPoints())])
        lower=vertices[np.abs(vertices[:,2]-vertices[:,2].min())<1e-9][:,:2]
        if len(lower)<3: continue
        center=lower.mean(axis=0)
        order=np.argsort(np.arctan2(lower[:,1]-center[1],lower[:,0]-center[0]))
        polys.append(lower[order]/mesh['height']);val.append(fields['u_x'][i]/mesh['u_ref'])
    fig,ax=plt.subplots(figsize=(12,4.7))
    obj=PolyCollection(polys,array=np.array(val),edgecolors='gray',linewidths=.15,cmap='coolwarm')
    ax.add_collection(obj);fig.colorbar(obj,ax=ax,label='U_x / U_in (CFDX natif)')
    ax.autoscale_view();ax.set(xlabel='x/H',ylabel='y/H',title=f'CFDX BFS152, iteration acceptee {iteration} — vraie carte VTU, non converge')
    ax.set_aspect('equal',adjustable='box');fig.tight_layout()
    fig.savefig(out/'01_cfdx_vtu_ux_map.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,5.2))
    d=[x for x in wall if x['patch']=='downstream_wall']
    px=np.array([x['x_over_h'] for x in d]);pcp=np.array([x['cp_owner_proxy_shifted'] for x in d])
    ax.plot(CP_EXP_BOTTOM[:,0],CP_EXP_BOTTOM[:,1],'o',ms=4,label='Exp. Driver–Seegmiller')
    ax.plot(CP_NASA[:,0],CP_NASA[:,1],'-',label='NASA CFL3D SSTm')
    ax.plot(px,pcp,'s--',label=f'CFDX iter {iteration}: Cp des cellules proches paroi, non converge')
    ax.set(xlabel='x/H',ylabel='Cp (reference aval approx.)',xlim=(-1,32),
           title='Cp — CFDX reel vs NASA / experience (comparaison DIAGNOSTIQUE seulement)')
    ax.grid(alpha=.25);ax.legend(fontsize=8);fig.tight_layout()
    fig.savefig(out/'02_cfdx_cp_owner_vs_nasa.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,5.2))
    ax.errorbar(CF_EXP[:,0],CF_EXP[:,1]*1000,yerr=CF_EXP[:,2]*1000,
                fmt='o',capsize=2,ms=3,label='Exp. Driver–Seegmiller')
    ax.plot(CF_NASA[:,0],CF_NASA[:,1]*1000,'-',label='NASA CFL3D SSTm')
    ax.plot(px,np.array([x['cf_one_point_viscous_proxy'] for x in d])*1000,
            's--',label='CFDX: proxy mu*U_owner/d (PAS Cf natif)')
    ax.axhline(0,linewidth=.7,color='gray')
    ax.set(xlim=(-1,32),xlabel='x/H',ylabel='1000 × Cf',
           title='Diagnostic de cisaillement estimé — CFDX non converge, non qualifiable')
    ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout()
    fig.savefig(out/'03_cfdx_cf_one_point_proxy_vs_nasa.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(14,5),sharey=True)
    for ax,station in zip(axes,(4,6,10)):
        rows=[x for x in profile if x['x_over_h']==station]
        e=np.array([[r['exp_u_u_ref'],r['y_over_h']] for r in rows])
        code=np.array([[r['cfdx_u_interp'],r['y_over_h']] for r in rows if r['cfdx_u_interp'] is not None])
        ax.plot(e[:,0],e[:,1],'o',label='Experience',ms=4)
        if len(code):ax.plot(code[:,0],code[:,1],'-s',ms=3,label=f'CFDX iter {iteration}')
        ax.axvline(0,color='gray',lw=.7)
        ax.set(xlabel='U_x/U_in',title=f'x/H={station}',ylim=(0,1.3))
        ax.grid(alpha=.2)
    axes[0].set_ylabel('y/H')
    axes[-1].legend(fontsize=8);fig.suptitle('CFDX vrai champ U vs vitesse experimentee — non convergé, entree differente')
    fig.tight_layout();fig.savefig(out/'04_cfdx_u_profiles_vs_experiment.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(14,5),sharey=True)
    for ax,station in zip(axes,(4,6,10)):
        rows=[x for x in profile if x['x_over_h']==station]
        e=np.array([[r['exp_uv_times1000'],r['y_over_h']] for r in rows])
        code=np.array([[r['cfdx_uv_wls_proxy'],r['y_over_h']] for r in rows if r['cfdx_uv_wls_proxy'] is not None])
        ax.plot(e[:,0],e[:,1],'o',label='Experience',ms=4)
        if len(code):ax.plot(code[:,0],code[:,1],'-s',ms=3,label='CFDX, -nut*Sxy via WLS externe')
        ax.axvline(0,color='gray',lw=.7)
        ax.set(xlabel="1000 × u'v'/Uref²",title=f'x/H={station}',ylim=(0,1.3))
        ax.grid(alpha=.2)
    axes[0].set_ylabel('y/H');axes[-1].legend(fontsize=8)
    fig.suptitle('CFDX eddy-viscosity stress proxy vs turbulence mesuree — non qualifie')
    fig.tight_layout();fig.savefig(out/'05_cfdx_uv_wls_proxy_vs_experiment.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.6))
    n=fields['nu_t']/ (mesh['mu']/mesh['rho'])
    axes[0].scatter(xyz[:,0]/mesh['height'],n,s=12)
    axes[0].set(yscale='log',xlabel='x/H',ylabel='nu_t/nu',title='CFDX natif : viscosite turbulente')
    axes[0].grid(alpha=.2)
    axes[1].scatter(xyz[:,0]/mesh['height'],fields['wall_distance']/mesh['height'],s=12)
    axes[1].set(xlabel='x/H',ylabel='d_paroi/H',title='CFDX natif : distance paroi')
    axes[1].grid(alpha=.2)
    fig.tight_layout();fig.savefig(out/'06_cfdx_nut_and_wall_distance.png',dpi=180);plt.close(fig)
    if mid_native is not None:
        it_mid, xyz_mid, f_mid=mid_native
        require(np.allclose(xyz_mid,xyz,rtol=0,atol=1e-12), 'inconsistent mesh between native VTUs')
        fig,ax=plt.subplots(figsize=(11,4.5))
        ax.plot(np.arange(len(xyz)),fields['u_x']-f_mid['u_x'],label='Delta U_x (m/s)')
        ax.plot(np.arange(len(xyz)),fields['p']-f_mid['p'],label='Delta p (Pa)')
        ax.set(xlabel='cell ID',ylabel='ecart brut (unites differentes)',
               title=f'Changement natif d etat CFDX iter {it_mid} vers {iteration} : PAS convergence')
        ax.grid(alpha=.2);ax.legend();fig.tight_layout()
        fig.savefig(out/'07_cfdx_field_change_accepted.png',dpi=180);plt.close(fig)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--case-h5',required=True)
    parser.add_argument('--vtu-dir',required=True)
    parser.add_argument('--exp-nearwall-csv',required=True)
    parser.add_argument('--output-dir',required=True)
    parser.add_argument('--allow-nonconverged-diagnostic',action='store_true')
    args=parser.parse_args()
    require(args.allow_nonconverged_diagnostic,
            'unconverged diagnostic requires explicit --allow-nonconverged-diagnostic')
    root=Path(args.vtu_dir)
    vtu=sorted([(int(m.group(1)),p) for p in root.glob('result_*.vtu')
                if (m:=re.fullmatch(r'result_(\d+)\.vtu',p.name))])
    require(vtu, 'No native accepted VTU output')
    mesh=get_mesh_and_case(args.case_h5)
    comparability=validate_one(read_case(args.case_h5))
    iteration,path=vtu[-1]
    grid,xyz,fields=read_native_vtu(path,mesh)
    require(np.allclose(mesh['points'],vtk_to_numpy(grid.GetPoints().GetData()),atol=1e-12),
            'VTU and original HDF5 geometry do not match')
    wall=wall_faces(mesh,xyz,fields,iteration)
    uv,cond=wls_uv_proxy(xyz,fields,mesh['u_ref'])
    samples=read_experimental_samples(args.exp_nearwall_csv)
    profiles=interpolated_profiles(xyz,fields,mesh['height'],mesh['u_ref'],samples,uv)
    out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
    write_csv(out/'cfdx_wall_approximation.csv',wall)
    write_csv(out/'cfdx_profiles_comparable_points.csv',profiles)
    other=None
    if len(vtu)>1:
        idx,path0=vtu[max(0,len(vtu)//2-1)]
        _,x0,f0=read_native_vtu(path0,mesh)
        other=(idx,x0,f0)
    plots(out,mesh,grid,xyz,fields,wall,profiles,iteration,other)
    code={'source':'CFDX actual accepted VTU + native HDF5 case',
          'native_output_file':str(path), 'iteration':iteration,
          'reference_geometry_status':comparability['geometry_status'],
          'benchmark_geometry_matches':comparability['geometry_checks']['pass_reference_domain'],
          'reference_comparability_gate':comparability['geometry_checks'],
          'model_status':'BFS152 pilot (NOT converged; no physical qualification)',
          'normalized_ref':'CFDX inlet U=44.2 m/s; NASA Uref near x/H=-4; geometry and BCs incompatible',
          'model':'CFDX legacy wall-resolved SST on insufficient mesh (no OpenFOAM triplet)',
          'fields_native':['u_x','u_y','p','k','omega','nu_t','wall_distance'],
          'cell_centers':'original HDF5 hex connectivity vertex average (exact for axis-aligned BFS152)',
          'vtkCellCenters_used':False,
          'fields_estimated':['wall owner Cp shifted at last downstream face',
                'viscous one-point Cf using owner Ux and normal distance',
                'modeled uv = -nu_t Sxy with external WLS gradients'],
          'unsupported':['native wall traction Cf','validated pressure-face Cp',
             'native Reynolds uv stress', 'qualified x_re/H'],
          'height_m':mesh['height'],'rho':mesh['rho'],'mu':mesh['mu'],
          'inlet_speed_m_s':mesh['u_ref'],'mesh_cells':mesh['ncell'],
          'wall_faces_downstream':int(len(mesh['patch']['downstream_wall'])),
          'pressure_zero_shift_x_over_h':wall[-1]['pressure_reference_x_over_h'],
          'model_uv_wls_valid_cells':int(np.isfinite(uv).sum()),
          'model_uv_wls_max_condition':float(np.nanmax(cond)),
          'metrics_diagnostic_only':compare_metrics(wall,profiles),
          'figures':sorted(p.name for p in out.glob('*.png')),
          'physical_validation_authorized':False,
          'absolute_comparison_warning':'NOT VALIDATION: wrong 7H/8H channel geometry vs 8H/9H experiment, 152 cells, nonconverged; previous vtkCellCenters positions incorrect.'}
    (out/'metrics_native_cfdx_v77.json').write_text(json.dumps(code,indent=2))
    print(json.dumps(code,indent=2))

if __name__=='__main__': main()
