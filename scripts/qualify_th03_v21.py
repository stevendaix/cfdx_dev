#!/usr/bin/env python3
"""Independent, fail-closed TH-03 VTK postprocessor (diagnostic, not certification).

Reads CFDX ASCII .vtu outputs WITHOUT CFDX solver internals or VTK bindings.
Surfaces the first-cell velocity sign-change as a *proxy* for reattachment,
not as a certified wall shear Cf=0. The 6.26+/-0.10 reference (Driver &
Seegmiller) is NOT directly comparable to this surrogate/coarse mesh and
surrogate initial/inlet fields; no comparison automatically validates SST.
"""
from __future__ import annotations
import argparse
import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

H = 0.0127
NU = 1.56e-5
U_NOMINAL = 44.2
RHO = 1.2
REFERENCE_EXPERIMENT = 6.26
REFERENCE_SIGMA = 0.10
REFERENCE_SSTM = 6.50


def read_ascii_vtu(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str,np.ndarray]]:
    root=ET.parse(path).getroot()
    piece=root.find('.//Piece')
    if piece is None:
        raise ValueError('VTU has no Piece')
    def load(node, dtype=float):
        if node is None or node.get('format') != 'ascii':
            raise ValueError('V21 requires ASCII VTU arrays (do not silently misparse binary)')
        data=np.fromstring(node.text or '',sep=' ',dtype=dtype)
        if not np.isfinite(data).all():
            raise ValueError('nonfinite VTU data')
        return data
    points_element=piece.find('./Points/DataArray')
    points=load(points_element).reshape(-1,3)
    cells=piece.find('./Cells')
    if cells is None: raise ValueError('VTU has no Cells')
    arrays={a.get('Name'):a for a in cells.findall('DataArray')}
    conn=load(arrays.get('connectivity'),np.int64)
    offsets=load(arrays.get('offsets'),np.int64)
    if len(offsets)!=int(piece.get('NumberOfCells','-1')) or (len(offsets)>0 and offsets[-1]!=len(conn)):
        raise ValueError('invalid VTU cell topology')
    if np.any(conn<0) or np.any(conn>=len(points)) or np.any(np.diff(offsets)<=0):
        raise ValueError('invalid connectivity')
    fields={}
    for a in piece.findall('./CellData/DataArray'):
        fields[a.get('Name')]=load(a)
    for name in ('p','u_x','u_y','u_z'):
        if name not in fields or len(fields[name])!=len(offsets):
            raise ValueError('missing/inconsistent VTU field '+name)
    centre=[];bounds=[];start=0
    for end in offsets:
        vertices=points[conn[start:end]]
        centre.append(vertices.mean(axis=0))
        bounds.append([vertices[:,0].min(),vertices[:,0].max(),vertices[:,1].min(),vertices[:,1].max()])
        start=end
    return np.asarray(centre),np.asarray(bounds),points,fields


def wall_sign_change(x_over_h: np.ndarray, ux: np.ndarray) -> list[float]:
    """One-cell reattachment proxy: linear zero-crossing from Ux<0 to Ux>=0."""
    if x_over_h.shape!=ux.shape or not np.isfinite(ux).all() or not np.isfinite(x_over_h).all():
        raise ValueError('nonfinite or inconsistent wall samples')
    if np.any(np.diff(x_over_h)<=0):
        raise ValueError('wall x coordinate must increase strictly')
    roots=[]
    for j in range(1,len(x_over_h)):
        if ux[j-1]<0<=ux[j]:
            roots.append(float(x_over_h[j-1] +
                (-ux[j-1])/(ux[j]-ux[j-1])*(x_over_h[j]-x_over_h[j-1])))
    return roots


def diagnostic(path: Path, *, require_sst: bool=False) -> dict:
    centres,bounds,points,fields=read_ascii_vtu(path)
    n=len(centres)
    wall=np.isclose(bounds[:,2],0,rtol=0,atol=1e-11)&(bounds[:,0]>=-1e-11)
    wall_idx=np.where(wall)[0]
    wall_idx=wall_idx[np.argsort(centres[wall_idx,0])]
    if len(wall_idx)==0:
        raise ValueError('cannot identify lower wall after the step')
    xx=centres[wall_idx,0]/H
    yy=centres[wall_idx,1]
    ux=fields['u_x'][wall_idx]
    if np.any(yy<=0): raise ValueError('invalid first-cell wall distance')
    roots=wall_sign_change(xx,ux)
    # One-point wall derivative, not a verified wall-traction reconstruction.
    cf_proxy=2.0*NU*ux/(U_NOMINAL**2*yy)
    # The experiment shifts Cp toward zero at x/H=40. The coarse mesh
    # does not include an exact wall sample at x/H=40: use its farthest
    # wall sample as a DIAGNOSTIC, not an experimental Cp comparison.
    downstream_p=float(fields['p'][wall_idx[-1]])
    cp_shift_proxy=(fields['p'][wall_idx]-downstream_p)/(0.5*RHO*U_NOMINAL**2)
    yplus_viscous_proxy=np.sqrt(np.abs(ux)*yy/NU)
    output={
        'file':str(path),'cells':n,
        'geometry_proxy':'first-cell vertex centroid; not finite-volume cell centroid on skewed cells',
        'reference':{'experiment_xr_h':REFERENCE_EXPERIMENT,'experiment_sigma':REFERENCE_SIGMA,
                     'NASA_SSTm_quasi_steady_xr_h':REFERENCE_SSTM,
                     'directly_comparable':False,
                     'reason':'surrogate coarse mesh, different inlet profile and wall treatment, first-order upwind'},
        'n_wall_cells':int(len(wall_idx)),'n_reversed_wall_cells':int((ux<0).sum()),
        'near_wall_ux_zero_crossings_x_h_proxy':roots,
        'distance_to_experiment_x_h_proxy':float(roots[-1]-REFERENCE_EXPERIMENT) if roots else None,
        'wall_cf_proxy_min':float(cf_proxy.min()),'wall_cf_proxy_max':float(cf_proxy.max()),
        'wall_cp_shift_proxy_min':float(cp_shift_proxy.min()),
        'wall_cp_shift_proxy_max':float(cp_shift_proxy.max()),
        'wall_cp_shift_proxy_min_at_x_h':float(xx[int(np.argmin(cp_shift_proxy))]),
        'wall_cp_shift_anchor_x_h':float(xx[-1]),
        'wall_cp_shift_anchor_p_pa':downstream_p,
        'wall_cp_note':'shifted by final lower-wall sample, not NASA experimental Cp alignment at x/H~40',
        'wall_yplus_viscous_proxy_min':float(yplus_viscous_proxy.min()),
        'wall_yplus_viscous_proxy_max':float(yplus_viscous_proxy.max()),
        'u_min_m_s':float(fields['u_x'].min()),'u_max_m_s':float(fields['u_x'].max()),
        'speed_max_m_s':float(np.hypot(fields['u_x'],fields['u_y']).max()),
        'pressure_span_pa':float(np.ptp(fields['p'])),
        'near_wall_table':[
            {'cell':int(i),'x_h':float(x),'y_h':float(y/H),
             'u_x':float(u),'cf_proxy':float(c),'yplus_proxy':float(yp),
             'p_pa':float(fields['p'][i]),'cp_shift_proxy':float(cp)}
            for i,x,y,u,c,yp,cp in zip(wall_idx,xx,yy,ux,cf_proxy,yplus_viscous_proxy,cp_shift_proxy)]
    }
    field_names=('k','omega','nu_t','wall_distance','sst_strain_rate')
    if require_sst and not all(x in fields for x in field_names):
        raise ValueError('SST VTU fields not exported; must run V21 production solver')
    if all(x in fields for x in field_names):
        k=fields['k'];omega=fields['omega'];nu_t=fields['nu_t']
        if np.any(k<0) or np.any(omega<=0) or np.any(nu_t<0):
            raise ValueError('inadmissible SST state')
        output['sst']={'k_min':float(k.min()),'k_max':float(k.max()),
                       'omega_min':float(omega.min()),'omega_max':float(omega.max()),
                       'nut_over_nu_max':float((nu_t/NU).max()),
                       'nut_over_nu_wall_max':float((nu_t[wall_idx]/NU).max()),
                       'k_near_wall_max':float(k[wall_idx].max())}
    summary_path=path.parent/'convergence.json'
    if summary_path.exists():
        cfg=json.loads(summary_path.read_text())
        output['convergence']={'converged':bool(cfg.get('converged',False)),
                               'iterations':int(cfg.get('iterations',0)),
                               'reason':cfg.get('convergence_reason'),
                               'last_momentum_rhs_relative':cfg.get('history',[-1])[-1].get('momentum_equation_residual_rhs_relative') if cfg.get('history') else None,
                               'last_continuity_normalized':cfg.get('history',[-1])[-1].get('continuity_normalized') if cfg.get('history') else None}
    output['physical_qualification']='NOT_QUALIFIED'
    return output


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('vtu',nargs='+',type=Path)
    ap.add_argument('--require-sst',action='store_true')
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    records=[diagnostic(p,require_sst=args.require_sst) for p in args.vtu]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps({'results':records},indent=2)+'\n')
    for row in records:
        print(row['file'], 'reattachment Ux proxy:',row['near_wall_ux_zero_crossings_x_h_proxy'],
              'momentum:',row.get('convergence',{}).get('last_momentum_rhs_relative'))
    return 0
if __name__=='__main__':raise SystemExit(main())
