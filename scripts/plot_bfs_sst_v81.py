@@ -0,0 +1,257 @@
#!/usr/bin/env python3
"""Read-only V81 evidence plots from actual CFDX logs, fail-closed on missing rows.

Do not interpret these as converged BFS physical results or NASA validation.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import re


def parse_line(line):
    items = line.strip().split()
    values = {}
    for token in items[1:]:
        if '=' not in token:
            continue
        key, value = token.split('=', 1)
        if key in values:
            raise ValueError('duplicate log field ' + key)
        values[key] = value
    return values


def get_number(row, name):
    if name not in row:
        raise ValueError('missing native field ' + name)
    value = float(row[name])
    if not math.isfinite(value):
        raise ValueError('nonfinite native field ' + name)
    return value


def read_native(path):
    text = Path(path).read_text()
    rows = text.splitlines()
    # Only the FIRST private SST transport is comparable across methods;
    # do not blend rejected retry histories into a spurious trajectory.
    first_end = next((i for i, s in enumerate(rows)
        if s.startswith('BFS_PHASE_V10') and 'phase=SST_TRANSPORT_END' in s), None)
    if first_end is None:
        raise ValueError('no complete native SST_TRANSPORT_END in ' + str(path))
    rows = rows[:first_end + 1]
    ends = parse_line(rows[-1]);
    if ends.get('committed') != '0' or ends.get('converged') != '0':
        raise ValueError('this diagnostic plot requires a rejected, nonconverged candidate')
    begin = [parse_line(s) for s in rows
             if s.startswith('BFS_PHASE_V10') and 'phase=SST_TRANSPORT_BEGIN' in s]
    if len(begin) != 1:
        raise ValueError('missing or repeated SST transport begin')
    target = get_number(begin[0], 'target')
    if not 0 < target <= 1e-8:
        raise ValueError('SST target unexpectedly changed')
    diagnostics = [parse_line(s) for s in rows if s.startswith('SST_TRANSPORT_DIAGNOSTIC ')]
    if not diagnostics:
        raise ValueError('no physical SST diagnostics in native log')
    intervals = [int(get_number(row, 'iteration')) for row in diagnostics]
    if sorted(set(intervals)) != intervals:
        raise ValueError('duplicate or out-of-order SST diagnostic iterations')
    if intervals[-1] != int(get_number(ends, 'iterations')):
        raise ValueError('last diagnostic does not match rejected SST end')
    for item in diagnostics:
        for key in ('k_returned_row_scaled','omega_returned_row_scaled',
                    'k_local_relative','omega_local_relative','dk_inf','domega_inf'):
            if get_number(item, key) < 0:
                raise ValueError('negative SST defect')
    budgets = [parse_line(s) for s in rows if s.startswith('SST_CELL_BUDGET ')]
    for b in budgets:
        signed=sum(get_number(b,f) for f in
                   ('convection','diffusion','physical_production',
                    'physical_destruction','physical_cross'))
        result=get_number(b,'residual')
        scale=sum(abs(get_number(b,f)) for f in
                  ('convection','diffusion','physical_production',
                   'physical_destruction','physical_cross'))
        if abs(signed-result)>max(1e-14,1e-9*scale):
            raise ValueError('SST physical-source decomposition inconsistent')
        if abs(get_number(b,'decomposition_error'))>max(1e-14,1e-9*scale):
            raise ValueError('SST assembled-residual closure inconsistent')
    return {'target':target,'diagnostics':diagnostics,'budgets':budgets,
            'end':ends,'log':str(path)}


def analyze_logs(logs):
    results = {key:read_native(path) for key,path in logs.items()}
    if 'picard800' not in results or 'baseline40' not in results:
        raise ValueError('both real Picard800 and baseline40 native logs required')
    return results


def export_csv(path, data):
    if not data:
        raise ValueError('empty CSV')
    keys=list(data[0].keys())
    with open(path,'w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys)
        writer.writeheader();writer.writerows(data)


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ap=argparse.ArgumentParser()
    ap.add_argument('--picard800',required=True)
    ap.add_argument('--baseline40',required=True)
    ap.add_argument('--molecular40',required=True)
    ap.add_argument('--omega-first40',required=True)
    ap.add_argument('--budget-k40',required=True)
    ap.add_argument('--budget-omega40',required=True)
    ap.add_argument('--mesh',default=None,help='Real BFS20,540 native HDF5 mesh')
    ap.add_argument('--output-dir',required=True)
    args=ap.parse_args()
    out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
    logs={'picard800':args.picard800,'baseline40':args.baseline40,
          'molecular40':args.molecular40,'omega_first40':args.omega_first40,
          'k5939':args.budget_k40,'w6028':args.budget_omega40}
    data=analyze_logs(logs)
    baseline=data['picard800'];d=baseline['diagnostics'];t=baseline['target']
    fig,ax=plt.subplots(figsize=(10.6,6.2))
    x=[int(get_number(row,'iteration')) for row in d]
    for field,label in [('k_returned_row_scaled','k: full reassembly'),
                        ('omega_returned_row_scaled','omega: full reassembly'),
                        ('dk_inf','delta k'),('domega_inf','delta omega')]:
        ax.semilogy(x,[get_number(row,field)/t for row in d],'-o',label=label,markersize=4)
    ax.axhline(1,ls='--',color='black',lw=1,label='existing acceptance threshold')
    ax.set(xlabel='SST private subcycle (NOT committed)',ylabel='Defect / unchanged 1e-8 target',
           title='CFDX V80, BFS 20,540: distance to SST acceptance')
    ax.grid(True,which='both',alpha=.25);ax.legend(fontsize=8);fig.tight_layout()
    fig.savefig(out/'01_distance_to_sst_acceptance.png',dpi=165);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,5.7))
    for field,label in [('k_local_relative','k local FV residual'),
                        ('omega_local_relative','omega local FV residual')]:
        ax.semilogy(x,[get_number(row,field) for row in d],'-o',label=label)
    ax.axhline(t,color='black',ls='--',lw=1,label='unchanged tolerance (reference only)')
    ax.set(xlabel='SST private subcycle',ylabel='Local normalized residual',
           title='Native reassembled near-wall/field residuals; no outer state accepted')
    ax.grid(True,which='both',alpha=.25);ax.legend();fig.tight_layout()
    fig.savefig(out/'02_local_residuals_v80.png',dpi=165);plt.close(fig)
    variants=['baseline40','molecular40','omega_first40']
    names=['V80 base','Molecular wall Gamma','Omega first']
    initial=data['baseline40']['diagnostics'][-1]
    fields=['k_returned_row_scaled','omega_returned_row_scaled']
    fig,ax=plt.subplots(figsize=(10.4,5.8))
    for j,field in enumerate(fields):
        for i,k in enumerate(variants):
            value=get_number(data[k]['diagnostics'][-1],field)/get_number(initial,field)
            ax.bar(j*4+i,value,width=.75,label=names[i] if j==0 else None)
            ax.text(j*4+i,value+0.009,f'{(value-1)*100:+.2f}%',ha='center',fontsize=9)
    ax.set(xticks=[1,5],xticklabels=['k physical row defect','omega physical row defect'],
           ylabel='At 40 private iterations / V80 baseline',ylim=(0,1.12),
           title='Native same-mesh controlled sensitivity: differences are small')
    ax.axhline(1,color='black',ls='--',lw=1);ax.legend();fig.tight_layout()
    fig.savefig(out/'03_sst_options_40_subcycles.png',dpi=165);plt.close(fig)
    labels=['convection','diffusion','physical_production','physical_destruction','physical_cross']
    for ident,field,source,name in [('k5939','k','k5939','04_k_near_wall_budget.png'),
                                    ('w6028','omega','w6028','05_omega_second_layer_budget.png')]:
        rows=[r for r in data[source]['budgets'] if r['field']==field]
        if len(rows)!=2 or [int(r['iteration']) for r in rows]!=[1,40]:
            raise ValueError('missing exact signed budget pair for '+ident)
        # Scale each iteration independently: for omega at cell 6028, the
        # first-iteration destruction is >8e3 whereas iteration 40 is ~142.
        # A shared axis would visually erase the physically important +/-142
        # cancellation. This is a budget comparison, NOT a convergence plot.
        fig,axes=plt.subplots(1,2,figsize=(13.2,5.8),sharey=False)
        for ax,row in zip(axes,rows):
            vals=[get_number(row,label) for label in labels]
            ax.bar(range(5),vals,width=.7)
            ax.set(xticks=range(5),xticklabels=['Convection','Diffusion','Production',
                                                'Destruction','Cross'])
            ax.tick_params(axis='x',labelrotation=23,labelsize=9)
            ax.set_ylabel('Signed FV contribution [native units]')
            residual=get_number(row,'residual')
            rel=get_number(row,'local_relative')
            ax.set_title(f'Private iteration {row["iteration"]}\n'
                         f'Signed residual={residual:.3e}, local ratio={rel:.3e}',fontsize=10)
            ax.axhline(0,color='black',lw=.8)
            ax.grid(True,axis='y',alpha=.2)
        fig.suptitle(f'CFDX {field} true signed sources: cell {ident[1:]}; '
                     'independent vertical scales, no accepted state')
        fig.tight_layout();fig.savefig(out/name,dpi=165);plt.close(fig)
    # A dimensionally grounded viscous-wall estimate; the omega max is NOT a
    # wall value certificate and is NEVER interpreted as complete wall equivalence.
    nu=1.56e-5;distance=5.808688154e-6
    omega_expected=6*nu/(.075*distance**2)
    omega_observed=get_number(d[0],'omega_reference_scale')
    fig,ax=plt.subplots(figsize=(8,5.5))
    ax.bar(['Viscous wall estimate (factor 6)','CFDX omega field scale (iter 1)'],
           [omega_expected,omega_observed]);ax.set(ylabel='omega [1/s]',
           title='Very high omega is consistent with fine wall spacing; not alone an instability')
    ax.tick_params(axis='x',labelrotation=8);fig.tight_layout()
    fig.savefig(out/'06_viscous_wall_scale.png',dpi=165);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10.2,5.6))
    ax.semilogy(x,[get_number(r,'k_reference_scale') for r in d],'-o',label='k: reported field scale')
    ax.axhline(.00109,ls='--',color='black',label='k inlet/initial = 0.00109 m2/s2')
    ax.set(xlabel='SST private subcycle (not a global steady solution)',
           ylabel='Reported field scale k [m2/s2]',
           title='Kinetic energy scale grows sharply before SST can converge')
    ax.grid(True,which='both',alpha=.25);ax.legend();fig.tight_layout()
    fig.savefig(out/'07_k_field_scale.png',dpi=165);plt.close(fig)
    if args.mesh:
        import h5py
        import numpy as np
        with h5py.File(args.mesh) as f:
            if int(f.attrs['n_cells'])!=20540:
                raise ValueError('hotspot map requires actual 20,540-cell BFS mesh')
            point=f['points'][:];co=f['cell_offsets'][:];cf=f['cell_faces'][:]
            fo=f['face_offsets'][:];fv=f['face_vertices'][:]
            centers={}
            for c in (5939,6028,20143,20407,20507):
                ff=cf[co[c]:co[c+1]]
                vv=np.unique(np.concatenate([fv[fo[i]:fo[i+1]] for i in ff]))
                centers[c]=point[vv].mean(axis=0)/.0127
        fig,axes=plt.subplots(2,1,figsize=(11.4,7.7))
        for ax in axes:
            ax.grid(True,alpha=.2)
            ax.set_xlabel('x / H')
        ax=axes[0];ax.axhline(9,ls='--',color='black',lw=1)
        for c in (20143,20407,20507):
            xy=centers[c];ax.scatter([xy[0]],[xy[1]],s=55)
            ax.annotate(f'cell {c}',(xy[0],xy[1]),xytext=(5,10),
                        textcoords='offset points',fontsize=9)
        ax.set(xlim=(-5,50),ylim=(8.96,9.006),ylabel='y / H',
               title='Early/late hotspots along the upper wall')
        ax=axes[1];ax.axhline(0,ls='--',color='black',lw=1)
        for c in (5939,6028):
            xy=centers[c];ax.scatter([xy[0]],[xy[1]],s=55)
            ax.annotate(f'cell {c}',(xy[0],xy[1]),xytext=(8,10),
                        textcoords='offset points',fontsize=9)
        ax.set(xlim=(4.5,10),ylim=(-.0004,.007),ylabel='y / H',
               title='Hotspots near the lower-wall reattachment region')
        fig.suptitle('CFDX BFS20,540 exact HDF5 cell centroids (not a converged field)')
        fig.tight_layout();fig.savefig(out/'08_native_hotspot_locations.png',dpi=165)
        plt.close(fig)
    table=[]
    for key,r in data.items():
        z=r['diagnostics'][-1]
        table.append({'run':key,'iteration':z['iteration'],
                      'k_row_scaled':z['k_returned_row_scaled'],
                      'omega_row_scaled':z['omega_returned_row_scaled'],
                      'k_local':z['k_local_relative'],'omega_local':z['omega_local_relative'],
                      'dk':z['dk_inf'],'domega':z['domega_inf'],
                      'transport_converged':r['end']['converged'],
                      'transport_committed':r['end']['committed']})
    export_csv(out/'native_sst_comparison_v81.csv',table)
    (out/'evidence_legend.json').write_text(json.dumps({
        'source_kind':'CFDX native stderr, first private SST transport only',
        'CFDX_solutions_converged':False,'CFDX_outer_accepted':0,
        'physical_comparison_valid':False,'target':t,
        'wall_distance_inferred_from_HDF5':distance,
        'logs':{k:str(v) for k,v in logs.items()}},indent=2))
    print(json.dumps({'plots':8 if args.mesh else 7,'rows':len(table),'target':t,'native_rejected':True}))

if __name__=='__main__':
    main()
