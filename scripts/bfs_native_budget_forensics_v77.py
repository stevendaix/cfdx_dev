#!/usr/bin/env python3
"""V77 native V75/V76 signed-budget forensics on the REAL production stderr.

Fail-closed: every physical term and normalization is traced to the input log.
The output is evidence of an *uncommitted, non-converged* BFS152 pilot, never
an algorithmic or turbulence-model qualification.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from bfs_comparison_plots_v76 import parse_budgets

TERMS=('convection','diffusion','production','destruction','cross')
NAMES=('Convection','Diffusion','Production','Destruction','Diffusion croisee')


def parse_kv_log(path, prefix):
    result=[]
    for line in Path(path).read_text(errors='replace').splitlines():
        if line.startswith(prefix+' '):
            fields={}
            for pair in line.split()[1:]:
                if '=' in pair:
                    k,v=pair.split('=',1)
                    if k in fields:
                        raise ValueError(f'duplicate {prefix} key: {k}')
                    fields[k]=v
            result.append(fields)
    return result


def check_native_log(path):
    path=Path(path)
    entries=parse_budgets(path)
    selected=[r for r in entries if r[0]==23 and r[1]==8 and r[2] in (1,2) and r[3] in (0,12)]
    expected={(cy,cell,eq) for cy in (1,2) for cell in (0,12) for eq in ('k','omega')}
    bykey={(cy,cell,eq):(vec,raw) for it,retry,cy,cell,eq,vec,raw in selected}
    if set(bykey)!=expected or len(selected)!=8:
        raise ValueError('missing/extra/duplicated native V75 budget rows for cycle1/2')
    statuses=parse_kv_log(path,'BFS_FIXED_DIRECTION_V74')
    statuses=[x for x in statuses if x.get('iteration')=='23' and x.get('retry')=='8']
    if len(statuses)!=2 or [(int(x['cycle']),x.get('selected')) for x in statuses]!=[(1,'1'),(2,'0')]:
        raise ValueError('V74 cycle 1 accepted privately and cycle 2 rejected is not proven')
    blocks=parse_kv_log(path,'BFS_SST_BLOCK_V76')
    blocks=[x for x in blocks if x.get('iteration')=='23']
    if len(blocks)<2 or any(x.get('selected')!='0' or x.get('reason')!='no_full_physical_descent' for x in blocks):
        raise ValueError('V76 local block rejected safeguard not proven')
    summary=parse_kv_log(path,'BFS_FIXED_DIRECTION_SUMMARY_V74')
    summary=[x for x in summary if x.get('iteration')=='23']
    if len(summary)!=1 or any(summary[0].get(k)!='0' for k in
                              ('native_momentum','native_sst','committed','qualified','mass_closed')):
        raise ValueError('native final certification missing or falsely claimed')
    argmax=parse_kv_log(path,'BFS_SST_ARGMAX_V75')
    argmax=[x for x in argmax if x.get('iteration')=='23' and x.get('cycle')=='2']
    if len(argmax)!=2 or any(x.get('after_cell')!='12' for x in argmax):
        raise ValueError('global SST residual maxima are not both at cell 12')
    records=[]
    for (cy,cell,eq), (vec,raw) in sorted(bykey.items()):
        assert len(vec)==5
        size=float(sum(abs(x) for x in vec))
        cancellation=size/max(abs(raw),1e-300)
        term_index=int(np.argmax(np.abs(vec)))
        # A largest absolute change is NOT a causal proof: budget is coupled.
        records.append(dict(cycle=cy,cell=cell,equation=eq,raw_residual_delta=raw,
              signed_terms_sum=float(sum(vec)),closure_error=float(sum(vec)-raw),
              cancellation_ratio=float(cancellation),
              largest_absolute_term=TERMS[term_index],
              **{f'{name}_delta':float(val) for name,val in zip(TERMS,vec)}))
    return {'source_log':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
          'status':'native pilot, uncommitted and not converged',
          'iteration':23,'retry':8,'private_cycle_1_selected':True,
          'private_cycle_2_rejected':True,'v76_local_block_rejected':True,
          'native_momentum_certificate_for_candidate':False,
          'native_sst_certificate_for_candidate':False,
          'global_argmax_cycle2':{r['equation']:int(r['after_cell']) for r in argmax},
          'records':records}


def plot_signed_budgets(out, audit):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    pathnames=[]
    for equation in ('k','omega'):
        fig,axes=plt.subplots(1,2,figsize=(13,5),sharey=True)
        for cy,ax in [(1,axes[0]),(2,axes[1])]:
            row=next(r for r in audit['records'] if r['cycle']==cy and r['cell']==12 and r['equation']==equation)
            x=[row[f'{name}_delta'] for name in TERMS]
            ax.barh(NAMES,x)
            ax.axvline(0,color='gray',linewidth=1)
            ax.axvline(row['raw_residual_delta'],color='black',linestyle=':',label='Delta R net')
            ax.set_title(f'Cycle {cy} — '+('accepte en prive' if cy==1 else 'rejete'))
            ax.set_xlabel(f'Delta terme signe dans R_{equation}')
            ax.text(.02,.025, f"Delta R={row['raw_residual_delta']:+.5e}\n"
                        f"annulation= {row['cancellation_ratio']:,.0f}x",
                    transform=ax.transAxes,ha='left',va='bottom',fontsize=9,
                    bbox={'facecolor':'white','edgecolor':'0.7','alpha':.85})
            ax.legend(loc='lower right',fontsize=8)
            ax.grid(axis='x',alpha=.2)
        fig.suptitle(f'CFDX V77 — bilans natifs de {equation} a la cellule 12 (iter 23)\n'
                     'Variations entre etats reassembles; ce ne sont pas des causes isolees')
        fig.tight_layout()
        name=f'09_bilan_signe_{equation}_c12_cycles1_2.png'
        fig.savefig(out/name,dpi=180,bbox_inches='tight');plt.close(fig)
        pathnames.append(name)
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),sharey=False)
    for ax,eq in zip(axes,['k','omega']):
        values=[next(r for r in audit['records'] if r['cycle']==cy and r['cell']==12 and r['equation']==eq)['raw_residual_delta'] for cy in (1,2)]
        ax.bar(['cycle 1','cycle 2'],values)
        ax.axhline(0,color='gray');ax.set_ylabel(f'Delta residu physique brut R_{eq} (cellule 12)')
        ax.set_title(f'R_{eq}: amelioration puis degradation')
        ax.grid(axis='y',alpha=.3)
    fig.suptitle('CFDX V77 — variation REELLE du residu non lineaire reassemble, cellule 12')
    fig.tight_layout();name='10_residu_reel_cellule12_cycles1_2.png'
    fig.savefig(out/name,dpi=180,bbox_inches='tight');plt.close(fig);pathnames.append(name)
    return pathnames



def plot_native_accepted_history(out, log):
    """Only accepted V65 production transactions, never rejected line-search tries."""
    accepted=parse_kv_log(log,'GLOBAL_SAME_OPERATOR_COMMITTED_V65')
    evidence=[]
    for a in accepted:
        if a.get('qualified')!='0' or a.get('physical_converged')!='0':
            raise ValueError('unexpected false qualification in native V65 history')
        evidence.append((int(a['iteration']),float(a['merit'])))
    if [i for i,v in evidence]!=list(range(1,23)) or any(not math.isfinite(v) or v<=0 for i,v in evidence):
        raise ValueError('non-contiguous / nonfinite accepted CFDX history')
    fig,ax=plt.subplots(figsize=(9,5))
    ax.semilogy([i for i,v in evidence],[v for i,v in evidence],'o-')
    ax.set_xlabel('Iteration CFDX reelle acceptee');ax.set_ylabel('Merite V65 historique')
    ax.set_title('CFDX natif — 22 transactions acceptees; iteration 23 rejetee')
    ax.set_xticks([1,2,4,6,8,10,12,14,16,18,20,22,23])
    ax.axvline(23,linestyle=':',color='black',label='Iteration 23 rejetee (aucun etat VTU)')
    ax.grid(alpha=.25);ax.legend(loc='best');fig.tight_layout()
    name='11_cfdx_merite_22_etats_acceptes.png'
    fig.savefig(Path(out)/name,dpi=180,bbox_inches='tight');plt.close(fig)
    return name,evidence

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--native-log',required=True)
    p.add_argument('--output-dir',required=True)
    args=p.parse_args()
    out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
    audit=check_native_log(args.native_log)
    audit['figures']=plot_signed_budgets(out,audit)
    history_name,history=plot_native_accepted_history(out,args.native_log)
    audit['figures'].append(history_name)
    audit['accepted_v65_merit_history']=history
    (out/'signed_budget_evidence_v77.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False)+'\n')
    with (out/'signed_budget_rows_v77.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(audit['records'][0]));w.writeheader();w.writerows(audit['records'])
    print(json.dumps({k:v for k,v in audit.items() if k!='records'},indent=2))


if __name__=='__main__':main()