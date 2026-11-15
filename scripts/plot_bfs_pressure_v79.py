#!/usr/bin/env python3
"""Reproducible plots from native V79 physical-face and CSR oracle logs.

Not a claim of accepted CFDX iterations or validation of Cf/Cp/uv.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np


def extract_native_oracle(log_path):
    observed=[]
    for line in Path(log_path).read_text(errors='replace').splitlines():
        if 'PRESSURE_SPLIT_ORACLE_V79 iter=' not in line:
            continue
        fields=dict(re.findall(r'([a-z_]+)=([^\s]+)',line))
        needed=('csr_relative','double_only_relative','physical_face_relative',
                'requested_tolerance','csr_closed','physical_face_closed',
                'production_commit_authorized','read_only')
        if any(key not in fields for key in needed):
            raise ValueError('truncated native pressure oracle certificate')
        vals={key:float(fields[key]) for key in needed[:4]}
        for key in ('csr_closed','physical_face_closed','production_commit_authorized','read_only'):
            vals[key]=fields[key] in ('1','true')
        if any(not math.isfinite(v) or v<0 for key,v in vals.items() if key in needed[:4]):
            raise ValueError('nonfinite or negative pressure residual')
        if vals['production_commit_authorized'] or not vals['read_only']:
            raise ValueError('V79 does not authorize a production state')
        if vals['csr_closed'] != (vals['csr_relative']<=vals['requested_tolerance']):
            raise ValueError('native CSR gate contradiction')
        if vals['physical_face_closed'] != (vals['physical_face_relative']<=vals['requested_tolerance']):
            raise ValueError('native face gate contradiction')
        observed.append(vals)
    if not observed:
        raise ValueError('no native V79 pressure lines')
    return observed


def run(csr_log, face_log, report_dir):
    before_rows=extract_native_oracle(csr_log)
    after_rows=extract_native_oracle(face_log)
    before=next((r for r in before_rows if r['csr_closed'] and not r['physical_face_closed']),None)
    after=next((r for r in after_rows if r['physical_face_closed'] and not r['csr_closed']),None)
    if before is None or after is None:
        raise ValueError('missing both independent closure modes')
    if before['requested_tolerance'] != after['requested_tolerance']:
        raise ValueError('changed physical pressure gate')
    out=Path(report_dir);out.mkdir(parents=True,exist_ok=True)
    data={'schema':'CFDX-PRESSURE-ORACLE-COMPARISON-V79/1',
          'state':'BFS20540_FIRST_PRESSURE_CORRECTOR_ONLY',
          'first_csr_oracle':before,'face_defect_corrected':after,
          'physical_face_improvement_factor':before['physical_face_relative']/after['physical_face_relative'],
          'native_linear_certificates':False,'production_commit_authorized':False,
          'physical_validation_authorized':False}
    (out/'native_pressure_oracle_comparison_v79.json').write_text(json.dumps(data,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(9.3,5.1))
    values=[[before['csr_relative'],before['physical_face_relative'],before['double_only_relative']],
            [after['csr_relative'],after['physical_face_relative'],after['double_only_relative']]]
    labels=['CSR high/low','Rhie-Chow face', 'CSR double only']
    x=np.arange(2);width=.22
    for i,label in enumerate(labels):
        ax.bar(x+(i-1)*width,[values[0][i],values[1][i]],width,label=label)
    ax.axhline(before['requested_tolerance'], ls='--', color='black',lw=1,label='Unchanged 1e-11 gate')
    ax.set_yscale('log');ax.set_ylim(3e-15,1e-8)
    ax.set_xticks(x,['After direct high/low solve','After one physical-face defect step'])
    ax.set_ylabel('Independent relative L2 residual (original operators)')
    ax.set_title('CFDX BFS20540: CSR closure is not physical-face closure')
    ax.legend(loc='upper left',fontsize=8)
    ax.grid(axis='y',alpha=.2)
    fig.tight_layout();fig.savefig(out/'01_csr_vs_rhie_chow_v79.png',dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.5,4.5))
    ax.plot([0,1],[before['physical_face_relative'],after['physical_face_relative']],marker='o',lw=2,
            label='True Rhie-Chow face continuity')
    ax.axhline(before['requested_tolerance'],ls='--',color='black',label='Unchanged gate')
    ax.set_yscale('log');ax.set_xticks([0,1],['CSR high/low oracle','One additional full-face Newton correction'])
    ax.set_ylabel('Relative cell-divergence L2')
    ax.set_title('Physical continuity reduction on the same 20,540-cell BFS')
    ax.legend();ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(out/'02_physical_face_contraction_v79.png',dpi=190);plt.close(fig)
    return data


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--csr-log',required=True)
    p.add_argument('--face-log',required=True)
    p.add_argument('--output-dir',required=True)
    a=p.parse_args()
    print(json.dumps(run(a.csr_log,a.face_log,a.output_dir),indent=2))

if __name__=='__main__':main()