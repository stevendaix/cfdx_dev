#!/usr/bin/env python3
"""V76 comparison plots. Never relabel reference or V74 results as CFDX-V76.

If NASA raw data files are available, read them; otherwise use the subset of
NASA reference points transcribed from the official public TMR source. No
unmeasured CFDX Cf/Cp/profile curves are invented.
"""
import argparse
import csv
import json
import math
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Public TMR numerical rows; source URLs in V76 report.
CF_EXP=np.array([
[-3.956,.00288,.000202],[-1.804,.00285,.000171],[-.804,.00311,.000311],
[.484,.00009,.000014],[1.804,-.00011,.000019],[2.804,-.00052,.000109],
[3.804,-.00102,.000163],[4.804,-.00082,.000123],[5.882,-.00022,.000062],
[7.09,.00045,.000167],[7.978,.00102,.000122],[9.196,.00118,.000106],
[11.196,.00154,.000123],[13.196,.00173,.000121],[16.196,.00191,.000134],
[20.196,.00202,.000141],[24.136,.00204,.000143],[28.13,.00205,.000205],
[32.196,.00212,.000148],[35.994,.00230,.000184]])
CP_EXP_BOTTOM=np.array([
[-8.5,-.1102],[-6.5,-.1141],[-4.5,-.1189],[-2.5,-.1372],
[-.5,-.1613],[.5,-.1777],[1,-.1806],[2,-.1883],[3,-.1903],
[3.5,-.1806],[4,-.1565],[4.5,-.1064],[5,-.0659],[5.5,-.0359],
[6,-.0012],[6.5,.0162],[7,.0248],[8,.0374],[9,.0394],
[11,.0336],[13,.0268],[15,.0201],[17,.0162],[19.5,.0144],
[21.5,.0105],[23.5,.0086],[25.5,.0086],[27.5,.0077],
[29.5,.0077],[31.5,.0058],[33.5,.0019],[35.5,.0048],[37.5,0.]])
CP_EXP_TOP=np.array([
[-5,-.1054],[3,-.0727],[5,-.0630],[7,-.0543],[9,-.0331],
[11,-.0148],[13,-.0023],[15,-.0003],[17,.0065],[19,.0075],[21,.0075],
[23,.0066],[25,.0095],[27,.0076],[29,.0067],[31,.0048],[33,.0019],
[35,.0038],[37,-.0010]])
# Sparse, original TMR CFL3D SSTm rows; not interpolated or digitized from an image.
CF_NASA=np.array([
[-4.011669922,.002944703447],[-1.013258934,.003105584066],
[-.1043625772,.003420332447],[.02083333395,-.000004398731107],
[.5,.00002829045116],[1,.00005834708645],[1.333333373,.000003125300054],
[1.354166627,-.00000935728076],[2,-.0005134993116],
[2.8125,-.001058889669],[3,-.001162969973],[3.5,-.001370295533],
[3.958333254,-.001435520942],[4,-.001434600912],
[4.5,-.001335233916],[4.8125,-.001198932994],
[5.25,-.0009376115049],[5.75,-.0005808479036],
[6,-.0003943655174],[6.5,-.00003034047404],
[6.541666508,-.000001252024276],[6.5625,.00001319786043],
[7,.0003000584256],[7.5,.0005834890762],[8.331223488,.0009460118599],
[9.949404716,.001359018148],[10.03099632,.001372555620],
[12.03322983,.001592228538],[14.00928783,.001697744941],
[16.25933647,.00176689832],[19.80643654,.001832234906],
[24.21634293,.001884140540],[28.57291985,.001924597891],
[32.10914230,.001953530591],[39.82146835,.001987882191]])
CP_NASA=np.array([
[.5,-.1961686313],[1.0625,-.1969980896],[2,-.2071116269],
[3,-.2002431452],[3.5,-.1836685240],[4,-.1579189003],
[4.5,-.1258537471],[5,-.09196963906],[5.5,-.06041380763],
[6,-.03384926915],[6.5,-.01328012161],[7,.001533718547],
[7.5,.01150176115],[8.0,.01860311814],[9.031744003,.02333978005],
[10.03099632,.02300504036],[11.04629993,.02067676373],
[12.03322983,.01791188307],[13.12444782,.01497969218],
[14.00928783,.01289150305],[16.25933647,.008889706805],
[20.28434181,.005098911002],[24.87245369,.003199906787],
[30.27099419,.001962991199],[37.39556503,.0003113465209]])

V74={
"stage":["V72 depart","Cycle 1 (prive)","Cycle 2 (rejete)"],
"mass":[6.248570238841e-4,5.892457607501e-4,5.556022008498e-4],
"merit":[.161860259313,.159609481254,.157455190171],
"Ux":[.058069822590,.058125430507,.059430025465],
"Uy":[.161860259313,.159609481254,.157455190171],
"cont":[5.969124869765e-4,5.628938125771e-4,5.307548427916e-4],
"k":[2.067910319789e-5,1.647446007467e-5,1.945520861587e-5],
"omega":[1.925824272711e-5,1.445226163835e-5,1.720517970796e-5],
"phi3":[3.428729133897e-4,3.184483145399e-4,2.954926949845e-4],
"rAU0":[1296.438930927,1297.942792364,1299.567376136],
"rAU12":[524.387175709,526.250674196,528.012739907],
"nut0":[.009027956139,.008999536084,.008971112137],
"nut12":[.006236538530,.006238132536,.006238802029]}

def load_tecplot(path,cols):
    zones={"default":[]}
    current="default"
    if path is None or not Path(path).exists():return {}
    for line in Path(path).read_text(errors="replace").splitlines():
        text=line.strip()
        if not text or text.startswith("#") or text.lower().startswith("variables"):
            continue
        if text.lower().startswith("zone"):
            found=re.search(r'\"([^\"]+)\"',text)
            current=found.group(1) if found else "zone"+str(len(zones))
            zones.setdefault(current,[])
            continue
        p=text.replace(","," ").split()
        if len(p)<cols:continue
        try: vals=[float(x.replace("D","E")) for x in p[:cols]]
        except ValueError:continue
        if all(math.isfinite(x) for x in vals):zones[current].append(vals)
    return {k:np.array(v) for k,v in zones.items() if v}

def bzero(points,near=(3.,12.)):
    p=points[(points[:,0]>=near[0])&(points[:,0]<=near[1])]
    for a,b in zip(p[:-1],p[1:]):
        if a[1]<0. and b[1]>=0.:
            return float(a[0]-a[1]*(b[0]-a[0])/(b[1]-a[1]))
    return None

def plot_bfs(out,cf_e,cf_m,cp_e,cp_m):
    exp_reat=bzero(cf_e[:,0:2])
    model_reat=bzero(cf_m)
    fig,ax=plt.subplots(figsize=(10,5.5))
    ax.errorbar(cf_e[:,0],cf_e[:,1]*1000,yerr=cf_e[:,2]*1000,fmt="o",
                ms=4,capsize=2,label="Driver-Seegmiller (NASA exp., avec erreurs)")
    ax.plot(cf_m[:,0],cf_m[:,1]*1000,"-s",ms=2,lw=1.5,
            label="NASA CFL3D SSTm (points officiels echantillonnes)")
    ax.axhline(0,color="gray",lw=.8);ax.axvspan(6.16,6.36,alpha=.13,label="x_r/H exp. 6.26 +/- 0.10")
    if model_reat:ax.axvline(model_reat,ls="--",alpha=.8,label=f"CFL3D x_r/H={model_reat:.3f}")
    ax.set(xlabel="x/H",ylabel="1000 x Cf",title="BFS - frottement parietal NASA SSTm vs experience",xlim=(-5,32))
    ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/"01_cf_reference.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,5))
    ax.plot(cp_e[:,0],cp_e[:,1],"o",label="Experience Cp paroi inferieure")
    ax.plot(cp_m[:,0],cp_m[:,1],"-s",ms=2,label="NASA CFL3D SSTm echantillonne")
    ax.axhline(0,color="gray",lw=.8)
    ax.set(xlabel="x/H",ylabel="Cp (decale, zero vers x/H=40)",
           title="BFS - recuperation de pression pres du recollement",xlim=(-1,30))
    ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(out/"02_cp_reference.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,4.7))
    ax.plot(CP_EXP_BOTTOM[:,0],CP_EXP_BOTTOM[:,1],"-o",ms=3,label="Cp experimental paroi inferieure")
    ax.plot(CP_EXP_TOP[:,0],CP_EXP_TOP[:,1],"-s",ms=3,label="Cp experimental paroi superieure")
    ax.axhline(0,color="gray",lw=.8)
    ax.set(xlabel="x/H",ylabel="Cp decale",xlim=(-5,38),
        title="Driver-Seegmiller - pression sur parois inferieure et superieure")
    ax.grid(alpha=.25);ax.legend();fig.tight_layout()
    fig.savefig(out/"12_exp_cp_top_bottom.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,4.7))
    inside=(cf_e[:,0]>=cf_m[:,0].min())&(cf_e[:,0]<=cf_m[:,0].max())
    diff=np.interp(cf_e[inside,0],cf_m[:,0],cf_m[:,1])-cf_e[inside,1]
    ax.errorbar(cf_e[inside,0],diff*1e3,yerr=cf_e[inside,2]*1e3,
                fmt="o-",capsize=2,label="CFL3D SSTm - experience")
    ax.axhline(0,color="black",lw=.8);ax.set(xlabel="x/H",ylabel="1000 x (Cf_SSTm - Cf_exp)",
        title="Ecart Cf au point de mesure (interpolation lineaire CFL3D)",xlim=(-4,30))
    ax.grid(alpha=.25);fig.tight_layout();fig.savefig(out/"03_delta_cf_reference.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,4.7))
    inside=(cp_e[:,0]>=cp_m[:,0].min())&(cp_e[:,0]<=cp_m[:,0].max())
    diff=np.interp(cp_e[inside,0],cp_m[:,0],cp_m[:,1])-cp_e[inside,1]
    ax.plot(cp_e[inside,0],diff,"o-");ax.axhline(0,color="black",lw=.8)
    ax.set(xlabel="x/H",ylabel="Delta Cp (CFL3D - exp)",
        title="Ecart de recuperation Cp a stations experimentales",xlim=(0,30))
    ax.grid(alpha=.25);fig.tight_layout();fig.savefig(out/"04_delta_cp_reference.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,5.2))
    keys=["Ux","Uy","cont","k","omega"]
    x=np.arange(3)
    for name in keys: ax.plot(x,V74[name],"-o",label=name)
    ax.set_yscale("log");ax.axvspan(1.5,2.5,alpha=.13,label="cycle 2 rejete")
    ax.set_xticks(x,V74["stage"]);ax.set(ylabel="residu V74 normalise (echelles communes)",
        title="CFDX BFS152 - six-groupes physiques: valeurs V74 documentees")
    ax.grid(alpha=.25,which="both");ax.legend(fontsize=9);fig.tight_layout();fig.savefig(out/"05_v74_residus.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,5.2))
    keys=["mass","merit","Ux","Uy","cont","k","omega"]
    for name in keys:
        y=np.array(V74[name]);ax.plot(np.arange(3),y/y[0],"-o",label=name)
    ax.axhline(1,color="gray",ls="--",lw=.7)
    ax.axvspan(1.5,2.5,alpha=.13);ax.set_xticks(np.arange(3),V74["stage"])
    ax.set(ylabel="rapport a l'etat V72 (=1)",title="CFDX BFS152 - echanges de contraction / croissance")
    ax.grid(alpha=.2);ax.legend(ncol=4,fontsize=8);fig.tight_layout();fig.savefig(out/"06_v74_ratios.png",dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    for ax,keys,title in [(axes[0],["phi3","rAU0","rAU12"],"Flux face 3 et rAU"),
                          (axes[1],["nut0","nut12"],"Viscosites turbulentes")]:
        for name in keys:
            values=np.array(V74[name]);ax.plot([0,1,2],values/values[0],"-o",label=name)
        ax.axhline(1,color="gray",ls="--",lw=.6);ax.axvspan(1.5,2.5,alpha=.12)
        ax.set_xticks([0,1,2],["depart","cycle 1","cycle 2"])
        ax.set(ylabel="valeur relative au depart",title=title);ax.legend();ax.grid(alpha=.2)
    fig.suptitle("CFDX V74 - coefficients changes a chaque co-iteration")
    fig.tight_layout();fig.savefig(out/"07_v74_face3_coefficients.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,4.7))
    nums=[0,1,11,22];labels=["y+ <= 1","1 < y+ < 30","30 <= y+ <= 300","y+ > 300"]
    ax.bar(labels,nums);ax.set(ylabel="Nombre de faces (N=34)",ylim=(0,25),
        title="CFDX BFS152 - regimes de paroi estimes V74 (solution non convergee)")
    for i,xv in enumerate(nums):ax.text(i,xv+.6,str(xv),ha="center")
    ax.tick_params(axis="x",labelrotation=15);ax.grid(axis="y",alpha=.25)
    fig.tight_layout();fig.savefig(out/"08_v74_yplus_bins.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,4.5))
    ax.errorbar([0],[6.26],yerr=[.1],fmt="o",capsize=5,label="Experience Driver-Seegmiller")
    if model_reat:ax.scatter([1],[model_reat],s=60,label="NASA CFL3D SSTm (Cf, interpolation)")
    ax.set_xticks([0,1],["Experience","NASA SSTm"])
    ax.set(ylabel="x_r/H",ylim=(5.9,6.9),title="Position de recollement - methode Cf=0")
    ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(out/"09_reattachement.png",dpi=180);plt.close(fig)
    return {"cf_exp_zero_linear":exp_reat,"cf_cfl3d_sstm_zero_linear":model_reat,
        "reattachment_nasa_minus_exp":model_reat-6.26 if model_reat else None,
        "cf_exp_min":float(cf_e[:,1].min()),"cf_sstm_sample_min":float(cf_m[:,1].min()),
        "cf_sstm_min_vs_exp_percent":100*(abs(cf_m[:,1].min())/abs(cf_e[:,1].min())-1),
        "cp_error_at_6_5":float(np.interp(6.5,cp_m[:,0],cp_m[:,1])-.0162),
        "cp_error_at_5":float(np.interp(5,cp_m[:,0],cp_m[:,1])-(-.0659)),
        "reference_only":True}

def plot_experimental_nearwall(out):
    """Experimental profiles only; never present them as CFDX or NASA SSTm."""
    data_file=Path(__file__).resolve().parents[1]/"data"/"nasa_tmr_profiles_exp_nearwall_v76.csv"
    if not data_file.is_file():return []
    measured=np.genfromtxt(data_file,delimiter=",",names=True)
    if measured.size != 52 or not np.isfinite(measured["uUref"]).all():
        raise ValueError("V76 incomplete nearwall transcription")
    plots=[]
    for ycol,label,title,name in [
        ("uUref","u / Uref", "Driver-Seegmiller - profils de vitesse pres de la marche",
         "10_experimental_nearwall_u.png"),
        ("uv_times1000","1000 x u'v' / Uref^2",
         "Driver-Seegmiller - contraintes de Reynolds pres de la marche",
         "11_experimental_nearwall_uv.png")]:
        fig,ax=plt.subplots(figsize=(8.5,5))
        for station in (1,4,6,10):
            a=measured[measured["xH"]==station]
            ax.plot(a[ycol],a["yH"],"-o",ms=3,label=f"x/H={station}")
        ax.axvline(0,color="gray",lw=.8)
        ax.set(xlabel=label,ylabel="y/H",title=title,ylim=(0,1.25))
        ax.grid(alpha=.2);ax.legend();fig.tight_layout()
        fig.savefig(out/name,dpi=180);plt.close(fig);plots.append(name)
    return plots

def load_code_curve_csv(path,kind):
    """Read REAL converged CFDX/OpenFOAM Cf/Cp samples; no normalization guesses."""
    if kind not in ("cf","cp"):
        raise ValueError("unknown comparison curve type")
    with Path(path).open(newline="") as f:
        reader=csv.DictReader(f)
        if not reader.fieldnames or "x_over_h" not in reader.fieldnames or kind not in reader.fieldnames:
            raise ValueError(f"{path}: expected x_over_h,{kind}")
        rows=[]
        for r in reader:
            rows.append([float(r["x_over_h"]),float(r[kind])])
    a=np.array(rows,dtype=float)
    if a.ndim!=2 or a.shape[0]<4 or a.shape[1]!=2 or not np.isfinite(a).all():
        raise ValueError("external CFD curve missing or nonfinite")
    if not (np.diff(a[:,0])>0).all():
        raise ValueError("external CFD x_over_h must increase strictly")
    return a


def plot_code_curves(out,args):
    """Compare only externally supplied converged curves with the same scales."""
    records=[]
    for label,kw in [("CFDX", "cfdx"),("OpenFOAM","openfoam")]:
        for kind in ("cf","cp"):
            path=getattr(args,f"{kw}_{kind}")
            if path:records.append((label,kind,path,load_code_curve_csv(path,kind)))
    if not records:return {}
    data={}
    generated=[]
    for kind, exp, nasa, file in (("cf", CF_EXP[:,:2],CF_NASA,"13_cf_multicode.png"),
                                  ("cp", CP_EXP_BOTTOM,CP_NASA,"14_cp_multicode.png")):
        series=[v for v in records if v[1]==kind]
        if not series:continue
        fig,ax=plt.subplots(figsize=(10,5.2))
        multiplier=1000 if kind=="cf" else 1
        ax.plot(exp[:,0],exp[:,1]*multiplier,"o",ms=4,label="Experience")
        ax.plot(nasa[:,0],nasa[:,1]*multiplier,"-",label="NASA CFL3D SSTm")
        for label,_,path,curve in series:
            ax.plot(curve[:,0],curve[:,1]*multiplier,"-",lw=2,label=label)
            mask=(exp[:,0]>=curve[0,0])&(exp[:,0]<=curve[-1,0])
            if mask.sum()<2:raise ValueError("external curve does not overlap >=2 measurements")
            err=np.interp(exp[mask,0],curve[:,0],curve[:,1])-exp[mask,1]
            info={"path":str(path),"measured_stations":int(mask.sum()),
                  "rmse_on_measured_stations":float(np.sqrt(np.mean(err**2))),
                  "max_abs_error":float(np.max(np.abs(err)))}
            if kind=="cf":info["zero_crossing_x_over_h"]=bzero(curve)
            data[f"{label.lower()}_{kind}"]=info
        ax.axhline(0,color="gray",lw=.8)
        ax.set(xlabel="x/H",ylabel="1000 x Cf" if kind=="cf" else "Cp (shifted)",
               title=f"BFS: {kind.upper()} NASA vs experience et courbes fournies")
        ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/file,dpi=180);plt.close(fig)
        generated.append(file)
    return {"user_supplied_comparison_plots":generated,"user_supplied_metrics":data}

def parse_budgets(log):
    entries=[]
    for line in Path(log).read_text(errors="replace").splitlines():
        if not line.startswith("BFS_SST_BUDGET_V75 "):continue
        d={}
        for kv in line.split()[1:]:
            if "=" in kv:
                k,v=kv.split("=",1)
                if k in d:raise ValueError("duplicated V75 key")
                d[k]=v
        mandatory=["iteration","retry","cycle","cell","equation","trace_selected", "committed",
                   "convection_delta","diffusion_delta","production_delta","destruction_delta",
                   "cross_delta","terms_delta_sum","raw_delta","closure_error"]
        if any(k not in d for k in mandatory):raise ValueError("incomplete V75 budget row")
        if d["committed"]!="0":raise ValueError("V75 fabricated commit")
        if d["trace_selected"] not in ("0","1"):raise ValueError("invalid private selection")
        if d["equation"] not in ["k","omega"]:raise ValueError("unknown SST variable")
        vals=np.array([float(d[k]) for k in ("convection_delta","diffusion_delta","production_delta",
                      "destruction_delta","cross_delta")]);raw=float(d["raw_delta"])
        if not np.isfinite(vals).all() or not np.isfinite(raw):raise ValueError("nonfinite V75 terms")
        signed_sum=float(sum(vals))
        reported_sum=float(d["terms_delta_sum"])
        reported_closure=float(d["closure_error"])
        if not np.isfinite([reported_sum,reported_closure]).all():
            raise ValueError("V75 nonfinite claimed certificate")
        eps=1e-12+1e-8*max(abs(raw),float(sum(abs(vals))))
        if abs(signed_sum-raw)>eps:
            raise ValueError("V75 budget nonconservative")
        if abs(reported_sum-signed_sum)>eps or abs(reported_closure-(signed_sum-raw))>eps:
            raise ValueError("V75 budget certificate inconsistent")
        if "raw_before" in d and "raw_after" in d:
            if abs((float(d["raw_after"])-float(d["raw_before"]))-raw)>eps:
                raise ValueError("V75 reassembled residual delta inconsistent")
        entries.append((int(d["iteration"]),int(d["retry"]),int(d["cycle"]),int(d["cell"]),d["equation"],vals,raw))
    if not entries:raise ValueError("no V75 data in log")
    unique={(x[0],x[1],x[2],x[3],x[4]) for x in entries}
    if len(unique)!=len(entries):raise ValueError("duplicate V75 budget")
    return entries

def plot_budgets(out,entries):
    named=["Convection","Diffusion","Production","Destruction","Cross-diff."]
    images=[]
    for iteration,retry,cycle,cell,eq,terms,raw in entries:
        fig,ax=plt.subplots(figsize=(9.5,4.5))
        ax.barh(named,terms)
        ax.axvline(0,color="gray",lw=.9)
        ax.set(xlabel="Delta contribution signed to physical residual",
               title=f"V75 SST budget: iteration {iteration}, retry {retry}, cycle {cycle}, cell {cell}, {eq}")
        ax.text(.99,.03,f"Sum={sum(terms):.4e}; Delta R={raw:.4e}",
                transform=ax.transAxes,ha="right",fontsize=9)
        ax.grid(axis="x",alpha=.2)
        fig.tight_layout()
        name=f"10_budget_it{iteration}_retry{retry}_cycle{cycle}_c{cell}_{eq}.png"
        fig.savefig(out/name,dpi=180);plt.close(fig);images.append(name)
    return images

def profile_plots(out,nasa_dir):
    # First user must provide raw TMR files; no synthetic velocity/RSS data.
    if not nasa_dir:return []
    dir=Path(nasa_dir)
    experimental=load_tecplot(dir/"profiles.exp.dat",11)
    velocity=load_tecplot(dir/"backstep_cfl3d_vel_sst.dat",3)
    shear=load_tecplot(dir/"backstep_cfl3d_uv_sst.dat",3)
    if not experimental or not velocity or not shear:
        return []
    result=[]
    for quantity,model,ecol,mcol,label in [("u",velocity,2,2,"u/Uref"),
                                           ("uv",shear,6,2,"1000 x u'v'/Uref^2 (CHECK scaling)")]:
        for station in (1,4,6,10):
            ez=[v for k,v in experimental.items() if f"x/H={station}" in k]
            mz=[v for k,v in model.items() if f"x/H={station}" in k]
            if not ez or not mz:continue
            ex,mod=ez[0],mz[0]
            fig,ax=plt.subplots(figsize=(6,6))
            ax.plot(ex[:,ecol],ex[:,1],"o",label="Experience")
            ax.plot(mod[:,mcol],mod[:,1],"-",label="NASA CFL3D SSTm")
            ax.set(xlabel=label,ylabel="y/H",title=f"x/H={station}: {quantity} - NASA vs exp")
            ax.grid(alpha=.25);ax.legend();fig.tight_layout()
            name=f"11_profile_{quantity}_x{station}.png";fig.savefig(out/name,dpi=180)
            plt.close(fig);result.append(name)
    return result

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-dir",default="plots_v76")
    parser.add_argument("--nasa-dir",default=None)
    parser.add_argument("--v75-log",default=None)
    parser.add_argument("--cfdx-cf",default=None,help="CSV x_over_h,cf ; actual converged CFDX data")
    parser.add_argument("--cfdx-cp",default=None,help="CSV x_over_h,cp ; actual converged CFDX data")
    parser.add_argument("--openfoam-cf",default=None,help="CSV x_over_h,cf ; actual OpenFOAM data")
    parser.add_argument("--openfoam-cp",default=None,help="CSV x_over_h,cp ; actual OpenFOAM data")
    args=parser.parse_args()
    out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
    cf_exp,cf_nasa,cp_exp,cp_nasa=CF_EXP,CF_NASA,CP_EXP_BOTTOM,CP_NASA
    raw_used=[]
    if args.nasa_dir:
        root=Path(args.nasa_dir)
        for name,ncols,attr in [("cf.exp.dat",3,"cf_exp"),
                                ("backstep_cfl3d_cf_sst.dat",2,"cf_nasa"),
                                ("cp.expnew.dat",3,"cp_exp"),
                                ("backstep_cfl3d_cp_sst.dat",2,"cp_nasa")]:
            f=load_tecplot(root/name,ncols)
            if not f:continue
            raw_used.append(name)
            if attr=="cf_exp":cf_exp=next(iter(f.values()))
            elif attr=="cf_nasa":cf_nasa=next(iter(f.values()))
            elif attr=="cp_exp":cp_exp=next(v[:,:2] for k,v in f.items() if "bottom" in k)
            elif attr=="cp_nasa":cp_nasa=next(iter(f.values()))
    values=plot_bfs(out,cf_exp,cf_nasa,cp_exp,cp_nasa)
    values["experimental_nearwall_plots"]=plot_experimental_nearwall(out)
    values.update(plot_code_curves(out,args))
    values["plot_source"]=("NASA TMR complete raw Cf/Cp" if len(raw_used)==4
        else "NASA TMR sparse numeric subset, with raw replacements where present")
    values["raw_reference_files_used"]=raw_used
    if args.v75_log:
        values["v75_budget_plots"]=plot_budgets(out,parse_budgets(args.v75_log))
    if args.nasa_dir:
        values["profile_plots"]=profile_plots(out,args.nasa_dir)
    values["plots"]=sorted(p.name for p in out.glob("*.png"))
    (out/"metrics_v76.json").write_text(json.dumps(values,indent=2))
    print(json.dumps(values,indent=2))

if __name__=="__main__":main()