#!/usr/bin/env python3
from pathlib import Path
import csv, sys

METHODS = [
("exact_geometric","Exact geometric","Reference"),
("search_based","NASA/SWIFT-style search","Near-wall exact; vertex-stage scalability issue"),
("mesh_wave","Mesh wave / Dijkstra","Graph-distance approximation"),
("directional_mesh_wave","Directional mesh wave","BUG: current direction factor is identically 1"),
("poisson","Poisson distance","Discrete PDE inconsistent with stated formulation"),
("eikonal","Eikonal","Not a genuine Eikonal solver"),
("hamilton_jacobi","Hamilton-Jacobi","Not a genuine H-J solver"),
("advection_diffusion","Advection-diffusion","Currently only Eikonal relaxation + smoothing"),
("hybrid_poisson_eikonal","Hybrid Poisson/Eikonal","Empirical blend, not literature-derived"),
]
EQS = {
"Exact geometric": r"d(x)=min_f dist(x,f), evaluated by exact point-to-triangle distance.",
"NASA/SWIFT-style search": r"d_v=min_v ||x-v||; within a threshold, replace it with exact face distance.",
"Mesh wave / Dijkstra": r"d_i=min_j(d_j+w_ij): shortest path on the grid graph.",
"Directional mesh wave": r"w=||dx||[1+0.15(1-A)], A=|dx·(-dx)|/||dx||²=1 in the current code, hence w=||dx||.",
"Poisson distance": r"∇²φ=-1, φ=0 at walls; d=-|∇φ|+sqrt(|∇φ|²+2φ).",
"Eikonal": r"|∇d|=1, d=0 on the wall; a proper discretization needs an upwind Eikonal update.",
"Hamilton-Jacobi": r"A modified wall-distance form is |∇d|=1+ε d ∇²d, or an equivalent advection-diffusion formulation.",
"Advection-diffusion": r"Tucker et al. use transport/advection-diffusion forms equivalent to Eikonal/H-J; the current code only smooths a graph distance.",
"Hybrid Poisson/Eikonal": r"Current implementation: d=0.35 d_P + 0.65 d_E; this is an empirical blend.",
}
REFS = [
"1. Tucker, Rumsey, Spalart, Bartels & Biedron, Computations of Wall Distances Based on Differential Equations, AIAA 2004-2232 / AIAA Journal 43(3), 539-549 (2005), DOI 10.2514/1.8626.",
"2. Zhang, Diskin, Walden & Nielsen, A Wall-Distance Method for Turbulence Modeling, ICCFD12 (2024), NASA NTRS 20240006379 (SWIFT/search-based method).",
"3. Zhang, Diskin, Walden & Nielsen, Improvements in Efficiency and Scalability of Wall-Distance Computations for Turbulence Modeling (2025), NASA NTRS 20250005253 (SWIFT-1).",
]

def load(path):
    with open(path,newline='') as f: return list(csv.DictReader(f))
def pct(x): return f"{100*x:.3g}%"

def main():
    src=Path(sys.argv[1]) if len(sys.argv)>1 else Path('wall_distance_complex_benchmark.csv')
    out=Path(sys.argv[2]) if len(sys.argv)>2 else Path('wall_distance_audit.md')
    rows=load(src); by={r['method']:r for r in rows}
    L=["# CFDX wall-distance: equation, code and result audit","",
       "Source: CI artifact wall-distance-benchmark, PR #467, commit f11228efe15c69994deda6a3583997eb1e7934b6.","",
       "## Benchmark","",
       "Complex wing/body/tail-like triangulated surface; Cartesian sampling grid 40×30×22 = 26,400 points; exact point-to-triangle distance is the reference. Near-wall metric uses 2h, with h=0.21.","",
       "| Method | L2 | L∞ | Near-wall L2 | Monotonicity violations | Time (ms) |","|---|---:|---:|---:|---:|---:|"]
    for r in rows: L.append(f"| {r['method']} | {pct(float(r['l2_relative']))} | {pct(float(r['linf_relative']))} | {pct(float(r['near_wall_l2_relative']))} | {int(float(r['monotonicity_violations']))} | {float(r['time_ms']):.3f} |")
    L += ["","## Equation/model audit",""]
    for key,label,status in METHODS:
        r=by[key]
        L += [f"### {label}",f"**Status:** {status}","",f"**Equation/model:** {EQS[label]}","",f"**Result:** L2={pct(float(r['l2_relative']))}; L∞={pct(float(r['linf_relative']))}; near-wall L2={pct(float(r['near_wall_l2_relative']))}; monotonicity violations={int(float(r['monotonicity_violations']))}.",""]
    L += ["## Critical code findings","",
    "1. Directional wave is exactly the same as mesh wave. Because g.points[u]-g.points[v] = -dv, the implemented alignment is identically 1. The 0.15 correction therefore never activates. The identical benchmark values confirm it.",
    "2. Eikonal and H-J are not PDE implementations. Both call the same eikonal_fast_sweep; only relaxation changes (1.0 versus 0.65). The kernel is neighbor relaxation d_i <- min(d_i,d_j+h), i.e. graph propagation, not an upwind Eikonal/H-J discretization.",
    "3. Poisson is discretized inconsistently. The source uses h=min(spacing) instead of directional metrics; the seed band is imposed as interior zero-Dirichlet data; the gradient reconstruction sums dp*q/q² over ± neighbors, giving twice the centered derivative on a uniform grid; outer boundary treatment is implicit through missing neighbors rather than an explicit Neumann condition.",
    "4. Advection-diffusion is only post-processing. It computes the same graph-based Eikonal field and applies local averaging. No advection-diffusion equation is solved.",
    "5. Hybrid is empirical. The 0.35/0.65 weights have no derivation or literature reference in the code. Its 9.59% near-wall error is therefore a benchmark result for this blend, not evidence for a validated hybrid equation.",
    "6. Search-based is the useful production direction. The BVH exact-face query fixes the previous large-face projection defect: the current complex benchmark gives 0% near-wall L2 error. However, the initial vertex search still loops over every wall vertex for every query point, so the complete algorithm is not yet spatially accelerated end-to-end.",
    "7. Benchmark scope: this validates numerical wall-distance algorithms against the triangulated benchmark surface. It is not yet a turbulence-model V&V case.","",
    "## What should be corrected next","",
    "- Fix and test the directional metric.",
    "- Implement a genuine anisotropic upwind fast-sweeping/fast-marching Eikonal solver.",
    "- Implement H-J separately, including its diffusion/modified-distance term and a convergence residual.",
    "- Rebuild Poisson with the correct anisotropic Laplacian, wall Dirichlet BC, outer Neumann BC and consistent gradient reconstruction.",
    "- Implement the actual transport/advection-diffusion formulation from Tucker et al.; remove the smoothing surrogate.",
    "- Keep search-based as the exact near-wall reference and accelerate its vertex stage with a spatial index; then measure MPI scalability.",
    "- Rerun refinement, skew/non-orthogonal, imported-surface and MPI campaigns after these corrections.","",
    "## References",""] + REFS + [""]
    out.write_text("\n".join(L),encoding='utf-8')
    try:
        import matplotlib.pyplot as plt
        labels=[r['method'] for r in rows]; x=list(range(len(rows)))
        fig,ax=plt.subplots(figsize=(12,6))
        ax.bar([i-.2 for i in x],[100*float(r['l2_relative']) for r in rows],.4,label='L2 (%)')
        ax.bar([i+.2 for i in x],[100*float(r['near_wall_l2_relative']) for r in rows],.4,label='near-wall L2 (%)')
        ax.set_xticks(x,labels,rotation=35,ha='right'); ax.set_ylabel('Relative error (%)'); ax.legend(); fig.tight_layout()
        fig.savefig(out.with_name('wall_distance_errors.png'),dpi=160); plt.close(fig)
        fig,ax=plt.subplots(figsize=(12,6)); ax.bar(labels,[float(r['time_ms']) for r in rows])
        ax.set_xticks(x,labels,rotation=35,ha='right'); ax.set_ylabel('Time (ms)'); fig.tight_layout()
        fig.savefig(out.with_name('wall_distance_timing.png'),dpi=160); plt.close(fig)
    except ImportError:
        pass

if __name__=='__main__':
    main()
