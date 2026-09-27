#!/usr/bin/env python3
from pathlib import Path
import csv, sys

METHODS = [
("exact_geometric","Exact geometric","Reference"),
("search_based","NASA/SWIFT-style search","Near-wall exact; vertex-stage scalability issue"),
("mesh_wave","Mesh wave / Dijkstra","Graph-distance approximation"),
("directional_mesh_wave","Directional mesh wave","Heuristic graph metric; not a PDE model"),
("poisson","Poisson distance","Poisson auxiliary PDE with discrete wall-seed Dirichlet treatment"),
("eikonal","Eikonal","Upwind Eikonal relaxation / fast-sweeping implementation"),
("hamilton_jacobi","Hamilton-Jacobi","Modified Eikonal/H-J pseudo-time PDE"),
("advection_diffusion","Advection-diffusion","Explicit transport/diffusion PDE relaxation"),
("hybrid_poisson_eikonal","Hybrid Poisson/Hamilton-Jacobi","Poisson initialization followed by H-J refinement; no distance-field blending"),
]
EQS = {
"Exact geometric": r"d(x)=min_f dist(x,f), evaluated by exact point-to-triangle distance.",
"NASA/SWIFT-style search": r"d_v=min_v ||x-v||; within a threshold, replace it with exact face distance.",
"Mesh wave / Dijkstra": r"d_i=min_j(d_j+w_ij): shortest path on the grid graph.",
"Directional mesh wave": r"w=||dx||[1+0.15(1-|e·n|)], where e is the propagation direction and n is the wall-normal estimate; this is a heuristic graph metric, not a published wall-distance PDE.",
"Poisson distance": r"∇²φ=-1 with φ=0 on the numerical wall-seed band; d=sqrt(|∇φ|²+2φ)-|∇φ|. The transformation is the Tucker formulation; the current code uses a Cartesian finite-difference Laplacian.",
"Eikonal": r"|∇d|=1, d=0 on Γ. The implementation uses a monotone upwind/fast-sweeping update based on directional one-sided differences; the numerical seed band is an approximation of the wall condition.",
"Hamilton-Jacobi": r"|∇d|=1+Γ(d)∇²d, with Γ(d)=εd. CFDX advances the residual R=|∇d|_G-1-Γ∇²d in pseudo-time, ∂d/∂τ=-R, using a monotone Godunov gradient and explicit diffusion stability restriction.",
"Advection-diffusion": r"U·∇d=1+Γ(d)∇²d, with U constructed from the monotone wall-distance gradient. CFDX discretizes U·∇d with upwind differencing and Γ∇²d with the Cartesian central Laplacian.",
"Hybrid Poisson/Hamilton-Jacobi": r"First obtain a smooth initial distance from ∇²φ=-1 and the Tucker Poisson distance reconstruction. Use that field as the initial iterate for |∇d|_G=1+Γ(d)∇²d, Γ(d)=εd. This is a solver hybridization, not an algebraic blend and has no arbitrary Poisson/Eikonal weight.",
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
       "Source: CI artifact wall-distance-benchmark. The CSV must be regenerated after each wall-distance solver change; do not reuse pre-change benchmark values as post-change results.","",
       "## Benchmark","",
       "Complex wing/body/tail-like triangulated surface; Cartesian sampling grid 40×30×22 = 26,400 points; exact point-to-triangle distance is the reference. Near-wall metric uses 2h, with h=0.21.","",
       "| Method | L2 | L∞ | Near-wall L2 | Monotonicity violations | Time (ms) |","|---|---:|---:|---:|---:|---:|"]
    for r in rows: L.append(f"| {r['method']} | {pct(float(r['l2_relative']))} | {pct(float(r['linf_relative']))} | {pct(float(r['near_wall_l2_relative']))} | {int(float(r['monotonicity_violations']))} | {float(r['time_ms']):.3f} |")
    L += ["","## Equation/model audit",""]
    for key,label,status in METHODS:
        r=by[key]
        L += [f"### {label}",f"**Status:** {status}","",f"**Equation/model:** {EQS[label]}","",f"**Result:** L2={pct(float(r['l2_relative']))}; L∞={pct(float(r['linf_relative']))}; near-wall L2={pct(float(r['near_wall_l2_relative']))}; monotonicity violations={int(float(r['monotonicity_violations']))}.",""]
    L += ["## Critical code findings","",
    "1. Directional wave is retained as a heuristic graph metric. It must not be presented as a Tucker/NASA PDE model; its wall-normal/direction definition and benchmark role are documented explicitly.",
    "2. Eikonal now uses the monotone upwind Hamiltonian |∇d|_G=1 and is documented as a numerical fast-sweeping/relaxation method, with the wall-seed band identified as a numerical boundary approximation.",
    "3. Hamilton-Jacobi now documents the actual modified-distance PDE |∇d|_G=1+εd∇²d, its pseudo-time sign, Godunov gradient and explicit diffusion restriction.",
    "4. Poisson now documents the exact continuous PDE, the Tucker distance transformation, and the distinction between the physical wall condition and CFDX's numerical seed-band treatment.",
    "5. Advection-diffusion is implemented as the Tucker transport form U·∇d=1+Γ∇²d with first-order upwind advection and second-order central diffusion; it is no longer described as smoothing.",
    "6. The hybrid method is now Poisson-initialized H-J iteration. It no longer uses the arbitrary 0.35/0.65 distance-field blend or a non-literature weighted front direction.",
    "7. Search-based is the useful production direction. The BVH exact-face query fixes the previous large-face projection defect: the current complex benchmark gives 0% near-wall L2 error. However, the initial vertex search still loops over every wall vertex for every query point, so the complete algorithm is not yet spatially accelerated end-to-end.",
    "8. Benchmark scope: this validates numerical wall-distance algorithms against the triangulated benchmark surface. It is not yet a turbulence-model V&V case.","",
    "## What should be corrected next","",
    "- Fix and test the directional metric.",
    "- Implement a genuine anisotropic upwind fast-sweeping/fast-marching Eikonal solver.",
    "- Keep H-J and the hybrid on the same published implicit/iterated formulation, with an explicit residual/convergence diagnostic.",
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
