"""Fail-closed geometry/BC comparability of CFDX TH03 against Driver-Seegmiller.

Two independent gates: matching benchmark GEOMETRY and full PHYSICAL VALIDATION.
A geometry pass does NOT assert exact OpenFOAM polyMesh connectivity, proper wall
functions, experimental inlet profile, or a converged field. Plot is schematic,
using the coordinates read from the actual HDF5 input, not a computed solution.
"""

import argparse
import json
import re
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

H_REF = 0.0127
NU_REF = 1.56e-5
U_REF = 44.2
K_REF = 1.09e-3
OMEGA_REF = K_REF / (NU_REF * 0.009)


def gate_reference_geometry(s):
    needed = ("h_m", "h_upstream", "h_downstream", "x_min_h", "x_max_h")
    if any(k not in s for k in needed):
        raise ValueError("incomplete geometry preflight")
    vals = np.array([s[k] for k in needed], float)
    if not np.isfinite(vals).all():
        raise ValueError("nonfinite benchmark geometry")
    result = {
        "step_height_match": bool(abs(s["h_m"] / H_REF - 1) < 1e-9),
        "upstream_height_8H": bool(abs(s["h_upstream"] - 8) < 1e-7),
        "downstream_height_9H": bool(abs(s["h_downstream"] - 9) < 1e-7),
        "long_inlet_for_uniform_profile": bool(s["x_min_h"] <= -110 + 1e-7),
        "outlet_at_least_40H": bool(s["x_max_h"] >= 40 - 1e-7),
    }
    result["pass_reference_domain"] = all(result.values())
    return result


def read_case(path):
    with h5py.File(path) as f:
        p = f["points"][:]
        fo = f["face_offsets"][:]
        fv = f["face_vertices"][:]
        cf = f["cell_faces"][:]
        co = f["cell_offsets"][:]
        owners = f["owner"][:]
        off = f["patch_face_offsets"][:]
        ids = f["patch_face_ids"][:]
        names = [s.split(":")[0] for s in f.attrs["boundary_patches"].split(";")]
        patch = {name: ids[off[i] : off[i + 1]] for i, name in enumerate(names)}
        cfg = json.loads(f.attrs["case_setup_json"])
        faces = {name: [] for name in names}
        for name, face_ids in patch.items():
            for fid in face_ids:
                v = p[fv[fo[fid] : fo[fid + 1]]]
                faces[name].append(v.mean(axis=0))
        # Locate horizontal lower walls on both sides of the step.
        low_candidates = [
            name
            for name in patch
            if name in ("upstream_wall", "downstream_wall", "lowerWall")
        ]
        lower = np.concatenate([np.asarray(faces[name]) for name in low_candidates])
        upstream = lower[(lower[:, 0] < -1e-10) & (np.abs(lower[:, 1] - H_REF) < 1e-5)]
        downstream = lower[(lower[:, 0] > 1e-10) & (np.abs(lower[:, 1]) < 1e-5)]
        if len(upstream) == 0 or len(downstream) == 0:
            raise ValueError("cannot locate both horizontal lower walls")
        y_up = float(np.median(upstream[:, 1]))
        y_down = float(np.median(downstream[:, 1]))
        h = y_up - y_down
        if not h > 0:
            raise ValueError("invalid step")
        top = float(p[:, 1].max())
        # Independent geometrical proxy only (never a converged wall y+).
        d_first = []
        for name in low_candidates:
            for fid in patch[name]:
                fp = p[fv[fo[fid] : fo[fid + 1]]].mean(axis=0)
                if fp[0] <= 0 or abs(fp[1]) > 1e-8:
                    continue
                owner = int(owners[fid])
                unique = np.unique(
                    np.concatenate(
                        [fv[fo[j] : fo[j + 1]] for j in cf[co[owner] : co[owner + 1]]]
                    )
                )
                cg = p[unique].mean(axis=0)
                d_first.append(abs(float(cg[1] - fp[1])) / h)
        if not d_first:
            raise ValueError("no downstream first-cell distances")
        inlet = next(
            b for b in cfg["boundary_conditions"] if b["patch_name"] == "inlet"
        )
        material = cfg["materials"][0]
        rho = float(material["density"])
        mu = float(material["dynamic_viscosity"])
        U = float(inlet["velocity_magnitude"])
        cond = {
            x["field"]: float(x.get("value", 0))
            for x in inlet.get("sst_scalar_conditions", [])
        }
        k = cond.get("k")
        w = cond.get("omega")
        nu = mu / rho
        if not (nu > 0 and U > 0):
            raise ValueError("invalid inlet material")
        nut_nu = (
            k / (w * nu)
            if k is not None and w is not None and k >= 0 and w > 0
            else None
        )
        return {
            "case_file": str(path),
            "case_name": cfg["source"].get("case_name"),
            "n_cells": int(f.attrs["n_cells"]),
            "n_faces": int(f.attrs["n_faces"]),
            "source_format": f.attrs.get("source_format", "unknown"),
            "x_min_h": float(p[:, 0].min() / h),
            "x_max_h": float(p[:, 0].max() / h),
            "h_m": h,
            "y_upper_h": float((top - y_down) / h),
            "downstream_first_cell_d_over_h": {
                "min": min(d_first),
                "median": float(np.median(d_first)),
                "max": max(d_first),
                "method": "HDF5 owner-cell vertex mean to horizontal face; geometry only",
            },
            "h_upstream": float((top - y_up) / h),
            "h_downstream": float((top - y_down) / h),
            "expansion_ratio": float((top - y_down) / (top - y_up)),
            "rho": rho,
            "mu": mu,
            "nu": nu,
            "u_in": U,
            "re_h": float(U * h / nu),
            "k_in": k,
            "omega_in": w,
            "nut_over_nu_in": nut_nu,
            "boundary_names": names,
            "mesh_fingerprint": f.attrs.get("mesh_hash", ""),
            "case_setup": cfg,
        }


def validate_one(case):
    check = gate_reference_geometry(case)
    out = {k: v for k, v in case.items() if k != "case_setup"}
    out["geometry_checks"] = check
    out["geometry_status"] = (
        "MATCHES_REFERENCE_DOMAIN"
        if check["pass_reference_domain"]
        else "NOT_REFERENCE_DOMAIN"
    )
    out["sst_wall_law_equivalence"] = "UNVERIFIED"
    out["native_full_linear_certificates"] = "MISSING"
    out["converged_physical_state"] = False
    out["physical_validation_authorized"] = False
    out["physical_validation_reason"] = (
        "Even correct geometry requires matched BCs, native wall laws, convergence and multi-grid V&V."
    )
    return out


def plot_geometry(cases, out):
    fig, ax = plt.subplots(figsize=(11.4, 5.5))
    # Actual measured HDF5 extents are represented as scaled channel profiles.
    # We use normalized step-aligned y=0 and upstream floor y=1 for both.
    colors = ["#b55e33", "#0a78a5"]
    for n, case in enumerate(cases):
        c = colors[n % len(colors)]
        x0 = max(case["x_min_h"], -25)
        x1 = min(case["x_max_h"], 52)
        top = case["h_downstream"]
        label = (
            f"{case['n_cells']:,} cellules: top={top:g}H, amont={case['h_upstream']:g}H"
        )
        ax.plot([x0, 0, 0, x1], [1, 1, 0, 0], color=c, lw=2, alpha=0.9)
        ax.plot([x0, x1], [top, top], color=c, lw=2, alpha=0.9, label=label)
        ax.scatter([x0, x1], [top, top], color=c, s=20)
    ax.hlines(
        9,
        -25,
        52,
        linestyles="dotted",
        colors="black",
        linewidth=1.1,
        label="Référence exp. (top=9H)",
    )
    ax.axvline(0, color="gray", lw=0.7)
    ax.set(
        xlim=(-25, 52),
        ylim=(-0.8, 10.5),
        xlabel="x/H (fenêtre près de la marche)",
        ylabel="y/H (origine paroi aval)",
        title="CFDX : le pilote 152 cellules n’a pas la géométrie 8H→9H du benchmark",
    )
    ax.set_aspect("auto")
    ax.grid(alpha=0.15)
    ax.legend(loc="lower left", fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "13_geometrie_BFS152_vs_20540.png", dpi=200)
    plt.close(fig)
    # Objective turbulence inlet contrast, two logarithmic panels.
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5))
    entries = [
        ("k entrée (m²/s²)", "k_in", K_REF),
        ("ω entrée (1/s)", "omega_in", OMEGA_REF),
        ("νt/ν entrée", "nut_over_nu_in", 0.009),
    ]
    for ax, (name, key, reference) in zip(axes, entries):
        labels = [f"CFDX {c['n_cells']}" for c in cases] + ["TH03 V4 / OpenFOAM"]
        vals = [c[key] or np.nan for c in cases] + [reference]
        ax.bar(labels, vals, color=[*colors[: len(cases)], "#777777"])
        ax.set_yscale("log")
        ax.set(title=name)
        ax.tick_params(axis="x", rotation=25)
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle(
        "Conditions SST réellement écrites dans les HDF5 — ne pas confondre turbulence et géométrie"
    )
    fig.tight_layout()
    fig.savefig(out / "14_contraste_entree_sst.png", dpi=200)
    plt.close(fig)


def pressure_failure_evidence(logfile, execution, relative_target=1e-11):
    """Read only genuine native pressure-rejection evidence, never certify a run."""
    source = Path(logfile).read_text(errors="replace")
    rejected_lines = re.findall(
        r"^CFDX nonlinear iteration \d+ rejected;[^\n]*$", source, re.MULTILINE
    )
    raw = []
    for line in rejected_lines:
        m = re.search(
            r"pressure-correction solve did not converge \(status=(\d+), iterations=(\d+), residual=([0-9eE+\-.]+), relative=([0-9eE+\-.]+)\)",
            line,
        )
        if m:
            raw.append(m.groups())
    if not raw:
        raise ValueError("no native pressure-rejection evidence")
    execution_data = json.loads(Path(execution).read_text())
    if (
        execution_data.get("converged") is not False
        or execution_data.get("process_exit_code") != 2
    ):
        raise ValueError("unexpected physical solution status for pressure rejection")
    trials = []
    for j, (status, iterations, residual, relative) in enumerate(raw, 1):
        row = {
            "trial": j,
            "status": int(status),
            "iterations": int(iterations),
            "residual_absolute": float(residual),
            "residual_relative": float(relative),
        }
        if not (
            np.isfinite([row["residual_absolute"], row["residual_relative"]]).all()
            and row["residual_relative"] > 0
            and row["residual_absolute"] >= 0
        ):
            raise ValueError("nonfinite native pressure residual")
        trials.append(row)
    return {
        "schema": "CFDX-PRESSURE-FAIL-REPLAY/V78",
        "target_relative": relative_target,
        "process_exit_code": 2,
        "outer_iterations_accepted": 0,
        "converged": False,
        "trials": trials,
        "min_relative": min(r["residual_relative"] for r in trials),
        "physical_validation_authorized": False,
    }


def plot_pressure_failures(evidence, out):
    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    trials = evidence["trials"]
    y = [r["residual_relative"] for r in trials]
    ax.bar([r["trial"] for r in trials], y, color="#555e72")
    ax.axhline(
        evidence["target_relative"],
        color="#b43732",
        linestyle="--",
        label="Tolérance native 1e-11, inchangée",
    )
    ax.set(
        xlabel="Tentative de correction pression (première itération extérieure)",
        ylabel="Résidu linéaire relatif vrai",
        yscale="log",
        title="CFDX sur 20 540 cellules : 3 corrections pression rejetées, 0 état accepté",
    )
    ax.set_xticks([r["trial"] for r in trials])
    ax.grid(axis="y", alpha=0.2)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "17_bfs20540_echecs_pression_natifs.png", dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("case_h5", nargs="+", type=Path)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--pressure-log", type=Path)
    ap.add_argument("--pressure-execution", type=Path)
    args = ap.parse_args()
    o = args.out_dir
    o.mkdir(parents=True, exist_ok=True)
    cases = [validate_one(read_case(p)) for p in args.case_h5]
    data = {
        "schema": "CFDX-BFS-REFERENCE-COMPARABILITY/V78",
        "cases": cases,
        "source_reference": "NASA TMR 2DBFS: experiment 8H upstream, 9H downstream, Uref near x/H=-4",
        "geometry_only_not_physics": True,
    }
    (o / "comparabilite_geometrie_v78.json").write_text(json.dumps(data, indent=2))
    plot_geometry(cases, o)
    if bool(args.pressure_log) != bool(args.pressure_execution):
        raise ValueError("provide both pressure log and execution record")
    if args.pressure_log:
        evidence = pressure_failure_evidence(args.pressure_log, args.pressure_execution)
        (o / "pression_native_20540_v78.json").write_text(
            json.dumps(evidence, indent=2)
        )
        plot_pressure_failures(evidence, o)
    for c in cases:
        print(
            f"{c['n_cells']}: {c['geometry_status']} "
            f"up={c['h_upstream']:.6g}H down={c['h_downstream']:.6g}H "
            f"k={c['k_in']:.5g} omega={c['omega_in']:.5g} "
            f"nut/nu={c['nut_over_nu_in']:.5g} VALIDATED=0"
        )


if __name__ == "__main__":
    main()
