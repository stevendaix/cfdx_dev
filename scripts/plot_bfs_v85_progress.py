"""CFDX V85 diagnostic: physical-residual ledger + identical-scale VTU maps.

Uses PyVista when installed, or explicitly reports VTK fallback. Never treats
provisional states as converged. No inference of wall-Cf or reattachment from
Ux alone. Every accepted iteration is certified by the native COMMITTED marker.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from PIL import Image
import itertools

V83_SCRIPT = Path(__file__).resolve().parent / "plot_bfs_v83_pyvista.py"
spec = importlib.util.spec_from_file_location("cfdx_v83_vtk_pyvista", V83_SCRIPT)
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
FIELDS = ("U", "U_global", "p_continuity", "k", "omega", "merit")


def scan_ledger(path):
    trials = []
    accepted = []

    def extract(line):
        return dict(re.findall(r"([A-Za-z_]+)=([^\s]+)", line))

    for raw in path.read_text(errors="replace").splitlines():
        if raw.startswith("GLOBAL_PROGRESS_TRIAL_V85 "):
            row = extract(raw)
            trials.append(row)
        if raw.startswith("GLOBAL_PROGRESS_COMMITTED_V85 "):
            accepted.append(extract(raw))
    numbers = []
    for row in accepted:
        entry = {"iteration": int(row["iteration"])}
        for key in FIELDS:
            entry[key] = float(row[key])
            if not math.isfinite(entry[key]) or entry[key] < 0:
                raise ValueError("Nonfinite/negative native physical residual")
        if row.get("qualified") != "0" or row.get("converged") != "0":
            raise ValueError(
                "V85 visualization accepts provisional states only, no final qualification"
            )
        numbers.append(entry)
    if any(b["iteration"] <= a["iteration"] for a, b in itertools.pairwise(numbers)):
        raise ValueError("Non-increasing committed iteration records")
    return trials, numbers


def mesh_centres(h5file):
    with h5py.File(h5file) as f:
        pts = f["points"][...]
        co = f["cell_offsets"][...]
        cf = f["cell_faces"][...]
        fo = f["face_offsets"][...]
        fv = f["face_vertices"][...]
        if len(co) - 1 != 20540:
            raise ValueError("Expected exact 20540-cell qualification geometry")
        centres = np.empty((len(co) - 1, 3), np.float64)
        for i in range(len(co) - 1):
            faces = cf[co[i] : co[i + 1]]
            vertex = np.unique(np.concatenate([fv[fo[j] : fo[j + 1]] for j in faces]))
            centres[i] = pts[vertex].mean(axis=0)
    return centres


def make_figures(log, cases, h5, out, backend="auto"):
    out.mkdir(parents=True, exist_ok=True)
    trials, rows = scan_ledger(log)
    if not rows:
        raise RuntimeError("No native V85 COMMITTED records: will not fabricate plots")
    accepted_ids = {r["iteration"] for r in rows}
    vtu_files = sorted(
        cases.glob("result_*.vtu"), key=lambda p: int(p.stem.split("_")[-1])
    )
    vtu_files = [p for p in vtu_files if int(p.stem.split("_")[-1]) in accepted_ids]
    if not vtu_files:
        raise RuntimeError("No VTU matches a native COMMITTED record")
    if backend == "auto":
        try:
            import pyvista
        except ImportError:
            backend = "vtk"
        else:
            backend = "pyvista"
    if backend == "pyvista":
        import pyvista
    from vtk.util.numpy_support import numpy_to_vtk

    states = []
    for file in vtu_files:
        grid, a = base.load_vtu(file)
        if grid.GetNumberOfCells() != 20540:
            raise ValueError("Wrong native mesh")
        states.append((int(file.stem.split("_")[-1]), file, grid, a))
    centres = mesh_centres(h5)
    x, y = centres[:, 0] / 0.0127, centres[:, 1] / 0.0127
    stats = []
    for iteration, file, _, arr in states:
        ux = arr["u_x"]
        uy = arr["u_y"]
        pressure = arr["p"]
        mask = (x >= 0) & (x <= 20) & (y >= 0) & (y <= 3)
        recirc = ux < 0
        on = recirc & mask
        stats.append(
            {
                "iteration": iteration,
                "vtu": str(file),
                "sha256": base.checksum(file),
                "n_cells": len(ux),
                "negative_ux_cells": int(np.count_nonzero(recirc)),
                "backflow_fraction_near_step": float(np.mean(recirc[mask])),
                "backflow_x_min_H": float(x[on].min()) if on.any() else None,
                "backflow_x_max_H": float(x[on].max()) if on.any() else None,
                "ux_min": float(ux.min()),
                "ux_max": float(ux.max()),
                "uy_min": float(uy.min()),
                "uy_max": float(uy.max()),
                "pressure_min_pa": float(pressure.min()),
                "pressure_max_pa": float(pressure.max()),
                "k_max": float(arr["k"].max()),
                "omega_max": float(arr["omega"].max()),
                "raw_velocity_pressures_only": "provisional, not physically converged",
            }
        )
    # Consistent domain, cell fields, and scales are mandatory for comparisons.
    lim = base.color_limits([s[3] for s in states])
    lim["u_x"] = (min(-10.0, lim["u_x"][0]), lim["u_x"][1])
    near = (-2, 20, -0.2, 3.2)

    def render(grid, field, limits, cmap, path):
        raw = out / (".tmp_" + path.name)
        if backend == "vtk":
            extent = base.render_vtk(
                base._vtk_surface(grid), field, limits, cmap, near, raw, (1600, 370)
            )
        else:
            extent = base.render_pyvista(
                grid, field, limits, cmap, near, raw, (1600, 370)
            )
        fig = plt.figure(figsize=(13.3, 4.25), facecolor="white")
        ax = fig.add_axes((0.055, 0.23, 0.90, 0.57))
        ax.imshow(Image.open(raw), extent=extent, origin="upper")
        ax.set_xlim(near[0], near[1])
        ax.set_ylim(near[2], near[3])
        ax.set_aspect("auto")
        ax.plot([-2, 0, 0, 20], [1, 1, 0, 0], color="white", linewidth=1.5)
        ax.set_xlabel("x / H")
        ax.set_ylabel("y / H")
        fig.text(0.06, 0.92, path.stem.replace("_", " "), fontsize=14, weight="bold")
        cax = fig.add_axes((0.25, 0.115, 0.52, 0.035))
        cb = fig.colorbar(
            ScalarMappable(norm=Normalize(*limits), cmap=cmap),
            cax=cax,
            orientation="horizontal",
        )
        cb.ax.tick_params(labelsize=9)
        fig.text(
            0.06,
            0.03,
            "CFDX natif · etat provisoire NON CONVERGE · echelle identique entre etats · moteur "
            + backend.upper(),
            fontsize=9,
        )
        fig.savefig(path, dpi=135)
        plt.close(fig)
        raw.unlink(missing_ok=True)

    maps = []
    for i, file, grid, arr in states:
        for field, cmap in [("u_x", "RdBu_r"), ("p_kpa", "coolwarm")]:
            dest = out / f"carte_iter{i}_{field}.png"
            render(grid, field, lim[field], cmap, dest)
            maps.append(str(dest))
    if len(states) >= 2:
        a = states[0]
        b = states[-1]
        if len(a[3]["u_x"]) != len(b[3]["u_x"]):
            raise ValueError("different mesh size")
        # Fail closed when node/cell topology differs; do not silently compare
        # arrays by index from unrelated grids.
        if a[2].GetNumberOfPoints() != b[2].GetNumberOfPoints():
            raise ValueError("different node count")
        for name, key in [("delta_u_x", "u_x"), ("delta_p_kpa", "p_kpa")]:
            dv = b[3][key] - a[3][key]
            vmax = max(1e-12, float(np.abs(dv).max()))
            arr = numpy_to_vtk(np.ascontiguousarray(dv), deep=True)
            arr.SetName(name)
            b[2].GetCellData().AddArray(arr)
            p = out / f"carte_delta_iter{a[0]}_iter{b[0]}_{name}.png"
            render(b[2], name, (-vmax, vmax), "RdBu_r", p)
            maps.append(str(p))
    # Only actual accepted, committed results (not any rejected probes).
    import csv

    with (out / "native_accepted_history.csv").open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=["iteration"] + list(FIELDS))
        writer.writeheader()
        writer.writerows(rows)
    with (out / "native_backflow_and_field_ranges.csv").open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(stats[0]))
        writer.writeheader()
        writer.writerows(stats)
    xx = [r["iteration"] for r in rows]
    fig, ax = plt.subplots(figsize=(10, 5))
    for field in ("U", "U_global", "k", "omega", "merit"):
        ax.plot(xx, [r[field] for r in rows], marker="o", label=field)
    ax.set_yscale("log")
    ax.set_xlabel("Iteration exterieure effectivement engagee")
    ax.set_ylabel("Residus physiques reels (echelles V8)")
    ax.set_title("BFS20540 · V85 · etats provisoires acceptes, pas de qualif.")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "01_residus_physiques_acceptes.png", dpi=155)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(
        [s["iteration"] for s in stats],
        [s["negative_ux_cells"] for s in stats],
        "o-",
        label="Ux < 0, cellules",
    )
    ax.set_xlabel("Iteration exterieure effectivement engagee")
    ax.set_ylabel("Nombre cellules Ux<0")
    ax.set_title("Recirculation : diagnostic de signe, pas une mesure de recollement")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "02_apparition_contre_courant.png", dpi=155)
    plt.close(fig)
    manifest = {
        "case": "BFS20540",
        "status": "provisional_not_converged",
        "actual_pyvista_execution": backend == "pyvista",
        "renderer": backend,
        "native_log": str(log),
        "native_log_sha256": base.checksum(log),
        "hdf5": str(h5),
        "hdf5_sha256": base.checksum(h5),
        "accepted_iterations": [r["iteration"] for r in rows],
        "rejected_trials": sum(t.get("accepted") == "0" for t in trials),
        "states": stats,
        "scales": {k: list(v) for k, v in lim.items()},
        "maps": maps,
        "limitations": [
            "No final V8 certificate",
            "No native wall friction Cf from these VTU fields",
            "Pressure exported is not verified static pressure/Cp",
            "Ux sign alone does not locate skin-friction reattachment",
            "No field comparison performed against an unaccepted candidate",
        ],
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--native-log", type=Path, required=True)
    ap.add_argument("--states", type=Path, required=True)
    ap.add_argument("--h5", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--backend", choices=("auto", "vtk", "pyvista"), default="auto")
    arg = ap.parse_args()
    print(
        json.dumps(
            make_figures(arg.native_log, arg.states, arg.h5, arg.output, arg.backend),
            indent=2,
        )
    )
