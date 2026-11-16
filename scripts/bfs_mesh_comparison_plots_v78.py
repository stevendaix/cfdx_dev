import argparse
import csv
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

parser = argparse.ArgumentParser(
    description="Mesh comparison diagnostic, no physical qualification"
)
parser.add_argument("--pilot-case", type=Path, required=True)
parser.add_argument("--reference-case", type=Path, required=True)
parser.add_argument("--before-csv", type=Path, required=True)
parser.add_argument("--after-csv", type=Path, required=True)
parser.add_argument("--out-dir", type=Path, required=True)
args = parser.parse_args()
p = args.out_dir
p.mkdir(parents=True, exist_ok=True)
cases = [
    (str(args.pilot_case), "CFDX BFS152 (pilote diagnostic)"),
    (str(args.reference_case), "CFDX TH03 20 540 (reconstruction six blocs)"),
]
fig, axs = plt.subplots(1, 2, figsize=(14, 4.8), sharex=True, sharey=True)
for (fn, label), ax in zip(cases, axs):
    with h5py.File(fn) as f:
        pts = f["points"][:]
        fv = f["face_vertices"][:]
        fo = f["face_offsets"][:]
        cf = f["cell_faces"][:]
        co = f["cell_offsets"][:]
        polys = []
        for i in range(len(co) - 1):
            allverts = np.unique(
                np.concatenate([fv[fo[j] : fo[j + 1]] for j in cf[co[i] : co[i + 1]]])
            )
            cell = pts[allverts]
            if (
                cell[:, 0].max() < -4 * 0.0127
                or cell[:, 0].min() > 14 * 0.0127
                or cell[:, 1].min() > 2.5 * 0.0127
            ):
                continue
            lower = cell[np.isclose(cell[:, 2], cell[:, 2].min(), atol=1e-12)][:, :2]
            if len(lower) < 3:
                continue
            cen = lower.mean(0)
            a = np.argsort(np.arctan2(lower[:, 1] - cen[1], lower[:, 0] - cen[0]))
            polys.append(lower[a] / 0.0127)
        pc = PolyCollection(
            polys,
            facecolors="none",
            edgecolors="black",
            linewidths=0.55 if len(polys) < 50 else 0.2,
        )
        ax.add_collection(pc)
        ax.text(
            0.02,
            0.96,
            f"{len(polys):,} cellules visibles",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
        )
        ax.set(xlim=(-4, 14), ylim=(-0.05, 2.5), xlabel="x/H", title=label)
        ax.grid(alpha=0.1)
axs[0].set_ylabel("y/H depuis la paroi aval")
fig.suptitle(
    "Mailles effectivement reconstruites près de la marche — même cadrage spatial, tailles très différentes"
)
fig.tight_layout()
fig.savefig(p / "15_densite_maillage_pres_marche.png", dpi=210)
plt.close(fig)

original = args.before_csv
correction = args.after_csv
rows = []
for file in (original, correction):
    with file.open(newline="") as stream:
        rows.append(list(csv.DictReader(stream)))
fig, axs = plt.subplots(1, 3, figsize=(13.7, 5), sharey=True)
for a, station in zip(axs, [4.0, 6.0, 10.0]):
    old = [r for r in rows[0] if float(r["x_over_h"]) == station]
    new = [r for r in rows[1] if float(r["x_over_h"]) == station]
    ex = np.array([(float(r["exp_u_u_ref"]), float(r["y_over_h"])) for r in new])
    uold = np.array(
        [
            (float(r["cfdx_u_interp"]), float(r["y_over_h"]))
            for r in old
            if r["cfdx_u_interp"]
        ]
    )
    unew = np.array(
        [
            (float(r["cfdx_u_interp"]), float(r["y_over_h"]))
            for r in new
            if r["cfdx_u_interp"]
        ]
    )
    a.plot(ex[:, 0], ex[:, 1], "o", ms=4, label="Expérience NASA")
    a.plot(
        uold[:, 0],
        uold[:, 1],
        "--",
        marker="x",
        ms=4,
        label="CFDX V77 (centres VTK biaisés)",
    )
    a.plot(unew[:, 0], unew[:, 1], "-s", ms=3, label="CFDX corrigé (centres HDF5)")
    a.axvline(0, color="gray", lw=0.7)
    a.set(
        xlabel="U/Uref (réf. non strictement alignées)",
        title=f"x/H={station}",
        ylim=(0, 1.3),
    )
    a.grid(alpha=0.18)
axs[0].set_ylabel("y/H")
axs[-1].legend(fontsize=7, loc="upper left")
fig.suptitle(
    "BFS152 non convergé : influence de la correction des centres d’interpolation — PAS une validation NASA"
)
fig.tight_layout()
fig.savefig(p / "16_v77_vs_v78_profils_vitesse.png", dpi=200)
plt.close(fig)
