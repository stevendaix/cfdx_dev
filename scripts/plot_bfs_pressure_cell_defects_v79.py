"""Compare native signed pressure face defects per cell before/after one step.

Mesh plotting uses the average of each cell's unique vertices as a nominal
location; it is NOT a certified FV polyhedron centroid. No velocity/force/Cf
or SST solution is inferred from this pre-commit pressure probe.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def load_defect(path, n):
    with Path(path).open() as f:
        if f.readline().strip() != "CFDX_PRESSURE_FACE_DEFECT_V79":
            raise ValueError("invalid native physical defect file")
        parts = f.readline().split()
        if len(parts) != 2 or int(parts[0]) != n or int(parts[1]) <= 0:
            raise ValueError("defect belongs to a different mesh/matrix")
        data = np.loadtxt(f, ndmin=1)
    if data.shape != (n,) or not np.all(np.isfinite(data)):
        raise ValueError("incomplete or nonfinite physical cell defect")
    return data


def mesh_nominal_cell_locations(path):
    import h5py

    with h5py.File(path) as f:
        pts = f["points"][:]
        cells = f["cell_faces"][:]
        offsets = f["cell_offsets"][:]
        fo = f["face_offsets"][:]
        face_vertices = f["face_vertices"][:]
        centroids = np.zeros((len(offsets) - 1, 3))
        for i in range(len(centroids)):
            faces = cells[offsets[i] : offsets[i + 1]]
            vertices = np.unique(
                np.concatenate([face_vertices[fo[k] : fo[k + 1]] for k in faces])
            )
            centroids[i] = pts[vertices].mean(axis=0)
    return centroids


def run(case, before, after, output_dir, h=0.0127):
    import matplotlib.pyplot as plt

    centers = mesh_nominal_cell_locations(case)
    n = len(centers)
    b = load_defect(before, n)
    a = load_defect(after, n)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rep = {
        "schema": "CFDX-PRESSURE-CELL-DEFECTS-V79/1",
        "cell_count": n,
        "before_l2": float(np.linalg.norm(b)),
        "after_l2": float(np.linalg.norm(a)),
        "before_linf": float(np.max(np.abs(b))),
        "after_linf": float(np.max(np.abs(a))),
        "before_global_sum": float(b.sum()),
        "after_global_sum": float(a.sum()),
        "max_before_cell": int(np.argmax(np.abs(b))),
        "max_after_cell": int(np.argmax(np.abs(a))),
        "native_linear_certificates": False,
        "production_commit_authorized": False,
    }
    (out / "cell_defect_evidence_v79.json").write_text(json.dumps(rep, indent=2) + "\n")
    with (out / "signed_pressure_face_defects_v79.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            (
                "cell",
                "x_over_h_nominal",
                "y_over_h_nominal",
                "before_kg_s",
                "after_kg_s",
            )
        )
        writer.writerows(
            (i, centers[i, 0] / h, centers[i, 1] / h, b[i], a[i]) for i in range(n)
        )
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 4.8), sharex=True, sharey=True)
    mask = (centers[:, 0] / h >= -2) & (centers[:, 0] / h <= 15)
    all_vals = np.log10(np.maximum(np.abs(np.concatenate((b, a))), 1e-20))
    vmin = float(np.percentile(all_vals, 2))
    vmax = float(np.max(all_vals))
    last = None
    for ax, vec, title in zip(
        axes, (b, a), ("Before face-defect correction", "After face-defect correction")
    ):
        last = ax.scatter(
            centers[mask, 0] / h,
            centers[mask, 1] / h,
            c=np.log10(np.maximum(np.abs(vec[mask]), 1e-20)),
            s=4,
            marker="s",
            vmin=vmin,
            vmax=vmax,
            rasterized=True,
        )
        ax.set_title(title)
        ax.set_xlabel("x / H")
        ax.set_xlim(-2, 15)
        ax.set_ylim(0, 9)
        ax.grid(alpha=0.15)
    axes[0].set_ylabel("y / H")
    fig.suptitle(
        "CFDX V79 - native first pressure corrector, 20,540 cells (not committed)"
    )
    fig.subplots_adjust(top=0.86, bottom=0.13, left=0.06, right=0.88, wspace=0.08)
    cax = fig.add_axes([0.90, 0.20, 0.016, 0.62])
    cb = fig.colorbar(last, cax=cax)
    cb.set_label("log10(abs(native face defect, kg/s))")
    fig.savefig(out / "03_native_cell_face_defect_map_v79.png", dpi=170)
    plt.close(fig)
    return rep


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", required=True)
    p.add_argument("--before", required=True)
    p.add_argument("--after", required=True)
    p.add_argument("--output-dir", required=True)
    args = p.parse_args()
    print(
        json.dumps(run(args.case, args.before, args.after, args.output_dir), indent=2)
    )
