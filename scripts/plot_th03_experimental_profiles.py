"""Plot CFDX TH-03 mean-velocity profiles against raw Driver--Seegmiller data.

The plots are diagnostic only. The script never promotes a run to a physical
qualification, and it refuses to plot a run whose nonlinear convergence
certificate is absent or false unless --include-unconverged is explicit.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from scripts.qualify_bfs_nasa_fields_v16 import read_nasa_profiles
except ImportError:
    from qualify_bfs_nasa_fields_v16 import read_nasa_profiles

STATIONS = (1.0, 4.0, 6.0, 10.0)
REQUIRED = {"x_m", "y_m", "u_x_m_s", "u_y_m_s"}


class ProfilePlotError(ValueError):
    pass


def read_cells(path: Path) -> np.ndarray:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or not REQUIRED.issubset(reader.fieldnames):
            raise ProfilePlotError(
                f"{path}: missing columns {sorted(REQUIRED - set(reader.fieldnames or []))}"
            )
        rows = []
        for row in reader:
            try:
                values = [
                    float(row[key]) for key in ("x_m", "y_m", "u_x_m_s", "u_y_m_s")
                ]
            except (TypeError, ValueError) as exc:
                raise ProfilePlotError(f"{path}: non-numeric cell field") from exc
            if not all(math.isfinite(value) for value in values):
                raise ProfilePlotError(f"{path}: non-finite cell field")
            rows.append(values)
    data = np.asarray(rows, dtype=float)
    if data.ndim != 2 or len(data) < 8:
        raise ProfilePlotError(f"{path}: insufficient cell samples")
    return data


def nearest_plane(
    data: np.ndarray, x_target: float, tolerance: float
) -> tuple[np.ndarray, float]:
    x_values = np.unique(data[:, 0])
    x_selected = float(x_values[np.argmin(np.abs(x_values - x_target))])
    offset = abs(x_selected - x_target)
    if offset > tolerance:
        raise ProfilePlotError(
            f"no cell-centre plane near x/H={x_target:g}; nearest offset={offset:g} m"
        )
    plane = data[
        np.isclose(
            data[:, 0], x_selected, rtol=0.0, atol=1e-12 * max(1.0, abs(x_selected))
        )
    ]
    order = np.argsort(plane[:, 1])
    plane = plane[order]
    if len(plane) < 4 or np.any(np.diff(plane[:, 1]) <= 0):
        raise ProfilePlotError(
            f"insufficient or duplicate y samples at x={x_selected:g}"
        )
    return plane, x_selected


def interp_profile(plane: np.ndarray, y: np.ndarray, component: int) -> np.ndarray:
    if y.min() < plane[:, 1].min() or y.max() > plane[:, 1].max():
        raise ProfilePlotError(
            "experimental y/H range lies outside CFDX sampled profile"
        )
    return np.interp(y, plane[:, 1], plane[:, component])


def load_run(run_dir: Path, label: str, *, include_unconverged: bool) -> dict:
    convergence_path = run_dir / "convergence.json"
    cells_path = run_dir / "benchmark_observables" / "cell_fields.csv"
    if not convergence_path.is_file() or not cells_path.is_file():
        raise ProfilePlotError(
            f"{label}: missing convergence.json or benchmark_observables/cell_fields.csv"
        )
    convergence = json.loads(convergence_path.read_text(encoding="utf-8"))
    converged = convergence.get("converged") is True
    if not converged and not include_unconverged:
        raise ProfilePlotError(
            f"{label}: run is not converged; use --include-unconverged for diagnostic-only plots"
        )
    return {"label": label, "converged": converged, "data": read_cells(cells_path)}


def make_plots(
    reference_dir: Path,
    runs: list[tuple[str, Path]],
    output: Path,
    *,
    height: float,
    step_x: float,
    plane_tolerance_h: float,
    include_unconverged: bool,
) -> dict:
    reference = read_nasa_profiles(reference_dir / "profiles.exp.dat")
    output.mkdir(parents=True, exist_ok=True)
    loaded = [
        load_run(path, label, include_unconverged=include_unconverged)
        for label, path in runs
    ]
    if not loaded:
        raise ProfilePlotError("no CFDX runs supplied")

    # Uref is evaluated at the upstream reference station and channel centre.
    urefs = {}
    for run in loaded:
        plane, _ = nearest_plane(
            run["data"], step_x - 4.0 * height, plane_tolerance_h * height
        )
        uref = float(np.interp(4.0 * height, plane[:, 1], plane[:, 2]))
        if not math.isfinite(uref) or uref <= 0:
            raise ProfilePlotError(f"{run['label']}: invalid upstream centreline Uref")
        urefs[run["label"]] = uref

    manifest = {
        "schema": "CFDX-TH03-PROFILE-PLOTS/v1",
        "physical_qualification": "NOT_QUALIFIED",
        "reference": "Driver & Seegmiller (1985), NASA TMR profiles.exp.dat",
        "stations_x_over_h": list(STATIONS),
        "runs": [],
        "plots": [],
        "limitations": [
            "The plot uses nearest cell-centre x planes; offsets are recorded.",
            "Only mean U and V are compared. Reynolds-stress tensor components are not exported by the current production observable CSV.",
            "No interpolation extrapolation is allowed in y.",
            "A plot is not a convergence or physical-validation certificate.",
        ],
    }

    for station in STATIONS:
        fig, axes = plt.subplots(1, 2, figsize=(11, 6), constrained_layout=True)
        ref = reference[station]
        y_ref = np.asarray([row["y_h"] for row in ref])
        axes[0].plot(
            [row["u_ref"] for row in ref], y_ref, "ko", ms=3, label="Driver–Seegmiller"
        )
        axes[1].plot(
            [row["v_ref"] for row in ref], y_ref, "ko", ms=3, label="Driver–Seegmiller"
        )
        station_runs = []
        for run in loaded:
            plane, x_actual = nearest_plane(
                run["data"], step_x + station * height, plane_tolerance_h * height
            )
            y_h = plane[:, 1] / height
            u = plane[:, 2] / urefs[run["label"]]
            v = plane[:, 3] / urefs[run["label"]]
            axes[0].plot(u, y_h, "-", lw=1.4, label=run["label"])
            axes[1].plot(v, y_h, "-", lw=1.4, label=run["label"])
            station_runs.append(
                {
                    "label": run["label"],
                    "converged": run["converged"],
                    "x_over_h_requested": station,
                    "x_over_h_sampled": (x_actual - step_x) / height,
                    "x_offset_over_h": abs((x_actual - step_x) / height - station),
                    "u_ref_m_s": urefs[run["label"]],
                    "n_cells_in_plane": len(plane),
                }
            )
        axes[0].set(
            xlabel="U / Uref",
            ylabel="y / H",
            title=f"Mean streamwise velocity, x/H={station:g}",
        )
        axes[1].set(
            xlabel="V / Uref",
            ylabel="y / H",
            title=f"Mean transverse velocity, x/H={station:g}",
        )
        for axis in axes:
            axis.grid(True, alpha=0.3)
            axis.legend(fontsize=7)
        image_path = output / f"profiles_xH_{station:g}.png"
        fig.savefig(image_path, dpi=180)
        plt.close(fig)
        manifest["plots"].append(str(image_path))
        manifest["runs"].append({"station_x_over_h": station, "samples": station_runs})

    (output / "profile-plot-manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument(
        "--run",
        action="append",
        nargs=2,
        metavar=("LABEL", "OUTPUT_DIR"),
        required=True,
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--height", type=float, default=0.0127)
    parser.add_argument("--step-x", type=float, default=0.0)
    parser.add_argument("--plane-tolerance-h", type=float, default=1.3)
    parser.add_argument("--include-unconverged", action="store_true")
    args = parser.parse_args()
    if not all(
        math.isfinite(x) and x > 0 for x in (args.height, args.plane_tolerance_h)
    ):
        parser.error("height and plane tolerance must be positive and finite")
    try:
        report = make_plots(
            args.reference_dir,
            [(label, Path(path)) for label, path in args.run],
            args.output_dir,
            height=args.height,
            step_x=args.step_x,
            plane_tolerance_h=args.plane_tolerance_h,
            include_unconverged=args.include_unconverged,
        )
    except (OSError, ValueError, ProfilePlotError, KeyError, TypeError) as exc:
        parser.exit(2, f"TH-03 profile plotting failed: {exc}\n")
    print(f"Wrote {len(report['plots'])} diagnostic profile plots to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
