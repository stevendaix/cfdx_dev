"""CFDX V82: compare real first private BFS20540 SST transactions (NOT Cf/Cp)."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def tokens(line):
    out = {}
    for p in line.split()[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            out[k] = v
    return out


def finite(row, key):
    v = float(row[key])
    if not math.isfinite(v):
        raise ValueError("nonfinite " + key)
    return v


def parse_private_cycle(path, *, implicit=False):
    lines = Path(path).read_text().splitlines()
    begin = [
        i
        for i, l in enumerate(lines)
        if l.startswith("BFS_PHASE_V10") and "phase=SST_TRANSPORT_BEGIN" in l
    ]
    end = [
        i
        for i, l in enumerate(lines)
        if l.startswith("BFS_PHASE_V10") and "phase=SST_TRANSPORT_END" in l
    ]
    if not begin or not end or end[0] <= begin[0]:
        raise ValueError("missing native first complete SST transaction")
    header = tokens(lines[begin[0]])
    finish = tokens(lines[end[0]])
    if int(header["budget"]) != 10 or finite(header, "target") != 1e-8:
        raise ValueError("BFS V82 comparison requires identical unrelaxed SST settings")
    if int(finish["committed"]) != 0 or int(finish["converged"]) != 0:
        raise ValueError(
            "expected a rejected private SST transaction; do not conceal accepted results"
        )
    if int(finish["iterations"]) > 10 or int(finish["iterations"]) < 1:
        raise ValueError("invalid SST iteration count")
    prefix = lines[: end[0] + 1]
    pressure = [
        tokens(l) for l in prefix if l.startswith("PRESSURE_HILO_PHYSICAL_FLUX")
    ]
    if not pressure or finite(pressure[0], "stored_fp64_relative") > 1e-11:
        raise ValueError("uncertified pressure/continuity on the same first corrector")
    momentum = [tokens(l) for l in prefix if l.startswith("INCOMPRESSIBLE_ITER")]
    if not momentum:
        raise ValueError("missing native momentum residual")
    samples = {}
    for l in prefix:
        if l.startswith("SST_TRANSPORT_DIAGNOSTIC"):
            d = tokens(l)
            i = int(d["iteration"])
            samples[i] = {
                k: finite(d, k)
                for k in (
                    "k_returned_row_scaled",
                    "omega_returned_row_scaled",
                    "k_local_relative",
                    "omega_local_relative",
                )
            }
    if not samples or 1 not in samples:
        raise ValueError("no independent reassembled SST sample")
    updates = []
    for l in prefix:
        if l.startswith("SST_IMPLICIT_V82"):
            d = tokens(l)
            if (
                d.get("fixed_rhie_chow") != "1"
                or d.get("outer_U_p_flux_qualified") != "0"
                or d.get("physical_qualified") != "0"
            ):
                raise ValueError("false certification from an inner candidate")
            u = dict(
                iteration=int(d["iteration"]),
                selected=d["selected"],
                accepted=int(d["accepted"]),
                legacy_picard_fallback=int(d.get("legacy_picard_fallback", "0")),
                evaluations=int(d["evaluations"]),
                k_before=finite(d, "k_before"),
                k_after=finite(d, "k_after"),
                omega_before=finite(d, "omega_before"),
                omega_after=finite(d, "omega_after"),
                merit_before=finite(d, "old_merit"),
                merit_after=finite(d, "new_merit"),
            )
            if (
                (u["accepted"] and not (u["merit_after"] < u["merit_before"]))
                or (not u["accepted"] and u["selected"] != "none")
                or (u["accepted"] and u["legacy_picard_fallback"])
            ):
                raise ValueError("incoherent V82 inner acceptance")
            updates.append(u)
    if implicit and not updates:
        raise ValueError("no V82 coupled source trace")
    if not implicit and updates:
        raise ValueError("baseline accidentally enabled V82")
    return {
        "case": str(path),
        "pressure_fp64": finite(pressure[0], "stored_fp64_relative"),
        "momentum_relative": finite(momentum[0], "momentum_residual_rhs_relative"),
        "sst_iterations": int(finish["iterations"]),
        "sst_committed": 0,
        "sst_converged": 0,
        "samples": samples,
        "updates": updates,
        "k_final_raw": finite(finish, "k_residual"),
        "omega_final_raw": finite(finish, "omega_residual"),
    }


def build_plots(base, implicit, directory):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for key, label, filename in [
        ("k_returned_row_scaled", "Rk FV row-scaled [m2/s2]", "01_k_reassembled.png"),
        (
            "omega_returned_row_scaled",
            "Romega FV row-scaled [1/s]",
            "02_omega_reassembled.png",
        ),
        ("k_local_relative", "Rk local relative [-]", "03_k_local.png"),
        ("omega_local_relative", "Romega local relative [-]", "04_omega_local.png"),
    ]:
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for data, style, name in [
            (base, "o-", "V80 Picard"),
            (implicit, "s-", "V82 implicit + Picard filter"),
        ]:
            xs = sorted(data["samples"])
            ys = [data["samples"][x][key] for x in xs]
            ax.plot(xs, ys, style, label=name)
        ax.set_yscale("symlog", linthresh=1e-10)
        ax.set_xlabel("Sous-cycle SST prive")
        ax.set_ylabel(label)
        ax.set_title("BFS20540 - premiere tentative uniquement (NON CONVERGEE)")
        ax.legend()
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(directory / filename, dpi=155)
        plt.close(fig)
    ups = implicit["updates"]
    if ups:
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(
            [u["iteration"] for u in ups],
            [u["merit_before"] for u in ups],
            "o--",
            label="avant",
        )
        ax.plot(
            [u["iteration"] for u in ups],
            [u["merit_after"] for u in ups],
            "o-",
            label="apres",
        )
        ax.set_yscale("log")
        ax.set_xlabel("Sous-cycle SST prive")
        ax.set_ylabel("Max(k, omega) local, reassemble")
        ax.set_title("V82 : recherche non lineaire sur l operateur FV reel")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(directory / "05_v82_merit_private.png", dpi=155)
        plt.close(fig)
    payload = {
        "baseline": base,
        "v82": implicit,
        "verdict": "unconverged_inner_and_outer; no qualifying Cf Cp U uv",
        "note": "V82 sample arrays may be sparse; do not interpolate unreported subcycles.",
    }
    (directory / "evidence_v82.json").write_text(json.dumps(payload, indent=2))
    with (directory / "reassembled_samples_v82.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "source",
                "iteration",
                "k_row_scaled",
                "omega_row_scaled",
                "k_local",
                "omega_local",
            ]
        )
        for name, data in [("picard", base), ("v82", implicit)]:
            for i, point in sorted(data["samples"].items()):
                w.writerow(
                    [
                        name,
                        i,
                        point["k_returned_row_scaled"],
                        point["omega_returned_row_scaled"],
                        point["k_local_relative"],
                        point["omega_local_relative"],
                    ]
                )
    return payload


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--picard", type=Path, required=True)
    ap.add_argument("--v82", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    payload = build_plots(
        parse_private_cycle(args.picard),
        parse_private_cycle(args.v82, implicit=True),
        args.out,
    )
    print(
        "BFS20540 comparison recorded, outer_accepted=0, samples",
        len(payload["baseline"]["samples"]),
        len(payload["v82"]["samples"]),
    )


if __name__ == "__main__":
    main()
