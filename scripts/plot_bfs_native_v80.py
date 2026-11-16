"""V80 evidence from native pressure & full U,p,k,omega logs, strictly fail closed.

A pressure certificate never stands in for an accepted or converged CFD state.
Only GLOBAL_NONLINEAR_COMMITTED_V13 marks an accepted global transaction.
"""

import argparse
import json
import math
import re
from pathlib import Path

TAG = re.compile(r"^([A-Z_0-9]+)\s+(.*)$")
KV = re.compile(r"(\w+)=([^\s]+)")
REQUIRED = (
    "PRESSURE_NATIVE_CSR_V80",
    "PRESSURE_NATIVE_FACE_V80",
    "PRESSURE_HILO_PHYSICAL_FLUX",
)


def parse(text):
    records = {}
    for line in text.splitlines():
        match = TAG.match(line.strip())
        if match:
            tag, rest = match.groups()
            records.setdefault(tag, []).append(dict(KV.findall(rest)))
    if any(not records.get(tag) for tag in REQUIRED):
        raise ValueError("incomplete V80 native pressure evidence")
    csr, face, transport = (records[tag][0] for tag in REQUIRED)
    for sample in (csr, face):
        if sample.get("production_commit_authorized") != "0":
            raise ValueError("a pressure-only row cannot authorize CFD commit")
    required_numeric = [
        (csr, "initial_relative"),
        (csr, "final_relative"),
        (face, "initial_relative"),
        (face, "final_relative"),
        (transport, "exact_relative"),
        (transport, "stored_fp64_relative"),
        (face, "requested_tolerance"),
    ]
    for row, key in required_numeric:
        try:
            value = float(row[key])
        except (KeyError, ValueError) as exc:
            raise ValueError("missing or invalid V80 numeric field: " + key) from exc
        if not math.isfinite(value) or value < 0:
            raise ValueError("nonfinite/negative evidence: " + key)
    tol = float(face["requested_tolerance"])
    if not 0 < tol <= 1e-11:
        raise ValueError("pressure tolerance weakened or missing")
    # CSR and physical face are deliberately DIFFERENT operators after the
    # final Newton correction. Require their own independently recomputed gates.
    pressure_closed = (
        float(csr["final_relative"]) <= tol
        and face.get("true_face_converged") == "1"
        and float(face["final_relative"]) <= tol
        and float(transport["exact_relative"]) <= tol
        and float(transport["stored_fp64_relative"]) <= tol
    )
    accepted = []
    for row in records.get("GLOBAL_NONLINEAR_COMMITTED_V13", []):
        if row.get("qualified") != "0" or row.get("iteration") is None:
            raise ValueError("invalid global accepted marker")
        accepted.append(int(row["iteration"]))
    accepted = sorted(set(accepted))
    if accepted and not records.get("GLOBAL_NONLINEAR_V8"):
        raise ValueError("no independent global physical residual for commit")
    has_sst = bool(records.get("BFS_PHASE_V10")) and any(
        row.get("phase") == "SST_TRANSPORT_BEGIN" for row in records["BFS_PHASE_V10"]
    )
    # SST must complete with committed native certificate before outer acceptance.
    sst_rows = [
        row
        for row in records.get("BFS_PHASE_V10", [])
        if row.get("phase") == "SST_TRANSPORT_END"
    ]
    if accepted and (
        not has_sst
        or not sst_rows
        or not any(
            x.get("committed") == "1" and x.get("converged") == "1" for x in sst_rows
        )
    ):
        raise ValueError("no valid SST native transport certificate")
    return {
        "schema": "CFDX-V80-NATIVE-EVIDENCE/1",
        "pressure_native_verified": bool(pressure_closed),
        "requested_tolerance": tol,
        "csr_initial_relative": float(csr["initial_relative"]),
        "csr_after_highlow_relative": float(csr["final_relative"]),
        "face_before_relative": float(face["initial_relative"]),
        "face_after_relative": float(face["final_relative"]),
        "full_face_after_relative": float(transport["exact_relative"]),
        "stored_fp64_after_relative": float(transport["stored_fp64_relative"]),
        "sst_stage_observed": has_sst,
        "sst_transport_complete": bool(sst_rows),
        "accepted_outer_iterations": accepted,
        "converged_cfd": False,
        "physical_comparison_qualified": False,
        "final_status": "accepted_nonconverged_state"
        if accepted
        else "pressure_only_no_global_commit",
    }


def plot(data, out_dir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    stages = [
        "CSR initial",
        "CSR high/low",
        "RC initial",
        "RC refined",
        "RC full state",
        "Flux FP64",
    ]
    vals = [
        data["csr_initial_relative"],
        data["csr_after_highlow_relative"],
        data["face_before_relative"],
        data["face_after_relative"],
        data["full_face_after_relative"],
        data["stored_fp64_after_relative"],
    ]
    fig, ax = plt.subplots(figsize=(10.6, 5.1))
    ax.bar(range(len(vals)), vals)
    ax.set_yscale("log")
    ax.axhline(
        data["requested_tolerance"], ls="--", lw=1.6, label="Original gate 1e-11"
    )
    ax.set_xticks(range(len(vals)), stages, rotation=20, ha="right")
    ax.set_ylabel("Relative residual (true / compensated)")
    ax.set_title("CFDX V80: real native BFS20540 first corrector — NO outer commit")
    ax.legend()
    ax.grid(axis="y", which="both", alpha=0.22)
    fig.tight_layout()
    fig.savefig(out_dir / "01_v80_native_csr_vs_face.png", dpi=180)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    data = parse(args.log.read_text(errors="replace"))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "evidence_v80.json").write_text(
        json.dumps(data, indent=2) + "\n"
    )
    plot(data, args.output_dir)
    print(json.dumps(data, indent=2))


if __name__ == "__main__":
    main()
