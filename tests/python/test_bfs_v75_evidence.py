"""V75 independent log countercheck; synthetic self-tests with --log omitted.

V75 observes the existing V74 private direction, not a native linear solve.
"""

import argparse
import json
import math
import pathlib
import sys
from collections import defaultdict

TAG = "BFS_SST_BUDGET_V75 "
REQUIRED = {
    "iteration",
    "cycle",
    "cell",
    "equation",
    "raw_before",
    "raw_after",
    "raw_delta",
    "terms_delta_sum",
    "closure_error",
    "convection_delta",
    "diffusion_delta",
    "production_delta",
    "destruction_delta",
    "cross_delta",
    "largest_abs_delta",
    "k_before",
    "k_after",
    "omega_before",
    "omega_after",
    "nut_before",
    "nut_after",
    "rAUx_before",
    "rAUx_after",
    "phi3_before",
    "phi3_after",
    "trace_selected",
    "independently_closed",
    "private_trial",
    "committed",
    "qualified",
    "limiter_branches_captured",
    "normalized_before",
    "normalized_after",
}
PIECES = ("convection", "diffusion", "production", "destruction", "cross")


def check_close(a, b, scale=0.0, atol=1e-12, rtol=1e-8):
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b), scale)


def records(log):
    out = []
    for line in log.splitlines():
        if not line.startswith(TAG):
            continue
        d = {}
        for pair in line[len(TAG) :].split():
            if "=" not in pair:
                raise ValueError("unparsed evidence item")
            k, v = pair.split("=", 1)
            if k in d:
                raise ValueError("duplicate evidence key")
            d[k] = v
        if not REQUIRED <= d.keys():
            raise ValueError("missing V75 evidence fields")
        for k in REQUIRED - {"equation", "largest_abs_delta"}:
            x = float(d[k])
            if not math.isfinite(x):
                raise ValueError("nonfinite V75 evidence")
        if d["equation"] not in ("k", "omega"):
            raise ValueError("invalid SST equation name")
        if d["largest_abs_delta"] not in (*PIECES, "none"):
            raise ValueError("invalid dominant term")
        for f in ("private_trial", "independently_closed"):
            if d[f] != "1":
                raise ValueError("V75 is not a closed private diagnostic")
        for f in ("committed", "qualified", "limiter_branches_captured"):
            if d[f] != "0":
                raise ValueError("V75 falsely claimed production certification")
        out.append(d)
    return out


def verify(log, historical_pilot=True):
    items = records(log)
    if not items:
        raise ValueError("no V75 evidence")
    groups = defaultdict(dict)
    for d in items:
        key = (int(d["iteration"]), int(d["cycle"]), int(d["cell"]))
        if d["equation"] in groups[key]:
            raise ValueError("duplicate V75 equation row")
        groups[key][d["equation"]] = d

        def get(k, row=d):
            return float(row[k])

        pieces = []
        for k in PIECES:
            previous = "cross_before" if k == "cross" else k + "_before"
            current = "cross_after" if k == "cross" else k + "_after"
            # When auditing a real log, before/after terms must both exist.
            if previous not in d or current not in d:
                raise ValueError("missing before/after SST source terms")
            delta = get(k + "_delta")
            if not check_close(get(current) - get(previous), delta):
                raise ValueError("term difference does not match raw before/after")
            pieces.append(delta)
        raw_delta = get("raw_after") - get("raw_before")
        if not check_close(raw_delta, get("raw_delta")):
            raise ValueError("physical raw delta mismatch")
        if not check_close(sum(pieces), get("terms_delta_sum"), sum(map(abs, pieces))):
            raise ValueError("SST budget sum corrupted")
        if not check_close(
            get("terms_delta_sum") - get("raw_delta"), get("closure_error")
        ):
            raise ValueError("SST budget signed closure error corrupted")
        if not check_close(
            get("terms_delta_sum"), get("raw_delta"), sum(map(abs, pieces))
        ):
            raise ValueError("SST physical budget did not close")
        dominant = max(range(5), key=lambda i: abs(pieces[i]))
        expected = PIECES[dominant] if abs(pieces[dominant]) > 0 else "none"
        if d["largest_abs_delta"] != expected:
            raise ValueError("dominant SST contribution field corrupted")
        for field in ("k", "omega", "nut", "rAUx", "phi3"):
            if field in ("k", "omega", "nut") and get(field + "_before") < 0:
                raise ValueError("negative physical quantity")
            if field == "omega" and get(field + "_before") <= 0:
                raise ValueError("nonpositive omega")
            if field == "rAUx" and get(field + "_before") <= 0:
                raise ValueError("nonpositive momentum inverse diagonal")
    for key, two in groups.items():
        if set(two) != {"k", "omega"}:
            raise ValueError("missing k or omega budget")
        d = two["k"]
        e = two["omega"]
        for k in (
            "trace_selected",
            "phi3_before",
            "phi3_after",
            "rAUx_before",
            "rAUx_after",
            "nut_before",
            "nut_after",
            "k_before",
            "k_after",
            "omega_before",
            "omega_after",
        ):
            if d[k] != e[k]:
                raise ValueError("k and omega do not share the same state")
    for cycle in set(key[1] for key in groups):
        cells = {key[2] for key in groups if key[1] == cycle}
        if cells != {0, 12}:
            raise ValueError("V75 must cover both sensitive cells")
    if historical_pilot:
        if set(key[0] for key in groups) != {23}:
            raise ValueError("V74 pilot iteration 23 missing")
        if set(key[1] for key in groups) != {1, 2}:
            raise ValueError("V74 pilot needs both cycles")
        for key, eqs in groups.items():
            chosen = eqs["k"]["trace_selected"]
            if chosen != ("1" if key[1] == 1 else "0"):
                raise ValueError("V74 accepted/rejected path was changed")
    # Complete-domain maxima distinguish a truly local cause (cells 0/12)
    # from a larger residual elsewhere in the BFS152 mesh.
    max_rows = {}
    for line in log.splitlines():
        if not line.startswith("BFS_SST_ARGMAX_V75 "):
            continue
        parts = dict(piece.split("=", 1) for piece in line.split()[1:] if "=" in piece)
        required = {
            "iteration",
            "cycle",
            "equation",
            "before_cell",
            "after_cell",
            "before_abs",
            "after_abs",
            "source_cells_probed",
            "committed",
            "qualified",
        }
        if not required <= parts.keys():
            raise ValueError("incomplete V75 global SST maximum")
        key = (int(parts["iteration"]), int(parts["cycle"]), parts["equation"])
        if key in max_rows:
            raise ValueError("duplicate global SST maximum")
        if (
            not 0 <= int(parts["before_cell"]) < 152
            or not 0 <= int(parts["after_cell"]) < 152
        ):
            raise ValueError("invalid SST argmax cell")
        if not all(
            math.isfinite(float(parts[k])) and float(parts[k]) >= 0
            for k in ("before_abs", "after_abs")
        ):
            raise ValueError("invalid SST argmax magnitude")
        if (
            parts["source_cells_probed"] != "0,12"
            or parts["committed"] != "0"
            or parts["qualified"] != "0"
        ):
            raise ValueError("V75 argmax certification falsified")
        max_rows[key] = parts
    if set(max_rows) != {
        (key[0], key[1], eq) for key in groups for eq in ("k", "omega")
    }:
        raise ValueError("missing global SST argmax for an observed cycle")
    # Independent cross-check against the already-existing V74 acceptance
    # ledger. A self-consistent forged SST ledger is not sufficient.
    v74 = {}
    for line in log.splitlines():
        if not line.startswith("BFS_FIXED_DIRECTION_V74 "):
            continue
        parts = dict(piece.split("=", 1) for piece in line.split()[1:] if "=" in piece)
        if (
            not {
                "iteration",
                "cycle",
                "selected",
                "phi3_before",
                "phi3_after",
                "rAU0_before",
                "rAU0_after",
                "rAU12_before",
                "rAU12_after",
                "nut0_before",
                "nut0_after",
                "nut12_before",
                "nut12_after",
                "k_before",
                "k_after",
                "omega_before",
                "omega_after",
            }
            <= parts.keys()
        ):
            raise ValueError("incomplete V74 reference ledger")
        key = (int(parts["iteration"]), int(parts["cycle"]))
        if key in v74:
            raise ValueError("duplicate V74 reference cycle")
        v74[key] = parts
    if not v74:
        raise ValueError("missing independent V74 stage ledger")
    for key, eqs in groups.items():
        stage = v74.get(key[:2])
        if stage is None:
            raise ValueError("V75 row has no corresponding V74 cycle")
        d = eqs["k"]
        e = eqs["omega"]
        if d["trace_selected"] != stage["selected"]:
            raise ValueError("V75 selection differs from V74")
        cell = key[2]
        for prop, stageprop in (
            ("phi3_before", "phi3_before"),
            ("phi3_after", "phi3_after"),
            ("rAUx_before", f"rAU{cell}_before"),
            ("rAUx_after", f"rAU{cell}_after"),
            ("nut_before", f"nut{cell}_before"),
            ("nut_after", f"nut{cell}_after"),
            ("normalized_before", "k_before"),
            ("normalized_after", "k_after"),
        ):
            if not check_close(float(d[prop]), float(stage[stageprop])):
                raise ValueError("V75 snapshot does not match independent V74 state")
        for prop, stageprop in (
            ("normalized_before", "omega_before"),
            ("normalized_after", "omega_after"),
        ):
            if not check_close(float(e[prop]), float(stage[stageprop])):
                raise ValueError("V75 omega norm differs from V74")
    if set(v74) != {(k[0], k[1]) for k in groups}:
        raise ValueError("V74/V75 observation set differs")
    summaries = [
        line
        for line in log.splitlines()
        if line.startswith("BFS_SST_BUDGET_SUMMARY_V75 ")
    ]
    if len(summaries) != 1:
        raise ValueError("exactly one V75 diagnostic summary required")
    summary = dict(
        pair.split("=", 1) for pair in summaries[0].split()[1:] if "=" in pair
    )
    if (
        summary.get("status") != "forensic_only"
        or summary.get("committed") != "0"
        or summary.get("qualified") != "0"
        or summary.get("native_momentum") != "0"
        or summary.get("native_sst") != "0"
    ):
        raise ValueError("V75 summary falsely claimed native certification")
    if int(summary.get("cycles", "-1")) != len({k[1] for k in groups}):
        raise ValueError("V75 summary cycle count differs from actual evidence")
    if int(summary.get("private_accepted", "-1")) != sum(
        eqs["k"]["trace_selected"] == "1" for k, eqs in groups.items() if k[2] == 0
    ):
        raise ValueError("V75 summary acceptance count differs")
    return {
        "ok": True,
        "equation_records": len(items),
        "observations": len(groups),
        "cycles": sorted({k[1] for k in groups}),
        "native_certified": False,
    }


def fixture():
    def entry(cycle, cell, equation):
        terms = {
            "convection": 2.0,
            "diffusion": -1.0,
            "production": -3.0,
            "destruction": 1.0,
            "cross": 0.5,
        }
        updates = {
            "convection": -0.5,
            "diffusion": 0.5,
            "production": 1.0,
            "destruction": -0.5,
            "cross": 0.2,
        }
        raw_before = sum(terms.values())
        raw_after = raw_before + sum(updates.values())
        d = {
            "iteration": 23,
            "cycle": cycle,
            "cell": cell,
            "equation": equation,
            "raw_before": raw_before,
            "raw_after": raw_after,
            "raw_delta": raw_after - raw_before,
            "terms_delta_sum": sum(updates.values()),
            "closure_error": 0.0,
            "largest_abs_delta": "production",
            "k_before": 0.2,
            "k_after": 0.21,
            "omega_before": 2.0,
            "omega_after": 2.01,
            "nut_before": 0.03,
            "nut_after": 0.028,
            "rAUx_before": 500.0,
            "rAUx_after": 510.0,
            "phi3_before": 3e-4,
            "phi3_after": 2.9e-4,
            "normalized_before": 1e-5,
            "normalized_after": 2e-5,
            "trace_selected": 1 if cycle == 1 else 0,
            "private_trial": 1,
            "independently_closed": 1,
            "committed": 0,
            "qualified": 0,
            "limiter_branches_captured": 0,
        }
        for key, val in terms.items():
            kk = "cross" if key == "cross" else key
            d[kk + "_before"] = val
            d[kk + "_after"] = val + updates[key]
            d[kk + "_delta"] = updates[key]
        return TAG + " ".join(f"{k}={v}" for k, v in d.items())

    lines = [
        entry(c, cell, e) for c in (1, 2) for cell in (0, 12) for e in ("k", "omega")
    ]
    for cycle in (1, 2):
        for equation in ("k", "omega"):
            parts = {
                "iteration": 23,
                "cycle": cycle,
                "equation": equation,
                "before_cell": 0,
                "after_cell": 12,
                "before_abs": 1e-5,
                "after_abs": 2e-5,
                "source_cells_probed": "0,12",
                "committed": 0,
                "qualified": 0,
            }
            lines.append(
                "BFS_SST_ARGMAX_V75 " + " ".join(f"{k}={v}" for k, v in parts.items())
            )
        v74 = {
            "iteration": 23,
            "cycle": cycle,
            "selected": 1 if cycle == 1 else 0,
            "phi3_before": 3e-4,
            "phi3_after": 2.9e-4,
            "rAU0_before": 500.0,
            "rAU0_after": 510.0,
            "rAU12_before": 500.0,
            "rAU12_after": 510.0,
            "nut0_before": 0.03,
            "nut0_after": 0.028,
            "nut12_before": 0.03,
            "nut12_after": 0.028,
            "k_before": 1e-5,
            "k_after": 2e-5,
            "omega_before": 1e-5,
            "omega_after": 2e-5,
        }
        lines.append(
            "BFS_FIXED_DIRECTION_V74 " + " ".join(f"{k}={v}" for k, v in v74.items())
        )
    lines.append(
        "BFS_SST_BUDGET_SUMMARY_V75 iteration=23 cycles=2 private_accepted=1 "
        "status=forensic_only native_momentum=0 native_sst=0 "
        "committed=0 qualified=0"
    )
    # Historical V65 accepted transactions legitimately contain committed=1.
    # They are NOT V75 commits and must never invalidate this witness.
    lines.append(
        "GLOBAL_SAME_OPERATOR_COMMITTED_V65 iteration=4 committed=1 qualified=1"
    )
    return "\n".join(lines)


def self_test():
    baseline = fixture()
    assert verify(baseline)["equation_records"] == 8
    mutations = [
        baseline.replace("raw_delta=0.7000000000000002", "raw_delta=0.8", 1),
        baseline.replace("production_delta=1.0", "production_delta=2.0", 1),
        baseline.replace("convection_after=1.5", "convection_after=10.0", 1),
        baseline.replace("committed=0", "committed=1", 1),
        baseline.replace("independently_closed=1", "independently_closed=0", 1),
        baseline.replace("cell=12", "cell=9"),
        baseline.replace("trace_selected=0", "trace_selected=1"),
        baseline.replace(
            "limiter_branches_captured=0", "limiter_branches_captured=1", 1
        ),
        baseline.replace("phi3_before=0.0003", "phi3_before=0.001", 1),
        baseline.replace("rAU12_before=500.0", "rAU12_before=550.0", 1),
        baseline.replace("after_cell=12", "after_cell=200", 1),
        baseline.replace("native_sst=0", "native_sst=1", 1),
        baseline.replace("private_accepted=1", "private_accepted=2", 1),
    ]
    # Note: the first mutation checks a numeric token exactly as formatted by
    # Python repr; use a guaranteed mutation fallback if absent.
    if mutations[0] == baseline:
        mutations[0] = baseline.replace("raw_delta=0.7", "raw_delta=0.8", 1)
    for i, bad in enumerate(mutations):
        if bad == baseline:
            raise AssertionError(f"mutation {i} was not applied")
        try:
            verify(bad)
        except (ValueError, AssertionError):
            continue
        raise AssertionError(f"evidence falsification {i} was accepted")
    print(f"V75_EVIDENCE_PASS cases=1 counterfeit_rejected={len(mutations)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--any-cycles", action="store_true")
    args = parser.parse_args()
    try:
        if args.log is None:
            self_test()
        else:
            r = verify(args.log.read_text(errors="replace"), not args.any_cycles)
            if args.output:
                args.output.write_text(json.dumps(r, indent=2) + "\n")
            print(json.dumps(r, sort_keys=True))
    except (ValueError, AssertionError) as e:
        print(f"V75_EVIDENCE_FAIL {e}", file=sys.stderr)
        raise SystemExit(1)
