"""V73 audit: exact V72 history and fail-closed independent physical trials.

Run with --baseline/--log to audit a real BFS152 campaign, or without options
for falsification tests. No output can certify a physical BFS convergence.
"""

import argparse
import math
import re
from pathlib import Path

HISTORY = (
    "GLOBAL_SAME_OPERATOR_TRIAL_V65",
    "GLOBAL_SAME_OPERATOR_COMMITTED_V65",
    "INCOMPRESSIBLE_ITER",
)
TRIAL = "BFS_POST_SST_DIRECTION_V73"
SUMMARY = "BFS_POST_SST_SUMMARY_V73"


def records(text, tag):
    return [line.strip() for line in text.splitlines() if line.startswith(tag + " ")]


def fields(line):
    return dict(re.findall(r"([A-Za-z_][A-Za-z_0-9]*)=([^ ]+)", line))


def audit(baseline, log):
    for tag in HISTORY:
        if records(baseline, tag) != records(log, tag):
            raise ValueError(f"production history changed: {tag}")
    sums = records(log, SUMMARY)
    if len(sums) != 1:
        raise ValueError("expected exactly one V73 summary")
    sm = fields(sums[0])
    trials = records(log, TRIAL)
    if not trials or int(sm.get("evaluated", "-1")) > len(trials):
        raise ValueError("no V73 physical trials")
    if int(sm["evaluated"]) != sum(int(fields(t)["evaluated"]) for t in trials):
        raise ValueError("inconsistent V73 evaluated count")
    if int(sm["original_guards_ok"]) != sum(
        int(fields(t)["original_guards"])
        for t in trials
        if fields(t)["evaluated"] == "1"
    ):
        raise ValueError("inconsistent original guard total")
    for t in trials:
        data = fields(t)
        if data.get("iteration") != "23" or data.get("retry") != "8":
            raise ValueError("V73 not confined to critical BFS152 retry")
        for tag in ("private_trial", "committed", "qualified", "commit_authorized"):
            if data.get(tag) != ("1" if tag == "private_trial" else "0"):
                raise ValueError(f"V73 false candidate assertion: {tag}")
        if data["evaluated"] == "1" and not data.get("rAUx0_after"):
            raise ValueError("missing after-SST rAUx audit")
        if data["evaluated"] == "1" and not data.get("phi3_after"):
            raise ValueError("missing face3 flux audit")
    for tag in (
        "committed",
        "qualified",
        "commit_authorized",
        "native_momentum_linear_recertified",
        "native_sst_transport_recertified",
    ):
        if sm.get(tag) != "0":
            raise ValueError(f"V73 false native certification: {tag}")
    sources = records(log, "BFS_POST_SST_SOURCE_V73")
    if not sources:
        raise ValueError("missing omega source/decomposition evidence")
    # Every successfully evaluated candidate MUST have two independently
    # closed omega budgets on cells 0 and 12. A single trace is not evidence.
    by_key = {}
    for line in sources:
        row = fields(line)
        if row.get("committed") != "0" or row.get("qualified") != "0":
            raise ValueError("false commit in V73 SST budget")
        key = (row["gamma"], row["pressure_fraction"])
        if row.get("cell") not in ("0", "12"):
            raise ValueError("wrong BFS152 SST cell")
        if row["cell"] in by_key.setdefault(key, {}):
            raise ValueError("duplicate cell source budget")
        parts = [
            float(row[t])
            for t in (
                "convection_delta",
                "diffusion_delta",
                "production_delta",
                "destruction_delta",
                "cross_delta",
            )
        ]
        expected = float(row["omega_raw_delta"])
        if not all(math.isfinite(q) for q in parts + [expected]):
            raise ValueError("nonfinite source term")
        if abs(sum(parts) - expected) > 1e-8 * max(1.0, abs(expected), abs(sum(parts))):
            raise ValueError("nonconservative signed omega source decomposition")
        by_key[key][row["cell"]] = row
    valid_keys = set()
    for line in trials:
        row = fields(line)
        if row["evaluated"] == "1":
            valid_keys.add((row["gamma"], row["pressure_fraction"]))
            key = (row["gamma"], row["pressure_fraction"])
            if set(by_key.get(key, {})) != {"0", "12"}:
                raise ValueError("missing two-cell signed SST budgets")
            for name in (
                "mass_after",
                "merit_ratio",
                "worst_growth",
                "phi3_before",
                "phi3_after",
                "rAUx0_before",
                "rAUx0_after",
                "rAUx12_before",
                "rAUx12_after",
            ):
                if not math.isfinite(float(row[name])):
                    raise ValueError("nonfinite direction measurement")
    if set(by_key) != valid_keys:
        raise ValueError("source terms refer to an unqualified or missing trial")
    if int(sm["simultaneous_descent"]) != sum(
        int(fields(line)["all_blocks_improved"])
        for line in trials
        if fields(line)["evaluated"] == "1"
    ):
        raise ValueError("inconsistent simultaneous block descent count")
    return {
        "v73_trials": len(trials),
        "evaluated": int(sm["evaluated"]),
        "guards_ok": int(sm["original_guards_ok"]),
        "simultaneous_descent": int(sm["simultaneous_descent"]),
        "committed": 0,
        "qualified": False,
    }


def falsify():
    old = "\n".join(
        (
            "GLOBAL_SAME_OPERATOR_TRIAL_V65 id=1",
            "GLOBAL_SAME_OPERATOR_COMMITTED_V65 id=1",
            "INCOMPRESSIBLE_ITER id=1",
        )
    )
    trial = (
        "BFS_POST_SST_DIRECTION_V73 iteration=23 retry=8 evaluated=1 "
        "gamma=0.1 pressure_fraction=0.25 original_guards=1 "
        "all_blocks_improved=0 mass_after=0.1 merit_ratio=0.9 "
        "worst_growth=0.8 phi3_before=0.0 rAUx0_after=0.2 "
        "rAUx0_before=0.1 rAUx12_after=0.2 rAUx12_before=0.1 "
        "phi3_after=0.3 private_trial=1 committed=0 qualified=0 "
        "commit_authorized=0"
    )
    source = (
        "BFS_POST_SST_SOURCE_V73 iteration=23 retry=8 gamma=0.1 "
        "pressure_fraction=0.25 cell=0 convection_delta=0 diffusion_delta=0 "
        "production_delta=0 destruction_delta=0 cross_delta=0 sum=0 "
        "omega_raw_delta=0 private_trial=1 committed=0 qualified=0"
    )
    source12 = source.replace("cell=0", "cell=12")
    summary = (
        "BFS_POST_SST_SUMMARY_V73 evaluated=1 original_guards_ok=1 "
        "simultaneous_descent=0 native_momentum_linear_recertified=0 "
        "native_sst_transport_recertified=0 commit_authorized=0 "
        "committed=0 qualified=0"
    )
    log = old + "\n" + trial + "\n" + source + "\n" + source12 + "\n" + summary
    assert audit(old, log)["committed"] == 0
    attacks = [
        (old, log.replace("INCOMPRESSIBLE_ITER id=1", "INCOMPRESSIBLE_ITER id=2")),
        (old, log.replace("committed=0", "committed=1", 1)),
        (old, log.replace("commit_authorized=0", "commit_authorized=1", 1)),
        (
            old,
            log.replace(
                "native_momentum_linear_recertified=0",
                "native_momentum_linear_recertified=1",
            ),
        ),
        (
            old,
            log.replace(
                "native_sst_transport_recertified=0",
                "native_sst_transport_recertified=1",
            ),
        ),
        (
            old,
            log.replace(
                "evaluated=1 original_guards_ok=1", "evaluated=2 original_guards_ok=1"
            ),
        ),
        (old, log.replace("original_guards_ok=1", "original_guards_ok=2")),
        (old, log.replace("rAUx0_after=0.2", "broken_rau=0.2")),
        (old, log.replace("phi3_after=0.3", "broken_phi=0.3")),
        (old, log.replace("iteration=23 retry=8", "iteration=24 retry=8", 1)),
        (old, log.replace(source, "")),
        (old, log.replace(source12, "")),
        (old, log.replace("omega_raw_delta=0", "omega_raw_delta=1", 1)),
        (old, log.replace("all_blocks_improved=0", "all_blocks_improved=1")),
        (old, log + "\n" + summary),
    ]
    for i, (b, l) in enumerate(attacks):
        try:
            audit(b, l)
        except (ValueError, KeyError):
            continue
        raise AssertionError(f"attack {i} passed")
    print(
        "V73_EVIDENCE_FALSIFICATION_PASS " + str(len(attacks)) + "/" + str(len(attacks))
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--baseline", type=Path)
    p.add_argument("--log", type=Path)
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    if a.baseline and a.log:
        import json

        result = audit(a.baseline.read_text(), a.log.read_text())
        if a.output:
            a.output.write_text(json.dumps(result, indent=2) + "\n")
        print(result)
    else:
        falsify()


if __name__ == "__main__":
    main()
