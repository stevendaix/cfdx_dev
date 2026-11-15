@@ -0,0 +1,89 @@
#!/usr/bin/env python3
"""Fail-closed V69 same-candidate SST/RC Picard experiment log audit.

This is evidence for an *uncommitted* two-cycle SST transport, NOT proof of
physical convergence. Compare 149 V65 trial decisions against V68 literally.
"""
import argparse
import json
import math
from pathlib import Path

PREFIXES = ("GLOBAL_SAME_OPERATOR_TRIAL_V65", "BFS_SST_FLUX_SWAP_V68")

def extract(lines, prefix):
    return [line.strip() for line in lines if line.startswith(prefix + " ")]

def fields(line):
    values = {}
    for item in line.split()[1:]:
        if "=" in item:
            key, value = item.split("=", 1)
            values[key] = value
    return values

def audit(log, baseline, require_complete=False):
    lines = Path(log).read_text(errors="replace").splitlines()
    old = Path(baseline).read_text(errors="replace").splitlines()
    matching = {prefix: extract(lines, prefix) == extract(old, prefix)
                for prefix in PREFIXES}
    counts = {prefix: len(extract(lines, prefix)) for prefix in PREFIXES}
    attempts = extract(lines, "BFS_SST_RC_PICARD_V69")
    records = []
    negative = []
    for line in attempts:
        d = fields(line)
        if d.get("closed") == "0":
            negative.append(d)
            continue
        if d.get("committed") != "0" or d.get("qualified") != "0" or d.get("private_transport") != "1":
            raise ValueError("V69 probe was incorrectly committed or certified")
        numbers = {key: float(d[key]) for key in (
            "beta", "omega_phys_before", "omega_phys_after", "omega_ratio",
            "k_phys_after", "merit_before", "merit_after", "merit_ratio",
            "continuity_before", "continuity_after", "rc_phi_change_linf",
            "nut_change_linf", "v65_merit_ratio", "v65_worst_growth")}
        if not all(math.isfinite(x) and x >= 0 for k, x in numbers.items() if k != "beta"):
            raise ValueError("V69 nonfinite/negative evidence")
        if d.get("rc_closed") != "1":
            raise ValueError("V69 non-closed Rhie-Chow flux used for evidence")
        flags={key:int(d[key]) for key in
               ("v65_joint_ok","v65_monotone_ok","v65_would_pass")}
        if any(value not in (0,1) for value in flags.values()) or \
           flags["v65_would_pass"] != (flags["v65_joint_ok"] & flags["v65_monotone_ok"]):
            raise ValueError("V69 V65 counterfactual decision inconsistent")
        worst=int(d["v65_worst_block"])
        if not 0 <= worst < 6:
            raise ValueError("V69 invalid V65 residual block index")
        records.append({"iteration": int(d["iteration"]), "retry": int(d["retry"]),
                        "cycle": int(d["cycle"]), **numbers,
                        "v65_worst_block":worst, **flags})
    if not all(matching.values()) or counts[PREFIXES[0]] != 149 or counts[PREFIXES[1]] != 3:
        raise ValueError("V69 altered or lost the historical V65/V68 trace")
    if len(records) + len(negative) == 0:
        raise ValueError("V69 experiment did not execute")
    if not any("iteration=23" in line for line in attempts):
        raise ValueError("V69 did not reach the blocked iteration")
    expected={(r,c) for r in (0,8,10) for c in (1,2)}
    observed={(r["retry"],r["cycle"]) for r in records}
    complete=(len(records)==6 and observed==expected and not negative)
    if require_complete and not complete:
        raise ValueError("V69 incomplete private SST/RC cycle evidence")
    return {"private_cycle_evidence_complete": complete,"replay_exact": matching, "historical_counts": counts,
            "probe_records": records, "rejected_probes": negative,
            "physical_bfs_converged": False, "qualified": False,
            "interpretation": "read-only SST transport with candidate-owned RC flux; no physical commit"}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", required=True)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = audit(args.log, args.baseline, require_complete=True)
    Path(args.output).write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print("V69_AUDIT_RECORDED candidate_cycles=", len(result["probe_records"]),
          "negative=", len(result["rejected_probes"]),
          "v65_trials=", result["historical_counts"][PREFIXES[0]])
if __name__ == "__main__":
    main()
