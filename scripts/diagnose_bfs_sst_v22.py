"""Non-qualifying SST V22 JFNK vs V21 Picard log comparison (fail closed)."""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

KEY_VALUE = re.compile(r"([a-zA-Z_][a-zA-Z_0-9]*)=([^\s]+)")


def parse_run(text: str) -> dict:
    runs: list[dict] = []
    current: dict | None = None
    for line in text.splitlines():
        if "SST_TRANSPORT_DIAGNOSTIC iteration=1 " in line:
            current = {"iteration_diagnostics": [], "jfnk_steps": [], "gates": []}
            runs.append(current)
        elif current is None and "SST_JFNK_V22 iteration=1 " in line:
            current = {"iteration_diagnostics": [], "jfnk_steps": [], "gates": []}
            runs.append(current)
        if current is None:
            continue
        values = dict(KEY_VALUE.findall(line))
        if line.startswith("SST_JFNK_V22 "):
            current["jfnk_steps"].append(values)
        elif line.startswith("SST_TRANSPORT_DIAGNOSTIC "):
            current["iteration_diagnostics"].append(values)
        elif line.startswith("SST_ACCEPTANCE_GATES "):
            current["gates"].append(values)
    complete = []
    for run in runs:
        seen: set[tuple[str, str]] = set()
        diagnostic = []
        for row in run["iteration_diagnostics"]:
            key = (row.get("iteration", ""), row.get("k_returned_nonlinear_raw", ""))
            if key not in seen:
                seen.add(key)
                diagnostic.append(row)
        if not diagnostic:
            continue
        samples = []
        for row in diagnostic:
            record = {"iteration": int(row["iteration"])}
            for name in ("k_returned_nonlinear_raw", "omega_returned_nonlinear_raw",
                         "k_returned_row_scaled", "omega_returned_row_scaled",
                         "k_local_relative", "omega_local_relative"):
                number = float(row[name])
                if not math.isfinite(number) or number < 0:
                    raise ValueError(f"nonfinite or negative SST defect: {name}")
                record[name] = number
            samples.append(record)
        complete.append({
            "samples": samples,
            # "accepted" is a globalization step, not necessarily Newton!
            # JFNK-mode tries Picard as a certified fallback, so distinguish
            # which direction was selected when the new V22 trace provides it.
            "accepted_steps_in_jfnk_mode": sum(
                step.get("accepted") == "1" for step in run["jfnk_steps"]),
            "picard_steps_selected": sum(
                step.get("accepted") == "1" and step.get("selected") == "picard"
                for step in run["jfnk_steps"]),
            "newton_steps_selected": sum(
                step.get("accepted") == "1" and step.get("selected") == "jfnk"
                for step in run["jfnk_steps"]),
            "accepted_steps_unknown_family": sum(
                step.get("accepted") == "1" and
                step.get("selected") not in {"picard", "jfnk"}
                for step in run["jfnk_steps"]),
            "rejected_steps_in_jfnk_mode": sum(
                step.get("accepted") == "0" for step in run["jfnk_steps"]),
            "all_reported_gates_closed": bool(run["gates"]) and all(
                gate.get("outer_update_ok") == "1" and
                gate.get("returned_equations_ok") == "1" and
                gate.get("finite") == "1" for gate in run["gates"]),
        })
    if not complete:
        raise ValueError("no complete SST physical-residual samples; refusing report")
    outer_converged = "Converged YES" in text
    qualified = False  # No local SST log can certify independent NASA/BFS V&V.
    return {"schema": "CFDX-BFS-SST-V22", "runs": complete,
            "outer_converged_as_reported": outer_converged,
            "qualified_physics": qualified,
            "reason": "Diagnostic only: full coupled U-p-k-omega and NASA gates not certified"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        report = parse_run(args.log.read_text(encoding="utf-8", errors="replace"))
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))
    document = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.write_text(document, encoding="utf-8")
    else:
        print(document, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
