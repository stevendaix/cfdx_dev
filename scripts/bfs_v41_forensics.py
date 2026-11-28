"""Fail-closed, read-only BFS V41 diagnostics from *actual* production logs.

This does not solve the equations or grant physical qualification; a logged
inner SST convergence and a globally committed iteration are distinct gates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import pathlib
import re
from typing import Any

LOG_TAGS = (
    "SST_TWO_CYCLE_RESCUE_V41",
    "SST_CERTIFIED_INNER",
    "SST_FP64_PRECISION_STALL_V14",
    "GLOBAL_NONLINEAR_COMMITTED_V13",
    "COUPLED_MONOLITHIC_ILUT_V16",
    "BFS_PHASE_V10",
    "POST_SST_MOMENTUM",
)
PAIR = re.compile(r"([A-Za-z_][A-Za-z_0-9]*)=([^ ]+)")


def decode_number(x: str) -> Any:
    if x in ("true", "false"):
        return x == "true"
    try:
        if re.fullmatch(r"[+-]?\d+", x):
            return int(x)
        value = float(x)
        return value if math.isfinite(value) else x
    except (ValueError, OverflowError):
        return x


def audit_post_sst_operator_cycle(
    events: list[dict[str, Any]], *, window: int = 64
) -> dict[str, Any]:
    """Read-only periodicity screen for the post-SST momentum nu_t operator.

    This observes the outer momentum/SST coupling seam, unlike V41's k/omega
    inner-transport state detector. A candidate is diagnostic only: it never
    changes damping or a convergence/qualification gate.
    """
    samples = [
        (int(item["iteration"]), float(item["relative_nut_change"]))
        for item in events
        if "iteration" in item
        and "relative_nut_change" in item
        and math.isfinite(float(item["relative_nut_change"]))
        and float(item["relative_nut_change"]) >= 0.0
    ]
    samples.sort(key=lambda item: item[0])
    if len(samples) < window:
        return {
            "status": "INSUFFICIENT_DATA",
            "sample_count": len(samples),
            "window": window,
            "candidate_period": None,
            "period_errors": {},
            "candidate_period_four": False,
        }

    values = [value for _, value in samples[-window:]]
    rms = math.sqrt(sum(value * value for value in values) / len(values))
    scale = max(rms, 1e-15)
    errors: dict[str, float] = {}
    for period in (2, 3, 4, 5, 6):
        differences = [
            values[index] - values[index - period]
            for index in range(period, len(values))
        ]
        error = math.sqrt(
            sum(value * value for value in differences) / len(differences)
        ) / scale
        errors[str(period)] = error

    best_period = min(errors, key=errors.get)
    candidate_four = errors["4"] < 1e-2 and errors["2"] > 1e-1
    return {
        "status": "ANALYZED",
        "sample_count": len(samples),
        "window": window,
        "first_iteration": samples[-window][0],
        "last_iteration": samples[-1][0],
        "candidate_period": int(best_period),
        "period_errors": errors,
        "candidate_period_four": candidate_four,
        "interpretation": (
            "diagnostic_candidate_only_no_solver_action"
            if candidate_four
            else "no_period_four_candidate_detected",
        )[0],
    }


def parse_production_log(source: pathlib.Path) -> dict[str, Any]:
    if not source.is_file() or not source.stat().st_size:
        raise ValueError(f"missing or empty actual solver log: {source}")
    data = source.read_bytes()
    try:
        lines = data.decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise ValueError(f"non-UTF-8 solver trace: {source}") from exc
    events: dict[str, list[dict[str, Any]]] = {tag: [] for tag in LOG_TAGS}
    for lineno, line in enumerate(lines, 1):
        for tag in LOG_TAGS:
            if line.startswith(tag + " "):
                values = {k: decode_number(v) for k, v in PAIR.findall(line)}
                values["source_line"] = lineno
                events[tag].append(values)
                break
    first_certified = next(iter(events["SST_CERTIFIED_INNER"]), None)
    first_rescue = next(iter(events["SST_TWO_CYCLE_RESCUE_V41"]), None)
    committed = events["GLOBAL_NONLINEAR_COMMITTED_V13"]
    stalled = events["SST_FP64_PRECISION_STALL_V14"]
    for item in committed:
        # A committed nonlinear ITERATE is not a global convergence certificate.
        if "qualified" not in item or "merit" not in item:
            raise ValueError(f"incomplete global quality gate in {source}")
    return {
        "source": str(source),
        "sha256": hashlib.sha256(data).hexdigest(),
        "bytes": len(data),
        "first_sst_inner_certificate": first_certified,
        "first_two_cycle_rescue": first_rescue,
        "committed_outer_states": len(committed),
        "last_committed": committed[-1] if committed else None,
        "first_precision_stall": stalled[0] if stalled else None,
        "precision_stall_count": len(stalled),
        "global_physically_qualified": bool(
            committed
            and committed[-1].get("qualified") == 1
            and committed[-1].get("converged") == 1
        ),
        "ilut_recovery": events["COUPLED_MONOLITHIC_ILUT_V16"][:1],
        "post_sst_momentum_operator_cycle": audit_post_sst_operator_cycle(
            events["POST_SST_MOMENTUM"]
        ),
        "phase_end": [
            e for e in events["BFS_PHASE_V10"] if e.get("phase") == "SST_TRANSPORT_END"
        ],
    }


def audit(actual: dict[str, pathlib.Path]) -> dict[str, Any]:
    if not actual or len(set(actual.values())) != len(actual):
        raise ValueError("at least one distinct physical log is required")
    runs = {name: parse_production_log(path) for name, path in actual.items()}
    return {
        "schema": "CFDX_BFS_V41_FORENSICS_1",
        "qualification": "diagnostics_only",
        "never_infer_physical_convergence_from_inner_sst": True,
        "runs": runs,
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--log", nargs=2, action="append", metavar=("LABEL", "PATH"), required=True
    )
    p.add_argument("--output", type=pathlib.Path, required=True)
    a = p.parse_args()
    if len({x[0] for x in a.log}) != len(a.log):
        p.error("duplicate run label")
    result = audit({name: pathlib.Path(path) for name, path in a.log})
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"V41_FORENSICS: {len(result['runs'])} logs, JSON={a.output}; does not certify BFS"
    )


if __name__ == "__main__":
    main()
