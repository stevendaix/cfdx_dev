#!/usr/bin/env python3
"""Triage NACA_PRESSURE_CORRECTION chronologies; the linear-system versus
pre-solve coupling distinction is a diagnostic observation, never a gate."""

from __future__ import annotations

import argparse
import itertools
import math
import re
from pathlib import Path

PREFIX = "NACA_PRESSURE_CORRECTION"
GROWTH_METRICS = (
    "rhs_linf",
    "correction_linf",
    "pressure_before_linf",
    "rhs_diagonal_ratio_linf",
)
TOKEN = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=([^\s]+)")

# A converged solve with a tiny relative residual that still produced an
# oversized correction means the RHS/coupling was already extreme before
# Krylov ran: the linear system itself is not the origin of the growth.
CONVERGED_STATUS = 1
CLEAN_RELATIVE_THRESHOLD = 1.0e-8


def parse_pressure_correction(text: str) -> list[dict[str, float | int]]:
    records: list[dict[str, float | int]] = []
    for line in text.splitlines():
        if PREFIX not in line:
            continue
        fields = dict(TOKEN.findall(line))
        record: dict[str, float | int] = {}
        for key in ("iter", "corrector", "status", "iterations", "rows", "nnz",
                    "rhs_cell", "correction_cell", "pressure_cell",
                    "diagonal_cell", "minimum_diagonal_cell",
                    "rhs_diagonal_ratio_cell", "minimum_dominance_cell"):
            if key in fields:
                try:
                    record[key] = int(fields[key])
                except ValueError:
                    continue
        for key in ("residual", "relative", "rhs_linf", "correction_linf",
                    "pressure_before_linf", "diagonal_min", "diagonal_abs_max",
                    "minimum_diagonal_rhs", "minimum_diagonal_correction",
                    "rhs_diagonal_ratio_linf", "diagonal_dominance_min",
                    "row_sum_linf", "row_abs_sum_linf"):
            if key in fields:
                try:
                    record[key] = float(fields[key])
                except ValueError:
                    continue
        if "iter" in record and "rows" in record and "correction_linf" in record:
            records.append(record)
    return records


def first_growth(
    records: list[dict[str, float | int]],
    metric: str,
    ratio: float,
) -> dict[str, float | int | str] | None:
    if not math.isfinite(ratio) or ratio <= 1.0:
        raise ValueError("ratio must be finite and greater than one")
    for previous, current in itertools.pairwise(records):
        before = previous.get(metric)
        after = current.get(metric)
        if (
            isinstance(before, (int, float))
            and isinstance(after, (int, float))
            and math.isfinite(before)
            and math.isfinite(after)
            and before > 0.0
            and after >= ratio * before
        ):
            return {
                "metric": metric,
                "previous_iter": previous["iter"],
                "current_iter": current["iter"],
                "previous_value": before,
                "current_value": after,
                "factor": after / before,
            }
    return None


def is_clean_solve(record: dict[str, float | int]) -> bool:
    status = record.get("status")
    relative = record.get("relative")
    return (
        status == CONVERGED_STATUS
        and isinstance(relative, (int, float))
        and math.isfinite(relative)
        and relative <= CLEAN_RELATIVE_THRESHOLD
    )


def largest_clean_solve(
    records: list[dict[str, float | int]],
) -> dict[str, float | int] | None:
    clean = [
        record
        for record in records
        if is_clean_solve(record)
        and math.isfinite(record.get("correction_linf", float("inf")))
    ]
    if not clean:
        return None
    return max(clean, key=lambda record: record["correction_linf"])


def group_by_rows(
    records: list[dict[str, float | int]],
) -> dict[int, list[dict[str, float | int]]]:
    groups: dict[int, list[dict[str, float | int]]] = {}
    for record in records:
        rows = int(record["rows"])
        groups.setdefault(rows, []).append(record)
    return groups


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--ratio", type=float, default=10.0)
    args = parser.parse_args()
    records = parse_pressure_correction(
        args.log.read_text(encoding="utf-8", errors="replace"))
    if not records:
        print(f"No {PREFIX} records found in {args.log}")
        return 2
    print(
        f"records={len(records)} rows={sorted(group_by_rows(records))} "
        f"ratio={args.ratio:g} (diagnostic only)"
    )
    for metric in GROWTH_METRICS:
        event = first_growth(records, metric, args.ratio)
        if event is None:
            print(f"NO_GROWTH metric={metric}")
        else:
            print(
                "FIRST_GROWTH "
                + " ".join(f"{key}={value}" for key, value in event.items())
            )
    clean = largest_clean_solve(records)
    if clean is None:
        print("NO_CLEAN_SOLVE")
    else:
        print(
            "LARGEST_CLEAN_SOLVE"
            f" iter={clean.get('iter')}"
            f" rows={clean.get('rows')}"
            f" iterations={clean.get('iterations')}"
            f" relative={clean.get('relative')}"
            f" rhs_linf={clean.get('rhs_linf')}"
            f" correction_linf={clean.get('correction_linf')}"
            f" rhs_diagonal_ratio_linf={clean.get('rhs_diagonal_ratio_linf')}"
        )
    for rows in sorted(group_by_rows(records)):
        group = group_by_rows(records)[rows]
        worst = max(
            group, key=lambda record: record.get("correction_linf", 0.0)
        )
        print(
            f"MESH rows={rows} records={len(group)}"
            f" max_correction_linf={worst.get('correction_linf')}"
            f" max_correction_iter={worst.get('iter')}"
            f" max_correction_status={worst.get('status')}"
            f" max_correction_relative={worst.get('relative')}"
            f" max_rhs_linf={max(r.get('rhs_linf', 0.0) for r in group)}"
            f" first_iter={group[0].get('iter')}"
            f" last_iter={group[-1].get('iter')}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
