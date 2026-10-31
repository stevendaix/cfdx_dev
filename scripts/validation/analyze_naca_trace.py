#!/usr/bin/env python3
"""Triage NACA momentum traces; growth ratios are diagnostic, never gates."""
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

PREFIX = "NACA_MOMENTUM_TRACE"
METRICS = ("max_p", "max_gradp", "max_U", "max_phi", "max_divphi", "max_mu_eff")
TOKEN = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=([^\\s]+)")


def parse_trace(text: str) -> list[dict[str, float | int]]:
    records: list[dict[str, float | int]] = []
    for line in text.splitlines():
        if PREFIX not in line:
            continue
        fields = dict(TOKEN.findall(line))
        record: dict[str, float | int] = {}
        for key in ("seq", "iter"):
            if key in fields:
                record[key] = int(fields[key])
        for key in METRICS + ("boundary_flux_sum", "boundary_flux_abs_sum", "nonfinite_total"):
            if key in fields:
                try:
                    record[key] = float(fields[key])
                except ValueError:
                    continue
        if "seq" in record and "iter" in record:
            records.append(record)
    return records


def first_growth(records: list[dict[str, float | int]], ratio: float):
    if not math.isfinite(ratio) or ratio <= 1.0:
        raise ValueError("ratio must be finite and greater than one")
    for previous, current in zip(records, records[1:]):
        for metric in METRICS:
            before = previous.get(metric)
            after = current.get(metric)
            if (isinstance(before, (int, float)) and
                    isinstance(after, (int, float)) and
                    math.isfinite(before) and math.isfinite(after) and
                    before > 0.0 and after >= ratio * before):
                return {
                    "metric": metric,
                    "previous_seq": previous["seq"],
                    "current_seq": current["seq"],
                    "previous_iter": previous["iter"],
                    "current_iter": current["iter"],
                    "previous_value": before,
                    "current_value": after,
                    "factor": after / before,
                }
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--ratio", type=float, default=10.0)
    args = parser.parse_args()
    records = parse_trace(args.log.read_text(encoding="utf-8", errors="replace"))
    if not records:
        print(f"No {PREFIX} records found in {args.log}")
        return 2
    print(f"records={len(records)} ratio={args.ratio:g} (diagnostic only)")
    event = first_growth(records, args.ratio)
    if event is None:
        print("No consecutive-record growth event found.")
    else:
        print("FIRST_GROWTH " + " ".join(f"{k}={v}" for k, v in event.items()))
    nonfinite = [r for r in records if r.get("nonfinite_total", 0) > 0]
    if nonfinite:
        print(f"NONFINITE_RECORDS={len(nonfinite)} first_seq={nonfinite[0]['seq']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
