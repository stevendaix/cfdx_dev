#!/usr/bin/env python3
"""Generic computed-vs-reference QoI comparison for CFDX validation campaigns."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read_csv(path: Path) -> dict[str, dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "case" not in reader.fieldnames:
            raise ValueError(f"{path}: CSV must contain a 'case' column")
        result: dict[str, dict[str, float]] = {}
        for row_number, row in enumerate(reader, start=2):
            case = row.get("case", "")
            if not case:
                raise ValueError(f"{path}:{row_number}: empty case identifier")
            if case in result:
                raise ValueError(f"{path}:{row_number}: duplicate case '{case}'")
            values: dict[str, float] = {}
            for name, value in row.items():
                if name == "case" or value in ("", None):
                    continue
                try:
                    number = float(value)
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"{path}:{row_number}: non-numeric QoI '{name}'={value!r}"
                    ) from exc
                if not math.isfinite(number):
                    raise ValueError(
                        f"{path}:{row_number}: QoI '{name}' is not finite"
                    )
                values[name] = number
            result[case] = values
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--computed", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--qoi", nargs="+", required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    computed = read_csv(args.computed)
    reference = read_csv(args.reference)
    computed_cases = set(computed)
    reference_cases = set(reference)
    missing_computed = sorted(reference_cases - computed_cases)
    missing_reference = sorted(computed_cases - reference_cases)
    if missing_computed or missing_reference:
        raise SystemExit(
            "case-set mismatch: "
            f"missing computed={missing_computed}; "
            f"missing reference={missing_reference}"
        )

    results = []
    for case in sorted(computed_cases):
        for qoi in args.qoi:
            if qoi not in computed[case] or qoi not in reference[case]:
                raise SystemExit(f"missing QoI '{qoi}' for case '{case}'")
            cv = computed[case][qoi]
            rv = reference[case][qoi]
            absolute_error = abs(cv - rv)
            results.append({
                "case": case,
                "qoi": qoi,
                "computed": cv,
                "reference": rv,
                "absolute_error": absolute_error,
                "relative_error": (
                    absolute_error / abs(rv) if rv != 0.0 else None
                ),
            })

    report = {
        "computed": str(args.computed),
        "reference": str(args.reference),
        "qoi": args.qoi,
        "results": results,
    }
    data = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(data, encoding="utf-8")
    else:
        print(data, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
