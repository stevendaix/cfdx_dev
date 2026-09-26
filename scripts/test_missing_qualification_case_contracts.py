#!/usr/bin/env python3
"""Contract test for the non-VMFL qualification cases.

This test intentionally checks the qualification contracts, not solver PASS
status. A case is only promoted after a real solver driver is implemented.
"""
from pathlib import Path

CASES = {
    "EXT-NACA0012": "docs/validation/EXT-NACA0012_QUALIFICATION_CASE.md",
    "INT-CHANNEL": "docs/validation/INT-CHANNEL_QUALIFICATION_CASE.md",
    "SEP-BFS": "docs/validation/SEP-BFS_QUALIFICATION_CASE.md",
    "TUR-FLATPLATE": "docs/validation/TUR-FLATPLATE_QUALIFICATION_CASE.md",
    "TUR-CHANNEL": "docs/validation/TUR-CHANNEL_QUALIFICATION_CASE.md",
    "TUR-NACA0012": "docs/validation/TUR-NACA0012_QUALIFICATION_CASE.md",
}

REQUIRED = ("Reference:", "Required evidence:", "**Current status:**")


def main() -> int:
    failures = []
    for ident, filename in CASES.items():
        path = Path(filename)
        if not path.is_file():
            failures.append(f"{ident}: missing case contract {filename}")
            continue
        text = path.read_text(encoding="utf-8")
        for marker in REQUIRED:
            if marker not in text:
                failures.append(f"{ident}: missing contract section {marker!r}")
        print(f"QUALIFICATION_CASE_CONTRACT: {ident}: PASS")
    if failures:
        print("QUALIFICATION_CASE_CONTRACT: FAIL")
        for failure in failures:
            print("  -", failure)
        return 1
    print("QUALIFICATION_CASE_CONTRACT: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
