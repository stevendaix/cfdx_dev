#!/usr/bin/env python3
"""Reject deferred CFDX harness registrations without a run_all() call.

This narrow static guard recognizes literal-name/lambda run_case registrations
from common/test_harness.h. Local numerical helpers named run_case, comments,
and string literals are not test registrations. It complements execution; it
does not prove that an arbitrary call is reachable or that an oracle is sound.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re


def cpp_tokens(source: str) -> list[tuple[str, str]]:
    tokens: list[tuple[str, str]] = []
    i = 0
    while i < len(source):
        if source[i].isspace():
            i += 1
        elif source.startswith("//", i):
            end = source.find("\n", i + 2)
            i = len(source) if end < 0 else end + 1
        elif source.startswith("/*", i):
            end = source.find("*/", i + 2)
            i = len(source) if end < 0 else end + 2
        elif source.startswith('R"', i):
            opening = source.find("(", i + 2, i + 19)
            if opening < 0:
                i += 1
                continue
            closing = ")" + source[i + 2 : opening] + '"'
            end = source.find(closing, opening + 1)
            end = len(source) if end < 0 else end + len(closing)
            tokens.append(("string", source[i:end]))
            i = end
        elif source[i] in "\"'":
            quote, start = source[i], i
            i += 1
            while i < len(source):
                if source[i] == "\\":
                    i += 2
                elif source[i] == quote:
                    i += 1
                    break
                else:
                    i += 1
            tokens.append(("string", source[start:i]))
        elif source[i].isalpha() or source[i] == "_":
            match = re.match(r"[A-Za-z_][A-Za-z_0-9]*", source[i:])
            assert match is not None
            value = match.group(0)
            tokens.append(("identifier", value))
            i += len(value)
        else:
            tokens.append(("punctuation", source[i]))
            i += 1
    return tokens


def inspect_source(source: str) -> tuple[int, bool]:
    tokens = cpp_tokens(source)
    values = [token[1] for token in tokens]
    uses_harness = any(
        values[i : i + 3] == ["#", "include", '"common/test_harness.h"']
        for i in range(len(tokens) - 2)
    )
    if not uses_harness:
        return 0, False
    registrations = sum(
        values[i] == "run_case"
        and values[i + 1] == "("
        and tokens[i + 2][0] == "string"
        and values[i + 3 : i + 5] == [",", "["]
        for i in range(len(tokens) - 4)
    )
    executes = any(
        values[i : i + 3] == ["run_all", "(", ")"]
        for i in range(len(tokens) - 2)
    )
    return registrations, registrations > 0 and not executes


def self_test() -> None:
    prefix = '#include "common/test_harness.h"\n'
    registration = 'int main() { run_case("example", [] { fail(); }); return 0; }'
    assert inspect_source(prefix + registration) == (1, True)
    assert inspect_source(prefix + registration.replace("return 0", "return run_all()")) == (1, False)
    assert inspect_source(prefix + '// run_case("ignored", [] {});\n') == (0, False)
    assert inspect_source(prefix + 'const char *x = R"tag(run_case("ignored", [] {});)tag";') == (0, False)
    assert inspect_source(prefix + 'double run_case(int n) { return n; } int main() { run_case(16); }') == (0, False)
    assert inspect_source(prefix + registration + '// run_all()\nconst char *x="run_all()";') == (1, True)
    print("Harness guard self-tests: 6 passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    files, cases, missing = 0, 0, []
    for path in sorted((args.root / "tests").rglob("*.cpp")):
        count, absent = inspect_source(path.read_text(encoding="utf-8"))
        if count:
            files += 1
            cases += count
        if absent:
            missing.append((path.relative_to(args.root), count))
    for path, count in missing:
        print(f"FAIL: {path}: {count} deferred case(s) but no run_all() call")
    print(f"Harness execution guard: {files} files, {cases} registered cases, {len(missing)} missing executions")
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
