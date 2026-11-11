"""Fail-closed V40 audit for six identical original coupled operators.

Checks FNV64 fingerprints *alongside*, not instead of, identical physical HDF5
inputs. A missing capture or any mismatch is NOT an algebraically comparable
six-Schur campaign. FNV64 is not a cryptographic SHA256 guarantee.
"""
from __future__ import annotations
import re

PATTERN = re.compile(r"^COUPLED_MATRIX_FINGERPRINT_V40 A=([0-9a-f]{16}) "
                     r"b=([0-9a-f]{16}) x0=([0-9a-f]{16}) "
                     r"rows=(\d+) nnz=(\d+) fnv64_not_sha256=1 "
                     r"physical_converged=0$", re.M)
FIELDS = ("A", "b", "x0")


def parse_fingerprint(log: str) -> dict | None:
    entries = PATTERN.findall(log)
    if not entries:
        return None
    if len(entries) != 1:
        raise ValueError("V40 expected one, not multiple, frozen operator fingerprints")
    a, b, x0, rows, nnz = entries[0]
    rows, nnz = int(rows), int(nnz)
    if not rows or not nnz:
        raise ValueError("V40 empty captured coupled operator")
    return {"A": a, "b": b, "x0": x0, "rows": rows, "nnz": nnz,
            "hash_algorithm": "FNV64_NOT_CRYPTOGRAPHIC"}


def compare_six_fingerprints(results: list[dict], schur_names: tuple) -> dict:
    observed = {r.get("schur") for r in results}
    exact_six = len(results) == len(schur_names) and observed == set(schur_names)
    captures = [r.get("operator_fingerprint_v40") for r in results]
    complete = exact_six and all(isinstance(x, dict) and
                                 all(k in x for k in (*FIELDS,"rows","nnz"))
                                 for x in captures)
    agreement = complete and all(len({c[k] for c in captures}) == 1
                                 for k in (*FIELDS,"rows","nnz"))
    return {"all_six_present": exact_six,
            "all_six_original_operators_fingerprinted": bool(complete),
            "all_six_same_A_b_x0_fingerprint": bool(agreement),
            "identity_is_cryptographic_proof": False,
            "nonlinear_convergence_implied": False}
