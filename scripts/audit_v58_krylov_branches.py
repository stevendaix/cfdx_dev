"""Audit fail-closed d'une campagne V58 de 6 directions Krylov sur BFS.

Le diagnostic des pentes de flux utilise phi +/- h dphi a rAU fige. Ce n'est
pas une reevaluation du flux Rhie-Chow couple a l'etat perturbe et il ne
constitue pas une qualification de Newton.
"""
import argparse
import json
import math
import re
from pathlib import Path

MARKERS = ("FPHYS_KRYLOV_BRANCH_V58 ", "FPHYS_KRYLOV_FACE_V58 ",
           "FPHYS_LIMITER_STATE_V58 ", "FPHYS_KRYLOV_FACE_GAP_V58 ",
           "FPHYS_GAP_STATE_V58 ", "FPHYS_KRYLOV_TRUTH_V57 ",
           "FPHYS_JFNK_V52 ")


def fields(line):
    return dict(re.findall(r"([a-z][a-z_0-9]*)=([^\s]+)", line))


def number(row, name):
    if name not in row:
        raise ValueError("missing " + name)
    v = float(row[name])
    if not math.isfinite(v):
        raise ValueError("nonfinite " + name)
    return v


def parse_lines(lines):
    branch, faces, limiter, gap, gap_state, krylov, jfnk = {}, {}, {}, {}, {}, {}, []
    for line in lines:
        prefix = next((p for p in MARKERS if line.startswith(p)), None)
        if not prefix:
            continue
        r = fields(line)
        if prefix.startswith("FPHYS_JFNK_V52 "):
            if any(r.get(k) != "0" for k in ("production_committed", "qualified")):
                raise ValueError("unsafe Newton qualification or commit")
            for k in ("linear_true_relative", "arnoldi_relative", "jv_linearity", "v57_hybrid_convection_gap"):
                number(r, k)
            if r.get("selected_schur") != "bfbt":
                raise ValueError("not the captured BFBt Krylov experiment")
            jfnk.append(r)
            continue
        if prefix.startswith("FPHYS_KRYLOV_TRUTH_V57 "):
            if r.get("production_committed") != "0":
                raise ValueError("unexpected Krylov state commit")
            idx = int(r.get("direction", "-1"))
            if idx in krylov or idx not in range(1, 7):
                raise ValueError("duplicate/missing Krylov direction")
            for k in ("arnoldi", "true", "nonadditivity", "gap", "consistency_limit"):
                number(r, k)
            krylov[idx] = r
            continue
        if any(r.get(k) != "0" for k in ("physical_qualified", "committed")):
            raise ValueError("unsafe face qualification or commit")
        idx = int(r.get("sample", "-1"))
        comp = int(r.get("component", "-1"))
        if idx not in range(18) or comp != idx % 3:
            raise ValueError("invalid or shuffled Krylov sample")
        key = (idx, comp)
        if prefix.startswith("FPHYS_KRYLOV_BRANCH_V58 "):
            if key in branch or r.get("native_rc_reclosure_at_offsets") != "0":
                raise ValueError("duplicate branch or false native Rhie-Chow claim")
            if int(r.get("faces", "0")) != 634:
                raise ValueError("not the captured BFS152 topology")
            for k in ("ambiguous", "effective_kinks", "unresolved", "inert", "roundoff_limited"):
                n = int(r.get(k, "-1"))
                if n < 0 or n > 634:
                    raise ValueError("invalid branch count " + k)
            for k in ("sided", "halving"):
                if number(r, k) < 0:
                    raise ValueError("negative branch metric")
            branch[key] = r
        else:
            face = int(r.get("face", "-1"))
            if not (0 <= face < 634):
                raise ValueError("invalid face ID")
            if prefix.startswith("FPHYS_KRYLOV_FACE_V58 "):
                for k in ("phi", "dphi", "right", "left", "sided", "halving"):
                    number(r, k)
                faces.setdefault(key, []).append(r)
            elif prefix.startswith("FPHYS_LIMITER_STATE_V58 "):
                for k in ("up", "down", "increment", "full", "lower", "upper", "reconstructed"):
                    number(r, k)
                limiter.setdefault(key, []).append(r)
            elif prefix.startswith("FPHYS_KRYLOV_FACE_GAP_V58 "):
                for k in ("phi", "right", "left", "tangent", "central", "absolute_gap"):
                    number(r, k)
                gap.setdefault(key, []).append(r)
            else:
                for k in ("up", "down", "increment", "full", "lower", "upper"):
                    number(r, k)
                gap_state.setdefault(key, []).append(r)

    if len(branch) != 18 or len(krylov) != 6 or len(jfnk) != 1:
        raise ValueError("truncated or duplicate Krylov proof")
    for key, record in branch.items():
        rows = faces.get(key, [])
        lim = limiter.get(key, [])
        if len(rows) != len(lim) or len(rows) > 2:
            raise ValueError("missing or extra face sample and limiter")
        if [int(x["face"]) for x in rows] != [int(x["face"]) for x in lim]:
            raise ValueError("misaligned limiter evidence")
        gs, st = gap.get(key, []), gap_state.get(key, [])
        if len(gs) != len(st) or len(gs) > 2:
            raise ValueError("missing gap state")
        if [int(x["face"]) for x in gs] != [int(x["face"]) for x in st]:
            raise ValueError("misaligned gap evidence")
    root = faces[(0, 0)]
    if not {5, 65} <= {int(x["face"]) for x in root}:
        raise ValueError("root cause faces 5 and 65 not evidenced")
    if not all(number(x, "phi") > 0 for x in root):
        raise ValueError("face mass-flux signs not preserved")
    if int(branch[(0, 0)]["effective_kinks"]) != 2:
        raise ValueError("wrong effective kink count")
    s = gap_state.get((0, 0), [])
    if not {10, 593} <= {int(x["face"]) for x in s}:
        raise ValueError("missing multiple limiter equality cases")
    for x in s:
        if abs(number(x, "up") - number(x, "down")) > 1e-8:
            raise ValueError("not a uniform-velocity limiter equality")
        if abs(number(x, "increment")) > 1e-10:
            raise ValueError("not a near-zero SOU increment")
    last = krylov[6]
    max_a = max(number(x, "nonadditivity") for x in krylov.values())
    linear_limit = number(last, "consistency_limit")
    true = number(last, "true")
    if not (max_a > linear_limit and true > .3 and jfnk[0].get("private_accepted") == "0"):
        raise ValueError("this measured failed Krylov result has been falsified")
    if abs(true - number(jfnk[0], "linear_true_relative")) > 1e-9:
        raise ValueError("inconsistent final linear residual")
    return {"status": "BLOCKED_KRYLOV_NONADDITIVE", "global_newton_qualified": False,
            "face_ids_root": [5, 65], "near_triple_tie_face_ids": [10, 593],
            "base_mass_flux_nonzero": True, "krylov_directions": 6,
            "max_nonadditivity": max_a, "linear_truth_final": true,
            "linear_target_experimental": .3,
            "native_rc_reclosure_at_face_offsets": False,
            "production_committed": False}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", type=Path, required=True)
    ap.add_argument("--json", type=Path)
    args = ap.parse_args()
    result = parse_lines(args.log.read_text(errors="replace").splitlines())
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.write_text(serialized)
    print(serialized, end="")


if __name__ == "__main__":
    main()
