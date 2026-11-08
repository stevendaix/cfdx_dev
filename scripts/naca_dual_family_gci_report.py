"""Independent pressure/viscous GCI and dual-family cross-mesh comparison.

Never produce extrapolated CFD values from <3 physically converged, compatible,
second-order simulations per family. Preserve incomplete data as observations.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path

from naca_force_mesh_qualification import (certify, read_one_certificate,
                                           verify_physical_history)
from naca_compare_mesh_families import compare
from qualify_naca_grid import qualify

LEVELS = (64, 128, 256)
FIELDS = ("cd", "cdp", "cdv")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1048576), b""):
            h.update(b)
    return h.hexdigest()


def physical_setup(exe: Path, pressure_boundary_policy: str) -> dict:
    return {
        "airfoil": "NACA0012_finite_thickness_TE",
        "reynolds": 1000.0,
        "angle_deg": 0.0,
        "chord": 1.0,
        "farfield_radius_chords": 20.0,
        "velocity_boundary": "Dirichlet_Uinf_farfield_no_slip_airfoil",
        "pressure_boundary": pressure_boundary_policy,
        "convection_scheme": "TVD_VANLEER",
        "wall_shear_method": "FACE_WLS",
        "solver_revision": "bin-sha256:" + sha256(exe),
        "linear_tolerance": 1e-12,
        "nonlinear_tolerance": 1e-8,
    }


def assumed_order_two_grid(coarse: dict, fine: dict, *, order: float = 2.0,
                           safety: float = 3.0) -> dict:
    """Explicitly non-certified two-grid indicator, never an observed-order GCI.

    NASA recommends F_s=3 when using only two grid levels, and order is
    assumed (not measured). This function NEVER reports a qualified limit.
    """
    if not all(x.get("converged") for x in (coarse, fine)):
        raise ValueError("only two physically verified cases may enter estimator")
    if (coarse["N"], fine["N"]) != (64, 128):
        raise ValueError("only the tested N64/N128 pair may enter estimator")
    if (not math.isfinite(order) or order <= 0
            or not math.isfinite(safety) or safety < 1):
        raise ValueError("invalid assumed Richardson parameters")
    denominator = math.expm1(order * math.log(2.0))
    result = {"status": "ASSUMED_ORDER_ONLY_NOT_QUALIFIED",
              "assumed_order": order, "safety_factor": safety, "components": {}}
    for key in FIELDS:
        coarse_v, fine_v = (float(r[key]) for r in (coarse, fine))
        delta = fine_v - coarse_v
        uncertainty = safety * abs(delta) / denominator
        result["components"][key] = {
            "fine_value": fine_v,
            "hypothetical_order_based_limit": fine_v + delta / denominator,
            "gci_absolute_conditional": uncertainty,
            "gci_relative_conditional": uncertainty / abs(fine_v) if fine_v else None,
        }
    return result


def check_wall_shape(rows: list[dict]) -> dict:
    """Distinguish sustained 2-D wall aspect from singular TE/cosine panels."""
    p95 = [float(r["wall_edge_aspect_ratio_p95"]) for r in rows]
    p50 = [float(r["wall_edge_aspect_ratio_p50"]) for r in rows]
    maxr = [float(r["wall_edge_aspect_ratio_max"]) for r in rows]
    for values in (p95, p50, maxr):
        if len(values) != 3 or not all(math.isfinite(v) and v >= 1.0
                                      for v in values):
            raise ValueError(
                "missing or invalid independent near-wall aspect diagnostics")
    delta95 = max(p95) / min(p95)
    delta50 = max(p50) / min(p50)
    worst = max(maxr) / min(maxr)
    return {
        "status": ("REPEATABLE_BULK_WALL_ASPECT"
                   if max(delta95, delta50) <= 1.2
                   else "VARIABLE_BULK_WALL_ASPECT"),
        "wall_aspect_p95": p95,
        "wall_aspect_p50": p50,
        "wall_aspect_max": maxr,
        "worst_edge_aspect_growth": worst,
        "warning": (
            "worst wall-cell aspect changes strongly under refinement; "
            "inspect leading and trailing-edge panel topology" if worst > 1.4 else ""
        ),
    }


def inspect_family(family: str, family_dir: Path, runs: Path,
                   exe: Path) -> dict:
    quality = json.loads((family_dir / "quality.json").read_text())
    result = {"name": family, "status": "BLOCKED", "observations": [],
              "reasons": [], "mesh_quality": quality}
    rows = {int(x["n_surface"]): x for x in quality}
    if set(rows) != set(LEVELS):
        result["reasons"].append("mesh family does not have N64/N128/N256")
        return result
    try:
        result["wall_shape"] = check_wall_shape([rows[n] for n in LEVELS])
        if result["wall_shape"]["status"] != "REPEATABLE_BULK_WALL_ASPECT":
            result["reasons"].append(
                "bulk near-wall aspect is not comparable across grid levels")
    except (ValueError, TypeError, KeyError) as exc:
        result["reasons"].append("geometry: " + str(exc))
    certs = {}
    histories = {}
    verified_pressure_policies = set()
    policies_by_level = {}
    for n in LEVELS:
        location = runs / family / f"n{n}"
        meta = location / "run_summary.json"
        row = {"N": n, "converged": False}
        if not meta.is_file():
            result["reasons"].append(f"N{n}: missing solver run summary")
        else:
            try:
                m = json.loads(meta.read_text())
                if (m.get("exit_code") != 0
                        or m.get("scheme") != "TVD_VANLEER"
                        or m.get("linear_tolerance") != 1e-12
                        or m.get("wall_traction") != "FACE_WLS"
                        or m.get("nonlinear_tolerance") != 1e-8):
                    raise ValueError(
                        "solver failed or changed numerical contract")
                policy = m.get("pressure_boundary_policy")
                if policy not in ("DIRICHLET_ZERO", "ZERO_NORMAL_GRADIENT"):
                    raise ValueError(
                        "missing/unknown explicit farfield pressure policy")
                if m.get("kinematic_face_flux_certified") is not True:
                    raise ValueError(
                        "missing strict per-face kinematic flux proof")
                verified_pressure_policies.add(policy)
                policies_by_level[n] = policy
                expected = family_dir / f"naca0012_{family}_n{n}.h5"
                if Path(m["mesh"]).resolve() != expected.resolve():
                    raise ValueError("solver mesh does not match family")
                cp = location / "force.csv"
                hp = location / "history.csv"
                if not cp.is_file() or not hp.is_file():
                    raise ValueError(
                        "missing converged force or nonlinear history")
                with cp.open(newline="") as f:
                    rec = list(csv.DictReader(f))
                if len(rec) != 1:
                    raise ValueError("one force row required")
                strict = read_one_certificate(cp, n)
                physical = verify_physical_history(hp, n, strict)
                if physical["convection_scheme"] != m["scheme"]:
                    raise ValueError(
                        "historical convection scheme disagrees with run metadata")
                row.update({k: float(rec[0][k]) for k in FIELDS})
                row["converged"] = True
                certs[n] = cp
                histories[n] = hp
            except (ValueError, OSError, KeyError, TypeError) as exc:
                result["reasons"].append(f"N{n}: {exc}")
        result["observations"].append(row)
    if len(verified_pressure_policies) > 1:
        result["reasons"].append(
            "mixed farfield pressure policies across N64/N128/N256")
    result["pressure_boundary_policy"] = (
        next(iter(verified_pressure_policies))
        if len(verified_pressure_policies) == 1 else "UNVERIFIED_OR_MIXED")
    if (all(row["converged"] for row in result["observations"][:2])
            and policies_by_level.get(64) == policies_by_level.get(128)
            and policies_by_level.get(64) in ("DIRICHLET_ZERO",
                                              "ZERO_NORMAL_GRADIENT")):
        result["two_grid_assumed_order_only"] = assumed_order_two_grid(
            result["observations"][0], result["observations"][1])
    if len(certs) != 3 or result["reasons"]:
        # No extrapolated coefficients or GCI may appear on incomplete runs.
        return result
    certificate = certify(quality, certs, histories)
    result["solver_certificate"] = certificate
    physical_cases = {x.get("N"): x for x in certificate.get("cases", [])}
    for row in result["observations"]:
        row["converged"] = (
            physical_cases.get(row["N"], {}).get("status") == "CONVERGED"
        )
    if certificate["status"] != "GCI_CANDIDATE":
        result["reasons"].extend(
            certificate.get("reasons") or ["not all three quantities asymptotic"])
        # Preserve independently valid component estimates even when the total
        # three-component certificate fails. Never estimate if a *physical*
        # nonlinear gate, grid similarity, or mesh-provenance gate failed.
        if (all(x["converged"] for x in result["observations"])
                and certificate.get("mesh_family", {}).get("status")
                == "COMPATIBLE_WITH_SINGLE_RATIO_GCI"
                and result["wall_shape"]["status"]
                == "REPEATABLE_BULK_WALL_ASPECT"):
            rows3 = [read_one_certificate(certs[n], n) for n in LEVELS]
            result["gci"] = qualify(rows3)["quantities"]
            result["status"] = "PARTIAL_COMPONENT_GCI_ONLY"
        return result
    result["status"] = "GCI_CANDIDATE"
    result["gci"] = certificate["gci"]["result"]["quantities"]
    result["physical_setup"] = physical_setup(
        exe, result["pressure_boundary_policy"])
    return result


def analyze(families: Path, runs: Path, exe: Path,
            first: str = "quasi_isotropic", second: str = "anisotropic") -> dict:
    result = {"status": "BLOCKED", "families": {}, "comparison": {}, "reasons": [],
              "reference": {"cd": 0.119,
                            "type": "published_numerical_not_experiment"},
              "causal_anisotropy_isolated": False}
    for name in dict.fromkeys((first, second, "sheared_isosurface",
                               "controlled_anisotropic",
                               "controlled_anisotropic_mild")):
        if name in result["families"] or not (families / name
                                              / "quality.json").is_file():
            continue
        try:
            result["families"][name] = inspect_family(
                name, families / name, runs, exe)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            result["families"][name] = {"status": "BLOCKED", "reasons": [str(exc)]}
    a, b = result["families"][first], result["families"][second]
    wall_nodes_identical = False  # Cannot assert from only aggregate quality metrics.
    geometry_a = {r.get("surface_parameterization")
                  for r in a.get("mesh_quality", [])}
    geometry_b = {r.get("surface_parameterization")
                  for r in b.get("mesh_quality", [])}
    result["comparison_scope"] = (
        "TWO_MESH_DISTRIBUTION_FAMILIES_AND_POLYGONAL_WALL_APPROXIMATIONS"
        if geometry_a != geometry_b
        else "SAME_WALL_PARAMETERIZATION_NON_IDENTICAL_NODES_NOT_VERIFIED"
    )
    result["causal_anisotropy_isolated"] = wall_nodes_identical
    if a.get("pressure_boundary_policy") != b.get("pressure_boundary_policy"):
        result["reasons"].append(
            "incompatible pressure BC policy across mesh families")
        return result
    if a["status"] != "GCI_CANDIDATE" or b["status"] != "GCI_CANDIDATE":
        result["reasons"].append(
            "two complete physically converged three-grid GCI candidates required")
        return result
    try:
        raw = compare(a["mesh_quality"], b["mesh_quality"],
                      a["solver_certificate"], b["solver_certificate"],
                      a["physical_setup"], b["physical_setup"])
        result["comparison"] = raw
        if raw["status"] == "CROSS_MESH_GCI_AGREEMENT_CANDIDATE":
            result["status"] = "CROSS_MESH_GCI_AGREEMENT_CANDIDATE"
        else:
            result["reasons"].extend(raw["reasons"])
    except (KeyError, ValueError, TypeError) as exc:
        result["reasons"].append(str(exc))
    return result


def markdown(result: dict) -> str:
    lines = ["# NACA0012 — GCI par contribution de traînée", "",
             "**Convention :** Cd = Cd_pressure + Cd_viscous. "
             "Trois solutions convergées avec le même schéma sont indispensables.",
             "",
             "| Famille | Niveau | Cd | Cd pression | Cd viscosité | Statut |",
             "|---|---:|---:|---:|---:|---|"]
    for name, fam in result["families"].items():
        for obs in fam.get("observations", []):
            fmt = lambda k: (f"{obs[k]:.9f}" if k in obs else "—")
            status = "exécuté" if obs["converged"] else "manquant/échec"
            lines.append(f"| {name} | N{obs['N']} | {fmt('cd')} | {fmt('cdp')} "
                         f"| {fmt('cdv')} | {status} |")
    lines.extend([
        "", "## GCI conditionnel à deux grilles (p=2 supposé, Fs=3)", "",
        "**Non certifié :** ces chiffres utilisent N64/N128 uniquement; "
        "leur ordre n’est pas mesuré.",
        "",
        "| Famille | Grandeur | GCI conditionnel (%) "
        "| Limite sous hypothèse p=2 (non validée) |",
        "|---|---|---:|---:|"])
    for name, fam in result["families"].items():
        provisional = fam.get("two_grid_assumed_order_only", {}).get(
            "components", {})
        for label, key in (("Cd", "cd"), ("Cd pression", "cdp"),
                           ("Cd viscosité", "cdv")):
            if key in provisional:
                c = provisional[key]
                lines.append(
                    f"| {name} | {label} "
                    f"| {100 * c['gci_relative_conditional']:.3f} "
                    f"| {c['hypothetical_order_based_limit']:.9f} |")
    lines.extend([
        "", "## GCI à ordre observé sur trois grilles", "",
        "| Famille | Grandeur | ordre observé | Extrapolation h→0 "
        "| GCI fin (%) | Statut |",
        "|---|---|---:|---:|---:|---|"])
    for name, fam in result["families"].items():
        for label, key in (("Cd", "cd"), ("Cd pression", "cdp"),
                           ("Cd viscosité", "cdv")):
            q = fam.get("gci", {}).get(key)
            if q and q.get("status") == "asymptotic_candidate":
                qualification = (
                    "candidat numérique"
                    if fam.get("status") == "GCI_CANDIDATE"
                    else "indicatif (autre composante non qualifiée)")
                lines.append(
                    f"| {name} | {label} | {q['order']:.4f} "
                    f"| {q['extrapolated']:.9f} "
                    f"| {100 * q['GCI_relative']:.3f} | {qualification} |")
            else:
                reason = (q or {}).get("status", "calcul non qualifié")
                lines.append(f"| {name} | {label} | — | — | — | {reason} |")
    lines.extend(["", "## Comparaison des limites extrapolées", ""])
    comparison = result.get("comparison", {}).get("gci_comparison", {})
    if comparison:
        lines.extend([
            "| Grandeur | Limite quasi-isotrope "
            "| Limite famille anisotrope (cosinus) | Écart relatif "
            "| GCI compatibles |",
            "|---|---:|---:|---:|---|"])
        for label, key in (("Cd", "cd"), ("Cd pression", "cdp"),
                           ("Cd viscosité", "cdv")):
            q = comparison[key]
            lines.append(
                f"| {label} | {q['extrapolated_a']:.9f} "
                f"| {q['extrapolated_b']:.9f} "
                f"| {100 * q['relative_difference']:.3f}% "
                f"| {q['gci_bands_overlap']} |")
    else:
        lines.append(
            "Comparaison **BLOQUÉE** tant que les deux familles ne présentent "
            "pas chacune trois solutions physiquement convergées et des "
            "tendances extrapolables séparément.")
    lines.extend([
        "", "## Réserves et contrôles", "",
        "- Les tendances et GCI sont des estimations de discrétisation, "
        "non une validation expérimentale.",
        "- La valeur Cd≈0,119 est une référence numérique publiée, sans "
        "décomposition Cp/Cf indépendante ici.",
        "- Les deux familles comparées utilisent des échantillonnages "
        "différents du profil analytique; leur accord ne prouve pas une "
        "invariance à anisotropie pure, à géométrie discrète strictement fixée.",
        "- L’arête de fermeture du bord de fuite (épaisseur NACA non nulle) "
        "reste un seul segment : ce défaut de raffinement géométrique doit "
        "être suivi.",
        "- La famille cosinus historique accroît fortement son aspect maximal "
        "au bord d’attaque.",
        "- Une compensation entre erreurs de pression et de frottement "
        "n’est jamais une preuve d’indépendance du Cd total.",
        "", "Statut global : **" + result["status"] + "**.", ""])
    for name, fam in result["families"].items():
        for reason in fam.get("reasons", []):
            lines.append(f"- {name}: {reason}")
    for reason in result["reasons"]:
        lines.append("- " + reason)
    lines.append("")
    return "\n".join(lines)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--families", type=Path, required=True)
    p.add_argument("--runs", type=Path, required=True)
    p.add_argument("--exe", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--first", default="quasi_isotropic")
    p.add_argument("--second", default="anisotropic")
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result = analyze(args.families, args.runs, args.exe, args.first, args.second)
    tmp = args.output_dir / "dual_family_gci.json.tmp"
    tmp.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, args.output_dir / "dual_family_gci.json")
    (args.output_dir / "dual_family_gci.md").write_text(markdown(result))
    print(json.dumps({
        "status": result["status"],
        "families": {k: v["status"] for k, v in result["families"].items()},
        "reasons": result["reasons"]}, indent=2))
    return 0 if result["status"] == "CROSS_MESH_GCI_AGREEMENT_CANDIDATE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
