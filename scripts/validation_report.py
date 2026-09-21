#!/usr/bin/env python3
"""Generate the CFDX Fluent/VMFL validation report.

The report is intentionally data-driven:
- the VMFL matrix is parsed from docs/validation/FLUENT_VMFL_MATRIX.md;
- executable validation tests are run when available;
- raw stdout is archived next to the report;
- Python computes comparison metrics and plots;
- LaTeX is the final PDF renderer, driven by latexmk.

This keeps numerical data generation separate from document formatting and makes
the report reproducible in CI.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


@dataclass
class Case:
    ident: str
    name: str
    status: str
    comparison: str


def parse_matrix(path: Path) -> list[Case]:
    cases: list[Case] = []
    pattern = re.compile(r"^\|\s*(VMFL\d{3})\s*\|\s*(.*?)\s*\|\s*(.*?)\s*\|\s*(.*?)\s*\|$")
    for line in path.read_text(encoding="utf-8").splitlines():
        m = pattern.match(line)
        if m:
            cases.append(Case(*m.groups()))
    if len(cases) != 78:
        raise RuntimeError(f"Expected 78 VMFL cases, found {len(cases)} in {path}")
    return cases


def run_test(executable: Path, log_path: Path) -> tuple[int, str]:
    proc = subprocess.run(
        [str(executable)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log_path.write_text(proc.stdout, encoding="utf-8")
    return proc.returncode, proc.stdout


def parse_ghia(output: str) -> list[dict[str, float | str]]:
    pattern = re.compile(
        r"GHIA Re=(?P<re>[0-9.]+) grid=(?P<nx>\d+)x(?P<ny>\d+)"
        r".*?continuity=(?P<continuity>[-+0-9.eE]+)"
        r".*?continuity_norm=(?P<continuity_norm>[-+0-9.eE]+)"
        r".*?momentum_eq=(?P<momentum>[-+0-9.eE]+)"
        r".*?U_RMS=(?P<u_rms>[-+0-9.eE]+) U_max=(?P<u_max>[-+0-9.eE]+)"
        r" V_RMS=(?P<v_rms>[-+0-9.eE]+) V_max=(?P<v_max>[-+0-9.eE]+)"
    )
    result = []
    for m in pattern.finditer(output):
        d = m.groupdict()
        result.append({
            "re": float(d["re"]),
            "grid": f"{d['nx']}x{d['ny']}",
            "continuity": float(d["continuity"]),
            "continuity_norm": float(d["continuity_norm"]),
            "momentum": float(d["momentum"]),
            "u_rms": float(d["u_rms"]),
            "u_max": float(d["u_max"]),
            "v_rms": float(d["v_rms"]),
            "v_max": float(d["v_max"]),
        })
    return result


def parse_model_results(output: str) -> list[dict[str, str]]:
    pattern = re.compile(r"MODEL (?P<name>[A-Z0-9_]+) (?P<metric>[A-Za-z0-9_]+)=(?P<value>[-+0-9.eE]+) reference=(?P<reference>.*)")
    return [m.groupdict() for m in pattern.finditer(output)]


def run_ghia(executable: Path, log_path: Path) -> tuple[int, str, list[dict]]:
    rc, output = run_test(executable, log_path)
    return rc, output, parse_ghia(output)


def latex_escape(value: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#",
        "_": r"\_", "{": r"\{", "}": r"\}",
    }
    for old, new in replacements.items():
        value = value.replace(old, new)
    return value


def fmt(value: float | None, digits: int = 5) -> str:
    if value is None or not math.isfinite(value):
        return "n/a"
    return f"{value:.{digits}g}"


def make_plot(out: Path, ghia: list[dict]) -> str | None:
    if not ghia:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    re100 = [x for x in ghia if x["re"] == 100.0]
    if not re100:
        return None

    grids = [x["grid"] for x in re100]
    u = [x["u_rms"] for x in re100]
    v = [x["v_rms"] for x in re100]
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    x = range(len(grids))
    ax.plot(x, u, marker="o", label="U centreline RMS")
    ax.plot(x, v, marker="s", label="V centreline RMS")
    ax.set_xticks(list(x), grids)
    ax.set_xlabel("Mesh")
    ax.set_ylabel("RMS error against Ghia")
    ax.set_title("CFDX vs Ghia: mesh refinement, Re=100")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out.name


def write_tex(path: Path, cases: list[Case], ghia: list[dict], model_results: list[dict], logs: dict, suite_status: dict[str, int], plot_name: str | None, generated: str) -> None:
    counts: dict[str, int] = {}
    for c in cases:
        counts[c.status] = counts.get(c.status, 0) + 1

    active = [c for c in cases if c.status in {"READY", "PARTIAL"}]
    reference_values = [
        ("VMFL003", "Pressure drop", "21,744 Pa", "Ansys/White Moody-chart target"),
        ("VMFL004", "Mean velocity", "2.50 m/s", "Exact Couette + pressure-gradient solution"),
        ("VMFL005", "Pressure drop", "10.24 Pa", "Hagen--Poiseuille"),
    ]

    lines = [
        r"\documentclass[10pt,a4paper]{article}",
        r"\usepackage[margin=20mm]{geometry}",
        r"\usepackage{booktabs,longtable,array,graphicx,xcolor,hyperref}",
        r"\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue}",
        r"\setlength{\parindent}{0pt}",
        r"\setlength{\parskip}{5pt}",
        r"\begin{document}",
        r"\begin{titlepage}",
        r"\centering",
        r"{\LARGE\bfseries CFDX Validation and Verification Report\par}",
        r"\vspace{8mm}",
        r"{\Large Ansys Fluent VMFL-aligned campaign\par}",
        r"\vspace{12mm}",
        f"Generated: {latex_escape(generated)}\\",
        r"Reference basis: Ansys Fluid Dynamics Verification Manual, Release 2026 R1.",
        r"\vfill",
        r"\textbf{Important:} a reference oracle is not counted as a CFDX solver PASS.",
        r"\end{titlepage}",
        r"\tableofcontents",
        r"\newpage",
        r"\section{Executive summary}",
        f"The campaign contains {len(cases)} Fluent VMFL cases. Current matrix status: "
        f"READY={counts.get('READY',0)}, PARTIAL={counts.get('PARTIAL',0)}, "
        f"BLOCKED={counts.get('BLOCKED',0)}, PASS={counts.get('PASS',0)}.",
        r"",
        r"The report distinguishes three layers: (1) reference values and analytical oracles, "
        r"(2) executable CFDX solver comparisons, and (3) capability gaps. This prevents a "
        r"closed-form calculation from being presented as a complete CFD validation.",
        r"\section{Reference values currently encoded}",
        r"\begin{tabular}{llll}",
        r"\toprule Case & Quantity & Reference & Basis\\",
        r"\midrule",
    ]
    for row in reference_values:
        lines.append(" & ".join(latex_escape(x) for x in row) + r"\\")
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"",
        r"VMFL003--005 parameters and targets are taken from the Ansys verification material; "
        r"VMFL004 and VMFL005 also have exact analytical solutions in the CFDX reference suite.",
        r"\section{Executable CFDX results}",
    ]

    if plot_name:
        lines += [
            r"\begin{figure}[h]",
            r"\centering",
            f"\includegraphics[width=0.92\textwidth]{{{latex_escape(plot_name)}}}",
            r"\caption{Ghia Re=100 centreline RMS error reported by the CFDX solver at three meshes.}",
            r"\end{figure}",
        ]

    lines += [
        r"\\section{Validation executable matrix}",
        r"\\begin{longtable}{p{70mm}r}",
        r"\\toprule Executable & return code \\\\",
        r"\\midrule",
    ]
    for name, rc in suite_status.items():
        lines.append(f"\\texttt{{{latex_escape(name)}}} & {rc}\\\\")
    lines += [r"\\bottomrule", r"\\end{longtable}"]

    if ghia:
        lines += [
            r"\subsection{Ghia lid-driven cavity}",
            r"Reference: Ghia et al. centreline velocity tables. The table below reports the "
            r"CFDX solver result directly; continuity and physical momentum residuals are "
            r"reported independently from the profile error.",
            r"\begin{longtable}{rrrrrrrr}",
            r"\toprule Re & Mesh & $U_{RMS}$ & $U_{max}$ & $V_{RMS}$ & $V_{max}$ & cont. & mom.\\",
            r"\midrule",
        ]
        for x in ghia:
            lines.append(
                f"{fmt(x['re'],4)} & {latex_escape(x['grid'])} & {fmt(x['u_rms'])} & "
                f"{fmt(x['u_max'])} & {fmt(x['v_rms'])} & {fmt(x['v_max'])} & "
                f"{fmt(x['continuity'])} & {fmt(x['momentum'])}\\"
            )
        lines += [r"\bottomrule", r"\end{longtable}"]
    else:
        lines.append("No Ghia solver output was available in this report run.")

    lines += [
        r"\subsection{Reference-oracle execution}",
        r"The Fluent VMFL reference executable checks the published VMFL003--005 setup and "
        r"analytical targets. It is deliberately labelled as an oracle until the corresponding "
        r"CFDX geometry, boundary conditions and coupled solver are executed.",
    ]
    for ident, (rc, _) in logs.items():
        lines.append(f"\texttt{{{latex_escape(ident)}}}: return code {rc}.")

    lines += [
        r"\section{Numerical model verification}",
        r"The table below is populated directly from the deterministic numerical-model executable.",
        r"\begin{longtable}{p{38mm}p{30mm}p{30mm}p{65mm}}",
        r"\toprule Model & Metric & Value/error & Reference\\",
        r"\midrule",
    ]
    for m in model_results:
        lines.append(
            f"\texttt{{{latex_escape(m['name'])}}} & {latex_escape(m['metric'])} & "
            f"{latex_escape(m['value'])} & {latex_escape(m['reference'])}\\"
        )
    lines += [
        r"\bottomrule",
        r"\end{longtable}",
        r"\section{Full VMFL capability matrix}",

        r"\small",
        r"\begin{longtable}{p{12mm}p{48mm}p{18mm}p{72mm}}",
        r"\toprule ID & Fluent case & CFDX status & Main comparison\\",
        r"\midrule",
    ]
    for c in cases:
        lines.append(
            f"\texttt{{{c.ident}}} & {latex_escape(c.name)} & {latex_escape(c.status)} & "
            f"{latex_escape(c.comparison)}\\"
        )
    lines += [
        r"\bottomrule",
        r"\end{longtable}",
        r"\normalsize",
        r"\section{Acceptance rules}",
        r"Every solver-level case must report geometry/mesh, physical properties, boundary "
        r"conditions, discretisation and solver controls, residuals, conservation error, QoI, "
        r"reference value or curve, absolute/relative error, and refinement where meaningful.",
        r"",
        r"An executable case is PASS only when CFDX actually solves the physical configuration "
        r"and satisfies a quantitative criterion. READY/PARTIAL/BLOCKED are capability states, "
        r"not numerical passes.",
        r"\section{Reproducibility}",
        f"The raw executable logs are stored beside this report in the report output directory. "
        f"Report generation timestamp: {latex_escape(generated)}.",
        r"",
        r"Source matrix: \texttt{docs/validation/FLUENT\_VMFL\_MATRIX.md}.\\",
        r"Ansys reference: \url{https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/pdf/Ansys_Fluid_Dynamics_Verification_Manual.pdf}.",
        r"\end{document}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


VALIDATION_EXECUTABLES = [
    "test_analytical_benchmarks",
    "test_benchmark_matrix",
    "test_cht_validation",
    "test_fluent_vmfl_reference",
    "test_ghia_cavity",
    "test_level_b_reference_benchmarks",
    "test_level_c_coupled_verification",
    "test_m1_m4_validation",
    "test_m2_m4_solver_validation",
    "test_numerical_model_verification",
    "test_steady_incompressible_solver",
]

def run_validation_suite(build_dir: Path, output_dir: Path) -> dict[str, int]:
    status: dict[str, int] = {}
    for name in VALIDATION_EXECUTABLES:
        exe = build_dir / name
        if not exe.exists():
            status[name] = -1
            (output_dir / f"{name}.log").write_text("executable not found", encoding="utf-8")
            continue
        rc, _ = run_test(exe, output_dir / f"{name}.log")
        status[name] = rc
    return status


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    parser.add_argument("--output-dir", type=Path, default=Path("build/validation-report"))
    parser.add_argument("--matrix", type=Path, default=Path("docs/validation/FLUENT_VMFL_MATRIX.md"))
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    cases = parse_matrix(args.matrix)
    generated = datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%d %H:%M %Z")

    logs: dict[str, tuple[int, str]] = {}
    suite_status = run_validation_suite(args.build_dir, args.output_dir)
    ref_exe = args.build_dir / "test_fluent_vmfl_reference"
    if ref_exe.exists():
        logs["test_fluent_vmfl_reference"] = run_test(ref_exe, args.output_dir / "test_fluent_vmfl_reference.log")
    else:
        logs["test_fluent_vmfl_reference"] = (-1, "executable not found")

    model_results: list[dict] = []
    model_exe = args.build_dir / "test_numerical_model_verification"
    if model_exe.exists():
        rc, output = run_test(model_exe, args.output_dir / "test_numerical_model_verification.log")
        model_results = parse_model_results(output)
        logs["test_numerical_model_verification"] = (rc, output)
    else:
        logs["test_numerical_model_verification"] = (-1, "executable not found")

    ghia: list[dict] = []
    ghia_exe = args.build_dir / "test_ghia_cavity"
    if ghia_exe.exists():
        rc, output, ghia = run_ghia(ghia_exe, args.output_dir / "test_ghia_cavity.log")
        logs["test_ghia_cavity"] = (rc, output)
    else:
        logs["test_ghia_cavity"] = (-1, "executable not found")

    (args.output_dir / "results.json").write_text(
        json.dumps({"generated": generated, "ghia": ghia, "model_results": model_results,
                    "test_status": {k: v[0] for k, v in logs.items()},
                    "validation_suite_status": suite_status}, indent=2),
        encoding="utf-8",
    )

    plot_name = make_plot(args.output_dir / "ghia_mesh_convergence.png", ghia)

    tex = args.output_dir / "cfdx_validation_report.tex"
    write_tex(tex, cases, ghia, model_results, logs, suite_status, plot_name, generated)

    import shutil
    latexmk = shutil.which("latexmk")
    if latexmk is None:
        raise RuntimeError("latexmk is required to build the PDF")

    subprocess.run(
        [latexmk, "-pdf", "-interaction=nonstopmode", "-halt-on-error", tex.name],
        cwd=args.output_dir,
        check=True,
    )
    pdf = args.output_dir / "cfdx_validation_report.pdf"
    if not pdf.exists():
        raise RuntimeError(f"LaTeX completed without creating {pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
