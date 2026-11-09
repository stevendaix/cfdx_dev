"""Reproducible full-budget NACA/N9 and independent SST/thermal/BFS campaign.

By default no wall-clock timeout: each solver runs to its configured numerical
iteration budget. Every run has a distinct log, source/binary/mesh SHA256, and
explicit exit code. A successful CTest is NOT a physical certificate.

Examples:
  python scripts/run_cfdx_coupled_campaign.py --build build --meshes build/validation_meshes --out runs --only naca --models pcd block_local
  python scripts/run_cfdx_coupled_campaign.py --build build --meshes build/validation_meshes --out runs --only all
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

SCHUR_MODELS = ("pcd", "block_local", "lsc", "bfbt", "simple", "simplec")
N9_MODELS = SCHUR_MODELS
# (pressure-velocity algorithm, Schur approximation, outer preconditioner)
# The six default COUPLED variants are first; auxiliary PC variants follow.
# All runs use the same debug geometry, IC/BC and SST experiment settings.
BFS_VARIANTS = {
    "simple": ("simple", None, None),
    "simplec": ("simplec", None, None),
    "pcd_auto": ("coupled", "pcd", "auto"),
    "block_local_amg": ("coupled", "block_local", "coupled_block_schur"),
    "lsc_amg": ("coupled", "lsc", "coupled_block_schur"),
    "bfbt_amg": ("coupled", "bfbt", "coupled_block_schur"),
    "simple_schur_amg": ("coupled", "simple", "coupled_block_schur"),
    "simplec_schur_amg": ("coupled", "simplec", "coupled_block_schur"),
    "pcd_pcd": ("coupled", "pcd", "pcd"),
    "pcd_amg": ("coupled", "pcd", "coupled_block_schur"),
    "block_local_mgr": ("coupled", "block_local", "mgr"),
    "block_local_braess": ("coupled", "block_local", "braess_sarazin"),
}
CTESTS = {
    "bfs": ("test_bfs_quick",),
    "sst_thermal": (
        "test_turbulence_qualification_equations",
        "test_production_sst_thermal_coupling",
        "test_thermal_vv",
        "test_radiation_vv",
        "test_s2s_radiation_vv",
    ),
}


def verify_n9_shard(content: str, schur_model: str, exit_code: int) -> dict:
    """Reuse the strict physical oracle; an EXIT_ZERO is never enough."""
    script = Path(__file__).resolve().with_name("qualify_coupled_matrix.py")
    spec = importlib.util.spec_from_file_location("cfdx_n9_physical_gate", script)
    if spec is None or spec.loader is None:
        raise RuntimeError("N9 physical certificate implementation missing")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.summarize_n9_physical_shard(content, schur_model, exit_code)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        while chunk := file.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def execute(
    tag: str,
    cmd: list[str],
    env: dict[str, str],
    out: Path,
    *,
    wall_timeout: float | None,
) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    logfile = out / (tag + ".log")
    executable = Path(shutil.which(cmd[0]) or cmd[0]).resolve()
    present = executable.is_file()
    result = {
        "case": tag,
        "command": cmd,
        "binary_sha256": sha256(executable) if present else None,
        "mesh_sha256": sha256(Path(cmd[2]))
        if len(cmd) >= 3 and cmd[1] == "--single" and Path(cmd[2]).is_file()
        else (
            sha256(Path(cmd[cmd.index("--mesh") + 1]))
            if "--mesh" in cmd and Path(cmd[cmd.index("--mesh") + 1]).is_file()
            else None
        ),
        "env": {
            key: value for key, value in sorted(env.items()) if key.startswith("CFDX_")
        },
        "log": str(logfile),
        "exit_code": None,
        "wall_seconds": 0.0,
        "execution": "NOT_RUN",
        "physical_certificate": "NOT_QUALIFIED",
    }
    if not present:
        result["reason"] = "executable_missing"
        return result
    if result["mesh_sha256"] is None and len(cmd) >= 3 and cmd[1] == "--single":
        result["reason"] = "mesh_missing"
        return result
    start = time.monotonic()
    with logfile.open("w") as f:
        try:
            run = subprocess.run(
                cmd,
                env=env,
                stdout=f,
                stderr=subprocess.STDOUT,
                timeout=wall_timeout,
                check=False,
            )
            result["exit_code"] = run.returncode
            result["execution"] = "EXIT_ZERO" if run.returncode == 0 else "FAILED"
        except subprocess.TimeoutExpired:
            result["execution"] = "TIMEOUT_INCOMPLETE"
            result["reason"] = "explicit_wall_timeout"
        except OSError as exc:
            result["execution"] = "FAILED_TO_START"
            result["reason"] = str(exc)
    result["wall_seconds"] = time.monotonic() - start
    if tag.startswith("bfs_") and tag.removeprefix("bfs_") in BFS_VARIANTS:
        variant = tag.removeprefix("bfs_")
        algorithm, schur, preconditioner = BFS_VARIANTS[variant]
        from qualify_bfs_schur_evidence import validate_bfs_evidence
        output_flag = cmd.index("--output-dir") + 1 if "--output-dir" in cmd else -1
        run_dir = Path(cmd[output_flag]) if output_flag > 0 else out / "missing_run_dir"
        result["bfs_v8_gate"] = validate_bfs_evidence(
            run_dir,
            result["exit_code"],
            algorithm=algorithm,
            schur=schur,
            preconditioner=preconditioner,
            case_sha256=result["mesh_sha256"],
            case_path=Path(cmd[cmd.index("--mesh") + 1]) if "--mesh" in cmd else None,
        )
        if (
            result["execution"] == "EXIT_ZERO"
            and not result["bfs_v8_gate"]["nonlinear_v8_proven"]
        ):
            result["execution"] = "EXIT_ZERO_BUT_BFS_V8_EVIDENCE_INVALID"
    # Independently parse only current execution output. A historical log and
    # EXIT_ZERO with missing stages must not become qualification evidence.
    if tag.startswith("n9_") and result["exit_code"] is not None:
        with logfile.open(encoding="utf-8", errors="replace") as stream:
            current_log = stream.read()
        selected_model = tag.removeprefix("n9_")
        physical = verify_n9_shard(current_log, selected_model, result["exit_code"])
        result["n9_independent_shard_gate"] = physical
        result["physical_certificate"] = (
            "N9_FOUR_CASE_STRICT_GATE_PASSED_EXTERNAL_QUALIFICATION_STILL_REQUIRED"
            if physical["physical_shard_gate_complete"]
            else "NOT_QUALIFIED"
        )
        if (
            result["execution"] == "EXIT_ZERO"
            and not physical["physical_shard_gate_complete"]
        ):
            result["execution"] = "EXIT_ZERO_BUT_N9_PHYSICAL_EVIDENCE_INVALID"
    return result


def build_commands(args: argparse.Namespace):
    build = args.build.resolve()
    if args.only in ("naca", "all"):
        for n in args.levels:
            mesh = args.meshes.resolve() / f"naca0012_laminar_n{n}.h5"
            for model in args.models:
                tag = f"naca_n{n}_{model}"
                env = os.environ.copy()
                # Clear options that could silently change the numerical setup
                # between models/runs through inherited shell variables.
                for key in list(env):
                    if key.startswith(("CFDX_NACA_", "CFDX_PCD_FP_CAPTURE_")):
                        env.pop(key)
                env.update(
                    {
                        "CFDX_NACA_ALGORITHM": "COUPLED",
                        "CFDX_NACA_COUPLED_SCHUR": model,
                        "CFDX_NACA_LINEAR_TOLERANCE": "1e-12",
                        "CFDX_NACA_FARFIELD_PRESSURE_NEUMANN": "1",
                        "CFDX_NACA_REQUIRE_KINEMATIC_FACE_FLUX": "1",
                        "CFDX_NACA_CONVECTION_CELL_LIMITED": "1",
                        "CFDX_NACA_WALL_TRACTION_FACE": "1",
                        "CFDX_NACA_WALL_GRADIENT_WLS": "1",
                        "CFDX_NACA_PHYSICAL_HISTORY_CSV": str(
                            args.out.resolve() / (tag + "_history.csv")
                        ),
                        "CFDX_NACA_CONVERGED_FORCE_CSV": str(
                            args.out.resolve() / (tag + "_force.csv")
                        ),
                        "CFDX_NACA_RAW_WALL_FACES_CSV": str(
                            args.out.resolve() / (tag + "_raw_faces.csv")
                        ),
                    }
                )
                if getattr(args, "schur_action_guard", False):
                    # Same option for all six; never silently switch Schur.
                    env["CFDX_NACA_SCHUR_ACTION_GUARD"] = "1"
                if model == "pcd":
                    env["CFDX_PCD_FP_CAPTURE_PREFIX"] = str(
                        args.out.resolve() / (tag + "_Fp")
                    )
                    env["CFDX_PCD_FP_CAPTURE_DIFFICULT_PREFIX"] = str(
                        args.out.resolve() / (tag + "_difficult_Fp")
                    )
                yield (
                    tag,
                    [
                        str(build / "test_naca0012_qualification"),
                        "--single",
                        str(mesh),
                        "0.6",
                        "0.25",
                        str(args.iterations),
                    ],
                    env,
                )
    if args.only in ("n9", "all"):
        for model in args.n9_models:
            tag = f"n9_{model}"
            env = os.environ.copy()
            for key in list(env):
                if key.startswith(("CFDX_NACA_", "CFDX_PCD_FP_CAPTURE_")):
                    env.pop(key)
            env["CFDX_PCD_FP_CAPTURE_PREFIX"] = str(args.out.resolve() / (tag + "_Fp"))
            env["CFDX_PCD_FP_CAPTURE_DIFFICULT_PREFIX"] = str(
                args.out.resolve() / (tag + "_difficult_Fp")
            )
            yield (
                tag,
                [
                    str(build / "test_n9_physical_matrix"),
                    "--schur-model",
                    model,
                    "--require-independent-face",
                ],
                env,
            )
    if args.only in ("bfs", "all"):
        from generate_th03_bfs_case import reconfigure_bfs_reference_case, write_case
        reference = None
        for variant in getattr(args, "bfs_variants", list(BFS_VARIANTS)):
            algorithm, schur, preconditioner = BFS_VARIANTS[variant]
            tag = f"bfs_{variant}"
            case_file = args.out.resolve() / "bfs_cases" / f"{tag}.cfdx.h5"
            reference_case = getattr(args, "bfs_reference_case", None)
            if reference_case:
                generated = reconfigure_bfs_reference_case(
                    Path(reference_case),
                    case_file,
                    algorithm=algorithm,
                    schur=schur,
                    preconditioner=preconditioner,
                    iterations=args.iterations,
                )
            else:
                generated = write_case(
                    case_file,
                    n_up=getattr(args, "bfs_upstream", 4),
                    n_down=getattr(args, "bfs_downstream", 16),
                    ny=getattr(args, "bfs_ny", 10),
                    step_cells=getattr(args, "bfs_step_cells", 2),
                    wall_stretch=getattr(args, "bfs_wall_stretch", 0.0),
                    algorithm=algorithm,
                    schur=schur,
                    preconditioner=preconditioner,
                    iterations=args.iterations,
                )
            fingerprint = (generated["topology_hash"], generated["geometry_hash"])
            if reference is None:
                reference = fingerprint
            elif fingerprint != reference:
                raise ValueError(
                    "BFS cases do not share identical mesh topology and geometry"
                )
            run_dir = args.out.resolve() / "bfs_runs" / tag
            if run_dir.exists():
                # Never allow a prior successful run to satisfy this run's gate.
                shutil.rmtree(run_dir)
            env = os.environ.copy()
            for key in list(env):
                if key.startswith(("CFDX_NACA_", "CFDX_PCD_FP_CAPTURE_")):
                    env.pop(key)
            yield (
                tag,
                [
                    str(build / "cfdx_production_solver"),
                    "--mesh",
                    str(case_file),
                    "--output-dir",
                    str(run_dir),
                    "--iterations",
                    str(args.iterations),
                    "--experimental-sst",
                    "--benchmark-export",
                ],
                env,
            )
    for group in ("bfs", "sst_thermal"):
        if args.only not in (group, "all"):
            continue
        for test in CTESTS[group]:
            tag = f"{group}_{test}"
            env = os.environ.copy()
            yield (
                tag,
                [
                    "ctest",
                    "--test-dir",
                    str(build),
                    "-R",
                    f"^{test}$",
                    "--output-on-failure",
                ],
                env,
            )


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--build", type=Path, required=True)
    p.add_argument("--meshes", type=Path, default=Path("build/validation_meshes"))
    p.add_argument("--out", type=Path, required=True)
    p.add_argument(
        "--only", choices=("naca", "n9", "bfs", "sst_thermal", "all"), default="all"
    )
    p.add_argument(
        "--models", choices=SCHUR_MODELS, nargs="+", default=list(SCHUR_MODELS)
    )
    p.add_argument("--n9-models", choices=N9_MODELS, nargs="+", default=list(N9_MODELS))
    p.add_argument(
        "--levels", type=int, choices=(64, 128, 256), nargs="+", default=[64, 128, 256]
    )
    p.add_argument("--iterations", type=int, default=1500)
    p.add_argument(
        "--bfs-reference-case",
        type=Path,
        help="existing HDF5 BFS case with real inlet SST and bounded order-2 convection; preserve mesh and BCs",
    )
    p.add_argument(
        "--bfs-variants",
        choices=tuple(BFS_VARIANTS),
        nargs="+",
        default=list(BFS_VARIANTS),
        help="BFS debug campaign: six Schur plus SIMPLE/SIMPLEC and alternative block PCs",
    )
    p.add_argument("--bfs-upstream", type=int, default=4)
    p.add_argument("--bfs-downstream", type=int, default=16)
    p.add_argument("--bfs-ny", type=int, default=10)
    p.add_argument("--bfs-step-cells", type=int, default=2)
    p.add_argument("--bfs-wall-stretch", type=float, default=0.0)
    p.add_argument(
        "--source-archive", type=Path, help="optional immutable source ZIP; hash only"
    )
    p.add_argument(
        "--source-revision",
        help="declared Git commit (not an independently verified build link)",
    )
    p.add_argument(
        "--schur-action-guard",
        action="store_true",
        help="experimental common original-A action bound; no fallback; default OFF",
    )
    p.add_argument(
        "--wall-timeout",
        type=float,
        default=None,
        help="optional explicit local resource cap; default NONE",
    )
    args = p.parse_args(argv)
    if (
        args.iterations < 1
        or args.iterations > 100000
        or (args.wall_timeout is not None and args.wall_timeout <= 0)
    ):
        p.error("invalid numerical budget or runtime limit")
    if len(set(args.bfs_variants)) != len(args.bfs_variants):
        p.error("duplicate BFS variant would reuse evidence folder")
    args.out.mkdir(parents=True, exist_ok=True)
    records = []
    for tag, cmd, env in build_commands(args):
        print(f"CFDX_CAMPAIGN_START: {tag}", flush=True)
        result = execute(tag, cmd, env, args.out, wall_timeout=args.wall_timeout)
        records.append(result)
        (args.out / "coupled_matrix_campaign.json").write_text(
            json.dumps(
                {
                    "schema": "cfdx-full-coupled-matrix/v1",
                    "iterations": args.iterations,
                    "source_archive_sha256": sha256(args.source_archive)
                    if args.source_archive and args.source_archive.is_file()
                    else None,
                    "declared_source_revision": args.source_revision,
                    "binary_source_correspondence_proven": False,
                    "schur_action_guard": args.schur_action_guard,
                    "wall_timeout": args.wall_timeout,
                    "cases": records,
                    "physical_validation": "UNQUALIFIED_UNTIL_INDEPENDENT_EVIDENCE",
                },
                indent=2,
                allow_nan=False,
            )
            + "\n"
        )
        print(f"CFDX_CAMPAIGN_DONE: {tag} {result['execution']}", flush=True)
    software_ok = bool(records) and all(r["execution"] == "EXIT_ZERO" for r in records)
    bfs_ok = all(
        r.get("bfs_v8_gate", {}).get("nonlinear_v8_proven", False)
        for r in records
        if r["case"].removeprefix("bfs_") in BFS_VARIANTS
        and r["case"].startswith("bfs_")
    )
    return 0 if software_ok and bfs_ok else 2
if __name__ == "__main__":
    raise SystemExit(main())