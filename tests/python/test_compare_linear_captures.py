"""Batch fail-closed certificates for captured linear solvers."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "tools/diagnostics/compare_linear_captures.py"
spec = importlib.util.spec_from_file_location("cfdx_replay_batch", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def write_fake_executable(path: Path, status: str):
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "print('preconditioner,krylov,status,iterations,total_ms,true_abs_l2,true_relative')\n"
        "for method in sys.argv[2:]:\n"
        f" print(method + ',fgmres,{status},5,1.0,1e-8,1e-10')\n",
        encoding="utf-8",
    )
    path.chmod(0o755)


def test_fail_closed_if_no_capture(tmp_path):
    executable = tmp_path / "fake.py"
    write_fake_executable(executable, "converged")
    rows, errors = module.compare(tmp_path, executable, ("jacobi",), 10)
    assert not rows and errors


def test_good_multi_capture_results_retained(tmp_path):
    executable = tmp_path / "fake.py"
    write_fake_executable(executable, "converged")
    (tmp_path / "n128.cfdxlin").write_text("mock", encoding="utf-8")
    (tmp_path / "n256.cfdxlin").write_text("mock", encoding="utf-8")
    rows, errors = module.compare(tmp_path, executable, ("jacobi", "native"), 10)
    assert not errors
    assert len(rows) == 4
    assert sorted({r["capture"] for r in rows}) == ["n128.cfdxlin", "n256.cfdxlin"]


def test_failed_candidate_not_dropped(tmp_path):
    executable = tmp_path / "fake.py"
    write_fake_executable(executable, "max_iter_reached")
    (tmp_path / "n256.cfdxlin").write_text("mock", encoding="utf-8")
    rows, errors = module.compare(tmp_path, executable, ("native",), 10)
    assert len(rows) == 1
    assert errors


def test_truncated_csv_rejected(tmp_path):
    executable = tmp_path / "fake.py"
    executable.write_text(
        "#!/usr/bin/env python3\nprint('invalid')\n", encoding="utf-8"
    )
    executable.chmod(0o755)
    (tmp_path / "n256.cfdxlin").write_text("mock", encoding="utf-8")
    rows, errors = module.compare(tmp_path, executable, ("native",), 10)
    assert not rows and errors


def test_preserve_diagnostic_cause_for_coarse_amg_failure(tmp_path):
    executable = tmp_path / "failure.py"
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "print('preconditioner,krylov,status,iterations,total_ms,true_abs_l2,true_relative')\n"
        "print('directional,fgmres,not_applicable,0,1.0,2e-2,1')\n"
        "sys.stderr.write('PRECONDITIONER_SETUP_ERROR directional: invalid coarse level 7\\n')\n"
        "sys.exit(1)\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    (tmp_path / "n256.cfdxlin").write_text("mock", encoding="utf-8")
    rows, errors = module.compare(tmp_path, executable, ("directional",), 10)
    assert errors
    assert len(rows) == 1
    assert rows[0]["status"] == "not_applicable"
    assert rows[0]["failure_reason"] == "invalid coarse level 7"
