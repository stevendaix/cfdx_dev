"""Replay exit-code contract: true-linear failure cannot be a green benchmark."""
from pathlib import Path
import subprocess
import sys
import tempfile


def write_capture(folder: Path, name: str, diagonal: bool) -> Path:
    matrix_name = f"{name}.mtx"
    vectors_name = f"{name}.vectors.csv"
    manifest = folder / f"{name}.manifest"
    entries = "1 1 2\n2 2 3\n" if diagonal else "1 2 -1\n2 1 -1\n"
    (folder / matrix_name).write_text(
        "%%MatrixMarket matrix coordinate real general\n2 2 2\n" + entries,
        encoding="utf-8",
    )
    (folder / vectors_name).write_text(
        "row,cell,rhs,initial_guess\n0,0,2,0\n1,1,3,0\n",
        encoding="utf-8",
    )
    manifest.write_text(
        "format=CFDX_PRESSURE_CAPTURE_V1\n"
        + f"matrix={matrix_name}\nvectors={vectors_name}\n"
        + "nonlinear_iteration=1\ncorrector=0\ncapture_sequence=1\n"
        + "mesh_cells=2\nrows=2\nnnz=2\nreference_cell=0\n"
        + "fixed_pressure_boundary=1\nprojected_constant_null_space=0\ngauge_reduced=0\n",
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Expected paths to sequence and single-matrix replay executables")
    with tempfile.TemporaryDirectory(prefix="cfdx_replay_gate_") as tmp:
        base = Path(tmp)
        valid = write_capture(base, "valid", diagonal=True)
        invalid = write_capture(base, "invalid", diagonal=False)
        cmd = [sys.argv[1], "fgmres_sa_fast", "fresh", "1"]
        success = subprocess.run(cmd + [str(valid)], capture_output=True, text=True)
        failure = subprocess.run(cmd + [str(invalid)], capture_output=True, text=True)
        if success.returncode != 0 or ",PASS," not in success.stdout:
            raise AssertionError(f"Valid replay did not pass: {success.stderr}\n{success.stdout}")
        if failure.returncode == 0 or ",FAIL," not in failure.stdout:
            raise AssertionError(f"Invalid replay incorrectly passed: {failure.stderr}\n{failure.stdout}")
        if "PRESSURE_SEQUENCE_REPLAY_FAIL" not in failure.stderr:
            raise AssertionError("Missing explicit failure diagnostic")
        single_cmd = [sys.argv[2], str(valid), "50", "1e-10", "1"]
        single_success = subprocess.run(single_cmd, capture_output=True, text=True)
        single_cmd[1] = str(invalid)
        single_failure = subprocess.run(single_cmd, capture_output=True, text=True)
        if single_success.returncode or ",PASS," not in single_success.stdout:
            raise AssertionError("Valid single-matrix replay did not pass")
        if single_failure.returncode == 0 or "PRESSURE_REPLAY_FAIL" not in single_failure.stderr:
            raise AssertionError("Invalid single-matrix replay incorrectly passed")
    print("NACA_PRESSURE_REPLAY_FAIL_CLOSED: PASS")


if __name__ == "__main__":
    main()
