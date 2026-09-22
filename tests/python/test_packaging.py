from pathlib import Path
import subprocess
import sys


def test_package_builds(tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, "-m", "pip", "wheel", "--no-deps", ".", "-w", str(tmp_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert list(tmp_path.glob("cfdx-*.whl"))
