"""Capture a deterministic, populated CFDX Workbench screenshot for CI artifacts.

The screenshot opens a small synthetic HDF5 case in the metadata inspector so
that the exported image demonstrates the actual case-browser workflow rather
than only an empty application shell. The fixture is deliberately tiny and is
never used as a solver qualification case.
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

# Resolve the repository's packaged Python tree explicitly. Executing a script
# from scripts/ puts that directory first on sys.path and can otherwise allow
# an unrelated installed module named "cfdx" to shadow the project package.
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PYTHON_PACKAGE_ROOT = REPOSITORY_ROOT / "src" / "cfdx" / "python"
sys.path.insert(0, str(PYTHON_PACKAGE_ROOT))

import h5py
import numpy as np
from cfdx.case_inspector import inspect_cfdx_file
from cfdx.workbench import CFDXWorkbenchWindow
from PySide6.QtCore import QEventLoop, Qt, QTimer
from PySide6.QtWidgets import QApplication


def _create_demo_case(path: Path) -> None:
    """Write a small, valid interchange case for a representative UI capture."""
    setup = {
        "name": "CI cavity demonstration",
        "physics": {"flow": {"enabled": True}},
        "numerics": {"algorithm": "SIMPLE"},
        "boundaries": {"wall": {}, "inlet": {}, "outlet": {}},
    }
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["schema_version"] = "1"
        h5.attrs["n_cells"] = 1
        h5.attrs["source_case_name"] = "CI cavity demonstration"
        h5.attrs["case_setup_json"] = json.dumps(setup)
        h5.attrs["units"] = "SI"

        h5.create_dataset(
            "points",
            data=np.asarray(
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
                dtype=np.float64,
            ),
        )
        h5.create_dataset("owner", data=np.asarray([0], dtype=np.int64))
        h5.create_dataset("neighbour", data=np.asarray([-1], dtype=np.int64))
        h5.create_dataset(
            "face_vertices", data=np.asarray([0, 1, 2, 3], dtype=np.int64)
        )
        h5.create_dataset("face_offsets", data=np.asarray([0, 4], dtype=np.int64))
        h5.create_dataset("patch_face_ids", data=np.asarray([0], dtype=np.int64))
        h5.create_dataset("patch_face_offsets", data=np.asarray([0, 1], dtype=np.int64))

        scalar = h5.create_group("fields/scalar")
        scalar.create_dataset(
            "pressure", data=np.asarray([101325.0, 101320.0, 101315.0, 101310.0])
        )
        vector = h5.create_group("fields/vector")
        vector.create_dataset(
            "velocity",
            data=np.asarray(
                [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.1, 0.1, 0.0], [0.0, 0.1, 0.0]]
            ),
        )
        h5.create_dataset("metadata/source", data=np.bytes_("CI synthetic fixture"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/cfdx-workbench.png"),
        help="PNG path to write (default: artifacts/cfdx-workbench.png)",
    )
    args = parser.parse_args()

    app = QApplication.instance() or QApplication(["cfdx-workbench-screenshot"])
    app.setApplicationName("CFDX Workbench Screenshot")
    window = CFDXWorkbenchWindow()
    window.setWindowTitle("CFDX Workbench — Case Inspector")
    screen = app.primaryScreen()
    if screen is None:
        print(
            "No Qt screen available; run under Xvfb or a desktop session.",
            file=sys.stderr,
        )
        window.close()
        return 2

    with tempfile.TemporaryDirectory(prefix="cfdx-workbench-demo-") as directory:
        demo_path = Path(directory) / "cavity-demo.cfdx.h5"
        _create_demo_case(demo_path)
        summary = inspect_cfdx_file(demo_path)
        window.case_file_panel.set_summary(summary)
        window.case_file_panel.tree.expandAll()
        pressure_items = window.case_file_panel.tree.findItems(
            "/fields/scalar/pressure",
            Qt.MatchFlag.MatchExactly | Qt.MatchFlag.MatchRecursive,
            0,
        )
        if not pressure_items:
            print(
                "Synthetic case did not populate the pressure dataset row.",
                file=sys.stderr,
            )
            window.close()
            app.quit()
            return 3
        window.case_file_panel.tree.setCurrentItem(pressure_items[0])
        if not window.case_file_panel.mesh_quality_button.isEnabled():
            raise RuntimeError("surface quality action was not enabled for demo mesh")
        window.case_file_panel.mesh_quality_button.click()
        quality_report = window.case_file_panel.preview.toPlainText()
        if "faces: 1 · valid: 1 · degenerate: 0" not in quality_report:
            raise RuntimeError("surface quality action did not produce expected report")

        geometry = screen.availableGeometry()
        window.resize(min(1280, geometry.width()), min(900, geometry.height()))
        window.move(0, 0)
        window.show()
        window.raise_()
        window.activateWindow()

        # Let Qt paint the populated metadata tree and bounded sample preview.
        app.processEvents()
        loop = QEventLoop()
        QTimer.singleShot(1200, loop.quit)
        loop.exec()
        app.processEvents()

        args.output.parent.mkdir(parents=True, exist_ok=True)
        screenshot = screen.grabWindow(window.winId())
        if screenshot.isNull():
            print("Qt returned an empty Workbench screenshot.", file=sys.stderr)
            window.close()
            app.quit()
            return 1
        if not screenshot.save(str(args.output), "PNG"):
            print(f"Failed to write screenshot to {args.output}", file=sys.stderr)
            window.close()
            app.quit()
            return 1

        print(
            f"Populated Workbench screenshot saved: {args.output} "
            f"({screenshot.width()}x{screenshot.height()}); "
            f"datasets={len(summary.datasets)}, cells={summary.n_cells}, "
            f"sample=/fields/scalar/pressure, surface_quality=passed"
        )

    window.close()
    app.quit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
