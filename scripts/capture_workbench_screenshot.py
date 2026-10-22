"""Capture a deterministic screenshot of the CFDX Workbench for CI artifacts."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

# Resolve the repository's packaged Python tree explicitly. Executing a script
# from scripts/ puts that directory first on sys.path and can otherwise allow
# an unrelated installed module named "cfdx" to shadow the project package.
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PYTHON_PACKAGE_ROOT = REPOSITORY_ROOT / "src" / "cfdx" / "python"
sys.path.insert(0, str(PYTHON_PACKAGE_ROOT))

from PySide6.QtCore import QEventLoop, QTimer
from PySide6.QtWidgets import QApplication

from cfdx.workbench import CFDXWorkbenchWindow


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
    window.setWindowTitle("CFDX Workbench")
    screen = app.primaryScreen()
    if screen is None:
        print("No Qt screen available; run under Xvfb or a desktop session.", file=sys.stderr)
        window.close()
        return 2

    geometry = screen.availableGeometry()
    window.resize(min(1280, geometry.width()), min(900, geometry.height()))
    window.move(0, 0)
    window.show()
    window.raise_()
    window.activateWindow()

    # Let Qt paint the dock layout and let the optional VTK viewport initialize.
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

    print(f"Workbench screenshot saved: {args.output} ({screenshot.width()}x{screenshot.height()})")
    window.close()
    app.quit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
