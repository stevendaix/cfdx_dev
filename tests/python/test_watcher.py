from pathlib import Path
import time

from cfdx.watcher import ResultWatcher


def test_result_watcher_detects_file_change(tmp_path: Path) -> None:
    events = []
    watcher = ResultWatcher(tmp_path, lambda: events.append("changed"), interval=0.01)
    watcher.start()
    try:
        (tmp_path / "result.vtu").write_text("data", encoding="utf-8")
        deadline = time.monotonic() + 1.0
        while not events and time.monotonic() < deadline:
            time.sleep(0.01)
        assert events == ["changed"]
    finally:
        watcher.stop()
