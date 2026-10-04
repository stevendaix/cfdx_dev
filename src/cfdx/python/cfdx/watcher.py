"""Headless-safe result directory change detection."""
from __future__ import annotations

from pathlib import Path
import threading
import time
from typing import Callable


class ResultWatcher:
    """Poll a result directory and notify when its signature changes."""

    def __init__(
        self,
        directory: Path,
        callback: Callable[[], None],
        interval: float = 0.5,
    ) -> None:
        if interval <= 0:
            raise ValueError("interval must be positive")
        self.directory = directory
        self.callback = callback
        self.interval = interval
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._signature: tuple[tuple[str, int], ...] = ()

    def _signature_for(self) -> tuple[tuple[str, int], ...]:
        if not self.directory.exists():
            return ()
        return tuple(
            sorted(
                (str(path.relative_to(self.directory)), path.stat().st_mtime_ns)
                for path in self.directory.rglob("*")
                if path.is_file()
            )
        )

    def start(self) -> None:
        if self._thread and self._thread.is_alive():
            return
        self._signature = self._signature_for()
        self._stop.clear()
        self._thread = threading.Thread(target=self._watch, name="cfdx-result-watcher", daemon=True)
        self._thread.start()

    def _watch(self) -> None:
        while not self._stop.wait(self.interval):
            signature = self._signature_for()
            if signature != self._signature:
                self._signature = signature
                self.callback()

    def stop(self, timeout: float = 2.0) -> None:
        self._stop.set()
        if self._thread and self._thread is not threading.current_thread():
            self._thread.join(timeout)
