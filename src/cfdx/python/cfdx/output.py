"""Bounded solver output delivery for interactive clients."""
from __future__ import annotations

from collections import deque
import threading
import time
from typing import Callable


class OutputThrottle:
    """Batch output callbacks so high-volume solver logs do not flood a UI."""

    def __init__(self, emit: Callable[[str, bool], None], interval: float = 0.05, max_lines: int = 100) -> None:
        if interval <= 0:
            raise ValueError("interval must be positive")
        if max_lines <= 0:
            raise ValueError("max_lines must be positive")
        self._emit = emit
        self._interval = interval
        self._max_lines = max_lines
        self._queue: deque[tuple[str, bool]] = deque()
        self._lock = threading.Lock()
        self._last_flush = 0.0

    def push(self, line: str, is_stderr: bool = False) -> None:
        with self._lock:
            self._queue.append((line, is_stderr))

    def flush(self, now: float | None = None) -> int:
        now = time.monotonic() if now is None else now
        with self._lock:
            if not self._queue or now - self._last_flush < self._interval:
                return 0
            batch = [self._queue.popleft() for _ in range(min(self._max_lines, len(self._queue)))]
            self._last_flush = now
        for line, is_stderr in batch:
            self._emit(line, is_stderr)
        return len(batch)

    def pending(self) -> int:
        with self._lock:
            return len(self._queue)
