"""Optional Qt controls for a ResultSeries."""
from __future__ import annotations
from PySide6.QtCore import QTimer, Qt, Signal
from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QPushButton, QSlider, QVBoxLayout, QWidget
from .results_series import ResultSeries

class ResultsSeriesPanel(QWidget):
    frame_changed = Signal(object)

    def __init__(self,parent: QWidget | None=None) -> None:
        super().__init__(parent); self.series: ResultSeries | None=None
        self.play=QPushButton("Play"); self.speed=QComboBox(); self.speed.addItems(["0.5x","1x","2x","4x"])
        self.slider=QSlider(); self.slider.setOrientation(Qt.Orientation.Horizontal); self.slider.valueChanged.connect(self._select)
        self.status=QLabel("No result series loaded")
        row=QHBoxLayout(); row.addWidget(self.play); row.addWidget(self.speed)
        layout=QVBoxLayout(self); layout.addWidget(self.status); layout.addLayout(row); layout.addWidget(self.slider)
        self.timer=QTimer(self); self.timer.timeout.connect(self._advance); self.play.clicked.connect(self._toggle)

    def set_series(self, series: ResultSeries) -> None:
        self.series=series; self.slider.setRange(0,max(0,len(series.frames)-1)); self.slider.setValue(0)
        self._update_status()

    def _select(self,index: int) -> None:
        if self.series is not None and 0 <= index < len(self.series.frames):
            self.status.setText(self._label(index)); self.frame_changed.emit(self.series.frame(index))

    def _label(self,index: int) -> str:
        frame=self.series.frame(index); time=f"t={frame.time:g}" if frame.time is not None else "time=unknown"
        return f"{index+1}/{len(self.series.frames)} — {time} — {frame.path.name}"

    def _update_status(self) -> None:
        if self.series is None: self.status.setText("No result series loaded")
        elif self.series.frames: self.status.setText(self._label(self.slider.value()))

    def _toggle(self) -> None:
        if self.timer.isActive(): self.timer.stop(); self.play.setText("Play")
        else: self.timer.start(self._interval()); self.play.setText("Pause")

    def _interval(self) -> int:
        speed=float(self.speed.currentText()[:-1]); return max(25,int(500/speed))

    def _advance(self) -> None:
        if self.series is None or not self.series.frames: return
        nxt=(self.slider.value()+1)%len(self.series.frames); self.slider.setValue(nxt)
