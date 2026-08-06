"""Compact paired numeric inputs."""
from PySide6.QtWidgets import QDoubleSpinBox, QHBoxLayout, QLabel, QWidget


class RangeSelector(QWidget):
    """A minimum-to-maximum pair of floating point controls."""

    def __init__(self, minimum=0.0, maximum=100.0, parent=None):
        super().__init__(parent)
        self.minimum = self._spin(minimum)
        self.maximum = self._spin(maximum)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.minimum)
        layout.addWidget(QLabel('–', self))
        layout.addWidget(self.maximum)

    @staticmethod
    def _spin(value):
        spin = QDoubleSpinBox()
        spin.setRange(-1_000_000.0, 1_000_000.0)
        spin.setDecimals(4)
        spin.setSingleStep(0.1)
        spin.setValue(value)
        return spin

    def values(self):
        """Return the selected ``(minimum, maximum)`` pair."""
        return self.minimum.value(), self.maximum.value()
