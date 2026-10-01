"""Compact paired numeric inputs."""
from PySide6.QtWidgets import QDoubleSpinBox, QHBoxLayout, QLabel, QWidget


class RangeSelector(QWidget):
    """A minimum-to-maximum pair of floating point controls."""

    def __init__(self, minimum: float = 0.0, maximum: float = 100.0, parent: QWidget | None = None) -> None:
        """Create both controls side by side.

        Args:
           minimum: initial minimum
           maximum: initial maximum
           parent: parent widget
        """
        super().__init__(parent)
        self.minimum, self.maximum = QDoubleSpinBox(), QDoubleSpinBox()
        for spin, value in ((self.minimum, minimum), (self.maximum, maximum)):
            spin.setRange(-1_000_000.0, 1_000_000.0)
            spin.setDecimals(4)
            spin.setSingleStep(0.1)
            spin.setValue(value)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.minimum)
        layout.addWidget(QLabel('–', self))
        layout.addWidget(self.maximum)


    def values(self) -> tuple[float, float]:
        """Return the selected ``(minimum, maximum)`` pair.

        Returns:
           minimum, maximum
        """
        return self.minimum.value(), self.maximum.value()
