"""Main window for guided EBSD plotting and reproducible code export."""
from __future__ import annotations

from pathlib import Path
from typing import Any

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox,
                               QDoubleSpinBox, QFileDialog, QFormLayout,
                               QHBoxLayout, QLabel, QLineEdit, QMainWindow,
                               QMessageBox, QPushButton, QSpinBox, QSplitter,
                               QTabWidget, QTextEdit, QVBoxLayout, QWidget)

from ..ebsd import EBSD
from ..fileIO import LOADERS
from .rangeSelector import RangeSelector


class EBSDGui(QMainWindow):
    """A compact UI for the most common ebsdlab plot workflows."""

    def __init__(self) -> None:
        """Create the window with an empty plot area."""
        super().__init__()
        self.filePath: Path | None                      = None
        self.ebsd:     EBSD | None                      = None
        self.figure:   Any                              = None
        self.canvas:   FigureCanvasQTAgg | None         = None
        self.toolbar:  NavigationToolbar2QT | None      = None
        self.placeholder: QLabel | None                 = None
        self.overlays: list[tuple[float, float, float]] = []
        self.setWindowTitle('ebsdlab')
        self.resize(1280, 800)
        self._buildUi()


    def setFile(self, filePath: str | Path) -> None:
        """Select an input file while displaying only its filename in the UI.

        Args:
           filePath: path of the EBSD file
        """
        self.filePath = Path(filePath).expanduser().resolve()
        self.fileName.setText(self.filePath.name)
        self.fileName.setToolTip(str(self.filePath))


    def chooseFile(self) -> None:
        """Prompt for an EBSD file and create the selected plot."""
        name, _ = QFileDialog.getOpenFileName(
            self, 'Open EBSD file', '', 'EBSD files (*.ang *.osc *.txt *.crc)')
        if name:
            self.setFile(name)
            self.createPlot()


    def createPlot(self) -> None:
        """Load the selected data and render the chosen plot style."""
        try:
            ebsd = self._loadEbsd()
            self._applyFilters(ebsd)
            self.overlays = []
            plotType = self.plotType.currentText()
            if plotType == 'CI map':
                self.figure = ebsd.plot(ebsd.ci, show=False)
            elif plotType == 'IPF map':
                self.figure = ebsd.plotIPF(
                    direction=self.subplotStyle.currentText(), show=False)
            else:
                self.figure = ebsd.plotPF(
                    axis=[int(value) for value in self.subplotStyle.currentText().strip('[]').split(',')],
                    show=False)
            axis = self.figure.axes[0]
            if plotType != 'Pole figure':
                self.figure.subplots_adjust(left=0.0125, right=0.99,
                                             bottom=0.0125, top=0.99)
            if self.hideAxes.isChecked() and plotType != 'Pole figure':
                axis.axis('off')
            if self.scaleBar.isChecked() and plotType != 'Pole figure':
                ebsd.addScaleBarOverlay(axis)
            self._setFigure(self.figure)
            self._refreshCode()
            self._updateDetailVisibility()
            self.statusBar().showMessage('Plot created')
        except (OSError, ValueError, IndexError) as error:
            QMessageBox.critical(self, 'Could not create plot', str(error))
            self.statusBar().showMessage('Plot creation failed')


    def addUnitCellOverlay(self, event: Any) -> None:
        """Add a unit-cell overlay at a click in an IPF plot.

        Args:
           event: matplotlib mouse event
        """
        if self.ebsd is None or self.canvas is None or self.plotType.currentText() != 'IPF map':
            return
        if event.inaxes is None or event.xdata is None or event.ydata is None:  # click outside the map
            return
        scale = self.overlayScale.value()
        self.ebsd.addUnitCellOverlay(event.inaxes, event.xdata, event.ydata, scale)
        self.overlays.append((event.xdata, event.ydata, scale))
        self.canvas.draw_idle()
        self._refreshCode()
        self._updateDetailVisibility()
        self.statusBar().showMessage('Unit-cell overlay added')


    def clearUnitCellOverlays(self) -> None:
        """Recreate the plot without its unit-cell overlays."""
        if self.figure is not None and self.overlays:
            self.createPlot()


    def copyPythonCode(self) -> None:
        """Copy the generated, reproducible plotting code to the clipboard."""
        if self.code.toPlainText():
            QApplication.clipboard().setText(self.code.toPlainText())
            self.statusBar().showMessage('Python code copied to clipboard')


    def _buildUi(self) -> None:
        """Create the sidebar (file, plot style, code and plot-detail tabs) and the plot area."""
        central = QWidget(self)
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        splitter = QSplitter(Qt.Orientation.Horizontal, central)
        root.addWidget(splitter)

        sidebar = QWidget(splitter)
        sidebarLayout = QVBoxLayout(sidebar)
        mainForm = QFormLayout()
        self.fileName = QLineEdit(sidebar)
        self.fileName.setReadOnly(True)
        browse = QPushButton('Browse…', sidebar)
        browse.clicked.connect(self.chooseFile)
        fileRow = QWidget(sidebar)
        fileLayout = QHBoxLayout(fileRow)
        fileLayout.setContentsMargins(0, 0, 0, 0)
        fileLayout.addWidget(self.fileName)
        fileLayout.addWidget(browse)
        mainForm.addRow('EBSD file', fileRow)
        self.plotType = QComboBox(sidebar)
        self.plotType.addItems(['CI map', 'IPF map', 'Pole figure'])
        self.plotType.currentTextChanged.connect(self._updateMainStyle)
        mainForm.addRow('Plot style', self.plotType)
        self.subplotStyleLabel = QLabel('IPF direction', sidebar)
        self.subplotStyle = QComboBox(sidebar)
        mainForm.addRow(self.subplotStyleLabel, self.subplotStyle)
        self.recreate = QPushButton('Recreate plot', sidebar)
        font = self.recreate.font()
        font.setBold(True)
        self.recreate.setFont(font)
        self.recreate.clicked.connect(self.createPlot)
        mainForm.addRow(self.recreate)
        sidebarLayout.addLayout(mainForm)

        self.tabs = QTabWidget(sidebar)
        self.code = QTextEdit(self.tabs)
        self.code.setReadOnly(True)
        self.code.setPlaceholderText('Generated code appears after a plot is created.')
        codePage = QWidget(self.tabs)
        codeLayout = QVBoxLayout(codePage)
        codeLayout.addWidget(self.code)
        self.copyCode = QPushButton('Copy Python code', codePage)
        self.copyCode.clicked.connect(self.copyPythonCode)
        codeLayout.addWidget(self.copyCode)
        self.tabs.addTab(codePage, 'Code')

        detailsPage = QWidget(self.tabs)
        self.detailsForm = QFormLayout(detailsPage)
        self.ciEnabled = QCheckBox('Mask points below CI', detailsPage)
        self.ciEnabled.toggled.connect(self._updateDetailVisibility)
        self.detailsForm.addRow(self.ciEnabled)
        self.ciThreshold = self._doubleSpin(0.1, 0.0, 1.0, 0.01)
        self.detailsForm.addRow('CI threshold', self.ciThreshold)
        self.cropEnabled = QCheckBox('Show crop only', detailsPage)
        self.cropEnabled.toggled.connect(self._updateDetailVisibility)
        self.detailsForm.addRow(self.cropEnabled)
        self.cropX = RangeSelector(parent=detailsPage)
        self.cropY = RangeSelector(parent=detailsPage)
        self.detailsForm.addRow('Crop X', self.cropX)
        self.detailsForm.addRow('Crop Y', self.cropY)
        self.downsample = QSpinBox(detailsPage)
        self.downsample.setRange(1, 1000)
        self.downsample.setValue(1)
        self.detailsForm.addRow('Preview every nth point', self.downsample)
        self.overlayScale = self._doubleSpin(1.0, 0.01, 1_000_000.0, 0.1)
        self.detailsForm.addRow('Unit-cell scale', self.overlayScale)
        self.clearOverlays = QPushButton('Clear unit-cell overlays', detailsPage)
        self.clearOverlays.clicked.connect(self.clearUnitCellOverlays)
        self.detailsForm.addRow(self.clearOverlays)
        self.scaleBar = QCheckBox('Add scale bar', detailsPage)
        self.detailsForm.addRow(self.scaleBar)
        self.hideAxes = QCheckBox('Hide plot axes', detailsPage)
        self.detailsForm.addRow(self.hideAxes)
        self.tabs.addTab(detailsPage, 'Plot details')
        sidebarLayout.addWidget(self.tabs)
        splitter.addWidget(sidebar)

        self.plotHolder = QWidget(splitter)
        self.plotLayout = QVBoxLayout(self.plotHolder)
        self.plotLayout.setContentsMargins(0, 0, 0, 0)
        self.placeholder = QLabel('Choose an EBSD file to begin.', self.plotHolder)
        self.plotLayout.addWidget(self.placeholder)
        splitter.addWidget(self.plotHolder)
        splitter.setSizes([350, 930])
        self.statusBar().showMessage('Ready')
        self._updateMainStyle()


    @staticmethod
    def _doubleSpin(value: float = 0.0, minimum: float = -1_000_000.0,
                     maximum: float = 1_000_000.0, step: float = 0.1) -> QDoubleSpinBox:
        """Create a spin box for floating point values.

        Args:
           value: initial value
           minimum: minimum value
           maximum: maximum value
           step: step size

        Returns:
           spin box
        """
        spin = QDoubleSpinBox()
        spin.setRange(minimum, maximum)
        spin.setDecimals(4)
        spin.setSingleStep(step)
        spin.setValue(value)
        return spin


    def _updateMainStyle(self) -> None:
        """Update the sub-style choices (IPF direction, pole axis) for the selected plot style."""
        plotType = self.plotType.currentText()
        choices = {'IPF map': ('IPF direction', ['ND', 'RD', 'TD']),
                   'Pole figure': ('Pole axis', ['[1, 0, 0]', '[1, 1, 0]', '[1, 1, 1]'])}
        label, values = choices.get(plotType, ('Map variable', ['CI']))
        self.subplotStyleLabel.setText(label)
        self.subplotStyle.clear()
        self.subplotStyle.addItems(values)
        self.subplotStyle.setVisible(plotType != 'CI map')
        self.subplotStyleLabel.setVisible(plotType != 'CI map')
        self._updateDetailVisibility()


    def _updateDetailVisibility(self) -> None:
        """Show only the plot-detail rows relevant to the selected plot style and options."""
        isMap = self.plotType.currentText() in ('CI map', 'IPF map')
        isIpf = self.plotType.currentText() == 'IPF map'
        self.detailsForm.setRowVisible(self.ciThreshold, self.ciEnabled.isChecked())
        self.detailsForm.setRowVisible(self.cropX, self.cropEnabled.isChecked())
        self.detailsForm.setRowVisible(self.cropY, self.cropEnabled.isChecked())
        self.detailsForm.setRowVisible(self.overlayScale, isIpf)
        self.detailsForm.setRowVisible(self.scaleBar, isMap)
        self.detailsForm.setRowVisible(self.hideAxes, isMap)
        self.clearOverlays.setVisible(isIpf and bool(self.overlays))


    def _loadEbsd(self) -> EBSD:
        """Load the selected file.

        Returns:
           loaded EBSD data
        """
        if self.filePath is None or not self.filePath.is_file():
            raise ValueError('Choose an existing EBSD data file.')
        if self.filePath.suffix.lower() not in LOADERS:
            raise ValueError('Supported formats are ' + ', '.join(LOADERS) + '.')
        self.ebsd = EBSD(str(self.filePath))
        return self.ebsd


    def _applyFilters(self, ebsd: EBSD) -> None:
        """Apply CI mask, preview and crop from the plot details.

        Args:
           ebsd: EBSD data to filter
        """
        ebsd.maskReset()
        ebsd.setVMask(self.downsample.value())
        if self.ciEnabled.isChecked():
            ebsd.maskCI(self.ciThreshold.value())
        if self.cropEnabled.isChecked():
            xmin, xmax = self.cropX.values()
            ymin, ymax = self.cropY.values()
            ebsd.cropVMask(xmin, ymin, xmax, ymax)


    def _setFigure(self, figure: Any) -> None:
        """Show a new figure with its toolbar instead of the previous one.

        Args:
           figure: matplotlib figure
        """
        if self.canvas is not None:
            self.plotLayout.removeWidget(self.canvas)
            self.canvas.setParent(None)
        if self.toolbar is not None:
            self.plotLayout.removeWidget(self.toolbar)
            self.toolbar.setParent(None)
        if self.placeholder is not None:
            self.plotLayout.removeWidget(self.placeholder)
            self.placeholder.setParent(None)
            self.placeholder = None
        self.canvas = FigureCanvasQTAgg(figure)
        self.canvas.mpl_connect('button_press_event', self.addUnitCellOverlay)
        self.toolbar = NavigationToolbar2QT(self.canvas, self.plotHolder)
        self.plotLayout.addWidget(self.toolbar)
        self.plotLayout.addWidget(self.canvas)


    def _refreshCode(self) -> None:
        """Write Python code that reproduces the current plot into the code tab."""
        if self.ebsd is None:
            return
        lines = ['from ebsdlab import EBSD', 'import matplotlib.pyplot as plt', '',
                 f'emap = EBSD({str(self.filePath)!r})',
                 'emap.maskReset()', f'emap.setVMask({self.downsample.value()})']
        if self.ciEnabled.isChecked():
            lines.append(f'emap.maskCI({self.ciThreshold.value()})')
        if self.cropEnabled.isChecked():
            xmin, xmax = self.cropX.values()
            ymin, ymax = self.cropY.values()
            lines.append(f'emap.cropVMask({xmin}, {ymin}, {xmax}, {ymax})')
        plotType = self.plotType.currentText()
        if plotType == 'CI map':
            lines.append('fig = emap.plot(emap.ci, show=False)')
        elif plotType == 'IPF map':
            lines.append(
                f'fig = emap.plotIPF(direction={self.subplotStyle.currentText()!r}, show=False)'
            )
            for x, y, scale in self.overlays:
                lines.append(
                    f'emap.addUnitCellOverlay(fig.axes[0], {x:.6g}, {y:.6g}, scale={scale:.6g})'
                )
        else:
            lines.append(f'fig = emap.plotPF(axis={self.subplotStyle.currentText()}, show=False)')
        if plotType != 'Pole figure':
            lines.append('fig.subplots_adjust(left=0.0125, right=0.99, bottom=0.0125, top=0.99)')
        if self.scaleBar.isChecked() and plotType != 'Pole figure':
            lines.append('emap.addScaleBarOverlay(fig.axes[0])')
        if self.hideAxes.isChecked() and plotType != 'Pole figure':
            lines.append("fig.axes[0].axis('off')")
        lines.extend(['plt.show()', ''])
        self.code.setPlainText('\n'.join(lines))


def main() -> int:
    """Start the optional PySide6 GUI.

    Returns:
       exit code of the application
    """
    app = QApplication.instance() or QApplication([])
    window = EBSDGui()
    window.show()
    return app.exec()
