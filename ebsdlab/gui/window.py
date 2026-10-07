"""Main window: load, analyse, clean and plot an EBSD map; the generated Python code reproduces the plot."""
from __future__ import annotations
import math
import re
from pathlib import Path
from typing import Any
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout,
                               QHBoxLayout, QLabel, QLineEdit, QMainWindow, QMessageBox, QPushButton, QSpinBox,
                               QSplitter, QTextEdit, QVBoxLayout, QWidget)
from ..ebsd import EBSD
from ..fileIO import LOADERS

PREVIEW_POINTS = 250_000  # preview: pole figures, the slowest plots, need about 2 µs per point
IMPORTS = {'np': 'import numpy as np', 'plt': 'import matplotlib.pyplot as plt',
           'colormaps': 'from matplotlib import colormaps',
           'ListedColormap': 'from matplotlib.colors import ListedColormap',
           'Patch': 'from matplotlib.patches import Patch', 'EBSD': 'from ebsdlab.ebsd import EBSD'}
CHOICES = {'IPF': ['ND', 'RD', 'TD'], 'PF': ['100', '110', '111'], 'KAM': ['2°', '5°', '10°'],
           'Others': ['grains', 'phase', 'IQ', 'CI']}
UNIT_CELLS = {'no unit cells': 0., 'small unit cells': 0.03, 'medium unit cells': 0.06,
              'large unit cells': 0.12}  # size as fraction of the map width


def imports(lines: list[str]) -> list[str]:
    """Import lines for the names that the code uses.

    Args:
       lines: code lines; comments are ignored

    Returns:
       import lines
    """
    code = '\n'.join(line for line in lines if not line.lstrip().startswith('#'))
    return [line for name, line in IMPORTS.items() if re.search(rf'\b{name}\b', code)]


class EBSDGui(QMainWindow):
    """Load, analyse, clean and plot: few choices per step, every change replots."""

    def __init__(self) -> None:
        """Create the window with an empty plot area."""
        super().__init__()
        self.filePath: Path | None                 = None
        self.ebsd:     EBSD | None                 = None
        self.grains:   tuple[float, int] | None    = None   # tolerance and minimum size of the calculated grains
        self.dilated                               = False
        self.cells:    list[tuple[float, float]]   = []     # unit-cell positions in [um]
        self.second                                = ''     # last selected second choice, e.g. unit-cell size
        self.home:     tuple[float, ...] | None    = None   # limits of the whole map: xmin, xmax, ymin, ymax
        self.figure:   Any                         = None
        self.canvas:   FigureCanvasQTAgg | None    = None
        self.toolbar:  NavigationToolbar2QT | None = None
        self.setWindowTitle('ebsdlab')
        self.resize(1280, 800)
        self._buildUi()


    def setFile(self, filePath: str | Path) -> None:
        """Load an EBSD file and plot it.

        Args:
           filePath: path of the EBSD file
        """
        self.filePath = Path(filePath).expanduser().resolve()
        self.fileName.setText(self.filePath.name)
        self.fileName.setToolTip(str(self.filePath))
        try:
            self._load()
        except (OSError, ValueError, IndexError, KeyError) as error:
            self.ebsd = None
            QMessageBox.critical(self, 'Could not load file', str(error))
            return
        self.cells, self.home = [], None
        self._updateChoices()


    def chooseFile(self) -> None:
        """Prompt for an EBSD file and plot it."""
        name, _ = QFileDialog.getOpenFileName(self, 'Open EBSD file', '',
                                              'EBSD files (' + ' '.join('*'+suffix for suffix in LOADERS) + ')')
        if name:
            self.setFile(name)


    def replot(self, full: bool = False) -> None:
        """Analyse and clean as needed, then run the plot code: a preview, or the full detail of the view.

        Args:
           full: all points, cropped to the view; otherwise a preview of the whole map
        """
        if self.ebsd is None:
            return
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            self._process()
            _, mask, crop, plot = self._code()
            every = 1 if full else max(1, math.ceil(self.ebsd.nPoints/PREVIEW_POINTS))
            self.ebsd.maskReset()
            self.ebsd.setVMask(every)
            namespace = {'emap': self.ebsd}
            code = mask + (crop if full else []) + plot
            exec('\n'.join(imports(code) + code), namespace)  # pylint: disable=exec-used
            self._setFigure(namespace['fig'], full)
            self._refreshCode()
            self.statusBar().showMessage('Full detail' if every == 1 else
                                         f'Preview with every {every}th point; render full detail for all')
        except (OSError, ValueError, IndexError) as error:
            self.statusBar().showMessage(f'Plot failed: {error}')
        finally:
            QApplication.restoreOverrideCursor()


    def addUnitCell(self, event: Any) -> None:
        """Add a unit cell at a click into an IPF or grain map.

        Args:
           event: matplotlib mouse event
        """
        scale = self._cellScale()  # 0 unless an IPF map with unit cells
        if not scale or self.ebsd is None or event.inaxes is None or (self.toolbar and self.toolbar.mode):
            return
        self.cells.append((event.xdata, event.ydata))
        self.ebsd.addUnitCellOverlay(event.inaxes, event.xdata, event.ydata, scale)
        event.canvas.draw_idle()
        self._refreshCode()


    def copyPythonCode(self) -> None:
        """Copy the generated code to the clipboard."""
        if self.code.toPlainText():
            QApplication.clipboard().setText(self.code.toPlainText())
            self.statusBar().showMessage('Python code copied to clipboard')


    def _buildUi(self) -> None:
        """Create the sidebar (load, analyse, clean, plot, code) and the plot area."""
        splitter = QSplitter(Qt.Orientation.Horizontal, self)
        self.setCentralWidget(splitter)
        sidebar = QWidget(splitter)
        sidebarLayout = QVBoxLayout(sidebar)
        form = QFormLayout()
        sidebarLayout.addLayout(form)

        form.addRow(self._heading('Load', top=0))
        self.fileName = QLineEdit(sidebar)
        self.fileName.setReadOnly(True)
        browse = QPushButton('Browse…', sidebar)
        browse.clicked.connect(self.chooseFile)
        fileRow = QWidget(sidebar)
        fileLayout = QHBoxLayout(fileRow)
        fileLayout.setContentsMargins(0, 0, 0, 0)
        fileLayout.addWidget(self.fileName)
        fileLayout.addWidget(browse)
        form.addRow('File', fileRow)

        form.addRow(self._heading('Analyse'))
        self.tolerance = QDoubleSpinBox(sidebar)
        self.tolerance.setRange(0.1, 60.)
        self.tolerance.setSingleStep(0.5)
        self.tolerance.setValue(5.)
        self.tolerance.setSuffix(' °')
        form.addRow('Grain tolerance', self.tolerance)
        self.minSize = QSpinBox(sidebar)
        self.minSize.setRange(1, 1_000_000)
        self.minSize.setValue(6)
        form.addRow('Minimum grain size', self.minSize)

        form.addRow(self._heading('Clean'))
        self.minCI = QDoubleSpinBox(sidebar)
        self.minCI.setRange(0., 1.)
        self.minCI.setSingleStep(0.05)
        self.minCI.setSpecialValueText('off')
        form.addRow('Minimum CI', self.minCI)
        self.dilation = QCheckBox('Grain dilation', sidebar)
        form.addRow('', self.dilation)

        form.addRow(self._heading('Plot'))
        self.showList = QComboBox(sidebar)
        self.showList.addItems(list(CHOICES))
        form.addRow('Show', self.showList)
        self.choice1, self.choice2 = QComboBox(sidebar), QComboBox(sidebar)
        form.addRow('', self.choice1)
        form.addRow('', self.choice2)
        self.form = form
        fullDetail = QPushButton('Render full detail', sidebar)
        fullDetail.clicked.connect(lambda: self.replot(full=True))
        form.addRow('', fullDetail)

        sidebarLayout.addSpacing(12)
        self.code = QTextEdit(sidebar)
        self.code.setReadOnly(True)
        self.code.setPlaceholderText('Generated code appears after a plot is created.')
        sidebarLayout.addWidget(self.code)
        copyCode = QPushButton('Copy Python code', sidebar)
        copyCode.clicked.connect(self.copyPythonCode)
        sidebarLayout.addWidget(copyCode)

        # every change replots; spin boxes only after editing, not at each key
        for spin in (self.tolerance, self.minSize, self.minCI):
            spin.setKeyboardTracking(False)
            spin.valueChanged.connect(lambda _: self.replot())
        self.dilation.toggled.connect(lambda _: self.replot())
        self.showList.currentTextChanged.connect(self._updateChoices)
        self.choice1.currentIndexChanged.connect(lambda _: self._updateSecond())
        self.choice2.currentIndexChanged.connect(lambda _: self.replot())

        self.plotHolder = QWidget(splitter)
        self.plotLayout = QVBoxLayout(self.plotHolder)
        self.plotLayout.setContentsMargins(0, 0, 0, 0)
        splitter.setSizes([380, 900])
        self.info = QLabel('Choose an EBSD file to begin.', self)  # file details, right of the messages
        self.statusBar().addPermanentWidget(self.info)
        self.statusBar().showMessage('Ready')
        self._updateChoices()


    @staticmethod
    def _heading(text: str, top: int = 12) -> QLabel:
        """Bold label for a workflow step.

        Args:
           text: name of the step
           top: space above in pixels

        Returns:
           label
        """
        label = QLabel(text)
        font = label.font()
        font.setBold(True)
        label.setFont(font)
        label.setContentsMargins(0, top, 0, 0)
        return label


    def _load(self) -> None:
        """Load the file; forget grains and dilation."""
        if self.filePath is None:
            return
        self.ebsd = EBSD(str(self.filePath))
        self.grains, self.dilated = None, False
        e = self.ebsd
        phases = ', '.join(sym.lattice or 'unknown' for sym in e.sym[1:]) or 'none'
        self.info.setText(f'{self.filePath.suffix[1:]}: {e.nPoints:,} points, {e.grid}, step {e.stepSizeX:g} µm, '
                          f'phases: {phases}')


    def _updateChoices(self) -> None:
        """Fill the first choice for the selected 'Show' entry, then the second one."""
        show = self.showList.currentText()
        self.choice1.blockSignals(True)
        self.choice1.clear()
        self.choice1.addItems(CHOICES[show])
        if show == 'KAM':
            self.choice1.setCurrentText('5°')
        self.choice1.blockSignals(False)
        self._updateSecond()


    def _updateSecond(self) -> None:
        """Fill the second choice for the 'Show' entry and first choice, keep its selection if possible; replot."""
        show = self.showList.currentText()
        second: list[tuple[str, Any]]
        if show == 'IPF' or self.choice1.currentText() == 'grains':
            second = list(UNIT_CELLS.items())
        elif show == 'PF':
            phases = [(k, sym.lattice) for k, sym in enumerate(self.ebsd.sym) if sym.lattice] if self.ebsd else []
            second = [(style, (style == 'points', None)) for style in ('points', 'density')] if len(phases) < 2 else \
                     [(f'{style} - phase {k} ({lattice})', (style == 'points', k))
                      for style in ('points', 'density') for k, lattice in phases]
        else:
            second = []
        self.second = self.choice2.currentText() or self.second
        self.choice2.blockSignals(True)
        self.choice2.clear()
        for text, data in second:
            self.choice2.addItem(text, data)
        self.choice2.setCurrentText(self.second)
        self.choice2.blockSignals(False)
        self.form.setRowVisible(self.choice2, bool(second))
        self.replot()


    def _cellScale(self) -> float:
        """Unit-cell size in [um] for IPF and grain maps; 0 without unit cells.

        Returns:
           scale of the unit cells
        """
        if self.ebsd is None or not (self.showList.currentText() == 'IPF' or self.choice1.currentText() == 'grains'):
            return 0.
        return (self.choice2.currentData() or 0.) * self.ebsd.width


    def _process(self) -> None:
        """Calculate grains, dilation and KAM if needed; keep them until the file or the settings change."""
        grains = (self.tolerance.value(), self.minSize.value())
        dilate = self.dilation.isChecked()
        if self.dilated and (not dilate or self.grains != grains):  # undo the dilation
            self._load()
        e = self.ebsd
        assert e is not None
        if (dilate or self.choice1.currentText() == 'grains') and self.grains != grains:
            e.calcGrains(*grains)
            self.grains = grains
        if dilate and not self.dilated:
            e.grainDilation()
            e.kam, self.dilated = np.empty(0), True
        if self.showList.currentText() == 'KAM' and not e.kam.size:
            e.calcKAM()


    def _view(self) -> tuple[float, ...] | None:
        """Visible area of a zoomed map.

        Returns:
           xmin, xmax, ymin, ymax in [um]; None for the whole map or a pole figure
        """
        if self.figure is None or self.home is None or self.showList.currentText() == 'PF':
            return None
        ax = self.figure.axes[0]
        view = (*sorted(ax.get_xlim()), *sorted(ax.get_ylim()))
        return None if np.allclose(view, self.home) else view


    def _code(self) -> tuple[list[str], list[str], list[str], list[str]]:
        """Python code for the current choices: analysis, mask, crop and plot lines.

        Returns:
           lists of lines
        """
        e = self.ebsd
        assert e is not None
        show, choice1, choice2 = self.showList.currentText(), self.choice1.currentText(), self.choice2.currentData()
        analysis, mask, crop = [], [], []
        if self.dilation.isChecked() or choice1 == 'grains':
            analysis += ['# grains: neighbors misoriented less than tolerance [°] are linked; at least minSize points',
                         f'emap.calcGrains(tolerance={self.tolerance.value():g}, minSize={self.minSize.value()})']
        if self.dilation.isChecked():
            analysis += ['# grain dilation: fill points without grain from their neighbors (in place)',
                         'emap.grainDilation()']
        if show == 'KAM':
            analysis += ['# kernel average misorientation [°]: mean misorientation to the nearest neighbors',
                         'emap.calcKAM()']
        if self.minCI.value() > 0:
            mask += [f'# hide points with a confidence index of {self.minCI.value():g} or less (black)',
                     f'emap.maskCI({self.minCI.value():g})']
        if show == 'PF' and choice2 and choice2[1] is not None:
            mask += [f'# only phase {choice2[1]}', f'emap.mask &= emap.phaseID == {choice2[1]}']
        if view := self._view():
            crop += ['# crop to the zoomed view: xmin, ymin, xmax, ymax [µm]',
                     f'emap.cropVMask({view[0]:.4g}, {view[2]:.4g}, {view[1]:.4g}, {view[3]:.4g})']
        if show == 'IPF':
            plot = [f'# inverse pole figure map: color shows the crystal direction along {choice1}',
                    f"fig = emap.plotIPF('{choice1}', show=False)"]
        elif show == 'PF':
            plot = [f'# {{{choice1}}} pole figure, oriented as the map: X right, Y down',
                    f'fig = emap.plotPF(axis={[int(i) for i in choice1]}, points={choice2[0]}, show=False)']
        elif show == 'KAM':
            plot = ['# KAM map [°]; unindexed points are black',
                    f"fig = emap.plot(emap.kam, vmax={choice1.rstrip('°')}, show=False)"]
        elif choice1 == 'grains':
            plot = ['# grain map: 20 colors repeat; points without grain are black',
                    "cmap = colormaps['tab20'].with_extremes(bad='k')",
                    'fig = emap.plot(np.where(emap.grainID > 0, emap.grainID % 20, np.nan), cmap=cmap, cbar=False, '
                    'show=False)']
        elif choice1 == 'phase':
            labels = [f'{k}: {sym.lattice or "unknown"}' for k, sym in enumerate(e.sym[1:], 1)] or ['1: unknown']
            plot = ['# phase map: one color per phase; not identified points are black',
                    f"cmap = ListedColormap(colormaps['tab10'].colors[:{len(labels)}]).with_extremes(bad='k')",
                    'fig = emap.plot(np.where(emap.phaseID > 0, emap.phaseID, np.nan), cmap=cmap, vmin=0.5, '
                    f'vmax={len(labels)}.5, cbar=False, show=False)',
                    'fig.axes[0].legend(handles=[Patch(color=cmap(i), label=label) for i, label in '
                    f'enumerate({labels!r})])']
        elif choice1 == 'IQ':
            vector = 'bc' if e.bc.size and not e.iq.any() else 'iq'  # Oxford files: band contrast
            plot = [f'# {"band contrast" if vector == "bc" else "image quality"} map',
                    f'fig = emap.plot(emap.{vector}, show=False)']
        else:
            plot = ['# confidence index map', 'fig = emap.plot(emap.ci, show=False)']
        if scale := self._cellScale():
            plot += ['# unit cells at the clicked points, seen from above the sample'] + \
                    [f'emap.addUnitCellOverlay(fig.axes[0], {x:.4g}, {y:.4g}, scale={scale:.3g})'
                     for x, y in self.cells]
        if show != 'PF':
            plot += ['# map without axes, with scale bar', "fig.axes[0].axis('off')",
                     'emap.addScaleBarOverlay(fig.axes[0])']
        return analysis, mask, crop, plot


    def _refreshCode(self) -> None:
        """Write the Python code that reproduces the plot in full detail."""
        if self.ebsd is None or self.filePath is None:
            return
        analysis, mask, crop, plot = self._code()
        name = '_'.join([self.filePath.stem, self.showList.currentText(), self.choice1.currentText().rstrip('°')])
        lines = ['# load the map; phase 0 = not identified, phases are numbered 1, 2, ...',
                 f'emap = EBSD({str(self.filePath)!r})'] + analysis + mask + crop + plot + \
                ['# save the figure', f"plt.savefig('{name}.png', dpi=300, bbox_inches='tight')", '']
        self.code.setPlainText('\n'.join(['# EBSD plot generated by ebsdlab-gui'] + imports(lines) + [''] + lines))


    def _setFigure(self, figure: Any, full: bool) -> None:
        """Show a new figure with its toolbar; a map preview keeps the zoomed view of the previous map.

        Args:
           figure: matplotlib figure
           full: full detail, cropped to the view
        """
        view = self._view()
        for widget in (self.canvas, self.toolbar):
            if widget is not None:
                self.plotLayout.removeWidget(widget)
                widget.setParent(None)
        if self.figure is not None:
            plt.close(self.figure)
        self.figure = figure
        self.canvas = FigureCanvasQTAgg(figure)
        self.canvas.mpl_connect('button_press_event', self.addUnitCell)
        self.toolbar = NavigationToolbar2QT(self.canvas, self.plotHolder)
        self.plotLayout.addWidget(self.toolbar)
        self.plotLayout.addWidget(self.canvas)
        if self.showList.currentText() == 'PF':
            self.home = None
            return
        ax = figure.axes[0]
        if not full:
            self.home = (*sorted(ax.get_xlim()), *sorted(ax.get_ylim()))
            if view:
                self.toolbar.push_current()  # home button: whole map
                ax.set_xlim(view[0], view[1])
                ax.set_ylim(view[3], view[2])  # y down
                self.toolbar.push_current()
        for name in ('xlim_changed', 'ylim_changed'):  # the crop in the code follows the view
            ax.callbacks.connect(name, lambda _: self._refreshCode())


def main() -> int:
    """Start the optional PySide6 GUI.

    Returns:
       exit code of the application
    """
    app = QApplication.instance() or QApplication([])
    window = EBSDGui()
    window.show()
    return app.exec()
