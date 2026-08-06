"""Main window for guided EBSD plotting and reproducible code export."""
from __future__ import annotations

from pathlib import Path

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox,
                               QDoubleSpinBox, QFileDialog, QFormLayout,
                               QHBoxLayout, QLabel, QLineEdit, QMainWindow,
                               QMessageBox, QPushButton, QSpinBox, QSplitter,
                               QTabWidget, QTextEdit, QVBoxLayout, QWidget)

from ..ebsd import EBSD, SUPPORTED_SUFFIXES
from .widgets import RangeSelector



class EBSDGui(QMainWindow):
    """A compact UI for the most common ebsdlab plot workflows."""

    def __init__(self):
        super().__init__()
        self.file_path = None
        self.ebsd = None
        self.figure = None
        self.canvas = None
        self.toolbar = None
        self.overlays = []
        self.setWindowTitle('ebsdlab')
        self.resize(1280, 800)
        self._build_ui()

    def _build_ui(self):
        central = QWidget(self)
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        splitter = QSplitter(Qt.Orientation.Horizontal, central)
        root.addWidget(splitter)

        sidebar = QWidget(splitter)
        sidebar_layout = QVBoxLayout(sidebar)
        main_form = QFormLayout()
        self.file_name = QLineEdit(sidebar)
        self.file_name.setReadOnly(True)
        browse = QPushButton('Browse…', sidebar)
        browse.clicked.connect(self.choose_file)
        file_row = QWidget(sidebar)
        file_layout = QHBoxLayout(file_row)
        file_layout.setContentsMargins(0, 0, 0, 0)
        file_layout.addWidget(self.file_name)
        file_layout.addWidget(browse)
        main_form.addRow('EBSD file', file_row)
        self.plot_type = QComboBox(sidebar)
        self.plot_type.addItems(['CI map', 'IPF map', 'Pole figure'])
        self.plot_type.currentTextChanged.connect(self._update_main_style)
        main_form.addRow('Plot style', self.plot_type)
        self.subplot_style_label = QLabel('IPF direction', sidebar)
        self.subplot_style = QComboBox(sidebar)
        main_form.addRow(self.subplot_style_label, self.subplot_style)
        self.recreate = QPushButton('Recreate plot', sidebar)
        font = self.recreate.font()
        font.setBold(True)
        self.recreate.setFont(font)
        self.recreate.clicked.connect(self.create_plot)
        main_form.addRow(self.recreate)
        sidebar_layout.addLayout(main_form)

        self.tabs = QTabWidget(sidebar)
        self.code = QTextEdit(self.tabs)
        self.code.setReadOnly(True)
        self.code.setPlaceholderText('Generated code appears after a plot is created.')
        code_page = QWidget(self.tabs)
        code_layout = QVBoxLayout(code_page)
        code_layout.addWidget(self.code)
        self.copy_code = QPushButton('Copy Python code', code_page)
        self.copy_code.clicked.connect(self.copy_python_code)
        code_layout.addWidget(self.copy_code)
        self.tabs.addTab(code_page, 'Code')

        details_page = QWidget(self.tabs)
        self.details_form = QFormLayout(details_page)
        self.ci_enabled = QCheckBox('Mask points below CI', details_page)
        self.ci_enabled.toggled.connect(self._update_detail_visibility)
        self.details_form.addRow(self.ci_enabled)
        self.ci_threshold = self._double_spin(0.1, 0.0, 1.0, 0.01)
        self.details_form.addRow('CI threshold', self.ci_threshold)
        self.crop_enabled = QCheckBox('Show crop only', details_page)
        self.crop_enabled.toggled.connect(self._update_detail_visibility)
        self.details_form.addRow(self.crop_enabled)
        self.crop_x = RangeSelector(parent=details_page)
        self.crop_y = RangeSelector(parent=details_page)
        self.details_form.addRow('Crop X', self.crop_x)
        self.details_form.addRow('Crop Y', self.crop_y)
        self.downsample = QSpinBox(details_page)
        self.downsample.setRange(1, 1000)
        self.downsample.setValue(1)
        self.details_form.addRow('Preview every nth point', self.downsample)
        self.overlay_scale = self._double_spin(1.0, 0.01, 1_000_000.0, 0.1)
        self.details_form.addRow('Unit-cell scale', self.overlay_scale)
        self.clear_overlays = QPushButton('Clear unit-cell overlays', details_page)
        self.clear_overlays.clicked.connect(self.clear_unit_cell_overlays)
        self.details_form.addRow(self.clear_overlays)
        self.scale_bar = QCheckBox('Add scale bar', details_page)
        self.details_form.addRow(self.scale_bar)
        self.hide_axes = QCheckBox('Hide plot axes', details_page)
        self.details_form.addRow(self.hide_axes)
        self.tabs.addTab(details_page, 'Plot details')
        sidebar_layout.addWidget(self.tabs)
        splitter.addWidget(sidebar)

        self.plot_holder = QWidget(splitter)
        self.plot_layout = QVBoxLayout(self.plot_holder)
        self.plot_layout.setContentsMargins(0, 0, 0, 0)
        self.placeholder = QLabel('Choose an EBSD file to begin.', self.plot_holder)
        self.plot_layout.addWidget(self.placeholder)
        splitter.addWidget(self.plot_holder)
        splitter.setSizes([350, 930])
        self.statusBar().showMessage('Ready')
        self._update_main_style()

    @staticmethod
    def _double_spin(value=0.0, minimum=-1_000_000.0,
                     maximum=1_000_000.0, step=0.1):
        spin = QDoubleSpinBox()
        spin.setRange(minimum, maximum)
        spin.setDecimals(4)
        spin.setSingleStep(step)
        spin.setValue(value)
        return spin

    def _update_main_style(self):
        plot_type = self.plot_type.currentText()
        choices = {'IPF map': ('IPF direction', ['ND', 'RD', 'TD']),
                   'Pole figure': ('Pole axis', ['[1, 0, 0]', '[1, 1, 0]', '[1, 1, 1]'])}
        label, values = choices.get(plot_type, ('Map variable', ['CI']))
        self.subplot_style_label.setText(label)
        self.subplot_style.clear()
        self.subplot_style.addItems(values)
        self.subplot_style.setVisible(plot_type != 'CI map')
        self.subplot_style_label.setVisible(plot_type != 'CI map')
        self._update_detail_visibility()

    def _update_detail_visibility(self):
        is_map = self.plot_type.currentText() in ('CI map', 'IPF map')
        is_ipf = self.plot_type.currentText() == 'IPF map'
        self.details_form.setRowVisible(self.ci_threshold, self.ci_enabled.isChecked())
        self.details_form.setRowVisible(self.crop_x, self.crop_enabled.isChecked())
        self.details_form.setRowVisible(self.crop_y, self.crop_enabled.isChecked())
        self.details_form.setRowVisible(self.overlay_scale, is_ipf)
        self.details_form.setRowVisible(self.scale_bar, is_map)
        self.details_form.setRowVisible(self.hide_axes, is_map)
        self.clear_overlays.setVisible(is_ipf and bool(self.overlays))

    def set_file(self, file_path):
        """Select an input file while displaying only its filename in the UI."""
        self.file_path = Path(file_path).expanduser().resolve()
        self.file_name.setText(self.file_path.name)
        self.file_name.setToolTip(str(self.file_path))

    def choose_file(self):
        name, _ = QFileDialog.getOpenFileName(
            self, 'Open EBSD file', '', 'EBSD files (*.ang *.osc *.txt *.crc)')
        if name:
            self.set_file(name)
            self.create_plot()

    def _load_ebsd(self):
        if self.file_path is None or not self.file_path.is_file():
            raise ValueError('Choose an existing EBSD data file.')
        if self.file_path.suffix.lower() not in SUPPORTED_SUFFIXES:
            raise ValueError('Supported formats are .ang, .osc, .txt, and .crc.')
        self.ebsd = EBSD(str(self.file_path))

    def _apply_filters(self):
        self.ebsd.maskReset()
        self.ebsd.setVMask(self.downsample.value())
        if self.ci_enabled.isChecked():
            self.ebsd.maskCI(self.ci_threshold.value())
        if self.crop_enabled.isChecked():
            xmin, xmax = self.crop_x.values()
            ymin, ymax = self.crop_y.values()
            self.ebsd.cropVMask(xmin, ymin, xmax, ymax)

    def create_plot(self):
        try:
            self._load_ebsd()
            self._apply_filters()
            self.overlays = []
            plot_type = self.plot_type.currentText()
            if plot_type == 'CI map':
                self.figure = self.ebsd.plot(self.ebsd.CI, show=False)
            elif plot_type == 'IPF map':
                self.figure = self.ebsd.plotIPF(
                    direction=self.subplot_style.currentText(), show=False)
            else:
                self.figure = self.ebsd.plotPF(
                    axis=self._selected_pole_axis(), show=False)
            axis = self.figure.axes[0]
            if plot_type != 'Pole figure':
                self.figure.subplots_adjust(left=0.0125, right=0.99,
                                             bottom=0.0125, top=0.99)
            if self.hide_axes.isChecked() and plot_type != 'Pole figure':
                axis.axis('off')
            if self.scale_bar.isChecked() and plot_type != 'Pole figure':
                self.ebsd.addScaleBarOverlay(axis)
            self._set_figure(self.figure)
            self._refresh_code()
            self._update_detail_visibility()
            self.statusBar().showMessage('Plot created')
        except (OSError, ValueError, IndexError) as error:
            QMessageBox.critical(self, 'Could not create plot', str(error))
            self.statusBar().showMessage('Plot creation failed')

    def _set_figure(self, figure):
        if self.canvas is not None:
            self.plot_layout.removeWidget(self.canvas)
            self.canvas.setParent(None)
        if self.toolbar is not None:
            self.plot_layout.removeWidget(self.toolbar)
            self.toolbar.setParent(None)
        if self.placeholder is not None:
            self.plot_layout.removeWidget(self.placeholder)
            self.placeholder.setParent(None)
            self.placeholder = None
        self.canvas = FigureCanvasQTAgg(figure)
        self.canvas.mpl_connect('button_press_event', self.add_unit_cell_overlay)
        self.toolbar = NavigationToolbar2QT(self.canvas, self.plot_holder)
        self.plot_layout.addWidget(self.toolbar)
        self.plot_layout.addWidget(self.canvas)

    def add_unit_cell_overlay(self, event):
        if (self.ebsd is None or self.plot_type.currentText() != 'IPF map'
                or event.inaxes is None or event.xdata is None or event.ydata is None):
            return
        scale = self.overlay_scale.value()
        self.ebsd.addUnitCellOverlay(event.inaxes, event.xdata, event.ydata, scale)
        self.overlays.append((event.xdata, event.ydata, scale))
        self.canvas.draw_idle()
        self._refresh_code()
        self._update_detail_visibility()
        self.statusBar().showMessage('Unit-cell overlay added')

    def clear_unit_cell_overlays(self):
        if self.figure is not None and self.overlays:
            self.create_plot()

    def _selected_pole_axis(self):
        return [int(value.strip()) for value in
                self.subplot_style.currentText().strip('[]').split(',')]

    def _refresh_code(self):
        if self.ebsd is None:
            return
        lines = ['from ebsdlab import EBSD', 'import matplotlib.pyplot as plt', '',
                 'emap = EBSD({!r})'.format(str(self.file_path)),
                 'emap.maskReset()', 'emap.setVMask({})'.format(self.downsample.value())]
        if self.ci_enabled.isChecked():
            lines.append('emap.maskCI({})'.format(self.ci_threshold.value()))
        if self.crop_enabled.isChecked():
            xmin, xmax = self.crop_x.values()
            ymin, ymax = self.crop_y.values()
            lines.append('emap.cropVMask({}, {}, {}, {})'.format(xmin, ymin, xmax, ymax))
        plot_type = self.plot_type.currentText()
        if plot_type == 'CI map':
            lines.append('fig = emap.plot(emap.CI, show=False)')
        elif plot_type == 'IPF map':
            lines.append("fig = emap.plotIPF(direction={!r}, show=False)".format(
                self.subplot_style.currentText()))
            for x, y, scale in self.overlays:
                lines.append('emap.addUnitCellOverlay(fig.axes[0], {:.6g}, {:.6g}, scale={:.6g})'.format(
                    x, y, scale))
        else:
            lines.append('fig = emap.plotPF(axis={}, show=False)'.format(self._selected_pole_axis()))
        if plot_type != 'Pole figure':
            lines.append('fig.subplots_adjust(left=0.0125, right=0.99, bottom=0.0125, top=0.99)')
        if self.scale_bar.isChecked() and plot_type != 'Pole figure':
            lines.append('emap.addScaleBarOverlay(fig.axes[0])')
        if self.hide_axes.isChecked() and plot_type != 'Pole figure':
            lines.append("fig.axes[0].axis('off')")
        lines.extend(['plt.show()', ''])
        self.code.setPlainText('\n'.join(lines))

    def copy_python_code(self):
        if self.code.toPlainText():
            QApplication.clipboard().setText(self.code.toPlainText())
            self.statusBar().showMessage('Python code copied to clipboard')


def main():
    """Start the optional PySide6 GUI."""
    app = QApplication.instance() or QApplication([])
    window = EBSDGui()
    window.show()
    return app.exec()
