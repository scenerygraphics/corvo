import pathlib

import pandas as pd
import scanpy as sc
from PyQt5.QtCore import Qt, pyqtSlot, pyqtSignal
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QPushButton, QComboBox, QTabBar, QHBoxLayout, QLabel

from corvolauncher.gui.job_runners.anndata_worker import AnndataWorker
from corvolauncher.gui.job_runners.generic_worker import GenericWorker
from corvolauncher.gui.job_runners.processing_worker import ProcessingWorker


class CheckableComboBox(QComboBox):
    def __init__(self, parent=None):
        super(CheckableComboBox, self).__init__(parent)
        self.view().pressed.connect(self.handleItemPressed)
        self._changed = False

    def handleItemPressed(self, index):
        item = self.model().itemFromIndex(index)
        if item.checkState() == Qt.Checked:
            item.setCheckState(Qt.Unchecked)
        else:
            item.setCheckState(Qt.Checked)
        self._changed = True

    def hidePopup(self):
        if not self._changed:
            super(CheckableComboBox, self).hidePopup()
        self._changed = False

    def itemChecked(self, index):
        item = self.model().item(index, self.modelColumn())
        return item.checkState() == Qt.Checked

    def setItemChecked(self, index, checked=True):
        item = self.model().item(index, self.modelColumn())
        if checked:
            item.setCheckState(Qt.Checked)
        else:
            item.setCheckState(Qt.Unchecked)


class ProcessMenu(QWidget):
    interrupt_signal = pyqtSignal()
    def __init__(self, parent, dataset: str, threadpool):
        super().__init__()

        self.stateTracker = False

        self.parent = parent
        self.dataset = dataset
        self.threadpool = threadpool
        # self.adata = sc.read_h5ad("../resources/datasets/" + self.dataset)
        # self.adata = sc.read_h5ad(os.path.join(str(Path.home()), ".corvo", "resources", "datasets", self.dataset))

        self.layout = QVBoxLayout()
        self.layout.setAlignment(Qt.AlignTop)
        # self.setFixedWidth(self.parent.width() * 2)
        self.setFixedWidth(self.parent.width())
        self.setLayout(self.layout)

        self.process_layout = QHBoxLayout()

        self.process_button = QPushButton("Pre-process")
        self.process_button.setFixedWidth(65)
        self.process_button.setFixedHeight(20)
        self.process_button.clicked.connect(self.handle_launch)
        self.process_layout.addWidget(self.process_button)

        self.interrupt_button = QPushButton("Cancel")
        self.interrupt_button.setFixedWidth(65)
        self.interrupt_button.setFixedHeight(20)
        self.interrupt_button.clicked.connect(self.handle_interrupt)
        self.interrupt_button.hide()  # hidden until worker is launched
        self.process_layout.addWidget(self.interrupt_button)

        # self.process_label = QLabel(self.dataset[:-14].replace("_", " "))
        self.process_label = QLabel(self.dataset[:-14].replace("_", " "))
        self.process_label.setMinimumHeight(30)
        self.process_label.setFixedWidth(self.width() - self.process_button.width())
        self.process_label.setWordWrap(True)
        self.process_layout.addWidget(self.process_label)

        self.layout.addLayout(self.process_layout)

        self.load_adata_info()

    def load_adata_info(self):
        @pyqtSlot()
        def on_running():
            self.layout.addWidget(QLabel("loading dataset info..."))

        @pyqtSlot(pd.DataFrame, pd.DataFrame)
        def load_comboboxes(obs, var):
            self.layout.removeWidget(self.children()[4])  # remove loading label
            self.obs_combo = CheckableComboBox(self)
            c = 0
            for observation in obs:
                self.obs_combo.addItem(observation)
                self.obs_combo.setItemChecked(c, False)
                c += 1
            self.layout.addWidget(self.obs_combo)

            self.var_combo = CheckableComboBox(self)
            c = 0
            for variable in var:
                self.var_combo.addItem(variable)
                self.var_combo.setItemChecked(c, False)
                c += 1
            self.layout.addWidget(self.var_combo)

        worker = AnndataWorker(self.dataset)
        worker.signals.adata_result.connect(load_comboboxes)
        worker.signals.running.connect(on_running)
        self.threadpool.start(worker)

    def handle_interrupt(self):
        self.interrupt_signal.emit()

    def handle_launch(self):
        @pyqtSlot()
        def on_running():
            self.process_button.hide()
            self.interrupt_button.show()

        @pyqtSlot()
        def on_finished():
            self.process_button.show()
            self.interrupt_button.hide()
            self.parent.update_directory_index()

        processing_worker = ProcessingWorker(self.dataset)
        processing_worker.signals.running.connect(on_running)
        processing_worker.signals.finished.connect(on_finished)
        self.interrupt_signal.connect(processing_worker.stop)

        self.threadpool.start(processing_worker)
