import pathlib

import scanpy as sc
from PyQt5.QtCore import Qt, pyqtSlot, pyqtSignal
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QPushButton, QComboBox, QTabBar

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
        self.adata = sc.read_h5ad("../resources/datasets/" + self.dataset)

        self.layout = QVBoxLayout()
        self.layout.setAlignment(Qt.AlignTop)
        # self.setFixedWidth(self.parent.width() * 2)
        self.setFixedWidth(self.parent.width())
        self.setLayout(self.layout)

        self.process_button = QPushButton("Pre-process " + self.dataset)
        self.process_button.clicked.connect(self.handle_launch)
        self.layout.addWidget(self.process_button)

        self.interrupt_button = QPushButton("Cancel pre-processing of " + self.dataset)
        self.interrupt_button.clicked.connect(self.handle_interrupt)
        self.interrupt_button.hide()  # hidden until worker is launched
        self.layout.addWidget(self.interrupt_button)

        self.obs_combo = CheckableComboBox(self)
        c = 0
        for observation in self.adata.obs:
            self.obs_combo.addItem(observation)
            self.obs_combo.setItemChecked(c, False)
            c += 1
        self.layout.addWidget(self.obs_combo)

        self.var_combo = CheckableComboBox(self)
        c = 0
        for variable in self.adata.var:
            self.var_combo.addItem(variable)
            self.var_combo.setItemChecked(c, False)
            if variable == "feature_name":
                print(self.adata.var[variable])
            c += 1

        self.layout.addWidget(self.var_combo)

    def handle_interrupt(self):
        self.interrupt_signal.emit()
        pass
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