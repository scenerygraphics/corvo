import threading
from functools import partial

from PyQt5.QtCore import pyqtSlot, Qt, pyqtSignal
from PyQt5.QtGui import QFont, QMouseEvent
from PyQt5.QtWidgets import QVBoxLayout, QLabel, QWidget, QComboBox, QMenuBar, QAction, QHBoxLayout, QPushButton

from corvolauncher.gui.job_runners.generic_worker import GenericWorker
from corvolauncher.utilities.cellxgene_JSON_requests import CellxGeneJSONRequests


class DatasetFetcher(QWidget):
    download_signal = pyqtSignal()
    def __init__(self, parent, threadpool):
        super().__init__(parent)

        self.parent = parent
        self.threadpool = threadpool

        self.fetcher_layout = QVBoxLayout()
        self.setLayout(self.fetcher_layout)
        self.fetcher_layout.setSpacing(5)
        self.fetcher_layout.setAlignment(Qt.AlignLeft)

        self.combobox_layout = QHBoxLayout()
        self.combobox_layout.setAlignment(Qt.AlignLeft)
        self.combobox_layout.setSpacing(2)

        #  give section titles a bold font
        self.title = QLabel("CZCellxGene Hosted Datasets")
        self.title_font = QFont()
        self.title_font.setBold(True)
        self.title.setFont(self.title_font)
        self.fetcher_layout.addWidget(self.title)

        self.data_fetcher = CellxGeneJSONRequests()
        self.collection_map = self.data_fetcher.collections

        self.collection_label = QLabel("Collections:")

        self.collection_container = QComboBox()
        self.collection_container.blockSignals(True)
        self.collection_container.currentIndexChanged.connect(self.load_collection_dsets)

        for c in self.collection_map.keys():
            self.collection_container.addItem(c)

        self.combobox_layout.addWidget(self.collection_label)
        self.combobox_layout.addWidget(self.collection_container)

        self.fetcher_layout.addLayout(self.combobox_layout)

        self.fetcher_layout.addWidget(QLabel("Double-click to download:"))

        # self.collection_bar = QMenuBar()
        # self.file = self.collection_bar.addMenu("File")
        # self.file.addAction("New")
        #
        # self.save = QAction("Save", self)
        # self.save.setShortcut("Ctrl+S")
        # self.file.addAction(self.save)
        #
        # self.edit = self.file.addMenu("Edit")
        # self.edit.addAction("copy")
        # self.edit.addAction("paste")
        #
        # self.quit = QAction("Quit", self)
        # self.file.addAction(self.quit)
        # self.file.triggered[QAction].connect(self.processtrigger)
        #
        # self.layout.addWidget(self.collection_bar)

    def processtrigger(self, q):
        print(q.text() + " is triggered")

    @pyqtSlot(int)
    def load_collection_dsets(self, index):
        dset_button_layout = QVBoxLayout()
        dset_button_layout.setAlignment(Qt.AlignLeft)
        dataset_map = self.data_fetcher.get_collection_datasets(self.collection_map[self.collection_container.itemText(index)])

        for dset in dataset_map.keys():
            dset_button = self.DoubleClickButton(self, dset, dataset_map, self.data_fetcher, self.threadpool)
            dset_button_layout.addWidget(dset_button)
        self.fetcher_layout.addLayout(dset_button_layout)

    @pyqtSlot(dict)
    def fetch_dataset(self, dset_dict):
        print(type(dset_dict))

    class DoubleClickButton(QPushButton):
        def __init__(self, parent, text, map, fetcher, threadpool):
            super().__init__()
            self.parent = parent
            self.setText(text)
            self.map = map
            self.fetcher = fetcher
            self.threadpool = threadpool

        def mouseDoubleClickEvent(self, event):
            print("pos: ", event.pos())
            print("downloading")
            print(self.text())
            worker = GenericWorker(self.fetcher.download_dataset, self.map, self.text())
            worker.signals.running.connect(self.on_running)
            worker.signals.finished.connect(self.on_finished)
            self.threadpool.start(worker)

        @pyqtSlot()
        def on_running(self):
            # self.parent.fetcher_layout.removeWidget(self.parent.fetcher_layout.itemAt(3).widget())
            pass

        @pyqtSlot()
        def on_finished(self):
            self.parent.parent.update_directory_index()

