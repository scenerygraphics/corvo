import threading
from functools import partial

from PyQt5.QtCore import pyqtSlot, Qt, pyqtSignal, QStringListModel, QSortFilterProxyModel
from PyQt5.QtGui import QFont, QMouseEvent, QStandardItemModel, QStandardItem
from PyQt5.QtWidgets import QVBoxLayout, QLabel, QWidget, QComboBox, QMenuBar, QAction, QHBoxLayout, QPushButton, \
    QListView, QCompleter, QLineEdit, QFrame
from requests import Response

from corvolauncher.gui.job_runners.generic_worker import GenericWorker
from corvolauncher.gui.qt_line_break import QHLineBreakWidget
from corvolauncher.utilities.cellxgene_JSON_requests import CellxGeneJSONRequests


class DatasetFetcher(QWidget):
    download_signal = pyqtSignal()

    def __init__(self, parent, threadpool):
        super().__init__(parent)

        self.parent = parent
        self.threadpool = threadpool

        self.fetcher_layout = QVBoxLayout()
        self.setLayout(self.fetcher_layout)
        self.fetcher_layout.setSpacing(10)
        self.fetcher_layout.setAlignment(Qt.AlignLeft)

        #  give section titles a bold font
        self.title = QLabel("CZCellxGene Hosted Datasets")
        self.title_font = QFont()
        self.title_font.setPointSize(10)
        self.title_font.setBold(True)
        self.title.setFont(self.title_font)
        self.fetcher_layout.addWidget(self.title)

        self.data_fetcher = CellxGeneJSONRequests(self)
        self.collection_map = self.data_fetcher.collections
        self.collection_keys = self.collection_map.keys()
        self.index_list = list(self.collection_keys)
        self.index_list.insert(0, "-- search for a collection --")

        # block changed signal emitted on creation to stop addition to layout that doesn't exist yet
        self.combo_font = QFont()
        self.combo_font.setPointSize(10)

        self.collection_search = QLineEdit("-- search for a collection --")
        self.collection_completer = QCompleter(self.collection_keys, self.collection_search)
        self.collection_completer.setCaseSensitivity(False)
        self.collection_search.setCompleter(self.collection_completer)
        self.collection_search.setFont(self.combo_font)
        self.collection_search.setFixedHeight(25)

        self.collection_search.editingFinished.connect(self.search_box_handler)
        self.fetcher_layout.addWidget(self.collection_search)

        self.collection_container = QComboBox()
        self.collection_container.blockSignals(True)
        self.collection_container.currentIndexChanged.connect(self.load_collection_dsets)
        self.collection_container.setFont(self.combo_font)
        self.collection_container.setFixedHeight(25)

        listView = QListView()
        listView.setWordWrap(True)
        listView.setSpacing(5)
        self.collection_container.setView(listView)
        self.collection_container.addItem("-- select a collection to load its datasets --")

        for c in self.collection_keys:
            self.collection_container.addItem(c)
        self.fetcher_layout.addWidget(self.collection_container)

    @pyqtSlot()
    def search_box_handler(self):
        self.load_collection_dsets(self.index_list.index(self.collection_search.text()))

    @pyqtSlot(int)
    def load_collection_dsets(self, index):
        print("trying to load index: " + str(index))
        # loading datasets still runs blocking and overruns the size of the window if too many datasets exist
        if self.fetcher_layout.count() > 3:
            self.fetcher_layout.removeWidget(self.children()[4])

        if index != 0:
            dataset_map = self.data_fetcher.get_collection_datasets(
                self.collection_map[self.collection_container.itemText(index)])

            dset_button_layout = QVBoxLayout()
            dset_button_frame = QFrame()
            dset_button_layout.setAlignment(Qt.AlignLeft)

            for dset in dataset_map.keys():
                dset_layout = QHBoxLayout()

                dset_download_button = QPushButton("Download")
                dset_download_button.setFixedWidth(65)
                dset_download_button.clicked.connect(partial(self.on_download_click, dset, dataset_map))
                dset_layout.addWidget(dset_download_button)

                dset_label = QLabel(dset)
                dset_label.setFixedWidth(self.collection_container.width() - dset_download_button.width())
                # dset_label.setFixedWidth(self.collection_container.width() - 125)
                # dset_label.setMaximumWidth(self.collection_container.width() - dset_download_button.width())
                dset_label.setWordWrap(True)
                dset_layout.addWidget(dset_label)
                dset_button_layout.addLayout(dset_layout)
            dset_button_frame.setLayout(dset_button_layout)
            self.fetcher_layout.addWidget(dset_button_frame)

    @pyqtSlot(str, dict)
    def on_download_click(self, name, ids):
        print(name)
        print(ids)

        @pyqtSlot()
        def on_finished():
            self.parent.update_directory_index()

        @pyqtSlot(Response)
        def on_result(resp):
            print(resp)
            if len(resp.content) == 0:
                print("An internal server error has occurred. Please try again later or choose a different dataset.")
            else:
                print("saving!")
                print(resp.content)
                # unique_name = self.parent.parent.moniker_recursively(
                #     name.replace(" ", "_") + "_corvo_RAW.h5ad", "datasets")
                # open("../resources/datasets/" + name, "wb").write(resp.content)

        worker = GenericWorker(self.data_fetcher.download_dataset, ids, name)
        # worker.signals.running.connect(on_running)
        worker.signals.finished.connect(on_finished)
        # worker.signals.resp_result.connect(on_result)
        self.threadpool.start(worker)
